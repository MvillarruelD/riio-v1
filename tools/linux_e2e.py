"""Required native-Linux integration with real engines and online RefSeq context.

The shipped K-12 fragment certifies discovery. It lacks zntR's distal regulatory
targets and homolog promoter data, so complete prediction uses the unmodified
full GCF_000005845.2 assembly downloaded by the normal genome resolver.
No mocked engine, manufactured operator, or unavailable-network skip is used.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # The checkout-only tools/ validator is not an installed package.
ACCESSION = "GCF_000005845.2"
FIXTURES = ROOT / "tests" / "data"


def cli(out: Path, timeout: int, *arguments: str) -> None:
    """Keep logs and bound the entire CLI process group, including child engines."""
    command = [sys.executable, "-m", "predictor.run_cli", *map(str, arguments)]
    print("Running:", " ".join(command), flush=True)
    logs = out / "integration_logs"
    logs.mkdir(parents=True, exist_ok=True)
    logfile = logs / f"{len(list(logs.glob('*.log'))) + 1:02d}_{arguments[0]}.log"
    with logfile.open("w", encoding="utf-8") as stream:
        proc = subprocess.Popen(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                start_new_session=True)
        try:
            code = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            raise RuntimeError(f"CLI timed out after {timeout}s; inspect {logfile}") from None
    if code:
        print(logfile.read_text(encoding="utf-8", errors="replace")[-16000:], flush=True)
        raise RuntimeError(f"CLI exited {code}; inspect {logfile}")


def table(path: Path, delimiter=",") -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def require(condition, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sequence(path: Path) -> str:
    from Bio import SeqIO
    return str(SeqIO.read(path, "fasta").seq)


def select_fixture(run: Path) -> str:
    """Select the known positive by exact supplied protein, never by search order."""
    query = sequence(FIXTURES / "mini_tf.faa")
    matches = [row["run_name"] for row in table(run / "manifest.csv")
               if sequence(run / row["fasta"]) == query]
    require(len(matches) == 1, f"expected exactly one source zntR protein in {run}, got {matches}")
    return matches[0]


def validate_discovery(run: Path) -> dict:
    rows = table(run / "manifest.csv")
    require(rows and all(row["discovery_engine"] == "bitacora" for row in rows),
            f"real BITACORA discovery required in {run}")
    require(all("pfam_fallback" not in row["discovery_profile"] for row in rows),
            f"unexpected profile build fallback in {run}")
    raw = list((run / "discovery").rglob("*.tblout"))
    require(raw, f"raw HMMER evidence missing from {run}")
    for filename in ("DISCOVERY_REPORT.html", "candidates.tsv", "candidates.faa"):
        require((run / "discovery" / filename).is_file(), f"missing discovery artifact: {filename}")
    return {"candidates": len(rows), "raw_tblout_files": len(raw)}


def validate_run(run: Path, expected: str) -> dict:
    from openpyxl import load_workbook
    from tools.validate_bundle import validate_bundle
    from predictor.structure.af3_export import validate_jobs
    from predictor.structure.af3_jobset import validate_local

    bundles = sorted((run / "jobs").glob("*/bundle_manifest.json"))
    require([path.parent.name for path in bundles] == [expected],
            f"expected a complete selected bundle {expected} in {run}")
    summary = validate_bundle(bundles[0].parent)
    for filename in ("RUN_REPORT.html", "Supplementary_Table_S1.xlsx", "Supplementary_Table_S1.tsv",
                     "discovery/DISCOVERY_REPORT.html", "af3/af3_manifest.tsv"):
        require((run / filename).is_file(), f"missing {run / filename}")
    # Inspect links from the HTML and spreadsheet, not just the filenames.
    for html_path in (run / "RUN_REPORT.html", run / "discovery/DISCOVERY_REPORT.html"):
        for ref in re.findall(r'''(?:src|href)=["']([^"']+)["']''', html_path.read_text(encoding="utf-8")):
            parsed = urlsplit(ref)
            if not parsed.scheme and parsed.path:
                require((html_path.parent / unquote(parsed.path)).is_file(), f"broken report link: {ref}")
    workbook = load_workbook(run / "Supplementary_Table_S1.xlsx")
    links = [cell.hyperlink.target for sheet in workbook for row in sheet for cell in row if cell.hyperlink]
    require(links, "supplementary spreadsheet has no per-TF hyperlinks")
    require(all((run / unquote(urlsplit(link).path)).is_file() for link in links), "broken XLSX link")
    workbook.close()
    af3 = run / "af3"
    exported = table(af3 / "af3_manifest.tsv", delimiter="\t")
    for dialect in ("server", "local"):
        kinds = {row["kind"] for row in exported if row["format"] == dialect and row["status"] == "written"}
        require({"apo", "dna_r1"} <= kinds, f"missing apo/operator AF3 exports for {dialect}: {kinds}")
    for row in exported:
        if row["status"] == "written":
            require((af3 / row["file"]).is_file(), f"missing AF3 file: {row['file']}")
    server = list((af3 / "server").glob("*.json"))
    local = list((af3 / "local").glob("*.json"))
    require(server and local, "both AF3 formats must contain real jobs")
    for path in server:
        validate_jobs(json.loads(path.read_text(encoding="utf-8")))
    for path in local:
        validate_local(json.loads(path.read_text(encoding="utf-8")))
    return {"bundle": expected, "validation": summary, "server_batches": len(server), "local_jobs": len(local)}


def run(out: Path, timeout: int) -> dict:
    from predictor.discovery import bitacora
    from predictor.genome_run import ensure_genome

    require(sys.platform.startswith("linux"), "this integration certifies native Linux only")
    require(not out.exists() or not any(out.iterdir()), "use a fresh output directory; integration must not resume")
    out.mkdir(parents=True, exist_ok=True)
    for check in (bitacora.engine_available, bitacora.patch_status):
        ok, why = check()
        require(ok, why)
    require(bitacora._wsl_exe() is None, "native Linux must not use wsl.exe")

    mini = out / "mini_discovery"
    cli(out, timeout, "genome", "--fasta", FIXTURES / "mini_genome.fna", "--gff", FIXTURES / "mini_genome.gff",
        "--out", mini, "--discover-only", "--offline")
    evidence = {"mini_discovery": validate_discovery(mini)}
    select_fixture(mini)

    accession = out / "accession"
    cli(out, timeout, "genome", ACCESSION, "--out", accession, "--discover-only")
    evidence["accession_discovery"] = validate_discovery(accession)
    chosen = select_fixture(accession)
    # This call only locates the genome already downloaded by the accession CLI.
    fna, gff = ensure_genome(ACCESSION, allow_download=False)

    single = out / "single"
    cli(out, timeout, "predict", FIXTURES / "mini_tf.faa", "--organism", ACCESSION,
        "--name", "zntR", "--out", single, "--no-metalnet")
    cli(out, timeout, "report", single)
    cli(out, timeout, "af3", single, "--format", "both")
    evidence["single"] = validate_run(single, "zntR")

    batch = out / "batch"
    cli(out, timeout, "batch", FIXTURES / "mini_tf.faa", "--organism", ACCESSION,
        "--out", batch, "--no-metalnet", "--timeout", str(timeout))
    batch_name = table(batch / "manifest.csv")[0]["run_name"]
    evidence["batch"] = validate_run(batch, batch_name)

    local = out / "local_genome"
    local_args = ("genome", "--fasta", fna, "--gff", gff, "--organism", ACCESSION, "--out", local)
    cli(out, timeout, *local_args, "--discover-only")
    evidence["local_discovery"] = validate_discovery(local)
    local_name = select_fixture(local)
    cli(out, timeout + 60, *local_args, "--select", local_name, "--no-metalnet", "--timeout", str(timeout))
    evidence["local_genome"] = validate_run(local, local_name)

    cli(out, timeout + 60, "genome", ACCESSION, "--out", accession, "--select", chosen,
        "--no-metalnet", "--timeout", str(timeout))
    evidence["accession"] = validate_run(accession, chosen)
    return evidence


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=1200, help="per-command/prediction deadline, seconds")
    args = parser.parse_args()
    if args.timeout < 1:
        parser.error("--timeout must be positive")
    out = args.out.expanduser().resolve()
    try:
        evidence = run(out, args.timeout)
    except Exception as exc:
        if out.is_dir():
            (out / "integration_result.json").write_text(json.dumps({"status": "failed", "error": str(exc)}, indent=2))
        raise
    (out / "integration_result.json").write_text(json.dumps({"status": "passed", "evidence": evidence}, indent=2))
    print(json.dumps(evidence, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
