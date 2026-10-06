#!/usr/bin/env python3
"""bitacora_runner.py -- drive BITACORA and normalise its output into one joinable table.

Ours, not upstream's. It exists because BITACORA's own outputs cannot be joined to anything:

* its headline file is `*_proteins_trimmed.fasta`, cut down to the matched domain. That is why the
  survey's manifest lists *M. tuberculosis* IdeR (Rv2711) at 154 aa against its real 230, and why its
  `start`/`end` describe the trimmed span rather than the CDS. Matching those by sequence identity
  measures BITACORA's trimmer, not a census -- an early comparison that did so scored 62/150 and was
  wrong, where matching on the accession scores 149/150;
* it writes one directory tree per genome x seed profile, so "which profiles hit this protein" has to
  be reassembled from ~15 subdirectories.

So the runner emits **full-length records keyed by accession**, with the trim carried as an offset
annotation, and one row per protein listing every profile that hit it.

Invoked by `predictor/discovery/bitacora.py`, which is what the pipeline calls. On Windows the engine
runs through WSL (HMMER has no Windows build); this script may run on either side and converts paths
accordingly.

    python env/tools/bitacora_runner.py --genome g.fna --gff g.gff --proteins g.faa \\
        --profiles results/bitacora_profiles --engine env/tools/bitacora --out results/bitacora/g
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

_DISTRO = os.environ.get("BITACORA_WSL_DISTRO", "Ubuntu")
_WSL_BIN = os.environ.get("BITACORA_BIN", "$HOME/miniforge3/envs/bitacora/bin")


def win2wsl(p) -> str:
    p = str(Path(p).resolve())
    if len(p) > 1 and p[1] == ":":
        return "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")
    return p.replace("\\", "/")


def read_fasta(path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def run_bitacora(*, engine: Path, profiles: Path, proteins: Path, genome: Path | None,
                 gff: Path | None, out: Path, name: str, threads: int, evalue: str,
                 timeout: int) -> tuple[int, str]:
    """Run BITACORA in protein mode. Returns (returncode, combined output).

    Protein mode is deliberate: we already hold an ANNOTATED proteome (the accession-preserving route
    measured at 149/150), so there is nothing for BITACORA's gene prediction to add and running `full`
    mode would drag in GeMoMa and re-call genes we already have -- reintroducing the very ID-invention
    problem the annotated route exists to avoid.
    """
    on_windows = os.name == "nt"
    script = engine / "runBITACORA_command_line.sh"
    args = ["-m", "protein",
            "-q", win2wsl(profiles) if on_windows else str(profiles),
            "-p", win2wsl(proteins) if on_windows else str(proteins),
            "-n", name,
            "-sp", (win2wsl(engine / "Scripts") if on_windows else str(engine / "Scripts")),
            "-t", str(threads), "-e", evalue, "-c", "F"]
    # The coverage fractions are read by BITACORA's own Perl (see env/tools/bitacora_patches/).
    # WSL does not inherit the Windows environment, so they have to be exported into the inner shell
    # explicitly -- a variable set on this side and simply assumed to arrive would silently leave
    # upstream's defaults in force, and the run would look configured while behaving as stock.
    cov = "".join(f"export {k}={os.environ[k]}; " for k in ("BITACORA_QCOV", "BITACORA_SCOV")
                  if os.environ.get(k))
    if on_windows:
        inner = (f'export PATH={_WSL_BIN}:$PATH; ' + cov +
                 f'mkdir -p "{win2wsl(out)}" && cd "{win2wsl(out)}" && '
                 f'bash "{win2wsl(script)}" ' + " ".join(f'"{a}"' for a in args))
        cmd = [os.environ.get("SystemRoot", r"C:\Windows") + r"\System32\wsl.exe",
               "-d", _DISTRO, "-e", "bash", "-lc", inner]
    else:
        out.mkdir(parents=True, exist_ok=True)
        cmd = ["bash", str(script), *args]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                          cwd=(None if on_windows else str(out)))
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def collect(out: Path, proteins: dict[str, str]) -> list[dict]:
    """Reassemble BITACORA's per-profile trees into one row per PROTEIN.

    The trimmed FASTA gives which region matched; the full-length sequence comes from the proteome we
    submitted, keyed by accession. Every profile that hit a protein is recorded, because a protein
    hitting two seed profiles is a real observation about it and collapsing that to one loses it.
    """
    by_acc: dict[str, dict] = {}
    for trimmed in sorted(Path(out).rglob("*_proteins_trimmed.fasta")):
        profile = trimmed.name.replace("_proteins_trimmed.fasta", "")
        for acc, tseq in read_fasta(trimmed).items():
            base = acc.split()[0]
            full = proteins.get(base)
            if not full:
                # the trimmed id sometimes carries a suffix; fall back to a prefix match
                cand = next((k for k in proteins if base.startswith(k) or k.startswith(base)), None)
                full = proteins.get(cand) if cand else None
                if full:
                    base = cand
            if not full:
                continue
            row = by_acc.setdefault(base, {
                "protein_id": base, "sequence": full, "profiles": [],
                "trim_start": None, "trim_end": None, "search": "bitacora_protein_mode",
            })
            if profile not in row["profiles"]:
                row["profiles"].append(profile)
            i = full.find(tseq.strip().upper())
            if i >= 0 and row["trim_start"] is None:
                row["trim_start"], row["trim_end"] = i, i + len(tseq.strip())
    return list(by_acc.values())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteins", required=True, help="annotated proteome FASTA (accession-keyed)")
    ap.add_argument("--profiles", required=True, help="dir of <name>_db.fasta (+ optional _db.hmm)")
    ap.add_argument("--engine", required=True, help="the BITACORA checkout")
    ap.add_argument("--out", required=True)
    ap.add_argument("--genome", default=None)
    ap.add_argument("--gff", default=None)
    ap.add_argument("--name", default="tfop")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--evalue", default="1e-3")
    ap.add_argument("--timeout", type=int, default=7200)
    a = ap.parse_args(argv)

    proteins_path, out = Path(a.proteins), Path(a.out)
    proteome = read_fasta(proteins_path)
    if not proteome:
        print(f"ERROR: no sequences in {proteins_path}", file=sys.stderr)
        return 2
    out.mkdir(parents=True, exist_ok=True)

    rc, log = run_bitacora(engine=Path(a.engine), profiles=Path(a.profiles), proteins=proteins_path,
                           genome=Path(a.genome) if a.genome else None,
                           gff=Path(a.gff) if a.gff else None, out=out, name=a.name,
                           threads=a.threads, evalue=a.evalue, timeout=a.timeout)
    (out / "bitacora.log").write_text(log, encoding="utf-8", errors="replace")
    if rc != 0:
        print(f"ERROR: BITACORA exited {rc}; see {out / 'bitacora.log'}", file=sys.stderr)
        print(log[-1500:], file=sys.stderr)
        return rc

    rows = collect(out, proteome)
    table = out / "candidates.tsv"
    cols = ["protein_id", "profile", "sequence", "trim_start", "trim_end", "search",
            "contig", "start", "end", "strand", "evalue", "bits"]
    with table.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            r = dict(r)
            r["profile"] = "|".join(r.pop("profiles", []) or [])
            w.writerow(r)
    print(f"{len(rows)} candidate protein(s) -> {table}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
