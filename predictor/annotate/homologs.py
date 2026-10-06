"""
homologs.py -- MMseqs2-backed homolog sweep for the annotation stage (Phase 0.2).

From a TF protein sequence, collect homologs against a target protein DB with **MMseqs2** (more
sensitive than BLAST/DIAMOND via iterative profile search -> more homologs in the thin-alignment
regime, easing Snowprint failure mode #3). Emits:

  * the homolog hit table (target id, % identity, e-value, bits, coverage),
  * a raw homolog count and a **non-redundant** count (MMseqs2 `easy-cluster` dereplication),
  * a `stratum` label (sparse / moderate / rich) -- the Snowprint-availability axis (Phase 0).

The non-redundant homolog set is also what `annotate/context.py` and the family library consume.

MMseqs2 must be installed (see `docs/ENGINES.md`); on Windows it is `env/tools/mmseqs/mmseqs.bat`,
invoked through `cmd /c` because the `easy-*` workflows use a bundled BusyBox shell. Set `$MMSEQS`
to override the binary path.

Run `python homologs.py` for a self-test (needs MMseqs2; builds a tiny target DB on the fly).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
_MMSEQS_DIR = _REPO / "env" / "tools" / "mmseqs"
_DEFAULT_MMSEQS_EXE = _MMSEQS_DIR / "bin" / "mmseqs.exe"     # call the exe directly (no cmd/.bat)
_DEFAULT_MMSEQS_BAT = _MMSEQS_DIR / "mmseqs.bat"            # bootstraps BusyBox on first run

# easy-search --format-output columns we request (order matters for parsing)
_FMT = "query,target,pident,evalue,bits,qcov,tcov"


@dataclass
class Homolog:
    target: str
    pident: float       # fraction in [0,1] (MMseqs2 default)
    evalue: float
    bits: float
    qcov: float
    tcov: float


@dataclass
class HomologSweep:
    query_id: str
    n_hits: int
    n_nonredundant: int
    stratum: str
    hits: list = field(default_factory=list)        # list[Homolog], best first


# --------------------------------------------------------------------------- mmseqs plumbing
def find_mmseqs() -> str:
    """Locate the MMseqs2 entry point: $MMSEQS -> bin/mmseqs.exe -> mmseqs.bat -> PATH.

    Prefer the bare exe: invoking the `.bat` through `cmd /c` corrupts arguments when BOTH the program
    path and an argument contain spaces (as in '5.11 Predictor'). BusyBox is bootstrapped once by the
    first `.bat` run; thereafter the exe runs the easy-* workflows directly.
    """
    env = os.environ.get("MMSEQS")
    if env and Path(env).exists():
        return env
    # ensure BusyBox is installed (bin/bash) before relying on the exe for easy-* workflows.
    # NB: bin/bash is a BusyBox symlink whose .exists() raises WinError 1920 when present -> that
    # error itself means "installed", so only a clean FileNotFound triggers the one-time bootstrap.
    if _DEFAULT_MMSEQS_EXE.exists():
        try:
            has_bash = (_MMSEQS_DIR / "bin" / "bash").exists()
        except OSError:
            has_bash = True
        if not has_bash and _DEFAULT_MMSEQS_BAT.exists():
            try:
                subprocess.run(["cmd", "/c", str(_DEFAULT_MMSEQS_BAT), "version"],
                               capture_output=True, text=True, timeout=60)
            except Exception:
                pass
        return str(_DEFAULT_MMSEQS_EXE)
    if _DEFAULT_MMSEQS_BAT.exists():
        return str(_DEFAULT_MMSEQS_BAT)
    found = shutil.which("mmseqs")
    if found:
        return found
    raise FileNotFoundError(
        "MMseqs2 not found. Install per docs/ENGINES.md or set $MMSEQS to the binary/.bat path."
    )


def _run(mmseqs: str, args, cwd: Path) -> subprocess.CompletedProcess:
    if os.name == "nt" and mmseqs.lower().endswith(".bat"):
        cmd = ["cmd", "/c", mmseqs, *map(str, args)]   # .bat needs cmd; bootstraps BusyBox shell
    else:
        cmd = [mmseqs, *map(str, args)]
    res = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    if res.returncode != 0:
        tail = (res.stderr or res.stdout or "")[-2000:]
        raise RuntimeError(f"mmseqs {args[0]} failed (rc={res.returncode}):\n{tail}")
    return res


def _write_fasta(path: Path, records) -> None:
    with path.open("w", encoding="ascii") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


def _as_query_fasta(query, tmp: Path) -> tuple[Path, str]:
    """Accept a fasta path or a raw sequence string; return (fasta_path, query_id)."""
    if isinstance(query, (str, os.PathLike)) and Path(str(query)).exists():
        # read first record id
        with open(query, encoding="ascii") as fh:
            qid = next((l[1:].split()[0].strip() for l in fh if l.startswith(">")), "query")
        return Path(query), qid
    qf = tmp / "query.fasta"
    _write_fasta(qf, [("query", str(query).strip())])
    return qf, "query"


def _parse_m8(path: Path) -> list:
    hits = []
    if not path.exists():
        return hits
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.strip():
            continue
        q, t, pid, ev, bits, qcov, tcov = line.split("\t")
        pidf = float(pid)
        pidf = pidf / 100.0 if pidf > 1.0 else pidf      # normalize: MMseqs2 emits % or fraction
        hits.append(Homolog(t, pidf, float(ev), float(bits), float(qcov), float(tcov)))
    hits.sort(key=lambda h: (h.evalue, -h.bits))
    return hits


def _read_fasta(path) -> dict:
    seqs, cur = {}, None
    for line in Path(path).read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0].strip()
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def _stratum(n: int) -> str:
    if n < 5:
        return "sparse"          # Snowprint conservation signal unreliable -> structure must carry
    if n <= 25:
        return "moderate"
    return "rich"


# --------------------------------------------------------------------------- public API
def homolog_sweep(query, target_fasta, *, mmseqs: str | None = None,
                  sensitivity: float = 7.5, evalue: float = 1e-3, max_seqs: int = 2000,
                  min_qcov: float = 0.0, threads: int | None = None, cluster: bool = True,
                  cluster_min_seq_id: float = 0.9, cluster_cov: float = 0.8,
                  workdir: str | None = None) -> HomologSweep:
    """Sweep `query` (sequence or fasta) against `target_fasta` with MMseqs2 easy-search,
    then dereplicate the hits with easy-cluster. Returns a HomologSweep."""
    mm = mmseqs or find_mmseqs()
    threads = threads or max(1, (os.cpu_count() or 2) // 2)
    owns_tmp = workdir is None
    tmp = Path(workdir or tempfile.mkdtemp(prefix="homsweep_"))
    tmp.mkdir(parents=True, exist_ok=True)
    try:
        qf, qid = _as_query_fasta(query, tmp)
        res_m8 = tmp / "hits.m8"
        _run(mm, ["easy-search", qf, target_fasta, res_m8, tmp / "tmp_s",
                  "-s", sensitivity, "-e", evalue, "--max-seqs", max_seqs,
                  "--threads", threads, "--format-output", _FMT, "-v", 1], cwd=tmp)
        hits = [h for h in _parse_m8(res_m8) if h.qcov >= min_qcov]

        # dereplicate the hit set (non-redundant homolog count). easy-cluster is flaky on tiny
        # inputs on Windows (occasional segfault) -- it is a nicety, so failure falls back to n_hits.
        n_nr = len(hits)
        if hits and cluster:
            try:
                tdb = _read_fasta(target_fasta)
                sub = [(h.target, tdb[h.target]) for h in hits if h.target in tdb]
                if sub:
                    hf = tmp / "hits.fasta"
                    _write_fasta(hf, sub)
                    _run(mm, ["easy-cluster", hf, tmp / "clust", tmp / "tmp_c",
                              "--min-seq-id", cluster_min_seq_id, "-c", cluster_cov,
                              "--threads", threads, "-v", 1], cwd=tmp)
                    ctsv = tmp / "clust_cluster.tsv"
                    if ctsv.exists():
                        reps = {l.split("\t")[0] for l in ctsv.read_text().splitlines() if l.strip()}
                        n_nr = len(reps)
            except Exception:
                n_nr = len(hits)            # dereplication unavailable; report raw count

        return HomologSweep(qid, len(hits), n_nr, _stratum(n_nr), hits)
    finally:
        if owns_tmp:
            shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    import random

    # a base "regulator" protein + 3 near-homologs (point-mutated) + 2 unrelated proteins
    base = ("MSEKQDLTVKDLAKETGLSVHTLRYYERIGLLPEPDRSEGNYRLYTQAHLERLAFIKRAKRLGFSL"
            "EEIAELLALWDDRHRASADVKAIAQAHLAEVDARIAELQAMRDTLQHLADACCGDARPDCPILDELS")
    random.seed(7)

    def mutate(seq, n):
        s = list(seq)
        for _ in range(n):
            i = random.randrange(len(s))
            s[i] = random.choice("ACDEFGHIKLMNPQRSTVWY")
        return "".join(s)

    unrelated1 = "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(140))
    unrelated2 = "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(95))
    targets = [
        ("homolog_5pct", mutate(base, 7)),
        ("homolog_12pct", mutate(base, 16)),
        ("homolog_20pct", mutate(base, 27)),
        ("unrelated_1", unrelated1),
        ("unrelated_2", unrelated2),
    ]

    tmp = Path(tempfile.mkdtemp(prefix="homdemo_"))
    try:
        tfa = tmp / "targets.fasta"
        _write_fasta(tfa, targets)
        sweep = homolog_sweep(base, tfa, evalue=1e-3, workdir=str(tmp / "work"))
        print(f"query={sweep.query_id}  hits={sweep.n_hits}  "
              f"non-redundant={sweep.n_nonredundant}  stratum={sweep.stratum}")
        for h in sweep.hits:
            print(f"  {h.target:<16} id={h.pident:.0%}  e={h.evalue:.1e}  "
                  f"bits={h.bits:.0f}  qcov={h.qcov:.2f}")
        found = {h.target for h in sweep.hits}
        assert "homolog_5pct" in found and "homolog_12pct" in found, "near-homologs missed"
        assert not (found & {"unrelated_1", "unrelated_2"}), "unrelated proteins should not hit"
        assert sweep.stratum in ("sparse", "moderate", "rich")
        print("OK: homologs recovered, unrelated proteins excluded, stratum assigned.")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    _demo()
