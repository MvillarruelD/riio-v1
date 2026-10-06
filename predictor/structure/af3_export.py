"""
af3_export.py -- emit AlphaFold-3 *server* input JSON for the top-ranked operators (the loop's hand-off).

The pipeline iterates cheaply on ESMFold2 + DeepPBS and spends the
premium AF3 quota on only the **top 5** operators. This module turns those picks into the
alphafoldserver.com job JSON the user uploads by hand. Per operator we emit one job:

  * the TF folded as a **homodimer** (`proteinChain` count=2 -- inter-subunit DNA/effector sites need it),
  * the operator as **double-stranded DNA** = two complementary `dnaSequence` chains (the server treats a
    duplex as two ssDNA chains), 60 bp centred on the dyad (NOT 30-40 -- short duplexes break dimer
    docking; ITERATION2),
  * an optional **ion** (the metal effector, mapped to an AF3-supported ion) and/or ligand.

`build_tf_dna_input` (the ESM pipeline) maps non-AF3 metals to the closest supported ion; we mirror that
mapping here (ZN for divalent/metalloids, CO for Ni2+, CU for Cu+/Ag+/Au+). The schema is validated
(required keys, [ACGT] strands, the two DNA chains are exact reverse-complements) before it is written, so
a malformed job never reaches the server.

Run `python af3_export.py` for a self-test (build + validate jobs offline; no network).
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

from predictor import resources

# alphafoldserver.com job names: letters/digits/underscore/hyphen only -- '.' (as in a RefSeq accession
# like "NC_000913.3"), ':', spaces etc. are REJECTED by the server. Every job name is sanitised through
# `sanitize_job_name` at construction (`af3_job`) so this can't silently reappear from a new name source.
_NAME_UNSAFE_RE = re.compile(r"[^A-Za-z0-9_-]+")
_NAME_SAFE_RE = re.compile(r"^[A-Za-z0-9_-]+$")


def sanitize_job_name(name: str, *, max_len: int = 80) -> str:
    """AF3-server-safe job name. Collapses any run of disallowed characters ('.', ':', space, ...) to a
    single underscore; truncates over-long names to `max_len`, keeping a short content hash suffix so two
    truncated names don't collide."""
    safe = _NAME_UNSAFE_RE.sub("_", name or "").strip("_") or "job"
    if len(safe) > max_len:
        h = hashlib.sha1(name.encode("utf-8")).hexdigest()[:6]
        safe = safe[:max_len - 7] + "_" + h
    return safe

# AlphaFold-3 server supported ions (alphafoldserver.com), `dialect: alphafoldserver`.
# This is EXACTLY the ten the server accepts: Ca2+, Co2+, Cu2+, Fe3+, K+, Mg2+, Mn2+, Na+, Zn2+, Cl-.
# Anything else is rejected at submission, so every other metal must be mapped to a surrogate below.
# (Verified empirically: a job carrying FE2 was refused while an otherwise identical FE job ran.)
# Do NOT re-add NI/RB/CS/SR/HG/CD/GD/BA/YB/FE2 -- they are valid CCD codes but not server-supported.
AF3_IONS = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}

# effector metal -> AF3-supported ion (mirrors _lib.esmfold2_inputs metal mapping in the ESM pipeline)
_EFFECTOR_ION = {
    "ZN": "ZN", "ZN2+": "ZN", "ZINC": "ZN", "CD": "ZN", "CD2+": "ZN", "PB": "ZN", "PB2+": "ZN",
    "CO": "CO", "CO2+": "CO", "NI": "CO", "NI2+": "CO",          # Ni2+ -> CO per the ESM mapping
    "CU": "CU", "CU+": "CU", "CU1+": "CU", "AG": "CU", "AG+": "CU", "AU": "CU", "AU+": "CU",
    # the server's only iron is FE (Fe3+); Fe2+ sensors get it as the nearest available surrogate
    "MN": "MN", "MN2+": "MN", "FE": "FE", "FE2+": "FE", "FE3+": "FE",
    "MG": "MG", "CA": "CA", "HG": "ZN", "HG2+": "ZN",           # Hg2+ -> ZN divalent surrogate
    "AS": "ZN", "AS3+": "ZN", "SB": "ZN", "BI": "ZN",            # metalloids -> ZN divalent surrogate
}

_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def revcomp(seq: str) -> str:
    return "".join(_COMP[b] for b in reversed(seq.upper()))


def map_effector_to_ion(effector) -> str | None:
    """Map a curated effector / metal label to an AF3-supported ion code, or None if not a mappable
    metal (organic ligands go through the `ligand` field instead).

    Accepts POLYSPECIFIC labels ("Cd2+/Pb2+"): SSN clusters are routinely lumped, and a slash label is
    how that is recorded. Each alternative is resolved and the first mappable ion wins -- without this a
    polyspecific metal cluster silently yields ion=None and the AF3 job goes out with no metal at all.
    (Non-metal slash labels like "persulfide/RSS" still resolve to None, which is correct -- they belong
    in the `ligand` field.)
    """
    if not effector:
        return None
    key = str(effector).strip().upper().replace(" ", "")
    if key in AF3_IONS:
        return key
    hit = _EFFECTOR_ION.get(key)
    if hit:
        return hit
    if "/" in key:
        for part in key.split("/"):
            ion = part if part in AF3_IONS else _EFFECTOR_ION.get(part)
            if ion:
                return ion
    return None


def center_construct(genome_seq: str, dyad_center: int, length: int = 60) -> str:
    """The `length`-bp top strand centred on the operator dyad, sliced from the genome (clamped to ends).
    Re-deriving flanks from the genome is what gives the dimer enough DNA to read both half-sites."""
    g = genome_seq.upper()
    half = length // 2
    lo = max(0, min(len(g) - length, dyad_center - half))
    return g[lo:lo + length]


def _clean_dna(seq: str, length: int) -> str:
    """Uppercase, keep only ACGT, and pad/centre-trim to `length` (N->A as a last resort for ragged ends)."""
    s = "".join(b if b in "ACGT" else "A" for b in seq.upper())
    if len(s) >= length:
        off = (len(s) - length) // 2
        return s[off:off + length]
    pad = length - len(s)
    return "A" * (pad // 2) + s + "A" * (pad - pad // 2)


@dataclass
class AF3Pick:
    """The minimum a pick needs to become an AF3 job. `dna_top` is the 60-bp construct top strand; if
    absent, `to_af3_jobs` derives it from a supplied genome sequence or the operator `seq`."""
    name: str
    rank_reason: str = ""
    accession: str = ""
    start: int = 0
    end: int = 0
    dyad_center: int = 0
    dna_top: str | None = None
    ion: str | None = None


def af3_job(name: str, tf_seq: str, dna_top: str, *, copies: int = 2, ion: str | None = None,
            ligand_ccd: str | None = None, seeds=(1,)) -> dict:
    """One alphafoldserver job dict: homodimer protein + dsDNA (two complementary chains) + optional
    ion/ligand."""
    dna_top = dna_top.upper()
    # every alphafoldserver entity REQUIRES `count` (the server rejects "failed to find count" otherwise),
    # not just the protein chain -- each dsDNA strand is its own chain with count 1.
    seqs = [{"proteinChain": {"sequence": tf_seq.upper(), "count": copies}},
            {"dnaSequence": {"sequence": dna_top, "count": 1}},
            {"dnaSequence": {"sequence": revcomp(dna_top), "count": 1}}]
    if ion:
        seqs.append({"ion": {"ion": ion, "count": 1}})
    if ligand_ccd:
        seqs.append({"ligand": {"ligand": ligand_ccd, "count": 1}})
    return {"name": sanitize_job_name(name), "modelSeeds": list(seeds), "sequences": seqs,
            "dialect": "alphafoldserver", "version": 1}


def to_af3_jobs(picks, tf_seq: str, *, genome_seq: str | None = None, length: int = 60,
                copies: int = 2, ion: str | None = None, ligand_ccd: str | None = None,
                seeds=(1,), name_prefix: str = "TF") -> list[dict]:
    """Turn ranked operator picks into AF3 jobs. Each pick may be an `AF3Pick` or any object exposing
    `genome_accession/start/end/dyad_center/seq` (e.g. `schema.CanonicalOperator`) plus an optional
    `rank_reason`/`generator`. DNA top strand resolution order: explicit `dna_top` -> centred slice of
    `genome_seq` -> the operator `seq` padded/trimmed. A per-pick `ion` overrides the global `ion`."""
    jobs = []
    for p in picks:
        acc = getattr(p, "accession", None) or getattr(p, "genome_accession", "") or ""
        start = int(getattr(p, "start", 0) or 0)
        end = int(getattr(p, "end", 0) or 0)
        dyad = int(getattr(p, "dyad_center", (start + end) // 2) or (start + end) // 2)
        reason = getattr(p, "rank_reason", "") or getattr(p, "generator", "") or "pick"
        dna = getattr(p, "dna_top", None)
        if not dna and genome_seq:
            dna = center_construct(genome_seq, dyad, length)
        if not dna:
            dna = getattr(p, "seq", "") or ""
        dna = _clean_dna(dna, length)
        name = getattr(p, "name", None) or f"{name_prefix}__{acc}_{start}-{end}_{reason}"
        pick_ion = getattr(p, "ion", None) or ion
        jobs.append(af3_job(name, tf_seq, dna, copies=copies, ion=pick_ion,
                            ligand_ccd=ligand_ccd, seeds=seeds))
    return jobs


def validate_jobs(jobs) -> None:
    """Raise ValueError on any schema problem the server would reject. Checks: dialect/version,
    exactly one dimeric proteinChain, exactly two dnaSequence chains that are reverse-complements,
    [ACGT]-only DNA, and a supported ion code."""
    problems = []
    for j in jobs:
        nm = j.get("name", "<unnamed>")
        if not nm or not _NAME_SAFE_RE.match(nm):
            problems.append(f"{nm!r}: job name has characters the AF3 server rejects "
                            f"(only letters/digits/underscore/hyphen allowed)")
        if j.get("dialect") != "alphafoldserver" or j.get("version") != 1:
            problems.append(f"{nm}: bad dialect/version")
        seqs = j.get("sequences", [])
        prots = [s["proteinChain"] for s in seqs if "proteinChain" in s]
        dnas = [s["dnaSequence"] for s in seqs if "dnaSequence" in s]
        ions = [s["ion"] for s in seqs if "ion" in s]
        # the server requires `count` on EVERY entity, not just the protein chain
        for s in seqs:
            for key, ent in s.items():
                if "count" not in ent:
                    problems.append(f"{nm}: {key} entity is missing required 'count'")
        if len(prots) != 1 or prots[0].get("count", 1) < 2:
            problems.append(f"{nm}: expected one homodimeric proteinChain (count>=2)")
        if len(dnas) != 2:
            problems.append(f"{nm}: expected two dnaSequence chains (dsDNA), got {len(dnas)}")
        else:
            top, bot = dnas[0]["sequence"].upper(), dnas[1]["sequence"].upper()
            if any(b not in "ACGT" for b in top + bot):
                problems.append(f"{nm}: DNA has non-ACGT characters")
            elif revcomp(top) != bot:
                problems.append(f"{nm}: the two DNA chains are not reverse-complements")
        for io in ions:
            if io.get("ion") not in AF3_IONS:
                problems.append(f"{nm}: unsupported ion {io.get('ion')!r}")
    if problems:
        raise ValueError("invalid AF3 jobs:\n  " + "\n  ".join(problems))


def manifest(jobs) -> list[dict]:
    """A flat, printable summary (name -> DNA length, ion, chains) for the hand-off."""
    rows = []
    for j in jobs:
        seqs = j["sequences"]
        dna = next((s["dnaSequence"]["sequence"] for s in seqs if "dnaSequence" in s), "")
        ion = next((s["ion"]["ion"] for s in seqs if "ion" in s), None)
        rows.append({"name": j["name"], "dna_len": len(dna), "ion": ion, "n_chains": len(seqs)})
    return rows


def write_af3_jobs(jobs, path) -> Path:
    """Validate, then write the jobs as a JSON array ready to upload to the AF3 server. Returns the path."""
    validate_jobs(jobs)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jobs, indent=2), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    tf = "MARSKDELTAEQVFKALSDPNRLRILSLLAKGELCVCDL217CNLALSQSTVSHHLKLLREAGLVTGERRGREVYYRLA"

    # effector -> ion mapping (Ni2+ -> CO; Cu+ -> CU; metalloid As -> ZN)
    assert map_effector_to_ion("Ni2+") == "CO" and map_effector_to_ion("Cu+") == "CU"
    assert map_effector_to_ion("Zn") == "ZN" and map_effector_to_ion("As3+") == "ZN"
    assert map_effector_to_ion("glucose") is None
    # polyspecific SSN labels must still yield an ion -- without this the AF3 job goes out metal-free.
    # Cd2+ and Pb2+ both map to the ZN surrogate, so the lumped MerR_2 (CadR/PbrR) label is unambiguous.
    assert map_effector_to_ion("Cd2+/Pb2+") == "ZN" == map_effector_to_ion("Cd2+"), "polyspecific -> ion"
    assert map_effector_to_ion("persulfide/RSS") is None, "non-metal slash label is not an ion"
    print(f"ion map: Ni2+->{map_effector_to_ion('Ni2+')} Cu+->{map_effector_to_ion('Cu+')} "
          f"As3+->{map_effector_to_ion('As3+')}")

    # centre a 60-bp construct on a dyad inside a synthetic genome
    genome = ("ACGT" * 100)
    con = center_construct(genome, dyad_center=200, length=60)
    assert len(con) == 60 and con == genome[170:230]
    print(f"centred construct: len={len(con)} ...{con[:12]}")

    # job names are sanitised for the AF3 server (no '.', ':', spaces, ...): a RefSeq-accession-style name
    # loses its dot, and the validator now flags anything that still slips through.
    assert sanitize_job_name("CueR_Ecoli__NC_000913.3_513913-513949_rescan_tier1") == \
        "CueR_Ecoli__NC_000913_3_513913-513949_rescan_tier1"
    assert _NAME_SAFE_RE.match(sanitize_job_name("weird name: with spaces & dots.here"))
    long_name = "x" * 120
    assert len(sanitize_job_name(long_name)) <= 80
    j_dotted = af3_job("CzrA__NC_007795.1_1000-1030_logo", tf, "ACGT" * 15, copies=2, ion="ZN")
    assert "." not in j_dotted["name"] and _NAME_SAFE_RE.match(j_dotted["name"])
    try:
        validate_jobs([{**j_dotted, "name": "bad.name"}])
    except ValueError as e:
        print(f"name validator rejects a dotted name as expected ({str(e).splitlines()[-1].strip()})")
    else:
        raise AssertionError("validator missed a job name with a disallowed character")
    print(f"sanitize: 'CzrA__NC_007795.1_1000-1030_logo' -> {j_dotted['name']!r}")

    # build jobs from three picks with different rank reasons + a per-pick ion
    picks = [
        AF3Pick("CzrA__NC_007795.1_1000-1030_logo", "logo", "NC_007795.1", 1000, 1030, 1015, ion="ZN"),
        AF3Pick("CzrA__NC_007795.1_4000-4030_quality", "complex_quality", "NC_007795.1",
                4000, 4030, 4015, ion="ZN"),
        AF3Pick("CzrA__NC_007795.1_7000-7030_conservation", "conservation", "NC_007795.1",
                7000, 7030, 7015, ion="ZN"),
    ]
    jobs = to_af3_jobs(picks, tf, genome_seq=genome * 40, length=60, ion="ZN", seeds=(1,))
    validate_jobs(jobs)                                   # raises if malformed
    assert len(jobs) == 3
    j0 = jobs[0]
    seqs = j0["sequences"]
    prot = next(s["proteinChain"] for s in seqs if "proteinChain" in s)
    dnas = [s["dnaSequence"]["sequence"] for s in seqs if "dnaSequence" in s]
    assert prot["count"] == 2, "TF must be a homodimer"
    assert len(dnas) == 2 and revcomp(dnas[0]) == dnas[1], "DNA chains must be reverse-complements"
    assert len(dnas[0]) == 60, "construct must be 60 bp"
    assert j0["dialect"] == "alphafoldserver" and j0["version"] == 1
    print(f"job[0]: {j0['name']}  chains={len(seqs)} dna={len(dnas[0])}bp dimer={prot['count']}")

    # a broken job (single-stranded DNA) must be rejected
    bad = af3_job("bad", tf, "ACGTACGTACGT", copies=2, ion="ZN")
    bad["sequences"] = [s for s in bad["sequences"] if "dnaSequence" not in s][:1] + \
                       [{"dnaSequence": {"sequence": "ACGTACGTACGT"}}]
    try:
        validate_jobs([bad])
    except ValueError as e:
        print(f"rejected malformed job as expected ({str(e).splitlines()[-1].strip()})")
    else:
        raise AssertionError("validator missed a single-stranded DNA job")

    # write + manifest round-trip
    out = resources.cache_path("af3_demo", "jobs.json")
    write_af3_jobs(jobs, out)
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert len(loaded) == 3 and loaded[0]["name"] == jobs[0]["name"]
    for row in manifest(jobs):
        print(f"  manifest: {row['name']:<48} {row['dna_len']}bp ion={row['ion']} chains={row['n_chains']}")
    print(f"OK: AF3 jobs built (dimer + dsDNA + ion), validated (revcomp/[ACGT]/ion), written -> {out}")


if __name__ == "__main__":
    _demo()
