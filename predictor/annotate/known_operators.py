"""
known_operators.py -- prior operator evidence: per-family consensus PWMs + known-operator lookup.

Cross-references the curated operator sets (arsr_merr_families + operators_db, ~5600 operators) to:
  * group operator sequences by TF family (using family_db's TF->family map), and
  * build a **family consensus PWM** (palindrome-aware, via motifs/em_finder) that seeds the motif
    finder and the genome rescan for a new family member -- the Phase 1.6B family-library prior, and
    what the first `predict` pass scans the genome with.
  * return any **known operator** already recorded for a specific TF id.

Built once into packaged `predictor/data/refs/operators_by_family.json`. Consensus PWMs are consumable by
`signals/motif_rescan.py` and `motifs/pwm_scan.py`.

Run `python known_operators.py` for a self-test (builds the index, makes a family consensus PWM).
"""
from __future__ import annotations

import glob
import hashlib
import json
import re
from pathlib import Path

from predictor import resources

from predictor.annotate import family_db
from predictor.motifs import em_finder

_REPO = Path(__file__).resolve().parents[2]
_OP_DIRS = [
    _REPO.parent / "5.8 Promoters" / "arsr_merr_families" / "operators",
    _REPO.parent / "5.8 Promoters" / "operators_db" / "operators",
]
_TF_DIRS = [
    _REPO.parent / "5.8 Promoters" / "arsr_merr_families" / "transcription_factors",
    _REPO.parent / "5.8 Promoters" / "operators_db" / "transcription_factors",
]
_INDEX = resources.REFS_DIR / "operators_by_family.json"

_DNA = set("ACGTNacgtn")


def _is_dna_seq(s):
    return s and sum(c in _DNA for c in s) / len(s) > 0.9


#: Spelling variants the source corpora use for one family. The index key MUST match
#: `annotate.ssn_clusters.FAMILIES`, because that is the string `motif_finder.family_seeded` and
#: `FAMILY_PRIORS` look up: a family filed as "ArsR" is simply invisible to a query for "ArsR/SmtB".
#: Canonical names map to themselves so a subgroup suffix can be stripped back onto them.
_FAMILY_ALIASES = {
    "arsr": "ArsR/SmtB", "smtb": "ArsR/SmtB", "arsr/smtb": "ArsR/SmtB",
    "merr": "MerR", "fur": "Fur", "gntr": "GntR", "copy": "CopY", "rrf2": "Rrf2",
    "csor": "CsoR/FrmR", "frmr": "CsoR/FrmR", "csor/frmr": "CsoR/FrmR",
    "nikr": "NikR", "copg/nikr": "NikR",
    "dtxr": "DtxR/MntR", "mntr": "DtxR/MntR", "ider": "DtxR/MntR", "dtxr/mntr": "DtxR/MntR",
    "marr": "MarR/SlyA", "slya": "MarR/SlyA", "marr/slya": "MarR/SlyA",
    "tetr": "TetR/AcrR", "acrr": "TetR/AcrR", "tetr/acrr": "TetR/AcrR",
    "lysr": "LysR-type (LTTR)", "lttr": "LysR-type (LTTR)",
}


def _canon_family(name):
    """Canonical family name for `name`, or None if it carries no family.

    Collapses the corpus's spelling variants ('ArsR', 'LysR', 'CsoR', 'CopG/NikR') and subgroup
    suffixes ('ArsR/SmtB (persulfide-sensing subgroup)') onto the pipeline's names. A family
    outside the twelve is cleaned but returned unchanged -- it is still indexed, just not renamed.
    """
    if not name:
        return None
    s = str(name).strip().strip('"')
    if not s or s.lower() in ("null", "none", "~"):
        return None
    if s.lower() in _FAMILY_ALIASES:
        return _FAMILY_ALIASES[s.lower()]
    base = re.sub(r"\s*\(.*?\)\s*$", "", s).strip()          # drop a trailing subgroup qualifier
    return _FAMILY_ALIASES.get(base.lower(), s)


def _family_of_record(d):
    """The family of one TF YAML, by PRECEDENCE -- and the precedence is the whole point.

    `classification.tf_family` is the specific curated call. `classification.protein_family` is a
    coarser roll-up: it buckets ten FNR/CRP regulators as "GntR" and every two-component response
    regulator as "LuxR-type / FixJ". Reading protein_family first (the behaviour until 2026-08-31)
    filed **740 FNR/CRP sites under GntR -- 90 % of that family's index set** -- and, because it
    never looked at tf_family at all, dropped every TF whose protein_family is null: CopY, NikR and
    Rrf2 therefore had NO curated operators whatsoever, and their geometry priors could not be
    checked against a single sequence.
    """
    cls = d.get("classification") or {}
    return _canon_family(cls.get("tf_family")) or _canon_family(cls.get("protein_family"))


# --------------------------------------------------------------------------- index build
def build_index(out=_INDEX) -> dict:
    """Group operator sequences by TF family and by TF id; cache to JSON."""
    import yaml
    out.parent.mkdir(parents=True, exist_ok=True)

    # TF id -> family. The TF YAMLs are the curated source of truth and OVERRIDE family_db's
    # all_tfs.json, which is kept only to cover TFs that have no YAML record.
    fam_of = {}
    if family_db._ALL_JSON.exists():
        fam_of = {k: _canon_family(v.get("family")) for k, v in
                  json.loads(family_db._ALL_JSON.read_text(encoding="utf-8")).items()}
    op_to_tf = {}
    for root in _TF_DIRS:
        for f in glob.glob(str(root / "**" / "*.yaml"), recursive=True):
            try:
                d = yaml.safe_load(open(f, encoding="utf-8"))
            except Exception:
                continue
            if not isinstance(d, dict):
                continue
            fam = _family_of_record(d)
            if fam:
                fam_of[d.get("id")] = fam
            for op in (d.get("operators") or []):
                op_to_tf[op] = d.get("id")

    by_family, by_tf = {}, {}
    for root in _OP_DIRS:
        for f in glob.glob(str(root / "**" / "*.yaml"), recursive=True):
            try:
                d = yaml.safe_load(open(f, encoding="utf-8"))
            except Exception:
                continue
            if not isinstance(d, dict):
                continue
            seq = (d.get("sequence") or "")
            if isinstance(seq, dict):
                seq = seq.get("dna") or seq.get("sequence") or ""
            seq = str(seq).replace("\n", "").strip().upper()
            if not _is_dna_seq(seq) or len(seq) < 8:
                continue
            tf = d.get("transcription_factor_id") or op_to_tf.get(d.get("id"))
            fam = fam_of.get(tf)
            by_tf.setdefault(tf, []).append(seq) if tf else None
            if fam:
                by_family.setdefault(fam, []).append(seq)

    index = {"by_family": by_family, "by_tf": by_tf}

    # NEVER overwrite a good shipped index with an empty rebuild. `_TF_DIRS` are the source YAMLs, which
    # live OUTSIDE this repo and so exist only on a machine that has the operators_db checkout next to it.
    # On a plain `pip install` they are absent, the rebuild yields {} -- and writing that destroyed the
    # 212 KB vendored index this package ships (measured: 212610 -> 30 bytes just by running the
    # self-test). The shipped artifact is the source of truth for an installed copy; a rebuild may only
    # replace it when it actually produced something.
    if not by_family and not by_tf:
        existing = json.loads(out.read_text(encoding="utf-8")) if out.exists() else index
        if existing.get("by_family") or existing.get("by_tf"):
            print(f"  known_operators: no source YAMLs found under {[str(r) for r in _TF_DIRS]}; "
                  f"keeping the vendored index ({len(existing.get('by_family', {}))} families)")
            return existing
        return index
    out.write_text(json.dumps(index), encoding="utf-8", newline="\n")
    return index


def load_index():
    """The operator index: the vendored artifact if present, else a rebuild from the source YAMLs."""
    if not _INDEX.exists():
        return build_index()
    return json.loads(_INDEX.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- consensus / lookup
def known_operator_for(tf_id: str):
    """Operator sequence(s) already recorded for a specific TF id, if any."""
    return load_index()["by_tf"].get(tf_id, [])


_BUNDLED_PWM_DIR = resources.REFS_DIR / "family_pwms"
_PWM_CACHE = resources.cache_path("family_pwms")


def family_operator_sequences(family: str, *, exclude_tf: str | None = None,
                              strict_exclusion: bool = False) -> list[str]:
    """Return curated family operators, optionally leaving one TF name/accession out.

    Exclusion uses the TF-level index plus TF-family metadata. ``CueR_bench`` therefore excludes
    all ``CueR_*`` records, a conservative leave-protein-out test of literature transfer.
    """
    from collections import Counter

    index = load_index()
    family_pool = list(index.get("by_family", {}).get(family, []))
    if not exclude_tf:
        return family_pool
    meta = {}
    if family_db._ALL_JSON.exists():
        meta = json.loads(family_db._ALL_JSON.read_text(encoding="utf-8"))

    query = re.sub(r"_bench(?:_v\d+)?$", "", exclude_tf or "", flags=re.IGNORECASE)
    query_norm = re.sub(r"[^a-z0-9]", "", query.lower())

    def excluded(tf_id: str) -> bool:
        if not query_norm:
            return False
        tf_name = re.split(r"[_|]", tf_id, maxsplit=1)[0]
        candidates = {
            re.sub(r"[^a-z0-9]", "", tf_id.lower()),
            re.sub(r"[^a-z0-9]", "", tf_name.lower()),
            re.sub(r"[^a-z0-9]", "", str((meta.get(tf_id) or {}).get("uniprot") or "").lower()),
        }
        return query_norm in candidates

    mapped = []
    excluded_ops = []
    had_mapped_family = False
    for tf_id, sequences in index.get("by_tf", {}).items():
        if (meta.get(tf_id) or {}).get("family") == family:
            had_mapped_family = True
            if excluded(tf_id):
                excluded_ops.extend(sequences)
            else:
                mapped.extend(sequences)
    if strict_exclusion:
        return mapped if had_mapped_family else []
    if not excluded_ops:
        return family_pool
    remove = Counter(excluded_ops)
    kept = []
    for sequence in family_pool:
        if remove[sequence]:
            remove[sequence] -= 1
        else:
            kept.append(sequence)
    return kept


def family_consensus(family: str, *, max_ops: int = 60, widths=None,
                     palindrome: bool = True, cache: bool = True,
                     exclude_tf: str | None = None, strict_exclusion: bool = False):
    """Build (and disk-cache) a palindrome-aware consensus PWM for a family's operators via em_finder.
    Returns (pwm 4xW, n_used) or (None, 0) if too few operators.

    Every family searches the same width range (`operator_logo.OPERATOR_WIDTHS`); two of them used to
    carry hand-written narrower ranges, so the EM could pick widths for one family that it was not allowed
    to consider for another."""
    import numpy as np
    if widths is None:
        from predictor.signals.operator_logo import OPERATOR_WIDTHS
        widths = OPERATOR_WIDTHS
    ops = family_operator_sequences(
        family, exclude_tf=exclude_tf, strict_exclusion=strict_exclusion,
    )
    ops = [o for o in ops if len(o) >= min(widths)][:max_ops]
    signature = hashlib.sha256(json.dumps({
        "version": 2, "family": family, "exclude_tf": exclude_tf,
        "strict_exclusion": strict_exclusion,
        "widths": list(widths), "palindrome": palindrome, "ops": ops,
    }, sort_keys=True).encode()).hexdigest()[:16]
    safe = "".join(c if c.isalnum() else "_" for c in family)
    cpath = _PWM_CACHE / f"{safe}_{signature}.npz"
    bundled = _BUNDLED_PWM_DIR / cpath.name
    existing = bundled if bundled.exists() else cpath
    if cache and existing.exists():
        cached = np.load(existing)
        return cached["pwm"], int(cached["n_used"])
    if len(ops) < 3:
        return None, len(ops)
    motif = em_finder.discover(ops, widths=widths, palindrome=palindrome, restarts=2)
    if cache:
        _PWM_CACHE.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cpath, pwm=motif.prob, n_used=len(ops))
    return motif.prob, len(ops)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    from predictor.motifs.pwm_scan import scan

    idx = load_index()                     # the vendored index; only rebuilds if it is genuinely absent
    fams = idx["by_family"]
    print(f"indexed operators across {len(fams)} families; "
          f"top: " + ", ".join(f"{k}={len(v)}" for k, v in
                                sorted(fams.items(), key=lambda kv: -len(kv[1]))[:6]))
    if not fams:
        print("  no operator index available (neither the vendored JSON nor the source YAMLs) "
              "-> nothing to self-test here")
        return

    # a family consensus PWM that actually scores its own operators (the best-populated family, whichever
    # that is -- no family is singled out here)
    fam = max(fams, key=lambda k: len(fams[k]))
    pwm, n = family_consensus(fam, widths=range(12, 22))
    assert pwm is not None and pwm.shape[0] == 4, "no consensus PWM"
    consensus = "".join("ACGT"[i] for i in pwm.argmax(axis=0))
    print(f"{fam}: consensus from {n} operators (W={pwm.shape[1]}): {consensus}")

    ops = [o for o in fams[fam] if len(o) >= pwm.shape[1]][:30]
    hits = scan({f"op{i}": o for i, o in enumerate(ops)}, prob=pwm, pvalue_thresh=1e-2)
    detected = len({h.seq_id for h in hits})
    print(f"family PWM re-detects the motif in {detected}/{len(ops)} of its own operators")
    assert detected >= max(3, len(ops) // 3), "family consensus too weak to find its own operators"

    print("OK: operator index built; family consensus PWM scores its own operators.")


if __name__ == "__main__":
    _demo()
