"""
structure_backend.py -- the REAL biohub ESMFold2 backend for structure_scan (the live-compute path).

Bridges `signals/structure_scan.StructureVerifier` to the mature scanner library in
`5.8 Promoters/arsr_merr_families/scripts/_lib` (`esmfold2_inputs` -> `esmfold2_backend` ->
`esmfold2_adapter` -> `operator_score`). Two entry points:

  score_existing_fold(folder, tf_id)         -- score one of the 239 folds on disk (ZERO new folds).
  fold_and_score(tf_id, prot, dna, ...)      -- LIVE: build TF-DNA input -> biohub fold -> score.
  make_backend(tf_id, family=...)            -- a `fold_backend(tf_seq, dna, copies, seed, plddt)`
                                                closure for StructureVerifier.verify (debits the cap).

All `_lib`/`esm` imports are lazy so this module imports offline; only the live calls need the
`BIOHUB_TOKEN` (in `arsr_merr_families/.env`) + network. Folds land in a space-free work dir
(`~/.predictor/structure_work`). Budget: biohub allows ~100 folds/day -- always dry-run/count first.

Run `python structure_backend.py` for a self-test (offline: scores an existing fold if present; the
live fold path is exercised by `benchmark/live_structural.py`, not here).
"""
from __future__ import annotations

import sys
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
WORK = resources.cache_path("structure_work")
FAST_MODEL = "esmfold2-fast-2026-05"
FULL_MODEL = "esmfold2-2026-05"

# Contact pLDDT_min is FOLDER-SPECIFIC and does NOT transfer (2026-06-15 AF3 experiment): ESMFold2 folds
# the DBD at pLDDT 30-60 (use 30/50 by family), but AF3 folds it at ~90, so the 30/50 filter passes every
# contact incl. spurious ones (it inverted the CueR control until raised). Score AF3 folds at ~80-85.
AF3_PLDDT_MIN = 80.0

#: The ESMFold2 server's cap on ONE COMPLEX, in residues. Observed, not documented -- the server
#: answers an over-cap request with:
#:     422: Input sequence length (868) exceeds maximum allowed sequence length (768)
#: It applies to the whole assembly, so the monomer budget depends on oligomeric state: dimer
#: <= 384 aa, tetramer (CsoR/RcnR, LysR, LacI/GalR) <= 192 aa.
#:
#: Distinct from the LOCAL per-monomer guard in env/tools/esmfold_lib (500 aa), which is about what
#: this project is willing to send. This is what the backend will accept. Raising the local guard
#: from 150 to 500 took foldable coverage 66 -> 147 of 154 and was necessary but not sufficient:
#: in run A, 7 candidates still failed here (five MocR-type GntR at 442-490 aa, one ArsR at 434,
#: one MerR at 390 -- all fine as monomers, all over cap as dimers), each after ~30 min and a quota
#: unit spent on a request that could not succeed.
MAX_COMPLEX_RESIDUES = 768


def _esmfold_scripts_dir():
    """Dir containing `_lib/` (the ESMFold2 backend), via the central resource registry so a shipped
    package can point at an explicitly installed copy (``PREDICTOR_ESMFOLD_LIB``)."""
    return resources.esmfold_lib_dir()


def _lib():
    """Lazy import of the scanner library (adds the scripts dir to sys.path once). Also loads the API
    tokens from the resolved .env so the biohub backend's token loader picks them up wherever it lives."""
    try:
        from predictor import resources
        resources.load_tokens()
    except Exception:
        pass
    scripts = _esmfold_scripts_dir()
    if scripts is None:
        raise RuntimeError("ESMFold2 adapter is not configured (set PREDICTOR_ESMFOLD_LIB)")
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    from _lib import esmfold2_inputs, esmfold2_backend, esmfold2_adapter, operator_score
    return esmfold2_inputs, esmfold2_backend, esmfold2_adapter, operator_score


def backend_available() -> bool:
    """True if the biohub ESMFold2 backend is usable (token present + SDK importable). Lets the structure
    step decide whether to fold or degrade to an AFDB-only / request-descriptor path WITHOUT folding."""
    try:
        _, backend, _, _ = _lib()
        backend._load_token()                 # raises if BIOHUB_TOKEN missing/placeholder
        return True
    except Exception:
        return False


def _features(fold_score) -> dict:
    """FoldScore -> the feature dict our composite consumes, plus the _lib composite for parity."""
    _, _, _, operator_score = _lib()
    fd = fold_score.feature_dict()
    fd["_lib_composite"] = operator_score.composite_score(fold_score)
    fd["_per_pos_base"] = list(getattr(fold_score, "per_pos_base", []) or [])
    fd["_contact_centroid"] = getattr(fold_score, "contact_centroid", None)
    fd["_mean_plddt_contacts"] = getattr(fold_score, "mean_plddt_contacts", None)
    return fd


def score_existing_fold(folder, tf_id: str, *, plddt_min: float = 50.0) -> dict | None:
    """Score an existing fold folder on disk (no network, no quota)."""
    _, _, _, operator_score = _lib()
    folder = Path(folder)
    if not folder.exists():
        return None
    fs = operator_score.score_fold(folder, tf_id, plddt_min=plddt_min)
    return _features(fs)


def fold_and_score(tf_id: str, prot_seq: str, dna_seq: str, *, plddt_min: float = 50.0,
                   model: str = FAST_MODEL, num_loops: int = 1, num_sampling: int = 16,
                   tag: str | None = None, cached_ok: bool = True) -> dict | None:
    """LIVE: build a TF-DNA input, fold it on biohub ESMFold2, score the result. Reuses a cached fold
    of the same tag when present. Returns the feature dict (with `_folded`/`_fold_seconds`) or None."""
    import re
    import time
    inputs, backend, adapter, operator_score = _lib()
    tag = tag or f"{tf_id.lower()}_live"
    if not re.fullmatch(r"[ACGT]+", dna_seq.upper()):
        return None
    folder = WORK / "af3_results" / f"{tag}_dna_v1"
    cif = folder / f"fold_{tag}_dna_v1_model_0.cif"
    folded, secs = False, None
    if not (cached_ok and cif.exists() and cif.stat().st_size > 0):
        built = inputs.build_tf_dna_input(tf_id, prot_seq, dna_seq.upper())
        t0 = time.time()
        result = backend.fold(built.structure_input, model=model, num_loops=num_loops,
                              num_sampling_steps=num_sampling)
        secs = round(time.time() - t0, 1)
        adapter.write_af3_layout(result, tf_id=tag, pred_type="dna", run_idx=1,
                                 chain_ids=["A", "B", "C", "D"], base_dir=WORK,
                                 chain_lengths={"A": len(prot_seq), "B": len(prot_seq),
                                                "C": len(dna_seq), "D": len(dna_seq)})
        folded = True
    fs = operator_score.score_fold(folder, tf_id, plddt_min=plddt_min)
    feat = _features(fs)
    feat["_folded"] = folded
    feat["_fold_seconds"] = secs
    return feat


def fold_apo_homomer(tf_id: str, prot_seq: str, *, copies: int = 2, model: str = FAST_MODEL,
                     num_loops: int = 1, num_sampling: int = 16, tag: str | None = None,
                     cached_ok: bool = True) -> dict | None:
    """LIVE: fold the selected apo homo-oligomer and write the AF3-style layout.

    Returns {cif, folder, plddt_mean, ptm, iptm, n_tokens, folded, fold_seconds} or None on failure /
    no backend. Quota: one biohub fold/day budget unit when not cached -- always probe one TF first."""
    import time
    inputs, backend, adapter, _ = _lib()
    if not isinstance(copies, int) or not 1 <= copies <= 26:
        raise ValueError(f"copies must be an integer in [1, 26], got {copies!r}")
    # Keep the existing dimer cache path; use a distinct key for every higher-order state so a
    # previously cached dimer can never be presented as the requested tetramer.
    tag = tag or (f"{tf_id.lower()}_apo" if copies == 2 else f"{tf_id.lower()}_apo{copies}")
    folder = WORK / "af3_results" / f"{tag}_apo_v1"
    cif = folder / f"fold_{tag}_apo_v1_model_0.cif"
    if cached_ok and cif.exists() and cif.stat().st_size > 0:
        # re-read the cached confidence: ptm/iptm from the summary, mean pLDDT from full_data atom_plddts
        import glob
        import json as _json
        plddt_mean = ptm = iptm = None
        summ = next(iter(glob.glob(str(folder / "*summary_confidences_0.json"))), None)
        if summ:
            try:
                d = _json.loads(Path(summ).read_text(encoding="utf-8"))
                ptm, iptm = d.get("ptm"), d.get("iptm")
            except Exception:
                pass
        full = next(iter(glob.glob(str(folder / "*full_data_0.json"))), None)
        if full:
            try:
                ap = _json.loads(Path(full).read_text(encoding="utf-8")).get("atom_plddts") or []
                if ap:
                    plddt_mean = float(sum(ap) / len(ap))
            except Exception:
                pass
        return {"cif": str(cif), "folder": str(folder), "plddt_mean": plddt_mean, "ptm": ptm,
                "iptm": iptm, "n_tokens": None, "folded": False, "fold_seconds": None,
                "copies": copies}
    built = inputs.build_apo_homomer_input(tf_id, prot_seq, copies=copies)
    t0 = time.time()
    result = backend.fold(built.structure_input, model=model, num_loops=num_loops,
                          num_sampling_steps=num_sampling)
    secs = round(time.time() - t0, 1)
    chain_ids = built.chain_ids
    rep = adapter.write_af3_layout(
        result, tf_id=tag, pred_type="apo_dimer", run_idx=1,
        chain_ids=chain_ids, base_dir=WORK,
        chain_lengths={chain_id: len(prot_seq) for chain_id in chain_ids},
    )
    return {"cif": str(rep.cif_path), "folder": str(rep.folder), "plddt_mean": rep.plddt_mean,
            "ptm": result.ptm, "iptm": result.iptm, "n_tokens": rep.n_tokens,
            "folded": True, "fold_seconds": secs, "copies": copies}


def fold_apo_dimer(tf_id: str, prot_seq: str, **kwargs) -> dict | None:
    """Backward-compatible two-chain wrapper."""
    return fold_apo_homomer(tf_id, prot_seq, copies=2, **kwargs)


def make_backend(tf_id: str, *, plddt_min: float = 50.0, model: str = FAST_MODEL,
                 num_loops: int = 1, num_sampling: int = 16):
    """A `fold_backend(tf_seq, dna_seq, copies, seed, plddt)` closure for StructureVerifier.verify.
    `copies` is informational (the input builder always uses the dimer); `seed` tags the cache."""
    def backend(tf_seq, dna_seq, copies, seed, plddt):
        return fold_and_score(tf_id, tf_seq, dna_seq, plddt_min=plddt or plddt_min, model=model,
                              num_loops=num_loops, num_sampling=num_sampling,
                              tag=f"{tf_id.lower()}_s{seed}_{abs(hash(dna_seq)) % 10**6}")
    return backend


_BASE3 = {"DA": "A", "DT": "T", "DG": "G", "DC": "C", "A": "A", "T": "T", "G": "G", "C": "C"}
_COMP1 = {"A": "T", "T": "A", "G": "C", "C": "G"}
_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def structure_readout_pwm(folder, tf_id: str, *, plddt_min: float = 50.0, family: str = "auto",
                          pseudocount: float = 0.3, max_dist: float = 3.5) -> dict | None:
    """A structure-based base-specificity readout (our DeepPBS stand-in, since torch_geometric/DeepPBS
    is not installed). Reads the recognition-helix→DNA direct-base contacts in a folded TF–DNA complex
    and builds a TOP-STRAND PWM: each contacted position is weighted toward the base the helix reads
    (bottom-strand contacts complemented). Contacted (specificity-determining) positions get peaked
    columns / high information; uncontacted positions stay ~uniform. The dimer's two half-sites both
    contribute, so the dyad is reinforced.

    Returns dict(pwm 4xL, core_pwm over the contacted span, per_base_importance, contacted positions,
    n_base_contacts, dyad_center, L). The real DeepPBS network plugs in via deeppbs `model=` later.
    """
    import numpy as np
    inputs, backend, adapter, operator_score = _lib()
    from _lib import af3_parse
    from _lib.contacts import find_contacts
    from Bio.PDB import MMCIFParser

    folder = Path(folder)
    models = af3_parse.parse_result_folder(folder, tf_id, "tf_dna")
    if not models:
        return None
    model = models[0]
    struct = MMCIFParser(QUIET=True).get_structure(tf_id, str(model["cif_path"]))
    dna_res = operator_score._DNA_RESIDUES
    protein_chains, dna_chains = [], []
    for chain in struct.get_chains():
        is_dna = any(r.resname.strip() in dna_res for r in chain.get_residues())
        (dna_chains if is_dna else protein_chains).append(chain.id)
    top_chain = "C" if "C" in dna_chains else (dna_chains[0] if dna_chains else "C")
    plddt_arr = model.get("plddt_per_residue")
    plddt_chains = model.get("plddt_chain_ids") or []
    L = plddt_chains.count(top_chain) or 80
    helix_start = operator_score.detect_helix_start(struct, protein_chains, family=family)
    contacts, _ = find_contacts(struct, plddt_arr, max_dist=max_dist, plddt_min=plddt_min,
                                operator_center_pos=None, helix_start_resi=helix_start)

    counts = np.full((4, L), pseudocount)
    n_base = 0
    for c in contacts:
        if c.interaction_type != "direct_base":
            continue
        n = int(c.dna_position)
        b = _BASE3.get(c.dna_base.strip())
        if b is None:
            continue
        if c.dna_chain == top_chain:
            x = n - 1
        else:                                            # bottom-strand contact -> top base is complement
            x, b = L - n, _COMP1.get(b)
        if b is None or not (0 <= x < L):
            continue
        counts[_IDX[b], x] += 1.0
        n_base += 1
    pwm = counts / counts.sum(axis=0, keepdims=True)
    info = 2.0 + (pwm * np.log2(np.clip(pwm, 1e-9, 1.0))).sum(axis=0)
    contacted = [int(x) for x in np.where(info > 0.3)[0]]
    lo, hi = (min(contacted), max(contacted) + 1) if contacted else (0, L)
    return dict(pwm=pwm, core_pwm=pwm[:, lo:hi], per_base_importance=info, contacted=contacted,
                n_base_contacts=n_base, core_span=(lo, hi), dyad_center=(lo + hi) / 2.0, L=L)


def make_readout_model(tf_id: str, *, plddt_min: float = 50.0, family: str = "auto", core: bool = True):
    """A `model(fold_folder, seed) -> 4xL PWM` for `deeppbs.pwm_from_complex` that returns the structure
    readout (the contacted-core PWM by default). Seed-ensembling over several folds softens tolerant
    positions. Drop-in replaceable by the real DeepPBS network."""
    def model(fold_folder, seed):
        r = structure_readout_pwm(fold_folder, tf_id, plddt_min=plddt_min, family=family)
        if r is None:
            return None
        return r["core_pwm"] if core else r["pwm"]
    return model


# The DeepPBS bridge that used to live here (`_deeppbs_python`, `deeppbs_real_pwm`,
# `make_deeppbs_model`) was removed on 2026-09-02. It had exactly one consumer, `report.outputs
# .finalize_af3`, and that consumer imported `signals.deeppbs` -- a module already deleted from the
# package -- so the whole path raised ModuleNotFoundError before reaching any of this. The DeepPBS and
# FoldX ENGINES may still be installed under `env/tools/` on a given machine; what is gone is the
# package-side glue, and it is gone because nothing could reach it. See CHANGELOG 2026-09-02.


def dinuc_shuffled_decoys(seq: str, n: int, *, seed: int = 0) -> list:
    """n dinucleotide-preserving shuffles of `seq` (the structural null), via _lib.operator_score."""
    import numpy as np
    _, _, _, operator_score = _lib()
    rng = np.random.default_rng(seed)
    return [operator_score.dinuc_shuffle(seq.upper(), rng) for _ in range(n)]


def zscore(observed: float, null_values) -> float:
    _, _, _, operator_score = _lib()
    return operator_score.zscore_vs_null(observed, list(null_values))


# --------------------------------------------------------------------------- self-test (offline)
def _demo() -> None:
    # This adapter intentionally has no private-checkout fallback. Its self-test is a readiness probe;
    # score_existing_fold remains available to callers with an explicit fold folder.
    if backend_available():
        print("OK: ESMFold2 adapter and credentials resolve")
    else:
        print("SKIP: ESMFold2 adapter/credentials not configured")


if __name__ == "__main__":
    _demo()
