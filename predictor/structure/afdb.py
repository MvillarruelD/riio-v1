"""
afdb.py -- AlphaFold DB client (the precomputed-structure API, cached).

AFDB is a "big" resource, so we use its API rather than re-folding -- a free, instant, high-quality
**monomer** model keyed on UniProt accession. Because it is monomer-only it serves the per-chain needs
(pocket detection, AlphaFill ligand/metal transplant); the biologically relevant multimer (dimer /
CsoR tetramer) still has to come from the folding backend (see `fold.py`). Models + metadata are cached
on disk so a TF is fetched at most once.

  api:    GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}  -> [ {cifUrl, pdbUrl, ...} ]
  cache:  $PREDICTOR_STRUCT_CACHE/afdb/  (default ~/.predictor/struct_cache/afdb)

Run `python afdb.py` for a self-test (offline metadata parse; live lookup SKIPs if no network).
"""
from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path

from predictor import resources

_CACHE = Path(os.environ.get("PREDICTOR_STRUCT_CACHE",
                             resources.cache_path("structure"))) / "afdb"
_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


def _meta_from_record(rec: dict) -> dict:
    return {
        "uniprot": rec.get("uniprotAccession") or rec.get("uniprotId"),
        "entry_id": rec.get("entryId"),
        "organism": rec.get("organismScientificName"),
        "cifUrl": rec.get("cifUrl"),
        "pdbUrl": rec.get("pdbUrl"),
        "paeUrl": rec.get("paeDocUrl") or rec.get("paeImageUrl"),
        "version": rec.get("latestVersion") or rec.get("modelCreatedDate"),
        "is_multimer": False,                 # AFDB is monomer-only
    }


def lookup(uniprot_acc: str, *, use_cache: bool = True, timeout: int = 30) -> dict | None:
    """Return AFDB model metadata for a UniProt accession (or None if not in AFDB / unreachable)."""
    acc = uniprot_acc.split(".")[0]
    cp = _CACHE / f"{acc}.meta.json"
    if use_cache and cp.exists():
        try:
            d = json.loads(cp.read_text(encoding="utf-8"))
            return d or None
        except Exception:
            pass
    try:
        with urllib.request.urlopen(_API.format(acc=acc), timeout=timeout) as h:
            data = json.loads(h.read().decode("utf-8"))
        meta = _meta_from_record(data[0]) if data else None
    except Exception:
        meta = None
    if use_cache:
        _CACHE.mkdir(parents=True, exist_ok=True)
        cp.write_text(json.dumps(meta or {}), encoding="utf-8")
    return meta


def fetch_model(uniprot_acc: str, *, fmt: str = "cif", use_cache: bool = True,
                timeout: int = 60) -> str | None:
    """Download the AFDB monomer model (cif|pdb) for a UniProt accession; return the cached path."""
    acc = uniprot_acc.split(".")[0]
    out = _CACHE / f"{acc}.{fmt}"
    if use_cache and out.exists() and out.stat().st_size > 0:
        return str(out)
    meta = lookup(uniprot_acc, use_cache=use_cache, timeout=timeout)
    if not meta:
        return None
    url = meta.get("cifUrl") if fmt == "cif" else meta.get("pdbUrl")
    if not url:
        return None
    try:
        _CACHE.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=timeout) as h:
            out.write_bytes(h.read())
    except Exception:
        return None
    return str(out)


def _demo() -> None:
    # offline: metadata parse from a synthetic AFDB API record
    rec = {"entryId": "AF-P0AE22-F1", "uniprotAccession": "P0AE22",
           "organismScientificName": "Escherichia coli",
           "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P0AE22-F1-model_v4.cif",
           "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P0AE22-F1-model_v4.pdb",
           "paeDocUrl": "https://x/pae.json", "latestVersion": 4}
    m = _meta_from_record(rec)
    print(f"parsed AFDB meta: {m['uniprot']} {m['entry_id']} multimer={m['is_multimer']}")
    assert m["uniprot"] == "P0AE22" and m["cifUrl"].endswith(".cif") and m["is_multimer"] is False
    print("OK: AFDB metadata parsed (monomer-only flagged).")

    # optional live lookup (ArsR E. coli, P0ACS5 is a known AFDB entry); SKIP on no network
    try:
        live = lookup("P0ACS5", use_cache=False, timeout=15)
        if live:
            print(f"OK (network): AFDB has P0ACS5 -> {live.get('entry_id')} ({live.get('organism')}).")
        else:
            print("SKIP (network): AFDB returned nothing for P0ACS5.")
    except Exception as ex:
        print(f"SKIP (network): AFDB not exercised ({type(ex).__name__}).")


if __name__ == "__main__":
    _demo()
