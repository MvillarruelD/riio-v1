"""
tf_record.py -- classify a TF from its sequence (Phase 0.2): family, nearest curated neighbour, and
(optionally) UniProt enrichment.

From a bare TF sequence (AA or DNA), assign the **protein family** by nearest-neighbour against the
curated 54-TF set (ArsR/SmtB + MerR) via MMseqs2 -- local, fast, and grounded in the project's own
labels rather than a heavy hmmscan/Pfam install.

There used to be a monomer length-sanity rule here (ArsR/SmtB, MerR and Fur monomers are typically
<=150 aa, so an oversized one suggested a wrong accession or a fusion). It was removed because
nothing read it: the flag it set was never consulted outside this module's own self-test, so it
blocked nothing and warned no one.

The curated reference (`predictor/data/refs/curated_tfs.faa` + `.json`) is built once from the
`arsr_merr_families/transcription_factors/*.yaml` records. UniProt enrichment is best-effort (network).

Run `python tf_record.py` for a self-test (builds the reference, leave-one-out classifies a real TF).
"""
from __future__ import annotations

import glob
import json
from dataclasses import dataclass, field
from pathlib import Path

from predictor import resources

from predictor.annotate.genome_resolver import is_dna
from predictor.annotate import family_db

_REPO = Path(__file__).resolve().parents[2]
_CURATED_TFS = _REPO.parent / "5.8 Promoters" / "arsr_merr_families" / "transcription_factors"
_REFDIR = resources.REFS_DIR
_FAA = _REFDIR / "curated_tfs.faa"
_JSON = _REFDIR / "curated_tfs.json"



@dataclass
class TFRecord:
    query_len: int
    is_dna: bool
    family: str | None
    family_group: str | None
    nearest_tf: str | None
    nearest_identity: float | None
    nearest_evalue: float | None
    uniprot: str | None
    nearest_uniprot: str | None
    pfam: list = field(default_factory=list)
    flags: list = field(default_factory=list)


# --------------------------------------------------------------------------- reference build
def build_reference(out_faa: Path = _FAA, out_json: Path = _JSON) -> int:
    """Build the curated TF reference (FASTA + metadata JSON) from the YAML records."""
    import yaml
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    meta, records = {}, []
    for f in sorted(glob.glob(str(_CURATED_TFS / "*" / "*.yaml"))):
        d = yaml.safe_load(open(f, encoding="utf-8"))
        seq = ((d.get("sequence") or {}).get("protein") or "").replace("\n", "").strip()
        if not seq:
            continue
        cls = d.get("classification") or {}
        ids = d.get("identifiers") or {}
        org = d.get("organism") or {}
        tid = d["id"]
        meta[tid] = {
            "family": cls.get("protein_family"),
            "family_group": cls.get("family_group"),
            "pfam": cls.get("pfam") or [],
            "uniprot": ids.get("uniprot"),
            "organism": org.get("species"),
            "taxid": org.get("taxid"),
            "length": len(seq),
        }
        records.append((tid, seq))

    # Never replace a populated vendored reference with an empty rebuild. `_CURATED_TFS` is the curated
    # YAML tree, which lives OUTSIDE this repo and so is absent on any plain install -- there the loop
    # above yields nothing, and writing that emptied the 53-sequence reference this package ships
    # (measured on a fresh `pip install`: 53 -> 0 sequences, from running the self-test alone). Same
    # failure mode as `known_operators.build_index`; the shipped artifact wins for an installed copy.
    if not records and out_faa.exists():
        have = sum(1 for ln in out_faa.read_text(encoding="utf-8", errors="replace").splitlines()
                   if ln.startswith(">"))
        if have:
            print(f"  tf_record: no curated YAMLs under {_CURATED_TFS}; "
                  f"keeping the vendored reference ({have} sequences)")
            return have

    # newline="\n": Python text mode translates to CRLF on Windows, but .gitattributes checks these
    # vendored files out as LF. Every rebuild (the self-test triggers one on a machine that has the
    # curated YAMLs) therefore left them permanently "modified" in git with ZERO content change --
    # noise that invites someone to clean it up and commit a mangled vendored reference by accident.
    with out_faa.open("w", encoding="ascii", newline="\n") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")
    out_json.write_text(json.dumps(meta, indent=0), encoding="utf-8", newline="\n")
    return len(records)


def default_reference(rebuild: bool = False):
    if rebuild or not _FAA.exists() or not _JSON.exists():
        build_reference()
    return str(_FAA), json.loads(_JSON.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- classify
def _to_protein(seq: str) -> str:
    if not is_dna(seq):
        return seq.replace("\n", "").strip()
    from Bio.Seq import Seq
    return str(Seq(seq.replace("\n", "").strip()).translate(to_stop=True))


def classify(seq: str, *, exclude_id: str | None = None, prefer: str = "hmm") -> TFRecord:
    """Holistic family classification (any family) via `family_db`: a confident Pfam HMM hit wins,
    else nearest-neighbour blastp against the 307-TF broad reference. `exclude_id` enables leave-one-out testing."""
    import json
    prot = _to_protein(seq)
    fc = family_db.classify_family(seq, exclude_id=exclude_id, prefer=prefer)
    fam = fc.family

    meta = {}
    if family_db._ALL_JSON.exists():
        meta = json.loads(family_db._ALL_JSON.read_text(encoding="utf-8"))
    m = meta.get(fc.nearest_tf or "", {})

    flags = list(fc.flags)
    pfam = [fc.pfam_id] if fc.pfam_id else m.get("pfam", [])
    # A bare sequence does not identify a UniProt record. The accession below belongs to the
    # nearest curated reference and must not enter query-identity paths such as AFDB QC or SSN lookup.
    nearest_uniprot = m.get("uniprot") if fc.nearest_tf else None
    return TFRecord(
        query_len=len(prot), is_dna=is_dna(seq), family=fam,
        family_group=m.get("family_group"), nearest_tf=fc.nearest_tf,
        nearest_identity=fc.nearest_identity, nearest_evalue=fc.nearest_evalue,
        uniprot=None, nearest_uniprot=nearest_uniprot, pfam=pfam,
        flags=flags + [f"family via {fc.method}: {fc.evidence}"],
    )


def uniprot_lookup(accession: str) -> dict:
    """Best-effort UniProt entry fetch by accession (network)."""
    import urllib.request
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    with urllib.request.urlopen(url, timeout=20) as r:
        d = json.loads(r.read().decode())
    return {
        "accession": d.get("primaryAccession"),
        "id": d.get("uniProtkbId"),
        "protein_name": (((d.get("proteinDescription") or {}).get("recommendedName") or {})
                         .get("fullName") or {}).get("value"),
        "organism": (d.get("organism") or {}).get("scientificName"),
    }


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    import yaml
    n = build_reference()
    print(f"built curated reference: {n} TFs -> {_FAA.name}, {_JSON.name}")

    # leave-one-out: classify a real MerR TF with itself excluded -> should still call MerR
    tf_yaml = _CURATED_TFS / "merr" / "AdhR_Bsubtilis.yaml"
    seq = yaml.safe_load(open(tf_yaml, encoding="utf-8"))["sequence"]["protein"]
    rec = classify(seq, exclude_id="AdhR_Bsubtilis")
    print(f"AdhR_Bsubtilis (LOO): family={rec.family} via {rec.nearest_tf} "
          f"id={rec.nearest_identity:.0%} len={rec.query_len}")
    assert rec.family == "MerR", f"expected MerR, got {rec.family}"
    assert rec.nearest_tf != "AdhR_Bsubtilis", "self should be excluded"
    assert rec.query_len <= 150

    # an ArsR/SmtB TF, leave-one-out -> ArsR/SmtB.
    # Uses Rv2358 (M. tuberculosis), a genuine ArsR/SmtB. This assertion previously used
    # "SmtB_Msmegmatis", which is NOT an ArsR/SmtB protein: it hits Pfam FUR at 122 bits vs 16 for the
    # ArsR-mapped HTH_20, carries the Fur DNA-binding motif, and is only 22% identical to Rv2358. That
    # record has been re-filed as Fur; the classifier was right and the expectation was wrong.
    af = _CURATED_TFS / "arsr_smtb" / "Rv2358_Mtuberculosis.yaml"
    s2 = yaml.safe_load(open(af, encoding="utf-8"))["sequence"]["protein"]
    r2 = classify(s2, exclude_id="Rv2358_Mtuberculosis")
    print(f"Rv2358_Mtuberculosis (LOO): family={r2.family} via {r2.nearest_tf} id={r2.nearest_identity:.0%}")
    assert r2.family == "ArsR/SmtB", f"expected ArsR/SmtB, got {r2.family}"

    print("OK: family classified by nearest curated neighbour.")


if __name__ == "__main__":
    _demo()
