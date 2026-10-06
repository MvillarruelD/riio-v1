"""
ligand.py -- enzyme protein -> its reaction metabolites (the candidate effector small molecules).

Reimplements Ligify's enzyme<->reaction<->chemical step. Ligify maps chemicals to enzymes via Rhea +
UniProt; for TF->ligand we run it FORWARD: a neighbouring enzyme's UniProtKB entry already lists its
catalysed reactions under `CATALYTIC ACTIVITY` comments, each carrying the reaction equation and its
**Rhea + ChEBI** cross-references. So one cached UniProt call per enzyme yields the participant
metabolites directly -- no extra Rhea/ChEBI round-trips (keeps API calls minimal, per project policy).

Currency metabolites (water, ATP/ADP, NAD(P), CoA, Pi/PPi, O2, CO2, H+, ...) are filtered: a TF almost
never senses these. What remains are the substrate/product candidates handed to `rank.py`.

Disk cache: `$PREDICTOR_EFFECTOR_CACHE` (default `~/.predictor/effector_cache`), one JSON per accession.

Run `python ligand.py` for a self-test (synthetic UniProt entry parsed offline; live fetch SKIPs).
"""
from __future__ import annotations

import json
import os
import re
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from predictor import resources

_CACHE_DIR = Path(os.environ.get("PREDICTOR_EFFECTOR_CACHE",
                                 resources.cache_path("effector")))
_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
# UniProt accession pattern (so we can tell a bare UniProt id from a RefSeq/locus-tag cross-reference)
_UNIPROT_ACC_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")

# common "currency"/cofactor metabolites a TF is essentially never a sensor for (lower-cased names)
CURRENCY = {
    "h2o", "water", "h(+)", "h+", "proton", "hydron", "atp", "adp", "amp", "gtp", "gdp", "pi",
    "phosphate", "diphosphate", "ppi", "pyrophosphate", "nad(+)", "nadh", "nadp(+)", "nadph",
    "nad+", "nadp+", "fad", "fadh2", "fmn", "coa", "coenzyme a", "co2", "co3(2-)", "carbon dioxide",
    "o2", "dioxygen", "oxygen", "nh4(+)", "nh3", "ammonia", "ammonium", "h2o2", "hydrogen peroxide",
    "so4(2-)", "cl(-)", "na(+)", "k(+)", "mg(2+)", "electron", "acetyl-coa", "udp", "ump", "cmp",
    "s-adenosyl-l-methionine", "s-adenosyl-l-homocysteine", "2-oxoglutarate", "l-glutamate",
}

# stoichiometric / article prefixes to strip from an equation participant ("2 H2O", "a fatty acid")
_PREFIX_RE = re.compile(r"^\s*(\d+|a|an|\(n\)|n)\s+", re.I)


@dataclass
class Reaction:
    equation: str
    rhea: str = ""
    ec: str = ""
    chebi: list = field(default_factory=list)         # ChEBI ids from UniProt reactionCrossReferences
    participants: list = field(default_factory=list)  # human-readable names (currency NOT yet removed)


# --------------------------------------------------------------------------- parsing
def _strip_participant(p: str) -> str:
    return _PREFIX_RE.sub("", p).strip()


_EQ_SEP_RE = re.compile(r"\s=\s|\s<=>\s|\s->\s")


def parse_equation(eq: str) -> list:
    """Split a Rhea/UniProt equation 'A + B = C + D' into participant names (stoichiometry stripped).

    Returns [] for a prose reaction description (no '='/'<=>'/'->' separator) -- some UniProt entries
    carry free-text activity ("Hydrolysis of terminal ... residues") rather than a balanced equation,
    and those must not be mistaken for a metabolite. Precision matters more than recall here."""
    if not eq or not _EQ_SEP_RE.search(eq):
        return []
    parts = []
    for side in _EQ_SEP_RE.split(eq):
        for token in side.split(" + "):
            name = _strip_participant(token)
            if name:
                parts.append(name)
    return parts


def reactions_from_entry(entry: dict) -> list:
    """Extract catalysed reactions from a UniProtKB entry JSON (`comments` of type CATALYTIC ACTIVITY)."""
    out = []
    for c in entry.get("comments", []):
        if (c.get("commentType") or "").upper() != "CATALYTIC ACTIVITY":
            continue
        rx = c.get("reaction") or {}
        eq = rx.get("name", "")
        rhea, chebi, ec = "", [], rx.get("ecNumber", "")
        for xr in rx.get("reactionCrossReferences", []):
            db, rid = (xr.get("database") or ""), (xr.get("id") or "")
            if db.lower() == "rhea" and not rhea and rid.upper().startswith("RHEA"):
                rhea = rid
            elif db.lower() == "chebi":
                chebi.append(rid)
        out.append(Reaction(eq, rhea, ec, chebi, parse_equation(eq)))
    return out


def _ec_numbers(entry: dict) -> list:
    pd = entry.get("proteinDescription") or {}
    names = [pd.get("recommendedName") or {}] + (pd.get("alternativeNames") or [])
    return [e.get("value") for n in names for e in (n.get("ecNumbers") or []) if e.get("value")]


def is_enzyme(entry: dict) -> bool:
    """An entry is an enzyme if it has a catalytic-activity reaction OR an EC number (even partial).
    Only the former yields reaction metabolites; the latter still marks the gene as enzymatic."""
    if any((c.get("commentType") or "").upper() == "CATALYTIC ACTIVITY"
           for c in entry.get("comments", [])):
        return True
    return bool(_ec_numbers(entry))


def candidate_ligands(entry: dict) -> list:
    """Non-currency reaction participants of an enzyme entry, de-duplicated, with provenance.
    Returns list[dict(name, chebi, equation, rhea)]."""
    seen, out = set(), []
    for rx in reactions_from_entry(entry):
        for name in rx.participants:
            key = name.lower()
            if key in CURRENCY or key in seen or len(name) <= 1:
                continue
            seen.add(key)
            out.append({"name": name, "chebi": rx.chebi, "equation": rx.equation, "rhea": rx.rhea})
    return out


# --------------------------------------------------------------------------- cached UniProt fetch
def _cache_path(accession: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", accession)
    return _CACHE_DIR / f"{safe}.json"


def uniprot_entry(accession: str, *, use_cache: bool = True, timeout: int = 30) -> dict | None:
    """Fetch the best UniProtKB entry for a protein accession (RefSeq WP_/NP_, EMBL, or UniProt id).
    Cached to disk to limit API calls. Returns the entry JSON or None."""
    cp = _cache_path(accession)
    if use_cache and cp.exists():
        try:
            return json.loads(cp.read_text(encoding="utf-8")) or None
        except Exception:
            pass
    acc = accession.split(".")[0]
    # RefSeq/EMBL accessions (NP_/WP_/YP_/XP_/AP_..., or a locus tag) are cross-references; a bare
    # UniProt accession is queried directly. A compound `... OR accession:NP_...` is a 400 (NP_ is not a
    # valid UniProt accession), which silently broke every RefSeq lookup -- so branch on the id shape.
    field = "xref" if ("_" in acc or not _UNIPROT_ACC_RE.match(acc)) else "accession"
    q = f"{field}:{acc}"
    url = (f"{_UNIPROT_SEARCH}?query={urllib.parse.quote(q)}&format=json&size=1"
           "&fields=accession,protein_name,cc_catalytic_activity,ec,xref_refseq")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as h:
            data = json.loads(h.read().decode("utf-8"))
        entry = (data.get("results") or [None])[0]
    except Exception:
        entry = None
    if use_cache:
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cp.write_text(json.dumps(entry or {}), encoding="utf-8")
    return entry


# --------------------------------------------------------------------------- self-test
def _synthetic_entry() -> dict:
    """A UniProtKB-shaped entry for an arabinose-metabolising enzyme (the AraC paradigm)."""
    return {
        "primaryAccession": "P00000",
        "proteinDescription": {"recommendedName": {"fullName": {"value": "L-arabinose isomerase"}}},
        "comments": [{
            "commentType": "CATALYTIC ACTIVITY",
            "reaction": {
                "name": "beta-L-arabinopyranose = L-ribulose",
                "reactionCrossReferences": [
                    {"database": "Rhea", "id": "RHEA:14821"},
                    {"database": "ChEBI", "id": "CHEBI:16401"},
                    {"database": "ChEBI", "id": "CHEBI:16880"},
                ],
                "ecNumber": "5.3.1.4",
            }}, {
            "commentType": "CATALYTIC ACTIVITY",
            "reaction": {"name": "ATP + D-ribulose = ADP + D-ribulose 5-phosphate",
                         "reactionCrossReferences": [{"database": "Rhea", "id": "RHEA:17601"}]},
        }],
    }


def _demo() -> None:
    entry = _synthetic_entry()
    assert is_enzyme(entry), "should detect catalytic activity"
    rxns = reactions_from_entry(entry)
    print(f"parsed {len(rxns)} reaction(s); first: {rxns[0].equation} (Rhea {rxns[0].rhea})")
    assert rxns[0].rhea == "RHEA:14821" and rxns[0].chebi == ["CHEBI:16401", "CHEBI:16880"]

    ligs = candidate_ligands(entry)
    names = [l["name"] for l in ligs]
    print(f"candidate ligands (currency filtered): {names}")
    assert "L-ribulose" in names and "beta-L-arabinopyranose" in names, "sugar participants missing"
    assert "ATP" not in names and "ADP" not in names, "currency metabolites not filtered"
    assert "D-ribulose" in names, "second-reaction substrate missing"
    print("OK: UniProt entry -> reactions -> non-currency candidate ligands (offline).")

    # optional live fetch (best-effort; SKIP on no network)
    try:
        e = uniprot_entry("P0A9E0", use_cache=False, timeout=15)   # E. coli AraC (a regulator, no rxn)
        if e is not None:
            print(f"OK (network): UniProt returned {e.get('primaryAccession')} "
                  f"(enzyme={is_enzyme(e)}).")
        else:
            print("SKIP (network): UniProt returned no entry.")
    except Exception as ex:
        print(f"SKIP (network): UniProt fetch not exercised ({type(ex).__name__}).")


if __name__ == "__main__":
    _demo()
