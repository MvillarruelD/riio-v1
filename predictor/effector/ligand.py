"""
ligand.py -- enzyme protein -> its reaction metabolites (the candidate effector small molecules).

Reimplements Ligify's enzyme<->reaction<->chemical step. Ligify maps chemicals to enzymes via Rhea +
UniProt; for TF->ligand we run it FORWARD: a neighbouring enzyme's UniProtKB entry already lists its
catalysed reactions under `CATALYTIC ACTIVITY` comments, each carrying the reaction equation and its
**Rhea + ChEBI** cross-references. So one cached UniProt call per enzyme yields the participant
metabolites directly -- no extra Rhea/ChEBI round-trips (keeps API calls minimal, per project policy).

Currency metabolites (water, ATP/ADP, NAD(P), CoA, Pi/PPi, O2, CO2, H+, ...) are filtered: a TF almost
never senses these. What remains are the substrate/product candidates handed to `rank.py`.

Disk cache: `$PREDICTOR_EFFECTOR_CACHE` (default platform resource cache), with new versioned
annotation records and read-only historical accession JSON. See effector.annotation for the NEW
provenance/cache-only/refresh/pinned-record API.

Run `python ligand.py` for a self-test (synthetic UniProt entry parsed offline; live fetch SKIPs).
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from predictor import resources

_CACHE_DIR = Path(os.environ.get("PREDICTOR_EFFECTOR_CACHE",
                                 resources.cache_path("effector")))

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
# EXP067: compartments do not change currency identity. Use only for the membership test, so real
# metabolites such as molybdate(out) retain their original stored/displayed spelling.
_COMPARTMENT_SUFFIX_RE = re.compile(r"\s*\(\s*(?:in|out)\s*\)\s*$", re.I)

# Macromolecular reaction participants a TF cannot sense as a small-molecule effector: nucleic acids,
# whole proteins and protein-bound residues. Currency (CURRENCY) removes cofactors; this removes the
# polymers that leaked through as "ligands" -- e.g. tRNA(Gly), RNA(n), protein L-histidine,
# N(6)-acetyl-L-lysyl-[protein] (EXP066: 8 of 26 prototype operon-chemistry calls). Deliberately narrow
# so it cannot catch a genuine metabolite: anchored on the "[protein]" tag, nucleic-acid words with
# boundaries, and a leading bare "protein".
_MACROMOLECULE_RE = re.compile(
    r"\[protein\]"                                  # ...-[protein] residue tags
    r"|\b[a-z]{0,4}rna\b|\brna\("                   # RNA, tRNA, rRNA, mRNA, ncRNA, tmRNA, RNA(n)
    r"|\bdna\b"
    r"|\bpoly(?:nucleotide|peptide|ribonucleotide|deoxyribonucleotide)\b"
    r"|\boligo(?:nucleotide|peptide)\b"
    r"|\bribosom(?:e|al)\b"
    r"|^\s*(?:an?\s+)?protein\b",                   # "protein L-histidine", "a protein"
    re.I)


def _is_macromolecule(name: str) -> bool:
    return bool(name and _MACROMOLECULE_RE.search(name))


@dataclass
class Reaction:
    equation: str
    rhea: str = ""
    ec: str = ""
    chebi: list = field(default_factory=list)         # legacy REACTION-wide IDs; never participant identity
    participants: list = field(default_factory=list)  # human-readable names (currency NOT yet removed)
    reaction_chebi: list = field(default_factory=list)
    occurrences: list = field(default_factory=list)
    source_record: dict = field(default_factory=dict)
    #: EXP072: UniProt's sourced direction for this reaction ("left-to-right" | "right-to-left"), or ""
    #: when the entry annotates none or both. The equation's written order is only a convention.
    physiological_direction: str = ""


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


def participant_occurrences(eq: str) -> list:
    """Lossless occurrence roles alongside the legacy stripped-name parser."""
    from predictor.effector.chemical_identity import compartment
    if not eq or not _EQ_SEP_RE.search(eq):
        return []
    out = []
    for side_index, side in enumerate(_EQ_SEP_RE.split(eq)):
        for token in side.split(" + "):
            name = _strip_participant(token)
            prefix = _PREFIX_RE.match(token)
            if name:
                out.append({"name": name, "raw_token": token, "side": "left" if side_index == 0 else "right",
                            "coefficient": prefix.group(1) if prefix else "1",
                            "occurrence": len(out), "compartment": compartment(name)})
    return out


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
        directions = {(p.get("directionType") or "").lower() for p in c.get("physiologicalReactions") or []}
        directions &= {"left-to-right", "right-to-left"}
        out.append(Reaction(eq, rhea, ec, chebi, parse_equation(eq), list(chebi),
                            participant_occurrences(eq), dict(rx),
                            next(iter(directions)) if len(directions) == 1 else ""))
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
    """Non-currency reaction participant occurrences of an enzyme entry, with provenance.
    `chebi` is now explicitly participant-scoped (empty without a source mapping).
    Reaction-wide IDs remain under `reaction_chebi`. Every occurrence is retained.
    Generic AH2 stays in raw evidence with headline_eligible=False."""
    from predictor.effector.chemical_identity import identity_key, transport_name
    out = []
    for rx in reactions_from_entry(entry):
        for occurrence in rx.occurrences:
            name = occurrence["name"]
            key = name.lower()
            currency_key = _COMPARTMENT_SUFFIX_RE.sub("", key).strip()
            if currency_key in CURRENCY or len(name) <= 1 or _is_macromolecule(name):
                continue
            # UniProt reactionCrossReferences do not map names to IDs. No positional inference.
            out.append({**occurrence, "chebi": [], "participant_identity_status": "unresolved",
                        "reaction_chebi": list(rx.reaction_chebi), "equation": rx.equation,
                        "rhea": rx.rhea, "source_record": rx.source_record,
                        "physiological_direction": rx.physiological_direction,
                        "identity_key": identity_key(name),
                        "headline_eligible": transport_name(name).casefold() != "ah2"})
    return out


# --------------------------------------------------------------------------- cached UniProt fetch
def _cache_path(accession: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", accession)
    return _CACHE_DIR / f"{safe}.json"


def uniprot_entry(accession: str, *, use_cache: bool = True, timeout: int = 30) -> dict | None:
    """Return raw UniProt JSON or None; exact mapping and provenance use annotation_record.

    Historical entries remain unverified and are never silently refetched. Ambiguous,
    failed or incompatible new acquisitions return None. ``use_cache=False`` bypasses
    both cache reads and writes. The provenance envelope never enters enzyme parsing.
    """
    from predictor.effector.annotation import annotation_record, raw_entry
    return raw_entry(annotation_record(accession, use_cache=use_cache, timeout=timeout, cache_dir=_CACHE_DIR))


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
