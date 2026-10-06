"""inducer_vocab.py -- one canonical spelling per inducer, and one place that knows what is not one.

The anchor KB stores whatever each source called the inducer, so the same species arrives under several
spellings: `Ni` / `Ni2+` / `Ni(II)`, `Cobalt` / `Co2+`, `Copper` / `Cu+`, `Zn(II)` / `Zn2+`. Two entries
are not inducers at all -- they are sentences meaning "none is known":

    'none identified (possibly NodD protein competition)'
    'none (glutamine synthetase feedback state)'

Both facts have consequences today, not merely for a future mixture design:

* `AnchorSet.conflict` counts DISTINCT inducer strings, and `effector.inducer` scales the
  `ssn_cluster` confidence down whenever it is true. `ArsR_c7` carried `persulfide`,
  `H2S / persulfide (RSS)` and an abstention sentence, was labelled `mixed-anchor (3 inducers)`, and
  had any assigned TF's confidence scaled by 1/n for what is one answer written three ways.
* Four sources emit inducer strings and `schema.InducerConsensus.from_calls` computes `agreement` by
  RAW equality across them, so two sources naming one species differently count as disagreeing. The
  SSN cluster table and `substrate_map` write arsenite `As3+`; the anchor KB writes `As(III)`.
  Full survey: `analysis/family_kb/vocab_audit.py`.

Scope is deliberately narrow. This normalises SPELLING; it does not decide chemistry. Deciding that
`Fe2+` and `Mn2+` are Irving-Williams neighbours, or that `Co2+` and formaldehyde are different
inducer classes, belongs to the mixture design (INDUCER_CALL_V2 Design 2) and not here. Two genuinely
different species stay two entries: `As(III)` and `MAs(III)` are arsenite and methylarsenite,
`antimonite` is not arsenite, and α-D-galacturonate is not α-D-glucuronate.

handoff §5.4 applies throughout: element symbols are matched as whole tokens, never as substrings.
`Cobalt` must not be found inside `steroid-CoA`, and `nitrogen` is not nickel.

    python -m predictor.effector.inducer_vocab --self-test
"""
from __future__ import annotations

import re

#: Exact spellings seen in the vendored anchor KB, mapped to one canonical token each. Built from the
#: real vocabulary (37 distinct strings over 2,046 anchors), not invented: extend it when the KB grows,
#: and keep the canonical form on the right-hand side stable, because it reaches reports.
_CANONICAL = {
    # --- metal ions -------------------------------------------------------------------------------
    "zn2+": "Zn2+", "zn(ii)": "Zn2+", "zinc": "Zn2+",
    "co2+": "Co2+", "co(ii)": "Co2+", "cobalt": "Co2+",
    "fe2+": "Fe2+", "fe": "Fe2+", "fe(ii)": "Fe2+", "iron": "Fe2+",
    "ni2+": "Ni2+", "ni(ii)": "Ni2+", "ni": "Ni2+", "nickel": "Ni2+",
    "hg2+": "Hg2+", "hg(ii)": "Hg2+", "mercury": "Hg2+",
    "cu+": "Cu+", "cu(i)": "Cu+", "copper": "Cu+",
    "mn2+": "Mn2+", "mn(ii)": "Mn2+", "manganese": "Mn2+",
    "cd2+": "Cd2+", "cd(ii)": "Cd2+", "cadmium": "Cd2+",
    "pb2+": "Pb2+", "pb(ii)": "Pb2+", "lead": "Pb2+",
    "au(i)": "Au+", "au+": "Au+", "gold": "Au+",
    "ag+": "Ag+", "ag(i)": "Ag+", "silver": "Ag+",
    "molybdate": "MoO42-", "moo42-": "MoO42-", "moo4(2-)": "MoO42-",
    "moo4": "MoO42-", "moo4-": "MoO42-",
    # Oxyanions of metals/metalloids are METAL-class inducers (what metalloregulators such as ModE,
    # the chromate/arsenate/vanadate sensors detect), even though the formula carries oxygen. They are
    # enumerated here as explicit tokens rather than pattern-matched on the element symbol.
    "tungstate": "WO42-", "wo42-": "WO42-", "wo4(2-)": "WO42-",
    "wo4": "WO42-", "wo4-": "WO42-",
    # `CrO4` is the exact spelling emitted by substrate_map; charge-omitted formulae are accepted
    # aliases, not separate species.
    "chromate": "CrO42-", "cro42-": "CrO42-", "cro4": "CrO42-", "cro4-": "CrO42-",
    "cr(vi)": "CrO42-",
    "arsenate": "AsO43-", "aso43-": "AsO43-", "aso4": "AsO43-", "aso4-": "AsO43-",
    "as(v)": "AsO43-", "asv": "AsO43-",
    "selenite": "SeO32-", "seo32-": "SeO32-", "seo3": "SeO32-", "seo3-": "SeO32-",
    "selenate": "SeO42-", "seo42-": "SeO42-", "seo4": "SeO42-", "seo4-": "SeO42-",
    "orthovanadate": "VO43-", "vanadate": "VO43-", "vo43-": "VO43-",
    "vo4": "VO43-", "vo4-": "VO43-",
    # Metavanadate and orthovanadate are distinct species; classify both as metal without merging.
    "metavanadate": "VO3-", "vo3-": "VO3-",
    "antimonate": "SbO43-", "sbo43-": "SbO43-", "sbo4": "SbO43-",
    "tellurite": "TeO32-", "teo32-": "TeO32-", "teo3": "TeO32-",
    "tellurate": "TeO42-", "teo42-": "TeO42-", "teo4": "TeO42-",
    # --- metalloid / oxyanion species. NOT collapsed into one another: different molecules. --------
    # `As3+` is how the SSN cluster table and `substrate_map` write arsenite while the anchor KB writes
    # `As(III)`. Found by `analysis/family_kb/vocab_audit.py`: it is the one species that three sources
    # spell two ways, so without this line the SSN source and the regulon source count as DISAGREEING
    # about arsenite and `InducerConsensus.agreement` drops for a unanimous call.
    "as(iii)": "As(III)", "arsenite": "As(III)", "as3+": "As(III)", "asiii": "As(III)",
    "aso2-": "As(III)", "aso33-": "As(III)", "aso3": "As(III)",
    "h2aso3-": "As(III)", "haso32-": "As(III)",
    "sb3+": "antimonite", "sb(iii)": "antimonite", "sbo33-": "antimonite", "sbo3": "antimonite",
    "mas(iii)": "MAs(III)", "mas(iii) (methylarsenite / monomethylarsenous acid)": "MAs(III)",
    "methylarsenite": "MAs(III)",
    "antimonite": "antimonite",
    # --- reactive species -------------------------------------------------------------------------
    "persulfide": "persulfide", "h2s / persulfide (rss)": "persulfide",
    "h2o2": "H2O2", "organic hydroperoxide": "organic hydroperoxide",
    # --- organic ligands --------------------------------------------------------------------------
    "formaldehyde": "formaldehyde",
    "salicylate": "salicylate", "cumate": "cumate", "citrate": "citrate", "co2": "CO2",
    "α-d-galacturonate": "α-D-galacturonate",
    "α-d-glucuronate": "α-D-glucuronate",
}

#: Strings that mean "no inducer is known". They must resolve to None, not to a ligand named "none".
_ABSTENTION_PREFIXES = ("none", "not identified", "unknown", "unclear", "n/a", "na", "-")


def is_abstention(raw: str | None) -> bool:
    """True when the stored value is a sentence meaning no inducer is known."""
    s = (raw or "").strip().lower()
    if not s:
        return True
    return any(s == p or s.startswith(p + " ") or s.startswith(p + "(") or s.startswith(p + ",")
               for p in _ABSTENTION_PREFIXES)


def normalise(raw: str | None) -> str | None:
    """Canonical inducer token, or None when the value is an abstention or empty.

    Unknown-but-real values are passed through with their case preserved, so a new ligand entering the
    KB is reported rather than silently dropped -- it simply does not get merged with anything.
    """
    if raw is None:
        return None
    s = raw.strip()
    if not s or is_abstention(s):
        return None
    low = s.lower()
    if low in _CANONICAL:
        return _CANONICAL[low]
    # strip a trailing parenthetical gloss: "MAs(III) (methylarsenite / ...)" -> "MAs(III)"
    stripped = re.sub(r"\s*\([^)]*\)\s*$", "", s).strip()
    if stripped and stripped.lower() in _CANONICAL:
        return _CANONICAL[stripped.lower()]
    # Last resort for an ion written in a notation the table has not seen: normalise the OXIDATION
    # STATE only (`Ni3+` / `Ni(III)` -> `ni`) and accept it if the bare element is known. This is
    # anchored on the whole token, so it can never fire inside a word -- `Cobalt` is not `Co`, and
    # `steroid-CoA` is not cobalt (handoff §5.4).
    # `Ni(2+)` as well as `Ni(III)` and `Ni3+`: ChEBI parenthesises the charge, and it reaches us that
    # way through Ligify. Tried on the gloss-stripped form too, so ChEBI's compartment suffix
    # (`Ni(2+)(in)` -> `Ni(2+)`) resolves instead of reading as a second, distinct species -- NikR
    # reported the candidate set ["Ni2+", "Ni(2+)(in)"], which is one ion displayed as an ambiguity.
    for cand in (s, stripped):
        m = re.fullmatch(r"([A-Za-z]{1,2})\s*(?:\(\s*[IVXivx]+\s*\)|\(\s*\d\s*[+-]\s*\)|\d\s*[+-])",
                         cand or "")
        if m and m.group(1).lower() in _CANONICAL:
            return _CANONICAL[m.group(1).lower()]
    return s


#: Separators that join two SPECIES in one stored label. A comma is deliberately absent: it appears
#: inside chemical names ("2,3-dihydroxybenzoate") far more often than it joins two of them.
_MIXTURE_SPLIT = re.compile(r"\s*/\s*")

#: Labels that read like a mixture but name ONE thing, mapped to the species they actually name.
#: Checked before splitting, because splitting them invents species that were never claimed. The value
#: is written out rather than derived, so the intent is auditable: each of these was read off the SSN
#: cluster table and decided individually.
_NOT_A_MIXTURE = {
    # one species written verbosely -- H2S and persulfide are the same RSS pool here
    "h2s / persulfide (rss)": "persulfide",
    "persulfide/rss": "persulfide",
    # a specific ligand plus the functional CLASS it belongs to ("salicylate, i.e. a multidrug
    # regulator"), not two co-inducers. The class is dropped; the species is the claim.
    "salicylate/multidrug": "salicylate",
    "steroid-coa/organic": "steroid-CoA",
    # cysteine and O-acetylserine are one biosynthetic signal in the CymR/CysB sense
    "cysteine/oas": "cysteine",
    # NO and RSS are both reactive species; the redox class carries this, not an ion name
    "no/rss": "NO",
}


def parse_mixture(raw: str | None) -> tuple[str, ...]:
    """Canonical COMPONENTS of a stored inducer label.

    The SSN cluster table stores co-sensing as a flattened string -- `Ni2+/Co2+`, `Cd2+/Pb2+`,
    `Cu+/persulfide`. A flattened string can never compare equal to a single-ion call, so
    `InducerConsensus.agreement` counts a cluster that co-senses Ni and Co as disagreeing with a source
    that says Ni. Returning the components lets a caller build a real mixture instead.

    A one-species label returns a 1-tuple, so callers need no special case. Three things it must not do,
    each of which would invent chemistry:

    * split a parenthetical gloss -- `MAs(III) (methylarsenite / monomethylarsenous acid)` is one species;
    * split a species-plus-class label -- `salicylate/multidrug` is salicylate, described;
    * split a verbose name for one pool -- `H2S / persulfide (RSS)`.
    """
    if raw is None:
        return ()
    s = raw.strip()
    if not s or is_abstention(s):
        return ()
    if s.lower() in _NOT_A_MIXTURE:
        return (_NOT_A_MIXTURE[s.lower()],)
    # a trailing parenthetical gloss may itself contain "/" -- remove it before considering a split
    body = re.sub(r"\s*\([^)]*\)\s*$", "", s).strip() or s
    if body.lower() in _NOT_A_MIXTURE:
        return (_NOT_A_MIXTURE[body.lower()],)
    if body.lower() in _CANONICAL:
        return (_CANONICAL[body.lower()],)
    parts = [p for p in _MIXTURE_SPLIT.split(body) if p.strip()]
    if len(parts) < 2:
        n = normalise(s)
        return (n,) if n else ()
    out = []
    for part in parts:
        n = normalise(part)
        if n and n not in out:
            out.append(n)
    return tuple(out)


def is_mixture(raw: str | None) -> bool:
    """True when the stored label names more than one species."""
    return len(parse_mixture(raw)) > 1


#: Canonical tokens whose inducer is a metal ION or metalloid oxyanion -- i.e. what the coordination
#: gate exists to detect. Written out as a SET of canonical tokens rather than pattern-matched, because
#: pattern-matching element symbols inside chemistry strings is the trap handoff §5.4 records twice.
#: As(III)/MAs(III)/antimonite are metalloids and are included: ArsR-family sensors of them carry a
#: thiolate site and the pipeline classes them `metal`, which is what this set is used to score.
_METAL_TOKENS = frozenset({
    "Zn2+", "Co2+", "Fe2+", "Ni2+", "Hg2+", "Cu+", "Mn2+", "Cd2+", "Pb2+", "Au+", "Ag+",
    # metal / metalloid oxyanions -- oxygen in the formula does not make them organic
    "MoO42-", "WO42-", "CrO42-", "AsO43-", "SeO32-", "SeO42-", "VO43-", "VO3-",
    "SbO43-", "TeO32-", "TeO42-", "As(III)", "MAs(III)", "antimonite",
})

#: Canonical tokens whose inducer is a reactive species -- a separate class from both metal and organic
#: (the pipeline's `redox`), and NOT a metal-site call even when the sensor ligates an [Fe-S] cluster.
_REDOX_TOKENS = frozenset({"persulfide", "H2O2", "organic hydroperoxide", "NO"})


#: The CLASS-LEVEL calls `_finalize` emits when the class is settled but the ion is not. They name a
#: class rather than a species, so they cannot be looked up as tokens -- and the fallthrough below
#: treats every unrecognised token as "organic", which turned the pipeline's own
#: "divalent metal (ion unresolved)" into an organic call. That is the answer inverted, on exactly
#: the labels the gate produces when it is working correctly.
_CLASS_LEVEL_CALLS = {
    "divalent metal (ion unresolved)": "metal",
    "non-metal effector (undetermined)": "organic",
    "reactive oxygen/nitrogen species (redox)": "redox",
}


def inducer_class_of(raw: str | None) -> str | None:
    """'metal' | 'redox' | 'organic' for one stored inducer label, or None when it names nothing.

    A MIXTURE resolves only when its components agree; `Cu+/persulfide` spans two classes and returns
    None rather than silently picking one. That is the honest answer for a cluster whose anchors
    genuinely disagree about the class, and it keeps this out of the business of deciding chemistry.
    """
    if raw and raw.strip().lower() in _CLASS_LEVEL_CALLS:
        return _CLASS_LEVEL_CALLS[raw.strip().lower()]
    parts = parse_mixture(raw)
    if not parts:
        return None
    classes = {("metal" if p in _METAL_TOKENS else
                "redox" if p in _REDOX_TOKENS else "organic") for p in parts}
    return classes.pop() if len(classes) == 1 else None


def is_metal(raw: str | None) -> bool:
    """True when the label names a metal/metalloid inducer and nothing else."""
    return inducer_class_of(raw) == "metal"


def display_label(top: str | None, candidates=()) -> str:
    """Return presentation text without changing the stored, class-bearing ``top`` label.

    A class-level placeholder with at least two named candidates represents a real ambiguity rather
    than a missing result (for example, the RcnR/FrmR clade's Co2+ versus formaldehyde split). Expose
    that shortlist to readers. All other calls, including genuine co-sensed mixtures, pass through
    verbatim. Callers must continue to pass ``top`` -- never this string -- to ``inducer_class_of``.
    """
    if not top:
        return top or ""
    if isinstance(candidates, str):
        candidates = (candidates,)
    named = []
    seen = set()
    for candidate in candidates or ():
        if not candidate or is_abstention(candidate):
            continue
        text = str(candidate).strip()
        if not text or text.lower() in _CLASS_LEVEL_CALLS:
            continue
        key = normalise(text) or text
        if key not in seen:
            seen.add(key)
            named.append(text)
    if top.strip().lower() in _CLASS_LEVEL_CALLS and len(named) >= 2:
        return (f"inconclusive: {' or '.join(named)} "
                "(a divergent clade lumps distinct sensors; not resolvable from sequence -- "
                "manual curation required)")
    return top


def distinct(values) -> list[str]:
    """Canonical, de-duplicated, order-stable inducers from raw strings; abstentions removed."""
    out = []
    for v in values:
        n = normalise(v)
        if n and n not in out:
            out.append(n)
    return out


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    # the spellings that are actually costing confidence today
    assert normalise("Ni") == normalise("Ni2+") == normalise("Ni(II)") == "Ni2+"
    assert normalise("Cobalt") == normalise("Co2+") == "Co2+"
    assert normalise("Copper") == normalise("Cu+") == "Cu+"
    assert normalise("Zn(II)") == normalise("Zn2+") == "Zn2+"
    assert normalise("Cd(II)") == normalise("Cd2+") == "Cd2+"
    assert normalise("Fe") == normalise("Fe2+") == "Fe2+"
    assert normalise("Manganese") == "Mn2+"
    assert normalise("Pb(II)") == "Pb2+"
    assert normalise("Au(I)") == "Au+"
    assert normalise("H2S / persulfide (RSS)") == normalise("persulfide") == "persulfide"
    assert normalise("MAs(III) (methylarsenite / monomethylarsenous acid)") == "MAs(III)"

    # the two abstention sentences in the KB
    assert normalise("none identified (possibly NodD protein competition)") is None
    assert normalise("none (glutamine synthetase feedback state)") is None
    assert is_abstention("none") and is_abstention("") and is_abstention(None)

    # genuinely different species must NOT be merged
    assert normalise("As(III)") != normalise("MAs(III)"), "arsenite is not methylarsenite"
    assert normalise("As(III)") != normalise("antimonite"), "arsenite is not antimonite"
    assert normalise("α-D-galacturonate") != normalise("α-D-glucuronate")
    assert normalise("Fe2+") != normalise("Mn2+"), "spelling only -- chemistry is the mixture design"

    # handoff §5.4: no substring matching on element names
    assert normalise("steroid-CoA") == "steroid-CoA", "'CoA' contains 'Co' but is not cobalt"
    assert normalise("nitrogen") == "nitrogen", "'nitrogen' contains 'ni' but is not nickel"
    assert normalise("non-metal effector") == "non-metal effector"

    # cross-source spellings found by analysis/family_kb/vocab_audit.py
    assert normalise("As3+") == normalise("As(III)") == normalise("arsenite") == "As(III)",         "the SSN table writes As3+ where the anchor KB writes As(III)"
    assert normalise("Sb3+") == normalise("antimonite") == "antimonite"
    assert normalise("Citrate") == normalise("citrate") == "citrate"
    assert normalise("Formaldehyde") == normalise("formaldehyde") == "formaldehyde"

    # the generic oxidation-state fallback, and its guard rails
    assert normalise("Ni3+") == "Ni2+", "an unusual oxidation state still identifies the element"
    assert normalise("Cobalt") == "Co2+"
    assert normalise("CoA") == "CoA", "a two-letter element match must not fire inside a word"
    assert normalise("Na+") == "Na+", "an element the table does not carry is passed through"

    # unknown-but-real values survive
    assert normalise("tetracycline") == "tetracycline"

    assert distinct(["Ni", "Ni2+", "none identified (x)", "Zn2+", ""]) == ["Ni2+", "Zn2+"]

    # --- mixtures: the SSN table stores co-sensing as one flattened string --------------------------
    assert parse_mixture("Ni2+/Co2+") == ("Ni2+", "Co2+")
    assert parse_mixture("Cd2+/Pb2+") == ("Cd2+", "Pb2+")
    assert parse_mixture("Cu+/persulfide") == ("Cu+", "persulfide")
    assert is_mixture("Ni2+/Co2+") and not is_mixture("Zn2+")

    # a single species is a 1-tuple, so callers need no special case
    assert parse_mixture("Zn2+") == ("Zn2+",)
    assert parse_mixture("none identified (x)") == ()

    # three things splitting must NOT invent
    assert parse_mixture("MAs(III) (methylarsenite / monomethylarsenous acid)") == ("MAs(III)",),         "the slash is inside a parenthetical gloss, not between two species"
    assert parse_mixture("salicylate/multidrug") == ("salicylate",),         "a ligand plus the functional class it belongs to is one claim, not two"
    assert parse_mixture("H2S / persulfide (RSS)") == ("persulfide",),         "one RSS pool written verbosely"
    assert parse_mixture("steroid-CoA/organic") == ("steroid-CoA",)

    # --- inducer CLASS, used to build metal/non-metal truth sets for the gate harness -------------
    assert inducer_class_of("Zn2+") == inducer_class_of("Ni2+/Co2+") == "metal"
    assert inducer_class_of("As(III)") == "metal", "a metalloid thiolate sensor is a metal-site call"
    # metal/metalloid OXYANIONS are metal-class, not organic -- oxygen in the formula is not organic
    assert inducer_class_of("MoO42-") == inducer_class_of("MoO4-") == \
        inducer_class_of("molybdate") == "metal", "molybdate is a metal call"
    assert inducer_class_of("WO42-") == inducer_class_of("tungstate") == "metal"
    assert inducer_class_of("CrO4") == inducer_class_of("chromate") == "metal"
    assert inducer_class_of("arsenate") == inducer_class_of("vanadate") == "metal"
    assert inducer_class_of("AsO3") == "metal", "arsenite as an oxyanion is still a metalloid call"
    assert normalise("VO3-") != normalise("VO43-"), "meta- and orthovanadate are distinct species"
    assert inducer_class_of("formaldehyde") == inducer_class_of("salicylate") == "organic"
    assert inducer_class_of("persulfide") == inducer_class_of("H2O2") == "redox"
    assert inducer_class_of("Cu+/persulfide") is None, "a mixture spanning classes must not pick one"
    assert inducer_class_of("none identified (x)") is None
    assert is_metal("Cd2+/Pb2+") and not is_metal("formaldehyde") and not is_metal("Cu+/persulfide")

    print("OK: spellings merge, abstentions vanish, distinct species stay distinct, no substrings, "
          "compound labels resolve to their components without inventing any, and the inducer class "
          "abstains on a mixture that spans classes.")


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        _demo()
    else:
        print(__doc__)
        sys.exit(2)
