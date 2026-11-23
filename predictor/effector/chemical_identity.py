"""Conservative identity metadata; legacy label scoring remains in inducer_vocab.

Keys are spelling comparisons, not inferred chemical equivalence. Unmapped names
stay inspectable and cannot validate a new participant's exact identity.
"""
from __future__ import annotations

import re

from predictor.effector import inducer_vocab as vocab

_TRANSPORT = re.compile(r"\s*\(\s*(in|out)\s*\)\s*$", re.I)
_ION = re.compile(r"([A-Z][a-z]?)\s*(?:\(\s*(\d*)([+-])\s*\)|(\d*)([+-])|\((I|II|III|IV|V|VI|VII|VIII)\))$", re.I)
_METAL_ELEMENTS = {"Li", "Na", "K", "Mg", "Ca", "Al", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                   "Zn", "Ag", "Cd", "Hg", "Pb", "Au", "Mo", "W", "As", "Sb", "Se", "Te", "V"}
_ROMAN = {v: i for i, v in enumerate(("I", "II", "III", "IV", "V", "VI", "VII", "VIII"), 1)}
_FORMULA = re.compile(r"(?:[A-Z][a-z]?\d*)+(?:\(\d*[+-]\)|[+-])?$")
_EXPLICIT_STATE = re.compile(r"(?:[+-]|\([^)]*(?:[+-]|[IVXivx])[^)]*\))$")
SENSOR_CLASSES = frozenset({"metal", "redox", "organic"})


def _monatomic_ion(text):
    """A two-letter element symbol cannot be two uppercase formula atoms (CO)."""
    match = _ION.fullmatch(text)
    if match and len(match.group(1)) > 1 and match.group(1).isupper():
        return None
    return match


def _bare_metal_symbol(text):
    """A case-exact bare element symbol from `_METAL_ELEMENTS` (`Fe`, `Co`), never `CO`, `fe` or `C`.

    Case matters for the same reason it matters in `_monatomic_ion`: `CO` is carbon monoxide."""
    return text if isinstance(text, str) and text in _METAL_ELEMENTS else None


def _ion_element(raw):
    ion = _monatomic_ion(transport_name(raw))
    return ion.group(1).title() if ion else None


def subsume_bare_elements(names, context=None):
    """Drop a bare element symbol when the set already holds an explicit ion of that element (EXP070).

    `Fe` says less than `Fe2+` and never contradicts it, so `{Fe, Fe2+}` is one species written two
    ways. Explicit ions are kept exactly as written -- `{Fe2+, Fe3+}` stays two and `{Fe, Fe3+}` becomes
    `Fe3+`, so no oxidation state is ever invented. `context` (default: `names`) supplies the explicit
    ions, so a sub-group can be reduced against the whole set it came from. Equality (`identity_key`)
    is deliberately unchanged: this is a set-assembly rule, not a new identity."""
    names = tuple(names or ())
    ion_elements = {e for e in map(_ion_element, names if context is None else context) if e}
    return tuple(n for n in names if _bare_metal_symbol(transport_name(n)) not in ion_elements)


def transport_name(raw):
    """Remove exactly one terminal transport tag, preserving all other chemistry."""
    return _TRANSPORT.sub("", (raw or "").strip()).strip()


def compartment(raw):
    match = _TRANSPORT.search(raw or "")
    return match.group(1).lower() if match else ""


def plain_metal_ion(raw):
    """`Cu(+)(in)` -> `Cu+`, `Zn(2+)(out)` -> `Zn2+`; None for anything but a monatomic metal ion.

    The ion as written by a transport equation, without its compartment tag and in the plain
    element-charge spelling the rest of the vocabulary uses. The charge is the one written; none is
    inferred, and a Roman-numeral state (`Ni(II)`) is returned untouched apart from the tag."""
    text = transport_name(raw)
    ion = _monatomic_ion(text)
    if not ion or chemistry(raw)[0] != "metal":
        return None
    element, n1, s1, n2, s2, roman = ion.groups()
    if roman:
        return text
    return f"{element}{n1 or n2 or ''}{s1 or s2}"


def identity_key(raw):
    text = transport_name(raw)
    if not text:
        return None
    ion = _monatomic_ion(text)
    # Formula case matters (CO versus Co); preserve the written formula/charge
    # before aliases that historically filled in or changed oxidation states.
    if _FORMULA.fullmatch(text) and (not ion or len(re.findall(r"[A-Z][a-z]?", text)) > 1):
        return f"formula:{text}"
    if ion:
        element, n1, s1, n2, s2, roman = ion.groups()
        charge = str(_ROMAN[roman.upper()]) if roman else n1 or n2 or "1"
        return f"{element}{charge}{s1 or s2 or '+'}".casefold()
    if vocab.is_abstention(text):
        return None
    if _EXPLICIT_STATE.search(text):
        return f"explicit:{text}"
    # Reuse explicit entries only; never the oxidation-erasing normalise fallback.
    mapped = vocab._CANONICAL.get(text.casefold(), text)
    return identity_key(mapped) if mapped != text else mapped.casefold()


def unique_names(names):
    """One representative per safe identity; retain raw spellings in source records."""
    out, seen = [], set()
    for name in names:
        key = identity_key(name)
        if key and key not in seen:
            out.append(name)
            seen.add(key)
    return tuple(out)


def exact_eligible(call):
    """One exact-identity rule for fusion, support and every fallback search."""
    if not call.ligand or (call.confidence or 0) <= 0 or call.mixture_kind == "unresolved":
        return False
    if call.identity_status and call.identity_status not in {"resolved", "co_sensed"}:
        return False
    if call.mixture_kind == "co_sensed":
        parts = unique_names(tuple(call.candidates or ()) or parse_components(call.ligand))
        if len(parts) > 1:
            return all(chemistry(p)[1] == "resolved" for p in parts)
    return chemistry(call.ligand)[1] == "resolved"


def chemistry(raw):
    """Return (canonical class or None, identity status, explicit evidence).

    No structural formula guessing, default organic, or substring classification.
    An explicit ion notation supplies charge identity even when absent in the old
    vocabulary. Bare element aliases remain legacy mappings, not new classifiers.
    """
    text = transport_name(raw)
    if not text or (vocab.is_abstention(text) and not _monatomic_ion(text)):
        return None, "unknown", {}
    if text.casefold() in vocab._CLASS_LEVEL_CALLS:
        return vocab._CLASS_LEVEL_CALLS[text.casefold()], "class_only", {"basis": "explicit_class_label"}
    from predictor.effector.ligand import CURRENCY, _is_macromolecule
    if text.casefold() == "ah2":
        return None, "generic", {"basis": "generic_reducing_equivalent"}
    ion = _monatomic_ion(text)
    mapped = vocab._CANONICAL.get(text.casefold())
    if mapped and identity_key(text) != identity_key(mapped):
        mapped = None  # a legacy alias cannot validate a differently written state/formula
    if text in vocab._REDOX_TOKENS or mapped in vocab._REDOX_TOKENS:
        return "redox", "resolved", {"basis": "existing_explicit_vocabulary"}
    currency_ion = bool(ion and identity_key(text) in {identity_key(x) for x in CURRENCY if _monatomic_ion(x)})
    if text.casefold() in CURRENCY or currency_ion or _is_macromolecule(text):
        return None, "excluded", {"basis": "currency_or_polymer"}
    if ion and ion.group(1).title() in _METAL_ELEMENTS:
        return "metal", "resolved", {"basis": "explicit_ion_notation", "identity_key": identity_key(text)}
    bare = _bare_metal_symbol(text)
    if bare:
        # EXP070: the class is certain, the oxidation state is not asserted. An element whose ion is a
        # currency metabolite (K, Na, Mg) is excluded exactly as its written ion is.
        if bare in {_ion_element(x) for x in CURRENCY}:
            return None, "excluded", {"basis": "currency_or_polymer"}
        return "metal", "resolved", {"basis": "bare_element_symbol", "element": bare}
    if text in vocab._METAL_TOKENS:
        return "metal", "resolved", {"basis": "existing_canonical_vocabulary"}
    if mapped:
        cls = ("metal" if mapped in vocab._METAL_TOKENS else
               "redox" if mapped in vocab._REDOX_TOKENS else
               "organic" if mapped not in {"CO2"} else None)
        return cls, "resolved" if cls else "unmapped", {"basis": "existing_explicit_vocabulary"}
    return None, "unmapped", {}


def parse_components(raw):
    """Preserve states/names while using existing explicit mixture exceptions."""
    if not raw or vocab.is_abstention(raw):
        return ()
    text = raw.strip()
    if text.casefold() in vocab._NOT_A_MIXTURE:
        return (vocab._NOT_A_MIXTURE[text.casefold()],)
    # Slashes inside parentheses are part of a name/gloss, never a separator.
    parts, start, depth = [], 0, 0
    for i, char in enumerate(text):
        depth += (char == "(") - (char == ")")
        if char == "/" and depth == 0:
            parts.append(text[start:i].strip())
            start = i + 1
    parts.append(text[start:].strip())
    return tuple(dict.fromkeys(p for p in parts if p and not vocab.is_abstention(p)))


def collection_class(names):
    records = [chemistry(n) for n in names]
    classes = {r[0] for r in records if r[0]}
    known = bool(records) and all(r[0] for r in records)
    return (next(iter(classes)) if known and len(classes) == 1 else None,
            "conflicting" if len(classes) > 1 else "supported" if known else "unknown",
            records)


def anchor_names(anchor_set):
    """Raw literature-stated identities without the legacy charge-erasing view."""
    out, seen = [], set()
    for anchor in anchor_set.anchors:
        if not anchor.is_literature:
            continue
        for name in parse_components(anchor.inducer):
            key = identity_key(name)
            if key and key not in seen:
                out.append(name)
                seen.add(key)
    return tuple(out)


def class_from_record(record, legacy_top=None):
    """Read explicit metadata; return (class, provenance), marking legacy fallback."""
    if "sensor_class" in record:
        return record.get("sensor_class"), "explicit"
    return vocab.inducer_class_of(record.get("top", legacy_top)), "legacy_label_fallback"


def display_label(top, candidates=()):
    """New structured display, preserving charge-distinct alternatives."""
    names, seen = [], set()
    for name in candidates or ():
        key = identity_key(name)
        if key and key not in seen:
            names.append(name)
            seen.add(key)
    if not top:
        return f"unresolved: {' or '.join(names)}" if names else ""
    top_class, top_status, _ = chemistry(top)
    if top_status == "class_only":
        # EXP071: a class placeholder with members of ITS OWN class is a shortlist -- say so and list
        # only those (a Cu+/Au(I) metal shortlist must not read as "Cu+ or Au(I) or <operon organics>").
        # A lone same-class member is an alternative the policy declined to make the headline, not a
        # shortlist, so it keeps the previous display.
        same = [n for n in names if chemistry(n)[0] == top_class]
        if len(same) > 1:
            label = {"metal": "metal-ion", "organic": "non-metal"}.get(top_class, top_class)
            rest = [n for n in names if n not in same]
            return (f"{label} shortlist: {' or '.join(same)}"
                    + (f"; other proposals: {', '.join(rest)}" if rest else ""))
    if top_status == "class_only" and len(names) > 1:
        return f"inconclusive: {' or '.join(names)} (identity unresolved)"
    return top
