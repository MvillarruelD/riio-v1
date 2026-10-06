"""The metal-signature catalogue and its generated reference.

Two things are pinned here, and they exist for the same reason: a coordination claim written by hand
drifts away from the data it describes, and nothing notices. `structure/metal_site.py` stated the NikR
H-x-H-x(4,7)-C site occurs in "<=3.2 % of every other family (Fur 3.2 %)"; against the member DB the
package ships it occurs in 60.6 % of Fur. That claim was heading into a Methods section.

So: the reference is GENERATED from the catalogue plus the shipped member DBs, and these tests fail if
the generated files are stale or if a motif's scope contradicts its own measured specificity.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from predictor import resources
from predictor.annotate import ssn_clusters as ssn
from predictor.structure import metal_motifs as MM

REPO = Path(__file__).resolve().parents[1]
JSON_REF = REPO / "predictor" / "data" / "refs" / "metal_signatures.json"
MD_REF = REPO / "docs" / "METAL_SIGNATURES.md"


@pytest.fixture(scope="module")
def ref() -> dict:
    assert JSON_REF.is_file(), f"{JSON_REF} is missing; run tools/build_metal_reference.py"
    return json.loads(JSON_REF.read_text(encoding="utf-8"))


def _members(tag: str) -> list[str]:
    p = resources.SSN_DATABASE / f"ssn_members_{tag}.fasta"
    out, buf, started = [], [], False
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if started:
                out.append("".join(buf).upper())
            started, buf = True, []
        else:
            buf.append(ln.strip())
    if started:
        out.append("".join(buf).upper())
    return out


@pytest.mark.released_data
def test_generated_reference_is_current():
    """`--check` reruns the measurement and compares. Stale files mean the shipped numbers are fiction."""
    r = subprocess.run([sys.executable, str(REPO / "tools" / "build_metal_reference.py"), "--check"],
                       capture_output=True, text=True, cwd=REPO)
    assert r.returncode == 0, (
        "the metal-signature reference is stale -- regenerate with "
        f"`python tools/build_metal_reference.py`\n{r.stdout}{r.stderr}")


@pytest.mark.released_data
def test_reference_covers_every_catalogue_entry(ref):
    assert {m["name"] for m in ref["motifs"]} == {m.name for m in MM.MOTIFS}
    assert ref["corpus"]["n_total"] > 100_000, ref["corpus"]


@pytest.mark.released_data
def test_nikr_site_is_family_scoped_because_it_fires_on_fur(ref):
    """The specific regression. This motif is only meaningful inside NikR, and the data says why."""
    m = next(x for x in ref["motifs"] if x["name"] == "NikR_His3Cys")
    p = m["prevalence"]
    assert m["scope"] == "family", "NikR_His3Cys must never be promoted to global scope"
    assert p["in_family"] > 0.80, p["in_family"]
    assert p["per_family"]["Fur"] > 0.50, (
        "the point of this test is that the NikR spacing is COMMON in Fur; if that changed, the "
        f"comment in metal_site.py must change with it (measured {p['per_family']['Fur']:.3f})")
    assert not p["global_safe"]

    # and the live gate must not be reachable from a Fur query, which is what makes it harmless
    from predictor.structure import metal_site
    fur_like = "MSEHHHDHSHHCEECGKTLHFDDAGQLKAHCEEHHHHHDHHHEHHHHDHHHEH"
    assert metal_site.coordination_gate(fur_like, "Fur")["gate"] == "Fur-His/Cys"
    assert metal_site.coordination_gate(fur_like, "NikR")["gate"] == "NikR-HxHxC"


@pytest.mark.released_data
def test_global_scope_is_justified_by_the_measurement(ref):
    """A motif may only claim `global` scope if it is measurably rare in every other family."""
    thr = ref["global_safe_threshold"]
    for m in ref["motifs"]:
        if not m["implemented"] or m["scope"] != "global":
            continue
        p = m["prevalence"]
        assert p["global_safe"], (
            f"{m['name']} declares global scope but reaches {p['max_other_family']:.3f} in "
            f"{p['max_other_family_name']} (limit {thr})")


@pytest.mark.released_data
def test_ion_naming_motifs_concentrate_in_one_clade(ref):
    """A motif may only NAME an ion if it is concentrated, not merely present.

    `ArsR_a3_CVC` and `CsoR_WXYZ_metal_site` were both written with `names_ion=True` from the review's
    wording and both had to be demoted once measured: C-(V/A)-C fires on 41 % of the ZINC clade ArsR_c2
    as well as the arsenic clade, and the CsoR W-X-Y-Z frame scores as highly on the Ni/Co sensor InrS
    as on the Cu(I) clades. This test is what caught them.
    """
    for m in ref["motifs"]:
        if not m["names_ion"]:
            continue
        p = m["prevalence"]
        assert m["implies"], f"{m['name']} names an ion but lists none"
        assert p["top_clade_rate"] and p["top_clade_rate"] >= 0.50, (
            f"{m['name']} names {m['implies']} but its best clade is only "
            f"{p['top_clade_rate']} ({p['top_clade']}) -- too diffuse to name an ion")


@pytest.mark.released_data
@pytest.mark.parametrize("name,tag,clade,floor", [
    # each ion-naming motif, against the clade it is supposed to identify
    ("ArsR_a5_DxHx10Hx2HE", "ArsR", "ArsR_c5", 0.90),   # alpha5 Zn/Ni site
    ("MerR_ZntR_CCx3Hx4C", "MerR", "MerR_9", 0.50),     # ecZntR, Zn
    ("CopY_CxC", "CopY", "CopY_c1", 0.90),              # CopY Cu(I), vs the BlaI/MecI branch
])
def test_ion_motifs_hit_their_clade_in_the_shipped_data(name, tag, clade, floor):
    """Recomputed from the FASTA, not read from the JSON -- so a wrong generator cannot hide a wrong motif."""
    motif = MM.by_name(name)
    assert motif is not None and motif.names_ion
    rx = re.compile(motif.pattern)
    p = resources.SSN_DATABASE / f"ssn_members_{tag}.fasta"
    hits = total = 0
    cid, buf = None, []

    def _flush():
        nonlocal hits, total
        if cid == clade and buf:
            total += 1
            if rx.search("".join(buf).upper()):
                hits += 1

    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            _flush()
            cid, buf = ln[1:].strip().split("__", 1)[0], []
        else:
            buf.append(ln.strip())
    _flush()
    assert total > 0, f"no members of {clade} in {p.name}"
    assert hits / total >= floor, f"{name} in {clade}: {hits}/{total} = {hits/total:.3f} < {floor}"


@pytest.mark.released_data
def test_copy_motif_separates_the_copper_clade_from_the_beta_lactam_branch():
    """CopY's `CxC` names Cu(I), and may only do so because the BlaI/MecI clades do not carry it."""
    motif = MM.by_name("CopY_CxC")
    rx = re.compile(motif.pattern)
    by_clade: dict[str, list[int]] = {}
    cid, buf = None, []

    def _flush():
        if cid and buf:
            c = by_clade.setdefault(cid, [0, 0])
            c[1] += 1
            if rx.search("".join(buf).upper()):
                c[0] += 1

    for ln in (resources.SSN_DATABASE / "ssn_members_CopY.fasta").read_text(
            encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            _flush()
            cid, buf = ln[1:].strip().split("__", 1)[0], []
        else:
            buf.append(ln.strip())
    _flush()

    cu = by_clade["CopY_c1"]
    assert cu[0] / cu[1] > 0.90, cu
    for clade in ("CopY_c2", "CopY_c3", "CopY_c4"):        # the beta-lactam (BlaI/MecI) branch
        h, n = by_clade[clade]
        assert h / n < 0.05, f"{clade}: {h}/{n} carry the CopY Cu motif; it can no longer name Cu+"


def test_families_and_ions_are_canonical():
    for m in MM.MOTIFS:
        if m.family:
            assert m.family in ssn.FAMILIES, f"{m.name}: {m.family!r} is not a canonical family name"
        for ion in m.implies:
            assert ion in MM.SELECTIVITY, f"{m.name}: {ion} has no selectivity entry"


@pytest.mark.released_data
def test_markdown_reference_states_the_measured_nikr_number():
    """The document a reviewer reads must carry the corrected figure, not the old claim."""
    md = MD_REF.read_text(encoding="utf-8")
    assert "60.6 %" in md or "60.6%" in md, "the NikR/Fur correction must appear in the reference"
    assert "do not edit" in md.lower()


# --------------------------------------------------------------------------- the inducer axis
def _coordination_call(seq: str, family: str):
    """Build the coordination InducerCall exactly as `infer_inducer` source 2 does."""
    from predictor.schema import InducerCall
    from predictor.structure import metal_site
    g = metal_site.coordination_gate(seq, family)
    ions = g.get("implied_ions") or {}
    return g, InducerCall(
        source="coordination", ligand=None, confidence=0.0,
        role=("metal" if g.get("has_site") else ("non-metal" if g.get("has_site") is False else "")),
        candidates=tuple(sorted(ions)), mixture_kind=("unresolved" if len(ions) > 1 else ""),
        evidence={"gate": g.get("gate"), "family": family, "motifs": list(g.get("motifs") or []),
                  "implied_ions": {k: list(v) for k, v in ions.items()},
                  "motif_ion": sorted(ions)[0] if len(ions) == 1 else None})


#: E. coli ZntR, verbatim from the shipped member DB record `MerR_9__P0ACS5`.
ZNTR = ("MYRIGELAKMAEVTPDTIRYYEKQQMMEHEVRTEGGFRLYTESDLQRLKFIRHARQLGFSLESIRELLSIRIDPEHHT"
        "CQESKGIVQERLQEVEARIAELQSMQRSLQRLNDACCGTAHSSVYCSILEALEQGASGVKSGC")


def test_coordination_gate_reports_motifs_without_changing_has_site():
    """The catalogue is ADDITIVE. `has_site` is the gate the pipeline has always applied."""
    from predictor.structure import metal_site
    g = metal_site.coordination_gate(ZNTR, "MerR")
    assert g["has_site"] is True and g["gate"] == "CX7-8C"        # unchanged behaviour
    assert "MerR_ZntR_CCx3Hx4C" in g["motifs"]
    assert g["implied_ions"]["Zn2+"] == ["MerR_ZntR_CCx3Hx4C"]

    # and a family with no branch still gets catalogue fields rather than a KeyError downstream
    g2 = metal_site.coordination_gate(ZNTR, "GntR")
    assert g2["has_site"] is None and "motifs" in g2 and "implied_ions" in g2


def test_motif_corroborates_an_agreeing_clade_without_taking_the_headline():
    from predictor.effector import inducer as IND
    from predictor.schema import InducerCall
    g, cc = _coordination_call(ZNTR, "MerR")
    ssn_call = InducerCall(source="ssn_cluster", ligand="Zn2+", confidence=0.8, role="metal",
                           evidence={"cluster": "MerR_9"})
    r = IND._finalize([ssn_call, cc], g["has_site"], gate_name=g["gate"], sensor_class="metal")
    assert r.top == "Zn2+"
    assert any("agreeing with SSN cluster MerR_9" in n for n in r.notes), r.notes


def test_motif_disagreeing_with_a_clade_is_reported_never_resolved():
    """A sequence motif may not overturn curated clade anchors -- but the reader must see the conflict."""
    from predictor.effector import inducer as IND
    from predictor.schema import InducerCall
    g, cc = _coordination_call(ZNTR, "MerR")
    ssn_call = InducerCall(source="ssn_cluster", ligand="Cd2+", confidence=0.8, role="metal",
                           evidence={"cluster": "MerR_1"})
    r = IND._finalize([ssn_call, cc], g["has_site"], gate_name=g["gate"], sensor_class="metal")
    assert r.top == "Cd2+", "the clade keeps the headline"
    assert set(r.candidates) == {"Cd2+", "Zn2+"}
    assert r.mixture_kind == "unresolved"
    assert any(n.startswith("DISAGREEMENT") for n in r.notes), r.notes


def test_motif_may_name_an_ion_only_when_no_clade_spoke():
    """The one case the axis exists for: 88 of 140 candidates on the last run had no clade."""
    from predictor.effector import inducer as IND
    g, cc = _coordination_call(ZNTR, "MerR")
    r = IND._finalize([cc], g["has_site"], gate_name=g["gate"], sensor_class="metal")
    assert r.top == "Zn2+"
    coord = next(c for c in r.calls if c.source == "coordination")
    assert coord.ligand == "Zn2+" and 0 < coord.confidence < 0.5, (
        "a motif-only call must be marked weaker than a clade call")
    assert any("rests on the coordination signature alone" in n for n in r.notes), r.notes


def test_a_non_diagnostic_signature_supports_the_class_but_names_nothing():
    """CsoR's W-X-Y-Z frame scores as highly on the Ni/Co sensor InrS as on the Cu(I) clades."""
    from predictor.effector import inducer as IND
    csor_like = "MSHSHDHDHSHEHSHCEGCGCSEHSHSHEHDHSHEHGHSHDHSHEHSH"
    g, cc = _coordination_call(csor_like, "CsoR/FrmR")
    assert not (g.get("implied_ions") or {}), "this frame must not name an ion"
    r = IND._finalize([cc], g["has_site"], gate_name=g["gate"], sensor_class="metal")
    assert r.top != "Cu+", "a non-diagnostic signature must not produce an ion headline"
