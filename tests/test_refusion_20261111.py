"""Four fusion defects found in the 2026-11-10 final rerun (RUN_DIFF.md), fixed on `refusion-20261111`.

Each test rebuilds the calls of a real bundle from that run and checks the final consensus.
"""
from predictor.effector import chemical_identity as ID
from predictor.effector import inducer as IND
from predictor.schema import InducerCall as Call

NEG_GATE = "metalloregulator family (MerR) but its coordination gate is NEGATIVE (CX7-8C) -> non-metal member"


def _abstain(source, role=""):
    return Call(source, None, 0.0, role)


def test_negative_gate_alone_is_a_non_metal_call():
    """E. coli BluR: MerR, gate negative, nothing else. Was: no headline and no class."""
    calls = [_abstain("ssn_cluster"), _abstain("coordination", "non-metal"), _abstain("metalnet", "non-metal")]
    cons = IND._finalize(calls, False, gate_name="CX7-8C", gate_negative_reason=NEG_GATE)
    assert cons.top == "non-metal effector (undetermined)"
    assert cons.sensor_class == "organic" and cons.class_status == "supported"
    # without the gate's reason (a caller with no family verdict) the old abstention stands
    assert IND._finalize(calls, False, gate_name="CX7-8C").top is None


def test_negative_gate_does_not_overrule_a_metal_clade():
    """M. tuberculosis KmtR: gate negative, clade set Zn/Ni/Cd with role metal. Stays metal."""
    clade = Call("ssn_cluster", None, 0.99, "metal", candidates=("Zn2+", "Ni(II)", "Cd2+"),
                 mixture_kind="unresolved", evidence={"cluster": "ArsR_x"})
    cons = IND._finalize([clade, _abstain("coordination", "non-metal"), _abstain("metalnet", "non-metal")],
                         False, gate_name="ArsR-metal", gate_negative_reason=NEG_GATE)
    assert cons.sensor_class == "metal"
    assert cons.top is None and set(cons.candidates) == {"Zn2+", "Ni(II)", "Cd2+"}


def test_non_metal_clade_outranks_the_family_prior():
    """M. tuberculosis BlaI: CopY family prior says metal, clade CopY_c3 says beta-lactam."""
    clade = Call("ssn_cluster", "beta-lactam", 1.0, "non-metal", evidence={"cluster": "CopY_c3"})
    calls = [clade, _abstain("coordination"), _abstain("metalnet", "non-metal")]
    cons = IND._finalize(calls, None, sensor_class="metal",
                         class_reason="metalloregulator family (CopY); coordination gate silent")
    assert cons.sensor_class == "organic" and cons.class_status == "supported"
    assert cons.top == "non-metal effector (undetermined)" and "beta-lactam" in cons.candidates
    # a MetalNet site is direct evidence about the protein: the prior is then not bare, and stands
    with_site = [clade, _abstain("coordination"), _abstain("metalnet", "metal")]
    assert IND._finalize(with_site, None, sensor_class="metal").sensor_class == "metal"


def test_negative_gate_drops_the_metal_member_of_an_unclassed_clade_set():
    """Rv2034: clade ArsR_c9 `SpdC or Ni(II)`, gate negative. Was: class metal."""
    clade = Call("ssn_cluster", None, 0.99, "metal", candidates=("SpdC", "Ni(II)"),
                 mixture_kind="unresolved", evidence={"cluster": "ArsR_c9"})
    calls = [clade, _abstain("coordination", "non-metal"), _abstain("metalnet", "non-metal")]
    cons = IND._finalize(calls, False, gate_name="ArsR-metal", gate_negative_reason=NEG_GATE)
    assert cons.sensor_class == "organic"
    assert cons.top == "non-metal effector (undetermined)" and tuple(cons.candidates) == ("SpdC",)
    assert cons.calls[0].candidates == ("SpdC", "Ni(II)"), "the raw clade set is retained"
    # a positive gate says nothing against an unclassed member: untouched
    kept, split = IND._gate_split_ssn(calls, True)
    assert split is None and kept[0].candidates == ("SpdC", "Ni(II)")


def test_transported_metal_ion_is_written_plainly():
    assert ID.plain_metal_ion("Cu(+)(in)") == "Cu+"
    assert ID.plain_metal_ion("Zn(2+)(out)") == "Zn2+"
    assert ID.identity_key(ID.plain_metal_ion("Cu(+)(in)")) == ID.identity_key("Cu(+)(in)")
    for not_an_ion in ("H(+)(in)", "L-serine(in)", "CO", "beta-lactam", None):
        assert ID.plain_metal_ion(not_an_ion) is None
