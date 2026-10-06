"""Sensor-class regressions: the metal / redox / organic call and the inducer fusion guards.

Every case here is a bug that actually shipped. The pattern behind most of them is the same: a WEAK signal
(a bare family prior, a loose Cys-density heuristic, a cluster's plurality label) was allowed to outrank a
STRONG one (a specific cofactor motif, a positive coordination site). `metal_site.inducer_class` and
`inducer._finalize` now order their checks by strength of evidence, and these tests pin that ordering.

Offline and deterministic -- no network, no external engines, no reference data beyond the repo.
"""
import pytest

from predictor.effector import inducer as IND
from predictor.schema import InducerCall
from predictor.structure import metal_site as MS

# --- real E. coli sequences (UniProt), the proteins whose calls regressed ---------------------------
RCNR = ("MSHTIRDKQKLKARASKIQGQVVALKKMLDEPHECAAVLQQIAAIRGAVNGLMREVIKGHLTEHIVHQGDELKREEDLDVVLKVLDSYIK")
FRMR = ("MPSTPEEKKKVLTRVRRIRGQIDALERSLEGDAECRAILQQIAAVRGAANGLMAEVLESHIRETFDRNDCYSREVSQSVDDTIELVRAYLK")

# --- synthetic constructs isolating one signature each ----------------------------------------------
ISCR_LIKE = "M" + "A" * 88 + "CAAAAACAAAAAC" + "A" * 60      # Rrf2 [Fe-S] triad C-x(5)-C-x(5)-C
NIKR_LIKE = "M" + "A" * 84 + "HVHINHDDC" + "A" * 38          # NikR His87/His89/Cys95 Ni site
SOXR_LIKE = "M" + "A" * 80 + "CAACACAAAAAC" + "A" * 20       # SoxR [2Fe-2S] ferredoxin motif
MERR_METAL = "M" + "A" * 70 + "C" + "X" * 7 + "C" + "A" * 10  # C-terminal CX7C
MERR_NONMETAL = "M" + "A" * 60 + "L" * 60                     # no cysteines at all
CYS_CLUSTER = "M" + "A" * 30 + "CAACAAC" + "A" * 40           # 3 Cys in 7 residues


class TestRedoxFamilies:
    """Rrf2 members ligate an [Fe-S] cluster or read thiol status -- they never sense a free metal ion.

    The 2026 12-family enrichment filed Rrf2/IscR/NsrR under SENSOR_FAMILIES, so the family prior fired
    before any redox evidence was consulted. Worse, the resulting sensor_class="metal" then tripped the
    metalloregulator guard in _finalize and actively demoted correct non-metal headlines.
    """

    @pytest.mark.parametrize("family", ["Rrf2", "IscR", "NsrR", "RsrR"])
    def test_rrf2_family_is_redox_never_metal(self, family):
        assert MS.inducer_class(ISCR_LIKE, family)[0] == "redox"

    def test_rrf2_has_no_metal_ion_site_to_gate_on(self):
        assert MS.coordination_gate(ISCR_LIKE, "Rrf2")["has_site"] is False

    def test_fes_triad_detected(self):
        assert MS.fes_cluster_triad(ISCR_LIKE)

    def test_rrf2_suppresses_ligify(self):
        # operon chemistry is not the effector of a cluster-occupancy sensor
        assert MS.ligify_weight(ISCR_LIKE, "Rrf2")[0] == pytest.approx(0.10)

    def test_soxr_2fe2s_motif_still_wins(self):
        # the most specific signature must be checked first, ahead of the MerR metal gate
        assert MS.iron_sulfur_cluster(SOXR_LIKE)
        assert MS.inducer_class(SOXR_LIKE, "MerR")[0] == "redox"


class TestNegativeGateOutranksFamilyPrior:
    """The coordination gate exists to catch non-metal members INSIDE a metalloregulator family.

    Its verdict was being discarded in favour of the family prior, so a gate-negative MerR was still
    called metal -- defeating the purpose of the gate (plan SS6: the TnrA/GlnR/NmlR and BmrR/MtaN clades).
    """

    def test_gate_negative_merr_is_not_metal(self):
        assert MS.coordination_gate(MERR_NONMETAL, "MerR")["has_site"] is False
        assert MS.inducer_class(MERR_NONMETAL, "MerR")[0] == "organic"

    def test_gate_positive_merr_is_metal(self):
        assert MS.inducer_class(MERR_METAL, "MerR")[0] == "metal"

    def test_silent_gate_keeps_the_metal_prior(self):
        # DtxR/MntR's Mn site is sequence-dispersed: the gate cannot see it and returns None, which is
        # NOT the same as a negative verdict. The family prior must still apply.
        assert MS.coordination_gate("M" + "A" * 120, "DtxR/MntR")["has_site"] is None
        assert MS.inducer_class("M" + "A" * 120, "DtxR/MntR")[0] == "metal"

    def test_silent_vs_negative_gate_differ_in_ligify_weight(self):
        assert MS.ligify_weight("M" + "A" * 120, "DtxR/MntR")[0] == pytest.approx(0.30)   # advisory
        assert MS.ligify_weight("M" + "A" * 120, "Fur")[0] == pytest.approx(0.85)         # non-metal member


class TestLigandFamilyOutranksCysHeuristic:
    """">=3 Cys in 16 residues" fires on 17.9% of MarR and 58.1% of ArsR members (measured on the
    vendored SSN member DBs), so it must not override a family assignment. E. coli MarR was called redox
    on the strength of Cys47/51/54, which are structural -- its effector is salicylate.
    """

    def test_organic_family_beats_cys_cluster(self):
        assert MS.reactive_cys_cluster(CYS_CLUSTER)                       # the heuristic does fire
        assert MS.inducer_class(CYS_CLUSTER, "MarR/SlyA")[0] == "organic"  # and is correctly outranked
        assert MS.ligify_weight(CYS_CLUSTER, "LysR-type (LTTR)")[0] == pytest.approx(1.00)

    def test_cys_cluster_still_used_for_unknown_family(self):
        assert MS.inducer_class(CYS_CLUSTER, "")[0] == "redox"
        assert MS.ligify_weight(CYS_CLUSTER, "")[0] == pytest.approx(0.10)

    def test_two_distant_cys_are_not_a_cluster(self):
        assert not MS.reactive_cys_cluster("M" + "C" + "A" * 50 + "C")


class TestCsoRMetalBoundary:
    """CsoR/RcnR/FrmR share one fold and split on which metal ligands survive: CsoR (Cu(I)
    Cys36/His61/Cys65) and RcnR (Ni/Co His3/Cys35/His60/His64) keep them; FrmR has lost His64 and the
    N-terminal His and senses formaldehyde via a Pro2/Cys35 methylene bridge (Osman 2016 Nat Chem Biol
    12:839 -- the gain-of-function FrmR(E64H) variant restores metal sensing).
    """

    def test_rcnr_keeps_the_ligand_set(self):
        assert MS.coordination_gate(RCNR, "CsoR/FrmR")["has_site"] is True
        assert MS.inducer_class(RCNR, "CsoR/FrmR")[0] == "metal"

    def test_frmr_has_lost_it(self):
        assert MS.coordination_gate(FRMR, "CsoR/FrmR")["has_site"] is False
        assert MS.inducer_class(FRMR, "CsoR/FrmR")[0] == "organic"


class TestNikRSite:
    """NikR's square-planar Ni(II) site is His87-x-His89-x(5)-Cys95. Present in 85.5% of the vendored
    NikR member DB and <=3.2% of every other family's -- so the call rests on evidence, not the prior."""

    def test_site_detected_and_used(self):
        g = MS.coordination_gate(NIKR_LIKE, "NikR")
        assert g["has_site"] and g["gate"] == "NikR-HxHxC"
        assert MS.inducer_class(NIKR_LIKE, "NikR")[1].startswith("metal-coordination site")

    def test_absent_site_fails_the_gate(self):
        assert MS.coordination_gate("M" + "A" * 130, "NikR")["has_site"] is False


class TestMetalBoundaryGuard:
    """A POSITIVE coordination gate outranks a non-metal SSN cluster label.

    E. coli RcnR shipped the headline "formaldehyde": it assigns to CsoR_ecFrmR (the rcnR-rich FrmR clade)
    and inherited that label even though its own CsoR-His/Cys site was present AND the regulon voted Co2+.
    The pre-existing guard could not catch it because it only fires for weak-source (regulon/Ligify) tops.
    """

    @staticmethod
    def _rcnr_calls():
        return [
            InducerCall(source="ssn_cluster", ligand="formaldehyde", confidence=1.0, role="non-metal",
                        evidence={"cluster": "CsoR_ecFrmR"}),
            InducerCall(source="coordination", ligand=None, confidence=0.0, role="metal",
                        evidence={"gate": "CsoR-His/Cys"}),
            InducerCall(source="ligify", ligand="Co2+", confidence=0.393, role="metal",
                        evidence={"weight": 1.0, "votes": {"Co2+": 1.0}}),
        ]

    def test_positive_gate_overrides_wrong_cluster_label(self):
        # The override still fires; what it can PRODUCE has changed. Naming Co2+ took a source
        # allowed to supply an ion, and the regulon -- which did that here -- no longer feeds the
        # inducer call. Ligify is deliberately excluded from supplying the ion, so the honest result
        # is the class-level metal call: the gate knows there is a metal site, not which metal.
        got = IND._finalize(self._rcnr_calls(), True, gate_name="CsoR-His/Cys", sensor_class="metal")
        assert got.top == "divalent metal (ion unresolved)"
        assert got.top_display.startswith("inconclusive: formaldehyde or Co2+")
        assert any("metal-boundary" in n for n in got.notes), "the override must be explained in the notes"

    def test_silent_gate_leaves_the_cluster_label_alone(self):
        calls = [c for c in self._rcnr_calls() if c.source != "coordination"]
        assert IND._finalize(calls, None, sensor_class=None).top == "formaldehyde"

    def test_negative_gate_leaves_the_cluster_label_alone(self):
        # FrmR: gate and cluster agree that there is no metal site, and formaldehyde is correct
        calls = [
            InducerCall(source="ssn_cluster", ligand="formaldehyde", confidence=1.0, role="non-metal",
                        evidence={"cluster": "CsoR_ecFrmR"}),
            InducerCall(source="coordination", ligand=None, confidence=0.0, role="non-metal",
                        evidence={"gate": "CsoR-His/Cys"}),
        ]
        assert IND._finalize(calls, False, sensor_class="organic").top == "formaldehyde"

    def test_weak_organic_vote_cannot_hijack_a_metal_sensor(self):
        # the original MntR regression: one "multidrug" regulon substrate vote beat the Mn(II) sensor
        calls = [InducerCall(source="ligify", ligand="multidrug", confidence=0.3, role="",
                             evidence={"weight": 1.0})]
        assert IND._finalize(calls, None, sensor_class="metal").top == "divalent metal (ion unresolved)"

    def test_redox_sensor_keeps_its_headline_through_refusion(self):
        # SoxR reverted from the redox string to the SSN cofactor label "Fe" when the consensus was
        # re-finalized during regulon augmentation without the sensor class
        calls = [InducerCall(source="ssn_cluster", ligand="Fe", confidence=1.0, role="redox",
                             evidence={"cluster": "MerR_4"})]
        assert "redox" in IND._finalize(calls, True, sensor_class=None).top

    def test_never_goes_silent(self):
        assert IND._finalize([], None, sensor_class=None).notes, "must always explain itself"
        assert IND._finalize([], True, sensor_class="metal").top == "divalent metal (ion unresolved)"


class TestAbstentionIsRecordedNotSilent:
    """A clade abstention must travel with the prediction, and must change nothing.

    `homolog_selection` names why an assignment failed (low identity, no hit, no SSN for the
    family). That reason used to be printed and discarded, so the run log knew why 88 of 140
    candidates had no clade and the bundles did not. It is now emitted as a call that votes for
    nothing -- which is only safe if it really is inert.
    """

    @staticmethod
    def _abstention():
        return IND.InducerCall(source="ssn_cluster", ligand=None, confidence=0.0, role="",
                               evidence={"status": "unassigned_low_identity",
                                         "top_identity": 0.42, "support": 0.0})

    def test_it_does_not_move_the_headline(self):
        base = [IND.InducerCall(source="coordination", ligand=None, confidence=0.0, role="metal",
                                evidence={"gate": "thiol_cluster"}),
                IND.InducerCall(source="ligify", ligand="salicylate", confidence=0.4, role="",
                                evidence={"weight": 1.0})]
        for gate, cls in ((True, "metal"), (None, "organic"), (False, "organic"), (None, None)):
            without = IND._finalize(list(base), gate, sensor_class=cls).top
            with_it = IND._finalize(list(base) + [self._abstention()], gate,
                                    sensor_class=cls).top
            assert without == with_it, f"abstention moved the headline at gate={gate}, class={cls}"

    def test_the_reason_survives_into_the_calls(self):
        got = IND._finalize([self._abstention()], None, sensor_class="organic")
        call = next(c for c in got.calls if c.source == "ssn_cluster")
        assert call.evidence["status"] == "unassigned_low_identity"
        assert call.ligand is None and call.confidence == 0.0, "it must not vote"
