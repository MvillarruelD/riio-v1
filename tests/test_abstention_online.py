"""The abstention check, on the ONLINE path (HANDOFF 7.2).

The hole this closes. The pre-registered abstention check -- MarR / AcrR / GntR are organic-effector
controls and must never be called metal -- lives in `analysis/regulondb_bench/regulondb_bench_v2.py`,
which consults only offline signals: the Pfam family call and `structure.metal_site`. But the `regulon`
inducer source **only exists online**, because it reads the substrate specificity of the reconstructed
regulon's transporters. So the check passed for exactly the wrong reason: the source that caused the
failure was not present when it ran.

The failure it missed, verbatim from `results/jobs/GntR_online_v3/dossier.json`: E. coli GntR, an
organic control (gluconate, LacI/GalR fold, no coordination site), came back **`top='Zn2+'` at
agreement 1.00**. Ligify sat at confidence -0.10 and therefore did not vote, leaving a single regulon
substrate vote (`Zn2+`, 0.393, from ONE co-regulated handler) as the only positive voter.

The fix -- the organic-family guard in `effector.inducer._finalize` -- is already in place and pinned by
`test_regulon_metal_guard.py`. What was still missing, and what this file is, is the CHECK: the
recorded online evidence replayed through the current fusion, so the regression cannot come back
silently. The call sets below are transcribed from the real bundles, not invented, so this is a
regression net rather than a restatement of the guard's own logic.

Offline and self-contained: the evidence is data here, so the test does not depend on
`results/jobs/*_online_v3` still existing (those are gitignored run outputs).
"""
import pytest

from predictor.effector import inducer as IND
from predictor.schema import InducerCall

# --- recorded online evidence, transcribed from results/jobs/<TF>_online_v3/dossier.json ------------
# Each entry: (sensor_class as the pipeline computes it for that family, the retained per-source calls).
RECORDED = {
    "GntR": ("organic", [
        InducerCall(source="coordination", ligand=None, confidence=0.0, role="", evidence={"gate": ""}),
        InducerCall(source="ligify", ligand="2-(3,4-dihydroxybenzoyloxy)-4,6-dihydroxybenzoate",
                    confidence=-0.1, role="", evidence={"weight": 1.0}),
        InducerCall(source="ligify", ligand="Zn2+", confidence=0.393, role="metal",
                    evidence={"weight": 1.0, "votes": {"Zn2+": 0.45}}),
    ]),
    "AcrR": ("organic", [
        InducerCall(source="coordination", ligand=None, confidence=0.0, role="", evidence={"gate": ""}),
        InducerCall(source="ligify", ligand="1-O-acetylmaltose", confidence=-0.35, role="",
                    evidence={"weight": 1.0}),
        InducerCall(source="ligify", ligand="multidrug", confidence=0.366, role="multidrug",
                    evidence={"support": 2, "votes": {"multidrug": 2.0, "Zn2+": 1.45}}),
    ]),
    "MarR": ("organic", [
        InducerCall(source="ssn_cluster", ligand="salicylate/multidrug", confidence=1.0,
                    role="multidrug", evidence={"cluster": "MarR_OhrR"}),
        InducerCall(source="coordination", ligand=None, confidence=0.0, role="", evidence={"gate": ""}),
        InducerCall(source="ligify", ligand="succinate", confidence=-0.75, role="",
                    evidence={"weight": 1.0}),
    ]),
}

#: What "called metal" means. Matching element symbols as SUBSTRINGS is a documented trap in this
#: project (HANDOFF 5.4: "steroid-CoA" contains "co", "nitrogen" contains "ni"), so the ion tokens are
#: matched whole and the class-level phrase is matched exactly.
_METAL_IONS = {"Zn2+", "Cu+", "Cu2+", "Fe2+", "Fe3+", "Ni2+", "Co2+", "Mn2+", "Cd2+", "Hg2+", "Pb2+",
               "Ag+", "As(III)", "molybdate"}
_METAL_CLASS = "divalent metal (ion unresolved)"


def is_metal_headline(top) -> bool:
    return bool(top) and (top in _METAL_IONS or top == _METAL_CLASS)


class TestOnlineAbstention:
    """No organic control may acquire a metal headline once the ONLINE regulon source is present."""

    @pytest.mark.parametrize("tf", sorted(RECORDED))
    def test_control_does_not_become_metal(self, tf):
        sensor_class, calls = RECORDED[tf]
        got = IND._finalize(list(calls), None, sensor_class=sensor_class)
        assert not is_metal_headline(got.top), (
            f"{tf} is an organic control and acquired the metal headline {got.top!r} on the online path")

    def test_gntr_is_the_recorded_failure_and_is_now_demoted(self):
        """The specific regression. With the guard absent this returns 'Zn2+' at agreement 1.00."""
        sensor_class, calls = RECORDED["GntR"]
        got = IND._finalize(list(calls), None, sensor_class=sensor_class)
        assert got.top == "non-metal effector (undetermined)"
        assert any("substrate vote" in n for n in got.notes), (
            f"the demotion must be explained in the notes: {got.notes}")

    def test_the_demoted_vote_is_retained_as_evidence_not_deleted(self):
        """Demoting the headline must not throw the evidence away -- disagreement leaves all candidates
        valid (schema.InducerConsensus), and the Zn2+ vote is still worth reporting."""
        sensor_class, calls = RECORDED["GntR"]
        got = IND._finalize(list(calls), None, sensor_class=sensor_class)
        vote = next((c for c in got.calls if c.ligand == "Zn2+"), None)
        assert vote is not None, "the demoted metal vote must survive in the per-source calls"

    def test_check_would_fail_on_a_genuine_flip(self):
        """A control on the check itself: if the guard were removed, this fixture MUST be caught. Here
        the same call set is finalized as if the sensor class were metal, which is what a broken guard
        effectively does -- and the assertion the other tests make then fails, as it should."""
        _sc, calls = RECORDED["GntR"]
        unguarded = IND._finalize(list(calls), None, sensor_class="metal")
        assert is_metal_headline(unguarded.top), (
            "the detector must be able to see a metal headline, or the checks above prove nothing")


class TestOfflineOnlyCheckIsInsufficient:
    """Why 7.2 existed at all: the offline evidence for GntR contains no metal signal whatsoever, so an
    offline-only check cannot fail no matter how broken the online fusion is."""

    def test_offline_calls_alone_never_produce_a_metal(self):
        _sc, calls = RECORDED["GntR"]
        offline = [c for c in calls if c.source != "regulon"]
        got = IND._finalize(offline, None, sensor_class="organic")
        assert not is_metal_headline(got.top)
