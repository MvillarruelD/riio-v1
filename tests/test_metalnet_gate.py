"""MetalNet2 as the fifth inducer source -- precedence, abstention safety, and the ion-naming ban.

Context. `metal_site.coordination_gate` implements metal-site chemistry for six of the twelve families
(MerR, ArsR/SmtB, NikR, CsoR/FrmR, Rrf2, Fur). For the other six it returns `has_site=None` on every
member: GntR (32 candidates), TetR (24), MarR (17), DtxR (4), CopY (3), LysR (2) -- 82 of the 150 run-5
candidates, 55 % of the run, with no structural opinion at all. MetalNet2 is family-agnostic and closes
that gap.

The three properties these tests pin, because each one is a way the change could do damage:

  1. **Scope.** MetalNet is consulted ONLY where our own chemistry is silent. A family-specific gate that
     ran -- positive or negative -- keeps the verdict, because those verdicts encode chemistry MetalNet
     cannot see (the non-metal MerR clades; the CsoR-family FrmR branch). A negative MetalNet also never
     demotes a metalloregulator family whose site is merely dispersed (MntR).
  2. **Abstention is not denial.** No engine, no MSA, or too shallow an MSA all mean "no opinion". Every
     verdict must then be bit-identical to what it was before this source existed -- which is also what
     makes the source safe to ship before the engine is installed everywhere.
  3. **MetalNet never names the ion.** Its second model does emit a metal type, but it is trained on
     structural sites of every kind and a Ca/Mg/structural-Zn site is not a regulator's sensed effector.
     The type is evidence only; the headline falls back to the CLASS.

Offline and deterministic: the engine is never invoked -- its verdict is injected as the dict the adapter
returns.
"""
import pytest

from predictor.effector import inducer as IND
from predictor.schema import InducerCall
from predictor.structure import metal_site as MS
from predictor.structure import metalnet as MN

# --- what structure.metalnet.predict_site returns, in its three states ------------------------------
SITE = {"available": True, "has_site": True, "n_site_residues": 4, "n_site_pairs": 5,
        "site_residues": ["C12", "C15", "H40", "C44"], "max_prob": 0.97,
        "metal_type_top": "ZN", "metal_type_votes": {"ZN": 5}, "msa_source": "colabfold"}
NO_SITE = {"available": True, "has_site": False, "n_site_residues": 0, "n_coevo_pairs": 12,
           "max_prob": 0.21, "msa_source": "colabfold"}
ABSTAIN = {"available": False, "has_site": None, "reason": "engine unavailable"}
SHALLOW = {"available": True, "has_site": None, "reason": "MSA too shallow (7 rows < 16)"}

# --- sequences ---------------------------------------------------------------------------------------
GNTR = "M" + "ARQLGID" * 12                                    # organic family; our gate has no branch
MNTR = "M" + "A" * 120                                          # DtxR/MntR: dispersed site, gate silent
MERR_METAL = "M" + "A" * 70 + "C" + "X" * 7 + "C" + "A" * 10    # C-terminal CX7C -> gate POSITIVE
MERR_NONMETAL = "M" + "A" * 60 + "L" * 60                       # no cysteines -> gate NEGATIVE
FRMR = "MPSTPEEKKKVLTRVRRIRGQIDALERSLEGDAECRAILQQIAAVRGAANGLMAEVLESHIRETFDRNDCYSREVSQSVDDTIELVRAYLK"
ISCR_LIKE = "M" + "A" * 88 + "CAAAAACAAAAAC" + "A" * 60


def _mn(role, **ev):
    return InducerCall(source="metalnet", ligand=None, confidence=0.0, role=role,
                       evidence={"metal_type_is_advisory": True, **ev})


class TestVerdictNormalisation:
    @pytest.mark.parametrize("payload,want", [(SITE, True), (NO_SITE, False), (ABSTAIN, None),
                                              (SHALLOW, None), (None, None), (True, True), (False, False)])
    def test_verdict(self, payload, want):
        assert MS.metalnet_verdict(payload) is want

    def test_unavailable_engine_is_never_a_negative(self):
        """The distinction the whole design rests on: `has_site=False` is a claim, `None` is silence."""
        assert MS.metalnet_verdict({"available": False, "has_site": False}) is None


class TestCoverageGap:
    """The six families the coordination gate is silent on -- what this change is FOR."""

    @pytest.mark.parametrize("family", ["GntR", "TetR/AcrR", "MarR", "DtxR/MntR", "CopY", "LysR-type (LTTR)"])
    def test_our_gate_is_silent(self, family):
        assert MS.coordination_gate(GNTR, family)["has_site"] is None

    def test_metalnet_site_establishes_the_metal_class(self):
        """Demonstrated on TetR, not GntR: GntR is now a structural-site family (see TestStructuralSite)
        and is deliberately exempt, so it can no longer stand for the general rule."""
        assert MS.inducer_class(GNTR, "TetR/AcrR")[0] == "organic"           # before
        assert MS.inducer_class(GNTR, "TetR/AcrR", None, SITE)[0] == "metal"  # after

    def test_the_reason_names_the_source(self):
        assert "MetalNet" in MS.inducer_class(GNTR, "TetR/AcrR", None, SITE)[1]

    def test_negative_metalnet_leaves_an_organic_family_organic(self):
        cls, why = MS.inducer_class(GNTR, "GntR", None, NO_SITE)
        assert cls == "organic"
        # HANDOFF 7.1: the organic-family guard's trigger becomes "neither method fired"
        assert "neither" in why

    @pytest.mark.parametrize("payload", [ABSTAIN, SHALLOW, None])
    def test_abstention_changes_nothing(self, payload):
        assert MS.inducer_class(GNTR, "GntR", None, payload) == MS.inducer_class(GNTR, "GntR")


class TestScope:
    """MetalNet is consulted only where our own chemistry has no opinion."""

    def test_positive_family_gate_outranks_a_negative_metalnet(self):
        assert MS.inducer_class(MERR_METAL, "MerR", None, NO_SITE)[0] == "metal"

    def test_negative_family_gate_outranks_a_positive_metalnet(self):
        assert MS.inducer_class(MERR_NONMETAL, "MerR", None, SITE)[0] == "organic"

    def test_frmr_stays_non_metal(self):
        """FrmR lost the CsoR His/Cys ligand set and senses formaldehyde (Osman 2016). MetalNet sees its
        remaining cysteines; the family-scoped gate knows what they are not."""
        assert MS.inducer_class(FRMR, "CsoR/FrmR", None, SITE)[0] == "organic"

    def test_redox_family_stays_redox(self):
        """An [Fe-S] cluster is not a sensed metal ion, whatever a metal-site detector says."""
        assert MS.inducer_class(ISCR_LIKE, "Rrf2", None, SITE)[0] == "redox"

    def test_negative_metalnet_does_not_demote_a_metalloregulator_family(self):
        """MntR's dispersed His/Asp carboxylate site is exactly what BOTH methods are weakest at, so a
        miss is weak evidence. The dissent is recorded, not acted on."""
        cls, why = MS.inducer_class(MNTR, "DtxR/MntR", None, NO_SITE)
        assert cls == "metal" and "dissent" in why


class TestLigifyWeighting:
    def test_confirmed_site_suppresses_operon_chemistry(self):
        assert MS.ligify_weight(GNTR, "TetR/AcrR")[0] == 1.00
        assert MS.ligify_weight(GNTR, "TetR/AcrR", metalnet=SITE)[0] == 0.10

    def test_metalnet_upgrades_advisory_to_suppressed(self):
        """A silent-gate metalloregulator was only 'advisory' (0.30) because nothing confirmed the site.
        MetalNet confirming it is direct evidence, so Ligify is suppressed as hard as for our own gate."""
        assert MS.ligify_weight(MNTR, "DtxR/MntR")[0] == 0.30
        assert MS.ligify_weight(MNTR, "DtxR/MntR", metalnet=SITE)[0] == 0.10


class TestFusionGuards:
    """The E. coli GntR regression, before and after -- HANDOFF 7.1's 'neither method fired'."""

    BASE = [InducerCall(source="coordination", ligand=None, confidence=0.0, role="", evidence={"gate": ""}),
            InducerCall(source="ligify", ligand="D-threo-isocitrate", confidence=0.0, role="",
                        evidence={"weight": 1.0}),
            InducerCall(source="ligify", ligand="Zn2+", confidence=0.393, role="metal",
                        evidence={"weight": 1.0, "votes": {"Zn2+": 1.0}})]

    def test_guard_still_demotes_a_lone_substrate_vote(self):
        got = IND._finalize(list(self.BASE), None, sensor_class="organic")
        assert got.top == "non-metal effector (undetermined)"

    def test_a_metalnet_site_lets_the_metal_call_stand(self):
        got = IND._finalize(list(self.BASE) + [_mn("metal", has_site=True, n_site_residues=4)],
                            None, sensor_class="metal")
        assert got.top == "Zn2+"
        assert any("MetalNet finds a metal site" in n for n in got.notes)

    def test_abstaining_engine_reproduces_the_old_verdict(self):
        before = IND._finalize(list(self.BASE), None, sensor_class="organic")
        after = IND._finalize(list(self.BASE) + [_mn("", has_site=None, reason="engine unavailable")],
                              None, sensor_class="organic")
        assert after.top == before.top

    def test_metalnet_never_names_the_ion(self):
        got = IND._finalize([InducerCall(source="coordination", ligand=None, confidence=0.0, role="",
                                         evidence={"gate": ""}),
                             _mn("metal", has_site=True, n_site_residues=4, metal_type_top="ZN")],
                            None, sensor_class="metal")
        assert got.top == "divalent metal (ion unresolved)"
        assert all(c.ligand is None for c in got.calls if c.source == "metalnet")

    def test_metal_type_is_flagged_advisory_in_the_notes(self):
        got = IND._finalize([InducerCall(source="coordination", ligand=None, confidence=0.0, role="",
                                         evidence={"gate": ""}),
                             _mn("metal", has_site=True, n_site_residues=4, metal_type_top="ZN")],
                            None, sensor_class="metal")
        assert any("ADVISORY" in n for n in got.notes)

    def test_family_gate_wins_a_conflict_and_says_so(self):
        got = IND._finalize([InducerCall(source="ssn_cluster", ligand="formaldehyde", confidence=1.0,
                                         role="non-metal", evidence={"cluster": "CsoR_ecFrmR"}),
                             InducerCall(source="coordination", ligand=None, confidence=0.0,
                                         role="non-metal", evidence={"gate": "CsoR-His/Cys"}),
                             _mn("metal", has_site=True, n_site_residues=3)],
                            False, gate_name="CsoR-His/Cys", sensor_class="organic")
        assert got.top == "formaldehyde"
        assert any("dissent recorded, not acted on" in n for n in got.notes)

    def test_sensor_class_is_recovered_on_a_refuse(self):
        """`augment_with_regulon` re-finalizes without a sequence, so it cannot recompute the class. A
        metal class established only by MetalNet must survive that, or the organic-family guard fires on
        a TF that does have a site."""
        got = IND._finalize(list(self.BASE) + [_mn("metal", has_site=True, n_site_residues=4)],
                            None, sensor_class=None)
        assert got.top == "Zn2+"


class TestAdapterDegradesGracefully:
    """No engine must mean no opinion -- never an exception and never a fabricated verdict."""

    def test_empty_sequence(self):
        got = MN.predict_site("")
        assert got["available"] is False and got["has_site"] is None

    def test_unavailable_shape(self):
        got = MN._unavailable("no engine")
        assert got["has_site"] is None and got["available"] is False and got["engine"] == "MetalNet2"

    def test_status_never_raises(self):
        st = MN.status()
        assert "ready" in st and isinstance(st["ready"], bool)

    def test_summary_handles_every_state(self):
        assert "unavailable" in MN.summary(ABSTAIN)
        assert "no opinion" in MN.summary(SHALLOW)
        assert "metal site" in MN.summary({**SITE, "n_site_pairs": 5})
        assert "no metal site" in MN.summary({**NO_SITE, "n_coevo_pairs": 12})

    def test_cache_key_is_sequence_scoped(self):
        assert MN._key("MKV") == MN._key("mkv") != MN._key("MKVA")

    def test_nul_bytes_are_stripped_from_an_a3m(self):
        """ColabFold emits NUL bytes; left in, they shift every column of a positionally-read alignment."""
        assert MN._clean_a3m(">q\nMK\x00V\n") == ">q\nMKV\n"


class TestStructuralSite:
    """A metal site in GntR is reported, but is not read as a metal-sensor call.

    Chem. Rev. metallostasis review 3.2.10: ~70 % of GntR carry a conserved Zn site buried in the
    effector cavity that may exist to position a carboxylate-bearing organic acid, and whether it is
    regulatory at all is stated there as an open question. GntR is also our largest blind block (32 of
    150 candidates), and MetalNet finds a site in 9 of them -- so treating "site" as "sensor" here
    would convert a fifth of the family on evidence the literature will not carry.
    """

    def test_gntr_does_not_become_a_metal_sensor_on_a_metalnet_site(self):
        assert MS.inducer_class(GNTR, "GntR", None, SITE)[0] == "organic"

    def test_but_the_site_is_STATED_not_silently_dropped(self):
        why = MS.inducer_class(GNTR, "GntR", None, SITE)[1]
        assert "MetalNet DOES find a metal site" in why
        assert "structural" in why, "the reader must be told why it is not a sensor call"

    def test_the_site_does_not_suppress_the_organic_evidence(self):
        """Measured on run 6: all 9 GntR candidates with a MetalNet site had Ligify down-weighted x0.1
        by a site that may exist to bind the very organic acid Ligify was proposing."""
        assert MS.ligify_weight(GNTR, "GntR", metalnet=SITE)[0] == 1.00

    def test_the_caveat_is_scoped_to_the_families_that_earned_it(self):
        assert MS.structural_site_caveat("GntR")
        for fam in ("TetR/AcrR", "MarR/SlyA", "MerR", "Fur", "LysR-type (LTTR)"):
            assert MS.structural_site_caveat(fam) is None, f"{fam} must keep the normal rule"

    def test_a_positive_family_gate_still_wins_everywhere(self):
        """The caveat only governs a SILENT gate. If our own chemistry fired, it decides."""
        assert MS.inducer_class(MERR_METAL, "MerR", None, SITE)[0] == "metal"


class TestCofactorIsNotAnInducer:
    """A metal the protein CONTAINS is not the metal it SENSES.

    Curated-binder ions are mined from UniProt metal-BINDING features, which do not distinguish the
    two. Every SoxR ortholog therefore arrives carrying `ion="Fe"` from its [2Fe-2S] cluster -- 7 such
    rows in the MerR set. Leave-one-sequence-out showed the identity route then names Fe2+ for SoxR at
    0.96 identity: a false metal for a superoxide sensor
    (`analysis/regulondb_bench/curated_binder_loo.py`). Production masked it only because the redox
    short-circuit happens to fire first, which is luck rather than a safeguard.
    """

    def test_a_cofactor_cluster_binder_names_no_inducer(self):
        from predictor.effector import inducer as IND
        soxr = next(b for b in IND._curated_binders("MerR") if b.get("acc") == "P0ACS2")
        assert not soxr.get("ion"), "an [2Fe-2S] cofactor must not be offered as a sensed ion"
        assert soxr.get("ion_not_inducer") == "Fe", "but the mined value is kept, not silently dropped"
        assert "cofactor" in (soxr.get("ion_source") or "")

    def test_genuine_metal_binders_are_untouched(self):
        from predictor.effector import inducer as IND
        cuer = next((b for b in IND._curated_binders("MerR") if b.get("acc") == "P0A9G4"), None)
        if cuer is not None:
            assert cuer.get("ion"), "CueR's Cu(I) is a sensed inducer and must survive the filter"
        assert sum(1 for b in IND._curated_binders("MerR") if b.get("ion")) >= 20


class TestCandidateSetIsCanonical:
    """One species must not appear twice in a candidate set because two sources spelled it differently.

    Found in the v7 run: E. coli ArsR came out ('As3+', 'antimonite', 'As(III)'), which reads as a
    three-way ambiguity when it is a two-way one -- the SSN cluster table writes arsenite `As3+` and
    the anchor KB writes `As(III)`.
    """

    @staticmethod
    def _call(source, ligand, conf=0.9, **kw):
        from predictor.schema import InducerCall
        return InducerCall(source=source, ligand=ligand, confidence=conf, role="metal", **kw)

    def test_two_spellings_of_one_species_collapse(self):
        from predictor.schema import InducerConsensus
        cons = InducerConsensus.from_calls([
            self._call("ssn_cluster", "As3+", 1.0),
            self._call("regulon", "As(III)", 0.6),
            self._call("ligify", "antimonite", 0.5),
        ])
        assert cons.top == "As3+", "the headline is still the highest-confidence RAW string"
        assert list(cons.candidates) == ["As3+", "antimonite"], \
            "As3+ and As(III) are one species; antimonite is genuinely a second"

    def test_chebi_charge_and_compartment_notation_is_one_species(self):
        """`Ni(2+)(in)` reaches us from ChEBI through Ligify and is plain Ni2+.

        Found in the v7 panel report: NikR's candidate set read ['Ni2+', 'Ni(2+)(in)'] -- one ion
        displayed as an ambiguity. `normalise` stripped the trailing `(in)` gloss, but its
        oxidation-state fallback accepted only `Ni(III)` or `Ni2+`, never a parenthesised charge.
        """
        from predictor.effector import inducer_vocab as vocab
        from predictor.schema import InducerConsensus

        assert vocab.normalise("Ni(2+)(in)") == vocab.normalise("Ni2+")
        assert vocab.normalise("Ni(2+)") == vocab.normalise("Ni(II)")
        cons = InducerConsensus.from_calls([
            self._call("ssn_cluster", "Ni2+", 1.0),
            self._call("ligify", "Ni(2+)(in)", 0.5),
        ])
        assert cons.top == "Ni2+"
        # An unambiguous call carries NO candidate set -- so the two sources agreeing must leave it
        # empty, not populate a one-entry "tie". Before the fix this was a two-entry ambiguity.
        assert list(cons.candidates) == []

    def test_the_pipelines_own_class_level_calls_classify_correctly(self):
        """`_finalize` emits these when the CLASS is settled but the ion is not.

        They name a class rather than a species, so a token lookup misses them -- and the fallthrough
        calls every unrecognised token "organic". That turned "divalent metal (ion unresolved)", the
        label the gate produces when it is working, into an ORGANIC call: the answer inverted on
        precisely the cases the gate exists to handle.
        """
        from predictor.effector import inducer_vocab as vocab
        assert vocab.inducer_class_of("divalent metal (ion unresolved)") == "metal"
        assert vocab.is_metal("divalent metal (ion unresolved)")
        assert vocab.inducer_class_of("non-metal effector (undetermined)") == "organic"
        assert vocab.inducer_class_of("reactive oxygen/nitrogen species (redox)") == "redox"
        # and a genuine cross-class mixture must still refuse to pick a side
        assert vocab.inducer_class_of("Cu+/persulfide") is None

    @pytest.mark.parametrize("label", [
        "MoO42-", "MoO4-", "molybdate", "WO4", "tungstate", "CrO4", "chromate",
        "AsO3", "arsenite", "AsO4", "arsenate", "SeO3", "selenite", "SeO4", "selenate",
        "VO3-", "metavanadate", "VO4", "orthovanadate", "SbO4", "antimonate", "TeO3", "tellurite",
    ])
    def test_metal_and_metalloid_oxyanions_are_metal_class(self, label):
        from predictor.effector import inducer_vocab as vocab
        assert vocab.inducer_class_of(label) == "metal"

    def test_every_metal_substrate_rule_emits_a_metal_class_token(self):
        from predictor.effector import inducer_vocab as vocab
        from predictor.effector import substrate_map
        missed = [r.effector for r in substrate_map.RULES
                  if r.role == "metal" and vocab.inducer_class_of(r.effector) != "metal"]
        assert missed == []

    def test_distinct_oxyanions_are_not_collapsed_to_one_species(self):
        from predictor.effector import inducer_vocab as vocab
        assert vocab.normalise("VO3-") != vocab.normalise("VO43-")
        assert vocab.normalise("arsenite") != vocab.normalise("arsenate")
        assert vocab.normalise("selenite") != vocab.normalise("selenate")

    def test_display_label_exposes_ambiguity_without_reclassifying_it(self):
        from predictor.effector import inducer_vocab as vocab
        top = "divalent metal (ion unresolved)"
        shown = vocab.display_label(top, ("Co2+", "formaldehyde"))
        assert shown.startswith("inconclusive: Co2+ or formaldehyde")
        assert vocab.inducer_class_of(top) == "metal"
        # The display prose is deliberately not a classification token; only `top` carries class.
        assert shown != top

    def test_display_label_only_rewrites_a_placeholder_with_two_named_candidates(self):
        from predictor.effector import inducer_vocab as vocab
        assert vocab.display_label("Ni2+/Co2+", ("Ni2+", "Co2+")) == "Ni2+/Co2+"
        assert vocab.display_label("divalent metal (ion unresolved)", ("Co2+",)) == \
            "divalent metal (ion unresolved)"
        assert vocab.display_label("divalent metal (ion unresolved)", "Co2+") == \
            "divalent metal (ion unresolved)"
        # Duplicate spellings and class-level placeholders do not manufacture a shortlist.
        assert vocab.display_label("divalent metal (ion unresolved)",
                                   ("Cobalt", "Co2+", "non-metal effector (undetermined)")) == \
            "divalent metal (ion unresolved)"

    def test_the_5_4_substring_traps_still_do_not_normalise_to_an_ion(self):
        """The widened charge regex must stay anchored on the whole token (HANDOFF 5.4)."""
        from predictor.effector import inducer_vocab as vocab
        for trap in ("steroid-CoA", "non-metal effector", "2,3-dihydroxybenzoate"):
            assert vocab.normalise(trap) == trap, f"{trap} must pass through untouched"

    def test_distinct_species_are_not_merged(self):
        from predictor.schema import InducerConsensus
        cons = InducerConsensus.from_calls([
            self._call("ssn_cluster", "Fe2+", 1.0),
            self._call("regulon", "Zn2+", 0.6),
        ])
        assert list(cons.candidates) == ["Fe2+", "Zn2+"]

    def test_raw_spelling_is_preserved_not_canonicalised(self):
        """The report prints these beside `top`, which is raw -- only the dedup KEY is canonical."""
        from predictor.schema import InducerConsensus
        cons = InducerConsensus.from_calls([
            self._call("ssn_cluster", "Cobalt", 1.0),
            self._call("regulon", "Fe2+", 0.6),
        ])
        assert cons.candidates[0] == "Cobalt", "not rewritten to Co2+"
