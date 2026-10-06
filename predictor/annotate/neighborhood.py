"""neighborhood.py -- is there metal-homeostasis machinery next to this regulator?

A metalloregulator very often sits beside the transporter, efflux pump or reductase it controls, so
the identity of a TF's immediate genomic neighbours is evidence about what it senses. The survey
(Rondon et al.) used exactly this as its second axis: a candidate is "effective" when it carries both
a predicted metal site AND a metal-homeostasis gene within 500 bp.

This module ports that classification onto our own gene model. The category tables come from the
collaborator's `neighborhood_analysis.py` (vendored verbatim under `predictor/discovery/vendor/`);
the matching, the coordinate handling and the abstention policy are ours.

Two rules this module exists to respect:

1. **It reports a metal CLASS, never an ion.** A `copA` neighbour says copper is handled nearby -- not
   that this regulator senses copper. Plenty of metalloregulators sit beside a transporter for a
   different metal, and plenty of non-metal regulators sit beside one by accident. `ligand` is `None`
   by design, exactly as for `structure.metalnet`.
2. **Element symbols and the word "metal" are never matched as substrings** (handoff §5.4).
   `steroid-CoA` contains "co", `nitrogen` contains "ni", `non-metal effector` contains "metal", and
   short gene symbols like `cora`/`moda`/`dps`/`bfr` are the ones most likely to appear inside an
   unrelated word. Symbols match on word boundaries; phrases are multi-word and matched literally.

    python -m predictor.annotate.neighborhood --self-test
"""
from __future__ import annotations

import re
from dataclasses import dataclass

#: Gene symbols and product phrases that mark metal homeostasis, by category. Vendored from the
#: collaborator's `neighborhood_analysis.py` CATS table -- do not silently extend it: the survey's
#: published counts (13/150 with a metal neighbour <= 500 bp, 4 "effective") are reproducible only
#: against this exact vocabulary.
CATEGORIES: dict[str, dict[str, tuple]] = {
    "metal_efflux_transport": {
        "sym": (
            "copa", "copb", "copc", "copd", "copz", "znta", "zntb", "zupt", "znua", "znub", "znuc",
            "nika", "nikb", "nikc", "nikd", "nike", "feob", "feoa", "arsb", "acr3", "czca", "czcb",
            "czcc", "czcd", "cnra", "cnrb", "cnrc", "cora", "mgta", "mgtb", "mgte", "mnth", "sita",
            "sitb", "sitc", "sitd", "troa", "moda", "modb", "modc", "adca", "adcb", "adcc", "dmef",
            "yiip", "zosa", "ctpv", "ctpc", "ctpg", "ctpj", "pbra", "pbrt", "pbrb", "nrsd",
        ),
        "phr": (
            "p-type atpase", "cation diffusion facilitator", "cation efflux",
            "heavy metal translocating", "zinc transporter", "zinc-transporting", "zinc abc",
            "nickel transporter", "nickel abc", "nickel/cobalt", "cobalt transporter",
            "cobalt-zinc-cadmium", "manganese transporter", "manganese abc", "copper-transporting",
            "copper-translocating", "copper abc", "copper resistance", "copper efflux",
            "copper-exporting", "molybdate", "tungstate", "cadmium", "mercuric transport",
            "heavy metal", "metal abc transporter", "metal transporter", "divalent metal",
            "divalent cation", "cation-transporting", "metal cation", "arsenical",
            "arsenite efflux", "magnesium/cobalt", "lead-translocating",
        ),
    },
    "metal_resistance_enzyme": {
        "sym": ("arsc", "arsi", "arsm", "arsh", "mera", "merb", "merc", "merp", "cadi", "cueo", "czcm"),
        "phr": (
            "arsenate reductase", "arsenite methyltransferase", "mercuric reductase",
            "alkylmercury", "tellurite", "copper oxidase", "multicopper oxidase", "cupric reductase",
        ),
    },
    "iron_siderophore": {
        "sym": (
            "feca", "fecb", "fecc", "fecd", "fece", "fhua", "fhub", "fhuc", "fhud", "feua", "feub",
            "iuca", "iucb", "iucc", "iucd", "irga", "tonb", "exbb", "exbd", "vcta", "vctc", "vctd",
            "vctg", "vctp", "viua", "viub", "hmut", "hmuu", "hmuv", "entc", "ente", "entf", "fepa",
            "fepb", "fepc", "fepd", "fepg", "irtb",
        ),
        "phr": (
            "siderophore", "ferric", "ferrous", "iron uptake", "iron chelate", "iron abc",
            "iron transport", "iron-regulated", "iron-siderophore", "enterobactin", "vibrioferrin",
            "vibriobactin", "mycobactin", "hemin", "heme transport", "heme abc", "heme uptake",
            "tonb-dependent", "ferrichrome", "ferric uptake", "iron complex", "ferrioxamine",
        ),
    },
    "metal_storage_chelation": {
        "sym": ("dps", "ftna", "ftnb", "bfr", "bfd", "mymt"),
        "phr": (
            "ferritin", "bacterioferritin", "metallothionein", "metallochaperone",
            "copper chaperone", "iron storage",
        ),
    },
}

#: The categories that count as metal homeostasis. The remaining labels below are descriptive only.
METAL_CATEGORIES = frozenset(CATEGORIES)

_OTHER_REGULATOR = ("transcriptional regulator", "repressor", "activator", "dna-binding",
                    "helix-turn-helix")
_HYPOTHETICAL = ("hypothetical", "duf", "unknown", "uncharacterized")
_GENERIC_TRANSPORT = ("transporter", "permease", "efflux", "abc transporter")

#: Compiled once. Symbols are anchored on word boundaries precisely because several of them (`cora`,
#: `moda`, `dps`, `bfr`, `entc`) are short enough to occur inside unrelated words.
_SYM_RE = {cat: re.compile(r"\b(?:" + "|".join(re.escape(s) for s in kw["sym"]) + r")\b")
           for cat, kw in CATEGORIES.items()}

#: How many genes either side to consider. The collaborator's script used 6; kept, so the published
#: neighbour counts stay reproducible.
FLANK_GENES = 6

#: The survey's own thresholds. 500 bp is what its "effective candidate" definition uses.
NEAR_BP = 500
VERY_NEAR_BP = 200


@dataclass(frozen=True)
class Neighbour:
    gene: str
    product: str
    category: str
    gap: int                 # bp between the two gene bodies; 0 when they overlap
    side: str                # 'up' | 'down' relative to the regulator on the + strand
    protein_id: str = ""

    @property
    def is_metal(self) -> bool:
        return self.category in METAL_CATEGORIES


@dataclass
class NeighbourhoodCall:
    """What the regulator's genomic neighbourhood says, with the reason when it says nothing."""
    status: str                              # 'ok' | 'tf_not_found' | 'no_context'
    neighbours: tuple = ()                   # nearest-first
    nearest_metal: Neighbour | None = None
    metal_within_500: bool = False
    metal_within_200: bool = False
    categories: tuple = ()                   # distinct metal categories seen within NEAR_BP
    notes: tuple = ()

    @property
    def fired(self) -> bool:
        return self.metal_within_500


def categorize(product: str | None, gene: str | None) -> str:
    """Classify one neighbouring gene from its symbol and free-text product."""
    text = f"{gene or ''} {product or ''}".lower()
    for cat, kw in CATEGORIES.items():
        if _SYM_RE[cat].search(text):
            return cat
        if any(phrase in text for phrase in kw["phr"]):
            return cat
    if any(k in text for k in _OTHER_REGULATOR):
        return "other_regulator"
    if any(k in text for k in _HYPOTHETICAL):
        return "hypothetical_uncharacterized"
    if any(k in text for k in _GENERIC_TRANSPORT):
        return "generic_transport_nonmetal"
    return "other"


def _gap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """bp between two gene bodies; 0 if they overlap. Coordinates are 0-based half-open."""
    if a_end >= b_start and b_end >= a_start:
        return 0
    return (b_start - a_end) if b_start > a_end else (a_start - b_end)


def analyse(ctx, tf_gene, *, flank: int = FLANK_GENES) -> NeighbourhoodCall:
    """Classify the `flank` genes either side of `tf_gene` in `ctx`.

    `tf_gene` is a `Gene` from the same context (see `context.find_tf_gene`).
    """
    if ctx is None or not getattr(ctx, "genes", None):
        return NeighbourhoodCall(status="no_context",
                                 notes=("no genome context; neighbourhood not evaluated",))
    genes = sorted(ctx.genes, key=lambda g: g.start)
    try:
        i = next(k for k, g in enumerate(genes)
                 if g.start == tf_gene.start and g.end == tf_gene.end)
    except (StopIteration, AttributeError):
        return NeighbourhoodCall(status="tf_not_found",
                                 notes=("the TF gene is not in this context",))

    out = []
    for j in range(max(0, i - flank), min(len(genes), i + flank + 1)):
        if j == i:
            continue
        g = genes[j]
        out.append(Neighbour(
            gene=g.name or "", product=g.product or "",
            category=categorize(g.product, g.name),
            gap=_gap(tf_gene.start, tf_gene.end, g.start, g.end),
            side="up" if g.end <= tf_gene.start else "down",
            protein_id=getattr(g, "protein_id", "") or "",
        ))
    out.sort(key=lambda n: n.gap)

    metals = [n for n in out if n.is_metal]
    near = [n for n in metals if n.gap <= NEAR_BP]
    return NeighbourhoodCall(
        status="ok",
        neighbours=tuple(out),
        nearest_metal=metals[0] if metals else None,
        metal_within_500=bool(near),
        metal_within_200=any(n.gap <= VERY_NEAR_BP for n in metals),
        categories=tuple(dict.fromkeys(n.category for n in near)),
    )


def inducer_call(call: NeighbourhoodCall):
    """`InducerCall` for the fused inducer, or None when the neighbourhood says nothing.

    Confidence is deliberately modest and distance-graded, and `ligand` is ALWAYS None: co-location
    argues that the regulator works in metal homeostasis, and says nothing about which ion it binds.
    """
    from predictor.schema import InducerCall

    if not call.fired or call.nearest_metal is None:
        return None
    n = call.nearest_metal
    confidence = 0.45 if call.metal_within_200 else 0.30
    return InducerCall(
        source="neighborhood",
        ligand=None,
        confidence=confidence,
        role="metal",
        evidence={
            "nearest_metal_gene": n.gene or n.product[:40],
            "nearest_metal_protein_id": n.protein_id,
            "gap_bp": n.gap,
            "category": n.category,
            "categories_within_500bp": list(call.categories),
            "note": "co-location with metal-homeostasis machinery; names no ion by design",
        },
    )


def survey_effective(consensus) -> dict:
    """The survey's own "effective candidate" definition, recomputed from a fused inducer call.

    Rondon et al. call a candidate **effective** when it carries BOTH a predicted metal site AND a
    metal-homeostasis gene within 500 bp -- the intersection of their two axes. Published on the same
    150 candidates: 41 with a MetalNet site, 13 with a metal neighbour <= 500 bp, **4 effective** at
    500 bp and 3 at 200 bp.

    Both inputs already ride in the fused `InducerConsensus`: the site from the `metalnet` source, the
    proximity from `neighborhood`. Reading them back out here is what lets one run of our pipeline
    reproduce both halves of the published analysis instead of only ours.

    Note this is deliberately THEIR definition and not ours: it asks about a MetalNet site
    specifically, not about our coordination gate, because that is what their number counts.
    """
    calls = list(getattr(consensus, "calls", None) or [])
    mn = next((c for c in calls if c.source == "metalnet"), None)
    nb = next((c for c in calls if c.source == "neighborhood"), None)
    has_site = None if mn is None else (mn.role == "metal")
    ev = (nb.evidence or {}) if nb is not None else {}
    gap = ev.get("gap_bp")
    near_500 = nb is not None
    near_200 = bool(near_500 and isinstance(gap, int) and gap <= VERY_NEAR_BP)
    return {
        "has_metalnet_site": has_site,
        "metal_neighbour_500bp": near_500,
        "metal_neighbour_200bp": near_200,
        "nearest_metal_gene": ev.get("nearest_metal_gene"),
        "nearest_metal_gap_bp": gap,
        "effective_500bp": bool(has_site) and near_500,
        "effective_200bp": bool(has_site) and near_200,
        "definition": "Rondon et al.: a predicted metal SITE and a metal-homeostasis gene within 500 bp",
    }


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    from predictor.signals.motif_rescan import Gene, GenomeContext

    # handoff §5.4: the substring traps that produced wrong numbers here twice before.
    assert categorize("steroid-CoA dehydrogenase", "scdA") != "metal_efflux_transport", \
        "'steroid-CoA' contains 'co' but is not a cobalt transporter"
    assert categorize("nitrogen regulatory protein", "glnB") != "metal_efflux_transport", \
        "'nitrogen' contains 'ni' but is not a nickel transporter"
    assert categorize("non-metal effector binding protein", "yfoo") != "metal_efflux_transport", \
        "'non-metal effector' contains 'metal' but is not metal machinery"
    assert categorize("decorated protein", "decorA") != "metal_efflux_transport", \
        "'decorA' contains 'cora' but is not a magnesium transporter"
    assert categorize("RNA polymerase sigma factor SigB", "sigB") not in METAL_CATEGORIES, \
        "a sigma factor is not metal machinery"

    # ... and the positives it must still catch
    assert categorize("copper-translocating P-type ATPase", "copA") == "metal_efflux_transport"
    assert categorize("mercuric reductase", "merA") == "metal_resistance_enzyme"
    assert categorize("TonB-dependent siderophore receptor", "fhuA") == "iron_siderophore"
    assert categorize("bacterioferritin", "bfr") == "metal_storage_chelation"
    assert categorize("transcriptional regulator, MarR family", "marR") == "other_regulator"

    ctx = GenomeContext(accession="TEST", sequence="A" * 6000, genes=[
        Gene(0, 900, "+", "copA", "WP_1", "copper-translocating P-type ATPase"),
        Gene(1000, 1400, "+", "tfR", "WP_2", "transcriptional regulator"),
        Gene(4000, 4600, "+", "gyrA", "WP_3", "DNA gyrase subunit A"),
    ])
    tf = ctx.genes[1]
    call = analyse(ctx, tf)
    assert call.status == "ok", call.status
    assert call.metal_within_500 and call.metal_within_200, "copA is 100 bp away"
    assert call.nearest_metal.gene == "copA"
    assert call.nearest_metal.gap == 100, call.nearest_metal.gap

    ic = inducer_call(call)
    assert ic is not None and ic.source == "neighborhood" and ic.role == "metal"
    assert ic.ligand is None, "the neighbourhood source must never name an ion"

    # a regulator with no metal machinery nearby must abstain, not guess
    lonely = GenomeContext(accession="TEST", sequence="A" * 6000, genes=[
        Gene(1000, 1400, "+", "tfR", "WP_2", "transcriptional regulator"),
        Gene(1500, 2000, "+", "gyrA", "WP_3", "DNA gyrase subunit A"),
    ])
    quiet = analyse(lonely, lonely.genes[0])
    assert quiet.status == "ok" and not quiet.fired
    assert inducer_call(quiet) is None, "no metal neighbour -> no call at all"

    # a distant metal gene is recorded but does not fire the 500 bp rule
    far = GenomeContext(accession="TEST", sequence="A" * 20000, genes=[
        Gene(1000, 1400, "+", "tfR", "WP_2", "transcriptional regulator"),
        Gene(9000, 9600, "+", "copA", "WP_1", "copper-translocating P-type ATPase"),
    ])
    fcall = analyse(far, far.genes[0])
    assert fcall.nearest_metal is not None and not fcall.metal_within_500
    assert inducer_call(fcall) is None

    # the survey's "effective" definition needs BOTH axes; neither alone is enough
    class _C:
        def __init__(self, calls):
            self.calls = calls

    class _Call:
        def __init__(self, source, role, evidence=None):
            self.source, self.role, self.evidence = source, role, (evidence or {})

    site = _Call("metalnet", "metal")
    no_site = _Call("metalnet", "non-metal")
    near = _Call("neighborhood", "metal", {"gap_bp": 100, "nearest_metal_gene": "copA"})
    far_ish = _Call("neighborhood", "metal", {"gap_bp": 400, "nearest_metal_gene": "copA"})

    assert survey_effective(_C([site, near]))["effective_500bp"] is True
    assert survey_effective(_C([site, near]))["effective_200bp"] is True
    assert survey_effective(_C([site, far_ish]))["effective_500bp"] is True
    assert survey_effective(_C([site, far_ish]))["effective_200bp"] is False, "400 bp is not <= 200"
    assert survey_effective(_C([no_site, near]))["effective_500bp"] is False, "a neighbour alone is not effective"
    assert survey_effective(_C([site]))["effective_500bp"] is False, "a site alone is not effective"
    assert survey_effective(_C([]))["has_metalnet_site"] is None, "no MetalNet call -> unknown, not False"

    print("OK: categories, substring traps, gap arithmetic, abstention, the no-ion rule and the "
          "survey's effective definition all hold.")


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        _demo()
    else:
        print(__doc__)
        sys.exit(2)
