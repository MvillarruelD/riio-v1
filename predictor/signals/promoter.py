"""
promoter.py -- promoter positioning (-35/-10/TSS) + mode-of-regulation classifier (Phase 1).

Runs the LaFleur PromoterCalculator on an operator's region to place the sigma70 promoter, then
classifies the **regulation mode** from the operator's geometry relative to that promoter. The label
states WHERE the operator sits, and nothing more:

  elongated-spacer-overlap  -- operator inside an over-long (>=18 bp) -35/-10 spacer. The described
                               mechanism for this geometry is DNA distortion: the TF stays bound and
                               bends the DNA to re-phase the boxes.
  promoter-core-overlap     -- operator covers the -35/-10 box or the promoter core, where a bound TF
                               would occlude RNAP.
  no-promoter-overlap       -- operator does not overlap the predicted promoter.

THE SAME GEOMETRY AND THE SAME LABELS ARE APPLIED TO EVERY FAMILY. An earlier revision emitted mechanism
archetypes ('MerR-distortion-activation', 'ArsR-occlusion-repression') for MerR and ArsR/SmtB and neutral
tags for everyone else, and let those two families fall back to a family prior when the geometry showed no
overlap at all. That asserted a mechanism from a family name rather than from evidence, and it made the
output of two families non-comparable with the other ten. Both are gone: the mechanism is now stated in
the docs above as an interpretation of the geometry, and the call itself reports only what was measured.

The decision logic (`_decide_mode`) is pure and unit-tested independently of the calculator's promoter
choices.

Run `python promoter.py` for a self-test.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field

_CANONICAL_SPACER = 17
_LONG_SPACER = 18           # >= this in an operator-spanning spacer => MerR-style


@dataclass
class ModeCall:
    mode: str
    spacer_len: int | None
    promoter: dict = field(default_factory=dict)
    confidence: float = 0.0
    evidence: str = ""
    flags: list = field(default_factory=list)
    # --- footprint: which promoter elements the WHOLE operator span physically covers (not just its
    # centre), + the operator's edges relative to the TSS. `overlaps` is bp of overlap per element;
    # `occludes` is the ordered list of elements the operator touches (['-35','-10','spacer','TSS']).
    overlaps: dict = field(default_factory=dict)
    occludes: list = field(default_factory=list)
    footprint: dict = field(default_factory=dict)


def _ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def _promoter_calculator_backend():
    """Import PromoterCalculator and bridge its Python 3.14 dataclass lookup bug.

    ``promotercalculator==1.2.4`` implements mapping access on ``PromoCalcResults`` through
    ``self.__annotations__``. Python 3.14 no longer exposes that special attribute through an instance,
    although the class annotation map is intact. The backend itself uses mapping access while choosing
    its minimum-energy promoter, so every call otherwise fails before returning a prediction.
    """
    import promoter_calculator as backend

    if sys.version_info >= (3, 14):
        from promoter_calculator.promoter_calculator import PromoCalcResults

        if not getattr(PromoCalcResults, "_tfop_py314_compat", False):
            def _getitem(self, key):
                if key in type(self).__annotations__:
                    return getattr(self, key)
                raise KeyError(key)

            def _keys(self):
                return type(self).__annotations__.keys()

            PromoCalcResults.__getitem__ = _getitem
            PromoCalcResults.keys = _keys
            PromoCalcResults._tfop_py314_compat = True
    return backend


def find_promoters(region_seq):
    """Return LaFleur promoters as dicts (sorted by Tx_rate desc)."""
    pc = _promoter_calculator_backend()
    # promoter_calculator.revcomp has no entry for N or other ambiguous bases;
    # replace anything that isn't ACGT with 'A' (safe for N-padded multi-contig genomes)
    clean = "".join(b if b in "ACGTacgt" else "A" for b in region_seq)
    out = []
    for p in pc.promoter_calculator(clean):
        h35 = sorted(p.hex35_position)
        h10 = sorted(p.hex10_position)
        sp = sorted(p.spacer_position)
        out.append(dict(tss=p.TSS, strand=p.strand, tx_rate=float(p.Tx_rate),
                        hex35=tuple(h35), hex10=tuple(h10), spacer=tuple(sp),
                        spacer_len=len(p.spacer), hex35_seq=p.hex35, hex10_seq=p.hex10))
    out.sort(key=lambda d: -d["tx_rate"])
    return out


def _decide_mode(op0, op1, prom):
    """Pure geometry decision, identical for every family. `prom` is a find_promoters() dict. Returns
    (mode, conf, evidence, overlaps) where `overlaps` is the bp of overlap between the operator span and
    each promoter element. The mode names the measured geometry; it never asserts a mechanism, and it never
    consults the TF's family -- two TFs whose operators sit in the same place get the same call."""
    h35, h10, sp = prom["hex35"], prom["hex10"], prom["spacer"]
    slen = prom["spacer_len"]
    in_spacer = _ov(op0, op1, sp[0], sp[1])
    in_35 = _ov(op0, op1, h35[0], h35[1])
    in_10 = _ov(op0, op1, h10[0], h10[1])
    core_lo, core_hi = min(h35[0], h10[0], sp[0]), max(h35[1], h10[1], sp[1])
    in_core = _ov(op0, op1, core_lo, core_hi)
    ov = {"-35": in_35, "-10": in_10, "spacer": in_spacer, "core": in_core}

    if in_spacer and slen >= _LONG_SPACER:
        conf = 0.6 + 0.1 * min(3, slen - _CANONICAL_SPACER)
        return ("elongated-spacer-overlap", min(0.95, conf),
                f"operator in a {slen}-bp spacer (>{_CANONICAL_SPACER} canonical)", ov)
    if in_35 or in_10:
        return ("promoter-core-overlap", 0.7,
                f"operator overlaps the {'-35' if in_35 else '-10'} box", ov)
    if in_core:
        if slen >= _LONG_SPACER:
            return ("elongated-spacer-overlap", 0.45, f"operator in promoter core, spacer {slen} bp", ov)
        return ("promoter-core-overlap", 0.45, "operator overlaps promoter core", ov)
    return ("no-promoter-overlap", 0.1, "operator does not overlap the predicted promoter", ov)


def _footprint(op0, op1, prom) -> tuple[list, dict]:
    """Which promoter elements the operator span covers + its edges relative to the TSS. Returns
    (occludes, footprint) where occludes is the ordered element list and footprint carries the TSS-relative
    edges (op_start_to_tss / op_end_to_tss) so a viewer can draw the full footprint, not just the centre."""
    h35, h10, sp = prom["hex35"], prom["hex10"], prom["spacer"]
    tss = prom.get("tss")
    occludes = []
    if _ov(op0, op1, h35[0], h35[1]):
        occludes.append("-35")
    if _ov(op0, op1, sp[0], sp[1]):
        occludes.append("spacer")
    if _ov(op0, op1, h10[0], h10[1]):
        occludes.append("-10")
    if tss is not None and op0 <= tss < op1:
        occludes.append("TSS")
    fp = {"op": [int(op0), int(op1)], "tss": (int(tss) if tss is not None else None),
          "op_start_to_tss": (int(op0 - tss) if tss is not None else None),
          "op_end_to_tss": (int(op1 - tss) if tss is not None else None),
          "op_center_to_tss": (int((op0 + op1) // 2 - tss) if tss is not None else None)}
    if tss is not None:                              # promoter element boxes relative to the TSS (for a
        fp["elements"] = {                           # zoom schematic: -35 / spacer / -10 vs the operator)
            "-35": [int(h35[0] - tss), int(h35[1] - tss)],
            "spacer": [int(sp[0] - tss), int(sp[1] - tss)],
            "-10": [int(h10[0] - tss), int(h10[1] - tss)]}
        fp["spacer_len"] = prom.get("spacer_len")
    return occludes, fp


def classify_mode(region_seq, op_start, op_end):
    """Place the promoter near the operator and classify the regulation mode from its geometry."""
    proms = find_promoters(region_seq)
    if not proms:
        return ModeCall("ambiguous", None, {}, 0.1, "no sigma70 promoter predicted in region",
                        ["no promoter"])
    opc = (op_start + op_end) / 2.0
    # choose the promoter whose core best overlaps (else is closest to) the operator
    def keyf(p):
        lo, hi = min(p["hex35"][0], p["spacer"][0]), max(p["hex10"][1], p["spacer"][1])
        ov = _ov(op_start, op_end, lo, hi)
        dist = abs(opc - (lo + hi) / 2.0)
        return (ov, -dist, p["tx_rate"])
    prom = max(proms, key=keyf)
    mode, conf, ev, overlaps = _decide_mode(op_start, op_end, prom)
    occludes, footprint = _footprint(op_start, op_end, prom)
    flags = [] if conf >= 0.4 else ["low-confidence mode (weak operator-promoter overlap)"]
    return ModeCall(mode, prom["spacer_len"], prom, round(conf, 2), ev, flags,
                    overlaps=overlaps, occludes=occludes, footprint=footprint)


# --------------------------------------------------------------------------- self-test
def _demo():
    # 1) pure geometry logic (independent of the calculator)
    prom_long = dict(hex35=(10, 16), hex10=(35, 41), spacer=(16, 35), spacer_len=19, tss=44)
    m, c, e, ov = _decide_mode(20, 31, prom_long)          # operator inside the 19-bp spacer
    occ, fp = _footprint(20, 31, prom_long)
    print(f"geometry (long spacer): {m} (conf {c:.2f}) -- {e}; occludes={occ}; overlaps={ov}")
    assert m == "elongated-spacer-overlap"
    assert ov["spacer"] > 0 and occ == ["spacer"], (ov, occ)     # footprint sits in the spacer only

    prom_core = dict(hex35=(10, 16), hex10=(33, 39), spacer=(16, 33), spacer_len=17, tss=42)
    m2, c2, e2, ov2 = _decide_mode(31, 42, prom_core)      # operator over the -10 box + TSS
    occ2, fp2 = _footprint(31, 42, prom_core)
    print(f"geometry (core overlap): {m2} (conf {c2:.2f}) -- {e2}; occludes={occ2}; "
          f"op_end-TSS={fp2['op_end_to_tss']}; elements={fp2.get('elements')}")
    assert m2 == "promoter-core-overlap"
    assert "-10" in occ2 and fp2["op_start_to_tss"] == 31 - 42, (occ2, fp2)
    assert fp2["elements"]["-10"] == [33 - 42, 39 - 42], fp2["elements"]     # -10 box relative to TSS

    m3, c3, e3, ov3 = _decide_mode(0, 1, prom_core)        # no overlap anywhere
    assert m3 == "no-promoter-overlap", m3

    # the classifier takes no family at all, so identical geometry cannot produce different calls
    import inspect as _insp
    assert "family" not in _insp.signature(_decide_mode).parameters
    assert "family" not in _insp.signature(classify_mode).parameters
    print("OK: every family gets the same geometric call -- no per-family mechanism labels.")

    # 2) integration: the calculator runs and a promoter is found on a real-ish sequence
    seq = ("GGGCGCGAACT" + "TTGACA" + "GCTAGCATCGATCGAT" + "TATAAT"
           + "GCATACTGGGCATGCATGCATGCGGGCCCAAATTTGGGACGT" * 2)
    proms = find_promoters(seq)
    print(f"calculator found {len(proms)} promoters; strongest Tx_rate={proms[0]['tx_rate']:.0f} "
          f"spacer={proms[0]['spacer_len']}bp")
    assert proms, "no promoter found"
    mc = classify_mode(seq, 18, 30)
    print(f"classify_mode: {mc.mode} (conf {mc.confidence}) spacer={mc.spacer_len} -- {mc.evidence}")
    assert mc.mode in ("elongated-spacer-overlap", "promoter-core-overlap", "no-promoter-overlap",
                       "ambiguous")
    print("OK: promoter positioning + geometry-first mode classification work.")


if __name__ == "__main__":
    _demo()
