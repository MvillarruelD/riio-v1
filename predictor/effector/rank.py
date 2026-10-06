"""
rank.py -- operon-context score for a predicted TF<->ligand interaction.

A faithful reimplementation of Ligify's `ligify/predict/rank.py` (d'Oelsnitz et al., ACS Synth. Biol.
2024 / Ligify 2.0 2026). The premise: a TF most likely senses the substrate/product of the enzyme it is
co-located with, and that inference is *stronger* when the enzyme sits right next to the regulator, the
operon is small, and there are no competing regulators. The score starts at 100 and deducts:

    -10 * (enzyme-regulator gene distance - 1)     enzymes adjacent to the TF are most trustworthy
    -15 *  (number of OTHER regulators in operon)   competing regulators muddy the assignment
    - 5 * (operon gene count - 2)                    large operons dilute the one-enzyme-one-ligand logic

Colour bands match Ligify (green >=70, yellow 50-69, orange 30-49, red <30) for figure parity.

Run `python rank.py` for a self-test against the formula's reference values.
"""
from __future__ import annotations

import re

_REGULATOR_RE = re.compile(r"regulator|repressor|activator", re.I)


def _is_regulator(text: str) -> bool:
    return bool(text and _REGULATOR_RE.search(text))


def calculate_rank(operon, reg_index: int, enz_index: int, *, descriptions=None) -> dict:
    """Score an interaction. `operon` is the ordered list of genes (any objects); `reg_index`/`enz_index`
    are the regulator and ligand-associated-enzyme positions within it. Regulator-like neighbours are
    detected from `descriptions` (parallel list) or each gene's `.product`/`.name` attribute."""
    total_genes = len(operon)
    enz_reg_distance = abs(enz_index - reg_index)

    if descriptions is None:
        descriptions = [f"{getattr(g, 'product', '')} {getattr(g, 'name', '')}" for g in operon]
    # start at -1 so the regulator itself (always regulator-like) is not counted as a competitor
    total_regs = -1 + sum(1 for d in descriptions if _is_regulator(d))
    total_regs = max(0, total_regs)

    rank = 100
    distance_deduction = -10 * (enz_reg_distance - 1)
    extra_reg_deduction = -15 * total_regs
    total_genes_deduction = -5 * (total_genes - 2)
    rank = rank + distance_deduction + extra_reg_deduction + total_genes_deduction

    if rank >= 70:
        color = "#02a602"
    elif rank >= 50:
        color = "#d4d400"
    elif rank >= 30:
        color = "#d48302"
    else:
        color = "#f50b02"

    return {"rank": rank, "color": color, "metrics": {
        "Genes within operon": {"Value": total_genes, "Deduction": total_genes_deduction},
        "Enzyme-regulator distance": {"Value": enz_reg_distance, "Deduction": distance_deduction},
        "Additional regulators": {"Value": total_regs, "Deduction": extra_reg_deduction},
    }}


def _demo() -> None:
    # adjacent enzyme, 2-gene operon, no extra regulators -> perfect 100 (the ideal Ligify case)
    descs = ["transcriptional regulator", "alcohol dehydrogenase"]
    r = calculate_rank([0, 1], reg_index=0, enz_index=1, descriptions=descs)
    print(f"ideal: rank={r['rank']} color={r['color']} {r['metrics']}")
    assert r["rank"] == 100 and r["color"] == "#02a602", "ideal case should be 100/green"

    # enzyme 2 genes away, 4-gene operon, one competing regulator
    descs2 = ["ArsR family repressor", "hypothetical protein", "metal transporter",
              "MerR family activator"]
    r2 = calculate_rank([0, 1, 2, 3], reg_index=0, enz_index=2, descriptions=descs2)
    # 100 -10*(2-1) -15*1 -5*(4-2) = 100 -10 -15 -10 = 65
    print(f"messier: rank={r2['rank']} {r2['metrics']}")
    assert r2["rank"] == 65, f"expected 65, got {r2['rank']}"
    assert r2["metrics"]["Additional regulators"]["Value"] == 1, "competing regulator not counted"
    print("OK: rank formula matches Ligify reference values (100/green ideal; 65 messier).")


if __name__ == "__main__":
    _demo()
