#!/usr/bin/env python
"""Measure every metal-coordination signature against the shipped SSN member DBs and write the reference.

    python tools/build_metal_reference.py            # regenerate both outputs
    python tools/build_metal_reference.py --check    # fail if they are stale (for CI)

Two outputs, both derived, neither hand-maintained:

  * `predictor/data/refs/metal_signatures.json` -- machine-readable; shipped with the package.
  * `docs/METAL_SIGNATURES.md`                  -- the Methods-ready table.

WHY THIS IS GENERATED. `predictor/structure/metal_site.py` states that the NikR H-x-H-x(4,7)-C site
occurs in "<=3.2 % of every other family (Fur 3.2 %)". Measured against the member DB the package
actually ships, it occurs in 60.6 % of Fur. The claim was written once, by hand, and nothing ever
re-checked it against the data it described -- so it drifted silently as the vendored DBs changed, and it
was on its way into a Methods section. Every prevalence in the outputs of this script is computed from
`predictor/data/ssn/database/ssn_members_<TAG>.fasta` at build time, so the same class of error cannot
recur: if a motif or a member DB changes and nobody regenerates, `--check` fails.

The per-CLADE column is the one that matters scientifically. A motif that fires on 3 % of MerR looks
useless until you see it fires on 60 % of MerR_9 and 0 % of every other clade -- at which point it is a
zinc-sensor detector, and the family-level number was the wrong denominator.
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))

from predictor.annotate import ssn_clusters as ssn        # noqa: E402
from predictor.structure import metal_motifs as MM        # noqa: E402
from predictor import resources                           # noqa: E402

JSON_OUT = _REPO / "predictor" / "data" / "refs" / "metal_signatures.json"
MD_OUT = _REPO / "docs" / "METAL_SIGNATURES.md"

#: A motif is safe to use without a family gate when it stays this rare in every OTHER family. 2 % is
#: not a tuned value: it is an order of magnitude below the weakest in-family signal we act on, and the
#: point of the column is to separate "specific" from "fires everywhere", which no borderline case does.
GLOBAL_SAFE_MAX_OTHER = 0.02


def _read_members(tag: str) -> list[tuple[str, str]]:
    """[(clade_id, sequence)] from a vendored member DB. Headers are `<cluster_id>__<accession>`."""
    p = resources.SSN_DATABASE / f"ssn_members_{tag}.fasta"
    if not p.is_file():
        return []
    out, cid, buf = [], None, []
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if cid is not None:
                out.append((cid, "".join(buf).upper()))
            cid = ln[1:].strip().split("__", 1)[0]
            buf = []
        else:
            buf.append(ln.strip())
    if cid is not None:
        out.append((cid, "".join(buf).upper()))
    return out


def measure() -> dict:
    """Prevalence of every implemented motif, per family and per clade."""
    corpus = {fam: _read_members(ssn.family_tag(fam)) for fam in ssn.FAMILIES}
    n_by_family = {f: len(v) for f, v in corpus.items()}

    entries = []
    for m in MM.MOTIFS:
        row = {
            "name": m.name, "family": m.family, "pattern": m.pattern,
            "implies": list(m.implies), "implied_class": m.implied_class,
            "geometry": m.geometry, "coordination_number": m.coordination_number,
            "ligand_set": m.ligand_set, "scope": m.scope, "names_ion": m.names_ion,
            "implemented": m.implemented, "source": m.source, "note": m.note,
        }
        if not m.implemented:
            row["prevalence"] = None
            entries.append(row)
            continue

        rx = re.compile(m.pattern)
        per_family, per_clade = {}, {}
        for fam, recs in corpus.items():
            if not recs:
                continue
            hits = sum(1 for _cid, s in recs if rx.search(s))
            per_family[fam] = round(hits / len(recs), 4)
            if fam == m.family:
                by_clade = collections.defaultdict(lambda: [0, 0])
                for cid, s in recs:
                    by_clade[cid][1] += 1
                    if rx.search(s):
                        by_clade[cid][0] += 1
                per_clade = {c: round(h / n, 4) for c, (h, n) in sorted(by_clade.items())
                             if n}
        own = per_family.get(m.family)
        others = {f: v for f, v in per_family.items() if f != m.family}
        max_other = max(others.values()) if others else 0.0
        max_other_family = max(others, key=others.get) if others else None
        top_clade = max(per_clade, key=per_clade.get) if per_clade else None
        row["prevalence"] = {
            "in_family": own,
            "per_family": per_family,
            "per_clade": per_clade,
            "top_clade": top_clade,
            "top_clade_rate": per_clade.get(top_clade) if top_clade else None,
            "max_other_family": round(max_other, 4),
            "max_other_family_name": max_other_family,
            "global_safe": bool(max_other <= GLOBAL_SAFE_MAX_OTHER),
        }
        entries.append(row)

    return {
        "generated_by": "tools/build_metal_reference.py",
        "catalogue": "predictor/structure/metal_motifs.py",
        "corpus": {"source": "predictor/data/ssn/database/ssn_members_<TAG>.fasta",
                   "n_by_family": n_by_family, "n_total": sum(n_by_family.values())},
        "global_safe_threshold": GLOBAL_SAFE_MAX_OTHER,
        "selectivity": MM.SELECTIVITY,
        "irving_williams": list(MM.IRVING_WILLIAMS),
        "motifs": entries,
    }


def _pct(x) -> str:
    return "—" if x is None else f"{100 * x:.1f} %"


def render_markdown(data: dict) -> str:
    L: list[str] = []
    A = L.append
    n_total = data["corpus"]["n_total"]
    impl = [m for m in data["motifs"] if m["implemented"]]
    doc = [m for m in data["motifs"] if not m["implemented"]]

    A("# Metal-coordination signatures")
    A("")
    A("**Generated file — do not edit.** Regenerate with `python tools/build_metal_reference.py`; "
      "`--check` fails when it is stale. The signatures themselves live in "
      "`predictor/structure/metal_motifs.py`, and every prevalence below is measured at build time "
      f"against the {n_total:,} sequences of the shipped SSN member databases "
      "(`predictor/data/ssn/database/ssn_members_<TAG>.fasta`).")
    A("")
    A("This file exists because a prevalence stated by hand drifts. `structure/metal_site.py` records "
      "that the NikR H-x-H-x(4,7)-C site occurs in \"≤3.2 % of every other family\"; measured against "
      "the database this package ships, it occurs in **60.6 % of Fur**. Nothing re-checked the claim "
      "after it was written. Everything here is recomputed from the data instead.")
    A("")
    A("## How to read it")
    A("")
    A("- **in-family** — the share of that family's members carrying the motif. On its own this is often "
      "the wrong denominator: a motif specific to one clade looks weak against the whole family.")
    A("- **top clade** — the clade where it concentrates. This is usually the real signal.")
    A("- **max other family** — the highest rate in any family it does *not* belong to. This decides "
      "whether the motif can be trusted without a family gate.")
    A(f"- **scope** — `global` means the motif is specific enough to use anywhere (measured "
      f"≤ {100 * data['global_safe_threshold']:.0f} % in every other family); `family` means it is only "
      "meaningful inside its own family and demonstrably fires elsewhere.")
    A("- **names ion** — whether this motif alone may name an ion in the inducer call. A motif marking "
      "a chemistry shared by several ions (a bare Cys pair) documents the site without naming the metal.")
    A("")
    A("## Detected signatures")
    A("")
    A("| motif | family | pattern | implies | class | in-family | top clade | max other family | scope | names ion |")
    A("|---|---|---|---|---|---:|---|---:|---|:-:|")
    for m in impl:
        p = m["prevalence"]
        top = (f"`{p['top_clade']}` {_pct(p['top_clade_rate'])}" if p["top_clade"] else "—")
        other = (f"{_pct(p['max_other_family'])}"
                 + (f" ({p['max_other_family_name']})" if p["max_other_family_name"] else ""))
        flag = "⚠ " if not p["global_safe"] and m["scope"] == "global" else ""
        A(f"| `{m['name']}` | {m['family'] or '—'} | `{m['pattern']}` | "
          f"{', '.join(m['implies']) or '—'} | {m['implied_class'] or '—'} | "
          f"{_pct(p['in_family'])} | {top} | {flag}{other} | {m['scope']} | "
          f"{'yes' if m['names_ion'] else 'no'} |")
    A("")
    A("## Documented but not detected")
    A("")
    A("These are real chemistry from the literature that no sequence pattern captures — the ligands come "
      "from different structural elements, or the discriminator is a *negative* design element, or the "
      "difference is a fold rather than a composition. They are listed so this table is a complete "
      "account of what each family senses, not only the automatable part.")
    A("")
    A("| motif | family | implies | ligand set | why not detected |")
    A("|---|---|---|---|---|")
    for m in doc:
        why = m["note"].split(". ", 1)[-1] if ". " in m["note"] else m["note"]
        A(f"| `{m['name']}` | {m['family'] or '—'} | {', '.join(m['implies']) or '—'} | "
          f"{m['ligand_set'] or '—'} | {why} |")
    A("")
    A("## Per-clade detail")
    A("")
    for m in impl:
        pc = m["prevalence"]["per_clade"]
        if not pc:
            continue
        A(f"**`{m['name']}`** — {m['family']}")
        A("")
        A("| clade | prevalence |")
        A("|---|---:|")
        for cid, v in sorted(pc.items(), key=lambda kv: -kv[1]):
            A(f"| `{cid}` | {_pct(v)} |")
        A("")
    A("## Selectivity priors")
    A("")
    A("From Dudev & Lim, *Chem. Rev.* 2014, 114, 538–556 (their Table 1 for coordination number and "
      "geometry; Fig. 10 and §6.3 for ligand preference).")
    A("")
    A("**Use these to sanity-check, never to rank.** A geometry that contradicts an ion is evidence the "
      "site is not that ion's. A geometry consistent with three ions does not order them: Dudev & Lim's "
      "own conclusion is that the protein matrix and the cellular free-ion concentration decide the "
      "winner, and this pipeline models neither. NikR binds Cu(II) more tightly than Ni(II) *in vitro* "
      "and is still a nickel sensor *in vivo*, because free cytosolic Cu is ~10⁻¹⁸ M.")
    A("")
    A("| ion | coordination number | geometry | preferred ligand |")
    A("|---|---|---|---|")
    for ion, v in data["selectivity"].items():
        A(f"| {ion} | {v['cn']} | {v['geometry']} | {v['prefers']} |")
    A("")
    A("Irving–Williams order of intrinsic affinity for any ligand set: "
      + " < ".join(data["irving_williams"]) + ". Every metalloregulator works *against* this ordering, "
      "which is why the cell holds free Zn at 10⁻¹²–10⁻¹⁵ M and free Cu at ~10⁻¹⁸ M.")
    A("")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if either output is stale")
    a = ap.parse_args(argv)

    data = measure()
    js = json.dumps(data, indent=2, sort_keys=False) + "\n"
    md = render_markdown(data)

    if a.check:
        stale = []
        for path, want in ((JSON_OUT, js), (MD_OUT, md)):
            have = path.read_text(encoding="utf-8") if path.is_file() else None
            if have != want:
                stale.append(path.relative_to(_REPO))
        if stale:
            print("STALE: " + ", ".join(str(p) for p in stale))
            print("regenerate with:  python tools/build_metal_reference.py")
            return 1
        print("metal signature reference is current")
        return 0

    JSON_OUT.parent.mkdir(parents=True, exist_ok=True)
    MD_OUT.parent.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(js, encoding="utf-8", newline="\n")
    MD_OUT.write_text(md, encoding="utf-8", newline="\n")
    n_impl = sum(1 for m in data["motifs"] if m["implemented"])
    n_unsafe = sum(1 for m in data["motifs"]
                   if m["implemented"] and m["scope"] == "global"
                   and not m["prevalence"]["global_safe"])
    print(f"wrote {JSON_OUT.relative_to(_REPO)}")
    print(f"wrote {MD_OUT.relative_to(_REPO)}")
    print(f"{len(data['motifs'])} signatures ({n_impl} measured) over "
          f"{data['corpus']['n_total']:,} member sequences")
    if n_unsafe:
        print(f"WARNING: {n_unsafe} motif(s) declared global scope but exceed the "
              f"{100 * GLOBAL_SAFE_MAX_OTHER:.0f} % out-of-family limit -- review their scope")
    return 0


if __name__ == "__main__":
    sys.exit(main())
