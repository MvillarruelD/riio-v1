"""Extract the vector drawings and the SSN clade layout that Figure 1 reuses.

Run once; ``reports/build_fig1_artwork.py`` only reads the files written here.

    python extract_assets.py --art TF_DNA_REFs.ai --session MerR_AS32.cys

Art: the author's Illustrator file is read through its PDF compatibility layer
(PyMuPDF) and each object is rewritten as plain SVG paths, so colours and
strokes stay editable. SSN: the MerR session's thresholded view supplies the
node positions and edges of the clade whose representative node holds E. coli
CueR (P0A9G4), i.e. cluster 7.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ASSETS = HERE

# PDF-point bounding boxes of the objects on the artwork page.
ART_OBJECTS = {
    "cuer_complex": (66, 20, 212, 135),
    "cuer_protein": (70, 22, 209, 126),       # the same drawing with the DNA removed
    "dna_long": (105, 194, 306, 213),
    "dna_medium": (105, 213, 252, 237),
    "dna_short": (296, 213, 380, 237),
}
DNA_FILLS = {"#d79583", "#ba8172", "#bd7960", "#dfab9d", "#a26853"}
CLADE_ANCHORS = {"ecCueR": "P0A9G4", "vcCueR": "P0C6D2", "ypCueR": "Q8ZCA8", "rmHmrR1": "P58378"}


def _hex(c):
    return "none" if c is None else "#%02x%02x%02x" % tuple(round(v * 255) for v in c[:3])


def _num(v):
    return ("%.2f" % v).rstrip("0").rstrip(".")


def _path_data(items, ox, oy):
    out, cur = [], None
    for it in items:
        kind = it[0]
        if kind == "re":
            r = it[1]
            out.append(f"M{_num(r.x0 - ox)} {_num(r.y0 - oy)}H{_num(r.x1 - ox)}"
                       f"V{_num(r.y1 - oy)}H{_num(r.x0 - ox)}Z")
            cur = None
            continue
        if kind == "qu":
            q = it[1]
            pts = [(p.x - ox, p.y - oy) for p in (q.ul, q.ur, q.lr, q.ll)]
            out.append("M" + "L".join(f"{_num(x)} {_num(y)}" for x, y in pts) + "Z")
            cur = None
            continue
        p0 = it[1]
        if cur is None or abs(cur.x - p0.x) > 1e-3 or abs(cur.y - p0.y) > 1e-3:
            out.append(f"M{_num(p0.x - ox)} {_num(p0.y - oy)}")
        if kind == "l":
            cur = it[2]
            out.append(f"L{_num(cur.x - ox)} {_num(cur.y - oy)}")
        elif kind == "c":
            c1, c2, cur = it[2], it[3], it[4]
            out.append(f"C{_num(c1.x - ox)} {_num(c1.y - oy)} {_num(c2.x - ox)} {_num(c2.y - oy)} "
                       f"{_num(cur.x - ox)} {_num(cur.y - oy)}")
    return "".join(out)


def _subpaths(items):
    """Split a compound path's items into subpaths (a new one starts wherever the pen jumps)."""
    subs, cur, last = [], [], None
    for it in items:
        if it[0] in ("re", "qu"):
            subs.append([it])
            last = None
            continue
        if last is None or abs(last.x - it[1].x) > 1e-3 or abs(last.y - it[1].y) > 1e-3:
            if cur:
                subs.append(cur)
            cur = []
        cur.append(it)
        last = it[-1]
    if cur:
        subs.append(cur)
    return subs


def _bbox(sub):
    pts = []
    for it in sub:
        if it[0] == "re":
            pts += [it[1].tl, it[1].br]
        elif it[0] == "qu":
            pts += [it[1].ul, it[1].ur, it[1].lr, it[1].ll]
        else:
            pts += [p for p in it[1:] if hasattr(p, "x")]
    return min(p.x for p in pts), min(p.y for p in pts), max(p.x for p in pts), max(p.y for p in pts)


def _drop_holes(items):
    """Keep the outer contours: drop subpaths lying inside another subpath of the same path."""
    subs = _subpaths(items)
    boxes = [_bbox(s) for s in subs]
    keep = [s for i, (s, b) in enumerate(zip(subs, boxes))
            if not any(j != i and o[0] <= b[0] and o[1] <= b[1] and o[2] >= b[2] and o[3] >= b[3]
                       and (o[2] - o[0]) * (o[3] - o[1]) > (b[2] - b[0]) * (b[3] - b[1])
                       for j, o in enumerate(boxes))]
    return [it for s in keep for it in s]


# Where the DNA crossed the dark monomer's lower lobe, the drawing cut its outer contour too.
# Bridge that cut with one curve (start, control 1, control 2, end in PDF points).
NOTCH_BRIDGE = ((182.5, 109.54), (184.6, 113.4), (189.2, 118.6), (192.76, 119.24))


def _bridge_notch(items):
    import pymupdf

    a, c1, c2, b = (pymupdf.Point(*p) for p in NOTCH_BRIDGE)
    near = [k for k, it in enumerate(items) if it[0] in "lc" and abs(it[-1] - a) < 0.3]
    far = [k for k, it in enumerate(items) if it[0] in "lc" and abs(it[1] - b) < 0.3]
    if not near or not far or far[0] <= near[0]:
        return items
    ia, ib = near[0], far[0]
    return items[:ia + 1] + [("c", items[ia][-1], c1, c2, items[ib][1])] + items[ib:]


def _art_paths(page, bbox, *, protein_only=False):
    import pymupdf

    box = pymupdf.Rect(bbox)
    paths = []
    for g in page.get_drawings():
        r = g["rect"]
        if not box.contains(pymupdf.Point((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)):
            continue
        fill, stroke = _hex(g.get("fill")), _hex(g.get("color"))
        width = g.get("width") or 0
        items = g["items"]
        if protein_only:
            # Drop the DNA fills, the DNA outlines and the combined outline that traces the DNA.
            if fill in DNA_FILLS or (stroke != "none" and (r.y0 > 70 or width > 0.5)):
                continue
            # The protein surfaces are cut where the DNA passed in front of them; close those cuts.
            if g.get("even_odd") or stroke != "none":
                items = _drop_holes(_bridge_notch(items))
        d = _path_data(items, bbox[0], bbox[1]) + ("Z" if g.get("closePath") else "")
        attrs = [f'd="{d}"', f'fill="{fill}"']
        if g.get("even_odd") and fill != "none" and not protein_only:
            attrs.append('fill-rule="evenodd"')
        if stroke != "none" and width:
            attrs += [f'stroke="{stroke}"', f'stroke-width="{_num(width)}"',
                      'stroke-linecap="round"', 'stroke-linejoin="round"']
        paths.append("<path " + " ".join(attrs) + "/>")
    return paths


def extract_art(art: Path):
    import pymupdf

    page = pymupdf.open(art)[0]
    for name, bbox in ART_OBJECTS.items():
        paths = _art_paths(page, bbox, protein_only=name == "cuer_protein")
        w, h = bbox[2] - bbox[0], bbox[3] - bbox[1]
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
               f'width="{w}pt" height="{h}pt">\n' + "\n".join(paths) + "\n</svg>\n")
        (ASSETS / f"{name}.svg").write_text(svg, encoding="utf-8")
        print(f"{name}: {len(paths)} paths, {w} x {h} pt")


def extract_clade(session: Path):
    with zipfile.ZipFile(session) as z:
        names = z.namelist()

        def read(pattern):
            return z.read(next(n for n in names if re.search(pattern, n))).decode("utf-8")

        net = read(r"/networks/87-MerR\.xgmml$")
        view = read(r"/views/128-[^/]*MerR\.xgmml$")
        sub_edges = read(r"/tables/128-MerR/LOCAL_ATTRS-[^/]*CyEdge")
        node_tab = read(r"/tables/87-MerR/SHARED_ATTRS-[^/]*CyNode")

    id2acc = dict(re.findall(r'<node id="(\d+)" label="([^"]+)"', net))
    all_edges = {i: (s, t) for i, s, t in
                 re.findall(r'<edge id="(\d+)" label="[^"]*" source="(\d+)" target="(\d+)"', net)}
    csv.field_size_limit(10 ** 9)
    shown = [row[0] for row in list(csv.reader(sub_edges.splitlines()))[5:]]
    edges = [all_edges[i] for i in shown if i in all_edges]
    adj = collections.defaultdict(set)
    for s, t in edges:
        adj[s].add(t)
        adj[t].add(s)
    rows = list(csv.reader(node_tab.splitlines()))[5:]
    anchors = {name: next(r[0] for r in rows if any(acc in c for c in r[1:]))
               for name, acc in CLADE_ANCHORS.items()}
    comp, stack = {anchors["ecCueR"]}, [anchors["ecCueR"]]
    while stack:
        for m in adj[stack.pop()]:
            if m not in comp:
                comp.add(m)
                stack.append(m)
    xy = {acc: (float(x), float(y)) for acc, x, y in
          re.findall(r'<node id="\d+" label="([^"]+)"[^>]*>\s*<graphics x="([-\d.eE]+)" z="[^"]*" y="([-\d.eE]+)"', view)}
    nodes = {n: xy[id2acc[n]] for n in comp if id2acc.get(n) in xy}
    out = {
        "source": "MerR SSN (alignment score 32), cluster 7; Cytoscape view coordinates",
        "nodes": {n: [round(x, 1), round(y, 1)] for n, (x, y) in nodes.items()},
        "edges": [[s, t] for s, t in edges if s in nodes and t in nodes],
        "anchors": {k: v for k, v in anchors.items() if v in nodes},
    }
    (ASSETS / "merr_clade7.json").write_text(json.dumps(out), encoding="utf-8")
    print(f"clade 7: {len(out['nodes'])} nodes, {len(out['edges'])} edges, anchors {sorted(out['anchors'])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--art", type=Path, help="Illustrator/PDF file holding the TF-DNA drawings")
    ap.add_argument("--session", type=Path, help="Cytoscape .cys session of the MerR SSN")
    args = ap.parse_args()
    ASSETS.mkdir(exist_ok=True)
    if args.art:
        extract_art(args.art)
    if args.session:
        extract_clade(args.session)


if __name__ == "__main__":
    main()
