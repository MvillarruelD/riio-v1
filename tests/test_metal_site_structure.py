"""Structure-based sensing-site typing in the coordination gate (two-site regulators)."""
from __future__ import annotations

from predictor.structure import metal_site as MS

_COLS = [
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id",
    "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
]


def _cif(atoms) -> str:
    """atoms: (group, element, atom_name, comp, chain, seq, x, y, z) -> minimal ModelCIF text."""
    lines = ["data_test", "loop_"] + [f"_atom_site.{c}" for c in _COLS]
    for i, (grp, el, an, comp, ch, seq, x, y, z) in enumerate(atoms, 1):
        lines.append(f"{grp} {i} {el} {an} {comp} {ch} {seq} {x:.3f} {y:.3f} {z:.3f}")
    lines.append("#")
    return "\n".join(lines)


def _cys4(cx, cy, cz, seqs):
    # a Zn thiolate cage: 4 CYS SG at ~2.3 A
    off = [(2.3, 0, 0), (-2.3, 0, 0), (0, 2.3, 0), (0, -2.3, 0)]
    return [("ATOM", "S", "SG", "CYS", "A", s, cx + dx, cy + dy, cz + dz)
            for s, (dx, dy, dz) in zip(seqs, off)]


def _no_site(cx, cy, cz, his_seqs, ox_seqs):
    # a His/carboxylate N/O site: 3 HIS NE2 + 2 ASP/GLU O at ~2.1 A
    a = [("ATOM", "N", "NE2", "HIS", "A", s, cx + dx, cy, cz)
         for s, dx in zip(his_seqs, (2.1, -2.1, 0))]
    a += [("ATOM", "O", "OD1", "ASP", "A", s, cx, cy + dy, cz)
          for s, dy in zip(ox_seqs, (2.1, -2.1))]
    return a


def test_gate_backward_compatible_without_fold():
    g = MS.coordination_gate("MDVSSHHHHCVCGCC", "Fur")
    assert "structure_sites" not in g and "structure_implied_ions" not in g


def test_two_site_fold_reports_regulatory_metal_only():
    # structural Zn(Cys4) cage (C-terminal seqs) + regulatory Co N/O site (bridging low seqs)
    atoms = [("HETATM", "ZN", "ZN", "ZN", "Z", 1, 0, 0, 0)] + _cys4(0, 0, 0, [90, 93, 130, 133])
    atoms += [("HETATM", "CO", "CO", "CO", "Z", 2, 15, 0, 0)] + _no_site(15, 0, 0, [30, 85, 88], [80, 100])
    res = MS.sensing_site_from_structure(_cif(atoms), "Fur")
    assert res["mode"] == "multi_site" and res["n_metal_sites"] == 2
    assert res["sensing_element"] == "CO"                       # the N/O regulatory site, not the cage
    assert res["structural_metals"] == ["Zn2+"]
    g = MS.coordination_gate("X", "Fur", fold_cif=_cif(atoms))
    # a surrogate names the class + Fe/Mn/Ni shortlist, never the structural Zn
    assert set(g["structure_implied_ions"]) == {"Fe2+", "Mn2+", "Ni2+"}
    assert "Zn2+" not in g["structure_implied_ions"]


def test_single_site_fold_is_passed_through():
    atoms = [("HETATM", "ZN", "ZN", "ZN", "Z", 1, 0, 0, 0)] + _cys4(0, 0, 0, [90, 93, 130, 133])
    res = MS.sensing_site_from_structure(_cif(atoms), "Fur")
    assert res["mode"] == "single_site" and res["n_metal_sites"] == 1
    assert res["sensing_metal"] == "Zn2+"


def test_malformed_fold_leaves_gate_unchanged():
    g = MS.coordination_gate("MDVSSHHHH", "Fur", fold_cif="not a cif")
    assert "structure_implied_ions" not in g
