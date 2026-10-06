"""Unchanged benchmark family constants and frozen paper membership."""
import json
from pathlib import Path

STUDY_FAMILY_PFAM = {
    "MerR": {"PF00376", "PF13411"},
    "ArsR/SmtB": {"PF01022", "PF12840"},
    "Fur": {"PF01475"},
    "CopY": {"PF03965"},
    "CsoR/FrmR": {"PF02583"},
    "DtxR/MntR": {"PF01325", "PF02742"},
    "GntR": {"PF00392"},
    "LysR-type (LTTR)": {"PF00126"},
    "MarR/SlyA": {"PF01047", "PF12802"},
    "NikR": {"PF08753"},
    "Rrf2": {"PF02082"},
    "TetR/AcrR": {"PF00440"},
}
PFAM_TO_STUDY = {pf: fam for fam, pfs in STUDY_FAMILY_PFAM.items() for pf in pfs}


def paper_family_members():
    path = Path(__file__).with_name("paper_family_members.json")
    return {key: set(value) for key, value in json.loads(path.read_text(encoding="utf-8")).items()}
