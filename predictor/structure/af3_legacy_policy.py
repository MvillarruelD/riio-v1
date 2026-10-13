"""Frozen Figure 4 recipe, retained exactly for historical reproduction.

This module records the existing paper hand-off policy. Live AF3 job exports use
af3_jobset and do not infer a physiological DNA-binding form from these labels.
"""

SURROGATE = {
    "Cd2+": ("ZN", "close",
             "Cd(II) and Zn(II) are both d10 group-12 ions that fill Cys-thiolate sites in "
             "trigonal/tetrahedral geometry; Zn(II) is the usual structural stand-in for Cd(II) in "
             "ArsR/SmtB sensors"),
    "Cu+": ("CU", "same element, wrong oxidation state",
            "the server's CU is Cu(II). Cu(I) is d10 and prefers linear or trigonal S-donor "
            "geometry, which Cu(II) (Jahn-Teller, square-planar) does not reproduce; kept as the "
            "same element because no d10 monovalent coin-metal ion is available"),
    "Fe2+": ("MN", "close",
             "Mn(II) is the standard redox-inert Fe(II) surrogate in Fur and DtxR/IdeR "
             "crystallography: same +2 charge, high-spin radius 0.83 vs 0.78 A, same octahedral "
             "N/O preference. The server's FE is Fe(III) -- same element but +3 and much smaller, "
             "so it matches the coordination less well"),
    "As(III)": ("ZN", "weak",
                "a metalloid, not an ion: As(III) sits in a trigonal-pyramidal As(Cys)3 site with a "
                "stereochemically active lone pair, which no server ion reproduces. Zn(II) is the "
                "nearest thiophilic surrogate, so treat this job as a Cys3/Cys4 site placeholder "
                "rather than a model of arsenite binding"),
}

DNA_BOUND_FORM = {"ArsR/SmtB": "apo", "CsoR/FrmR": "apo", "MerR": "apo",
                  "Fur": "holo", "DtxR/MntR": "holo"}
