"""effector -- predict a TF's cognate ligand/effector (a local reimplementation of Ligify).

`ligify.predict_effector(seq)` is the entry point: TF sequence -> genome neighborhood (offline mirror
first) -> operon enzymes -> their Rhea/ChEBI reaction metabolites -> ranked candidate effectors.
"""
