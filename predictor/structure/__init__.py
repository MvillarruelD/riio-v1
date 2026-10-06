"""Structure: an OPTIONAL branch, plus two sequence-based modules that only live here.

**The optional branch** is `early_structure` (ESMFold2 apo homo-oligomer) and `afdb` (AlphaFold-DB
monomer retrieval and QC), with `oligomer` supplying the biological oligomeric state -- dimer by
default, tetramer for CsoR-type -- because effector and DNA contacts sit at the monomer-monomer
interface and are lost in a monomer fold.

It is OFF by default, it enriches the report, and it changes NO prediction: the family call, the
inducer call, the operators and the regulon are identical with it on or off. Its only trace in the
operator list is a `pending_af3` placeholder that is never ranked. `folddisco` structural retrieval
used to be part of it and was deleted -- measured over 139 candidates, it answered none of them.

Turn it on with `predict(..., fold=True)`. The default used to be "fold if a backend token happens to
exist", which made the contents of a bundle depend on the machine and on run history: 95 of the last
run's 140 candidates folded and 45 did not, successes cached and failures not.

**Not part of that branch, despite living here:**

* `metal_site` -- the family coordination gate. It reads the SEQUENCE, not a structure, and it is
  tier-1 and the single most load-bearing evidence source in the pipeline: removing it changes 74 of
  139 headline inducer calls. Turning "structure" off must never turn this off.
* `metalnet` -- the MetalNet2 site detector. Also sequence-based (ESM-2 embeddings), also tier-1, and
  deliberately never names an ion.

The package name is historical. Read it as "things about the protein's fold, real or inferred", and
check which group a module is in before assuming a switch reaches it.
"""
