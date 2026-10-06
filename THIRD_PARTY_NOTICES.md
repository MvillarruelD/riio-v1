# Third-party credit and license scope

The root MIT license covers original RIIO software. It does not replace the licenses of third-party
components, reference databases, model weights, executables, or service outputs.

| Component | Use and credit | Distribution and terms |
|---|---|---|
| Metalloregulator-mining / SSN source project | Johnma J. Rondón, Giuliano T. Antelo, David P. Giedroc, Daiana A. Capdevila, as recorded by the source citation | Runtime reference material and upstream citation retained; MIT notice in `LICENSES/Metalloregulator-mining-MIT.txt`. Reference-only mining scripts excluded. |
| [Ligify](https://github.com/groov-bio/ligify-ui) | d'Oelsnitz et al.; published sensor records and locally implemented operon-context rank method | Compressed reference subset retains source hash/citation metadata. Full upstream MIT notice: `LICENSES/Ligify-MIT.txt`. |
| [BITACORA](https://github.com/molevol-ub/bitacora) | Vizueta, Sánchez-Gracia and Rozas; separate discovery engine | Engine not bundled. The included modified `get_blastp_parsed_newv2.pl` is an upstream-derived GPL-3.0-only component. Full license: `LICENSES/BITACORA-GPL-3.0.txt`; modification and upstream commit documented in `docs/ENGINES.md`. |
| [MetalNet2](https://github.com/wangchulab/MetalNet2) | Wang laboratory / upstream project authors; metal-site prediction | External checkout and weights not bundled or relicensed. Owned subprocess runner is included; upstream setup and model-use terms apply. |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Steinegger, Söding and upstream contributors; sequence search | External executable, not bundled. |
| [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/), [Datasets](https://www.ncbi.nlm.nih.gov/datasets/) | NCBI authors; sequence search and genome acquisition | External executables and remotely retrieved records, not bundled. |
| [HMMER / pyhmmer](https://pyhmmer.readthedocs.io/) | Eddy / HMMER team; Larralde / pyhmmer contributors | Installed as a dependency; retains upstream terms. |
| [Pfam / InterPro](https://www.ebi.ac.uk/interpro/), [UniProt](https://www.uniprot.org/) | Upstream database contributors; family models and protein reference annotations | Bundled selected reference assets and accession membership retain source identifiers; database terms apply, independently of the RIIO code license. |
| [RegulonDB](https://regulondb.ccg.unam.mx/) | RegulonDB contributors; external benchmark reference | Selected benchmark tables retain `KB_SOURCES.json`. Full source-mirror acquisition/building is outside the runtime export. Upstream data terms apply. |
| [Rhea](https://www.rhea-db.org/), [ChEBI](https://www.ebi.ac.uk/chebi/) | Upstream database contributors; reaction and chemical annotations | Online lookup and cached evidence; upstream terms apply. |
| ESMFold2 / AlphaFold services | Upstream software/model/service authors | Optional adapter retained; no weights, engine binaries, service credentials, or manuscript structures distributed. |

The software relies on methods described in its module comments (including local motif search and
conservation workflows). Method inspiration is distinct from copied code. Python dependencies retain
their own package notices. No claim of newly inventing those methods or owning upstream data is made.

Full license notices are included in wheels and source distributions through `license-files` metadata.
The package metadata uses `MIT AND GPL-3.0-only` to accurately identify the separately licensed patch.
