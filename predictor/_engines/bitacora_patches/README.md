# BITACORA patches

`env/tools/bitacora/` is the upstream checkout: gitignored and user-installed, so anything changed
there is lost on reinstall and invisible to anyone reproducing the run. The patched files live here,
tracked, and are copied over the install by `predictor.discovery.bitacora.ensure_patched()`.

## `get_blastp_parsed_newv2.pl`

Upstream keeps a BLASTP HSP if it covers **2/3 of the query OR 80 % of the subject**. That rule, not
the E-value, is what determines our candidate set: of the 16 candidates the survey lists and our run
did not recover, 11 were matched at our own E-value and removed afterwards by this filter.

Measured over all four genomes, holding the E-value at 1e-5 and moving only this rule
(`analysis/benchmarking/discovery/adjudicate_survey_gap.py`):

| qcov / scov | blastp candidates | of the survey's 150 | not listed by the survey |
|---|---:|---:|---:|
| 0.67 / 0.80 (upstream) | 121 | 116 | 5 |
| 0.60 / 0.70 | 130 | 123 | 7 |
| 0.50 / 0.60 | 140 | 129 | 11 |
| **0.50 / 0.50** | **144** | **133** | **11** |
| 0.40 / 0.40 | 150 | 138 | 12 |

Loosening the E-value to 1e-3 instead buys a survey candidate for **6.9** proteins the survey does
not list; relaxing coverage to 0.50/0.50 buys one for about **0.35**. Roughly a twentyfold better
exchange rate, which is why this is the knob we turn.

The patch does not hardcode the new values. It reads `BITACORA_QCOV` and `BITACORA_SCOV`, defaulting
to upstream's 2/3 and 0.8, so an unset environment reproduces stock BITACORA exactly and the change
is visible in the run configuration rather than buried in a vendored script.
