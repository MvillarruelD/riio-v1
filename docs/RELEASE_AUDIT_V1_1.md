# RIIO v1.1 release preparation audit

These checks cover the local source/export snapshots listed below. A passing software
check does not establish experimental biological accuracy or historical prediction neutrality.

## Source provenance

- RIO: `5ab06c7d19207293944ddf383d9f6cc14e02efb1`
- jjr-ssn-analysis: `174d1a5bafdba8e7ae04bd946d27834c3d9a821b`
- publication_checkout: `1840765caa94d5861183651e7347126a2c138799`

Complete file checksums are in [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json).
The [v1.0 release audit](RELEASE_AUDIT.md) and `EXPORT_AUDIT.json` are historical records.

## Completed checks

- Export provenance and metadata: Clean predictor 5ab06c7 and analysis 174d1a5, preserving publication donor 1840765; all registered file checksums verified, 480 reference data files byte-identical to frozen source, package and citation version 1.1.0, README local links valid.
- Python 3.10 generic staging suite: 861 passed, 8 explicitly scoped skips, 51 released-data tests deselected.
- Python 3.10 released-data staging suite: 48 passed, 3 scoped skips (two unavailable MMseqs fixtures and a non-Git staging check), 869 generic tests deselected.
- Static checks: Both predictor and analysis Ruff scopes and compileall passed; data MANIFEST current.
- Isolated package installation: Wheel and sdist built; clean non-editable wheel install and outside-checkout CLI passed, 46 installed module selftests passed with 2 missing engine/network skips, no timeouts/failures, outputs outside the read-only installed package.
- Source archive installation closure: 1039 members, Linux installer/environment, Streamlit configuration, environment example and all supplied smoke input files present; no private environment or env/tools files.
- Latest source software verification: Predictor 5ab06c7: Windows Python 3.12 full suite 718 passed/4 skipped including 45 released-data tests; native Ubuntu Python 3.12 report/GUI/contracts 233 passed/0 skipped, six changed paths synchronized exactly, local environment preserved; software/report verification only.
- Browser reports and downloads: Stored prior real K12 fixtures rebuilt with latest software on Windows/native: each full-run ZIP 466 members, 486 HTML/XLSX links, zero broken links/CRC errors, identical member bytes; keyboard numeric sorting and saved per-TF report figure loaded; AF3 ZIP 9 files with 6 jobs per local/server dialect and manifest; no new browser prediction claimed.
- Compact Figure 1 review artifact: Approved editable compact Figure 1 and restored drawings included in the separately retained review DOCX; 12 embedded images, exact final Figure 1 PNG bytes and OOXML validation passed; original manuscript DOCX/PDF preserved.
- Current user Windows environment: Python 3.14.3 at frozen predictor 5ab06c7: 718 passed and 4 skipped in the full software suite including 45 released-data tests; Ruff, compileall and data manifest passed; recorded source, reference and dependency state unchanged; compatibility evidence only, historical dependency identity and prediction neutrality remain unestablished.
- Final canonical share: 238 saved candidates and 18,084 artifacts with zero HTML/XLSX link issues; all 16,946 original bundle files byte-identical before/after final 5ab06c7 renderer refresh; final compact review DOCX copied into a fresh sibling, no new biological predictions.
- Git checksum portability: Effective checked-in Git attributes normalized 17 text files before release hashing; original source SHA256 retained, all 480 reference files and 238 canonical FASTAs preserved; real commit/clone regression passed for CRLF text, LF shebangs and immutable or binary bytes.

## Outstanding gates

- Exact historical 12-TF prediction-neutrality comparison.
- Remote Linux/Windows CI after an approved push.
- Manual AF3 server upload and acceptance.
- PDF conversion of the new review DOCX (converter unavailable).

No push, tag, release or main-branch merge is represented by this audit.
