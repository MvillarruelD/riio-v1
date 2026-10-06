# Report file contract

`python -m predictor.report.bundle_contract BUNDLE` performs an installed,
filesystem-only audit. Its JSON result reports validity, manifest version, assurance
level, number of hashes checked, whether `REPORT.html` was inspected, and findings.
Use `--no-report` only for an intentionally machine-readable export whose manifest
sets `report_generated` to false and which contains no `REPORT.html`. The flag never
suppresses inspection of an HTML file that is present. This check does not establish
biological validity, original-input completeness, or computational reproducibility.

New prediction bundles separately include a run-level traceability record. The
normalized TF query, exact genome sequence scanned, parsed gene model, discovery
windows, retained homolog promoters, source alignment when available, effective
configuration, motif-width trial audit, stage ledger, evidence index, structure
files and complete output tables live below the individual TF directory. Paths in
`dossier.json` and `run/run_record.json` are relative to that directory; publication
fails if either record retains an absolute machine-local path. The bundle identifies
the package sources and packaged data release by hash. It does not duplicate external
binaries or the complete package reference database into every TF directory.

## Portable paths

Manifest paths and decoded local report references use the same canonical relative
path rules. Paths use `/`, contain no empty, `.` or `..` components, and have no
repeated or trailing separator. Components may contain ordinary Unicode, but not
control characters, Windows-invalid characters, trailing dots/spaces, or Windows
device basenames such as `CON`, `AUX`, `CONIN$`, `CONOUT$`, `COM1`, or `LPT1`
(with or without an extension). The Python 3.10-compatible fallback also rejects
the Windows superscript-digit forms `COM¹`/`²`/`³` and `LPT¹`/`²`/`³`, and a
device stem followed by spaces before its extension. The runtime
`os.path.isreserved` check is used when available, with the deterministic fallback
retained in every supported Python version.

Each component is NFC-normalized and case-folded for collision detection, while its
original ordinary-Unicode spelling remains permitted. NFC/case-fold collisions,
linked roots/descendants, and filesystem escapes are rejected. Percent-encoded UTF-8
filenames are decoded before these checks; encoded separators and malformed escapes
are rejected.

The control manifest is the exact, case-sensitive root filename
`bundle_manifest.json`. A case variant such as `BUNDLE_MANIFEST.JSON`, or a second
name that NFC/case-folds to the control name, is rejected during the path inventory
before any JSON is parsed. The manifest is excluded from checksum inventory only
under that exact canonical spelling.

## HTML resources

The generated report format deliberately supports only its emitted semantic tag
set, with local `a[href]` and `img[src]` references. HTTP(S) and `mailto:` anchors,
fragment/query-only anchors, and embedded `data:image/...` images are also allowed.
Global `class`, `dir`, `hidden`, `id`, `lang`, `role`, `style`, `tabindex`, and
`title` attributes, `aria-*`/`data-*` attributes, and the small documented
tag-specific metadata/table attributes are harmless supported metadata. All other
tags or attributes fail closed; adding one to the generator requires an explicit
contract and test update. This rejects link preload `href`/`imagesrcset` resources,
inline SVG, legacy image `lowsrc`, and unknown future resource-bearing attributes
rather than trying to approximate every browser surface.

Script, iframe, embed, applet and object elements; event handler and `srcdoc`
attributes; `srcset`; document bases; and meta refresh are outside the format and
fail the audit. Duplicate attributes are rejected as ambiguous. Attribute values
containing CSS resource functions or resource at-rules are rejected even when the
attribute name is otherwise allowed; inline SVG presentation forms such as
`fill="url(...)"`, `filter="url(...)"`, and `cursor="url(...)"` therefore cannot
receive assurance.

This checker implements a generated-report subset, not browser-complete HTML or CSS
security semantics. CSS comments and backslash escapes are forbidden. The only CSS
at-rule accepted is `@media`; the only accepted functions are `blur`, `clamp`,
`minmax`, `nth-child`, `repeat`, `rgba`, and `var`, matching the current generated
stylesheet. Consequently external-resource and active functions such as `url`,
`image-set`, `src`, and `expression`, plus `@import` and `@font-face`, fail
conservatively. Expanding generated CSS requires an explicit contract and test update.

HTML comments, CDATA, processing instructions, and declarations other than the
optional simple HTML doctype are also outside the generated subset. A lexical check
runs before `HTMLParser`, so malformed browser-tokenization sequences such as
`<!--><script ...>-->` and `<![CDATA[><script ...>]]>` fail closed instead of being
treated as harmless parser data. The optional doctype is constrained to `html` with
only ASCII HTML whitespace and case variation; public/system identifiers and other
declaration forms require an explicit contract change.

## Manifest assurance

Version 1 or an omitted version is a legacy existence-level check: required files
must exist, but their bytes are not authenticated; an optional checksum-shaped map
does not silently upgrade v1 assurance. Version 2 requires a valid size
and lowercase SHA-256 record for every bundle file other than the manifest itself;
the inventory must cover exactly those files. Manifest shape, canonical names,
duplicate keys, and artifact records are validated before declared artifacts are
opened. `status` must be `complete`, and `report_generated` must be a Boolean (and
true when report inspection is required) before those reads. The Boolean is also
bound to the actual file inventory: true requires `REPORT.html`, including under
`--no-report`, while false forbids it. Any present report is always audited, so an
inconsistent false flag cannot bypass active-resource checks. Strict JSON rejects
`NaN`/infinities (including exponent overflow), invalid UTF-8, malformed JSON and
nesting beyond 64 container levels; decode and recursion failures are returned as
machine-readable issues. The full tree is checked for links before file hashing or
HTML reads.

`validate_files(...) -> list[str]` remains the compatibility API. Use
`audit_files(...)` when the assurance metadata is required. Assurance remains
`unavailable` unless the selected contract succeeds completely, including filesystem
inventory, every required size/hash operation and the requested HTML check. A
checksum mismatch, unreadable artifact or rejected HTML therefore never receives the
`v2-complete-checksums` label. `hashes_checked` counts successful hash operations;
`html_checks_performed` becomes true only after HTML decoding, parsing and constrained
checks complete (including a completed check that reports unsupported constructs).
An unsupported integer version is still returned in `manifest_version` to make the
incompatibility explicit.
