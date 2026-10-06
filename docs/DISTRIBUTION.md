# Distribution contract

Build with `python -m build --wheel --sdist`. The wheel contains the `predictor` package, its
immutable reference assets, owned engine adapters, and upstream license notices. The source archive
also contains the analysis companion, its required candidate and scoring inputs, and release tests.

Run `python tools/release_check.py --selftest-timeout 120` before creating a version tag. This builds
both artifacts, installs the wheel non-editably into a new environment, runs the CLI outside the
checkout, checks packaged references/adapters, and verifies outputs are outside site-packages.

For the complete source suite, first install `.[dev,genome,gui]` and `analysis/requirements.txt`.
Engine weights and upstream executables remain separate prerequisites. An optional self-test skip
does not prove those engines work. Record skip/failure details and never describe an offline smoke
check as a canonical online prediction rerun.

Use the v1 GitHub release for downloads. Citation author metadata is provisional; replace it with
the agreed author list before final paper submission. Add a DOI only after a real archival deposit.

`requirements-tested-python312.txt` records the Python packages resolved in the local Windows
verification environment. It is an environment record, not a cross-platform lockfile or a replacement
for the pinned external engine instructions.
