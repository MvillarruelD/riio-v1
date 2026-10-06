# Contributing

Open an issue with the software version, operating system, command, and a minimal reproducible input.
Remove credentials and private sequences from logs before posting them. For pull requests, explain
the resulting behavior and run the source and released-data test suites in the README.

Keep reference data immutable at runtime; caches and output belong outside the installed package.
Preserve third-party attribution. Changes to prediction behavior should be explicit and reviewed by
the project authors. Database recuration and manuscript authorship updates require their own review.

Run `python tools/release_check.py --selftest-timeout 120` for packaging changes. It tests a wheel
installed outside the checkout. The research workspace and generated predictions are not part of
this repository. Submit new source files only when the supported pipeline or its verification uses them.
