# predictor -- TF -> operator/biosensor design pipeline (production package).
#
# Every intra-package import is absolute (`from predictor.annotate import context`). The package no
# longer puts its own directory on sys.path: doing so let one module load twice under two names
# (`annotate.context` and `predictor.annotate.context`, with separate caches) and exposed generic
# top-level names such as `schema` and `signals` that could shadow unrelated installed packages.

__version__ = "1.0.0"
