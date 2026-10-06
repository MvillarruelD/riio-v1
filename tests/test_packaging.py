"""Distribution smoke tests: does the package import, expose its entry points, and ship its data?

These are the checks that fail on a *user's* machine rather than a developer's. The import test in
particular catches the bare-sibling-import bug class (`from tf_record import ...` instead of
`from .tf_record import ...`), which works only when a specific subdirectory happens to be on sys.path
and therefore passes in-repo while breaking an installed package.
"""
import importlib
import json
import os
import subprocess
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 uses the dev-only compatibility dependency.
    import tomli as tomllib

import pytest

from predictor import resources

REPO = Path(__file__).resolve().parents[1]

#: Modules that legitimately need something not present in a bare checkout/CI: an external engine
#: binary, a network call at import, or an optional heavy dependency. Import is still attempted, but a
#: failure is reported as a skip rather than a hard failure.
OPTIONAL_AT_IMPORT: set[str] = {"predictor.gui_app"}  # imported only by the optional ``gui`` extra
RELEASED_DATA_AT_IMPORT: set[str] = {"predictor.gui_app"}  # its app startup audits packaged records


def _python_sources():
    package = REPO / "predictor"
    paths = []
    for root, dirs, files in os.walk(package):
        root = Path(root)
        dirs[:] = sorted(d for d in dirs if d != "__pycache__")
        if root == package and "data" in dirs:
            dirs.remove("data")  # package records are not Python import targets; never descend there
        for name in sorted(files):
            if name.endswith(".py"):
                paths.append(root / name)
    for p in sorted(paths):
        yield p


def _all_modules():
    for p in _python_sources():
        if p.name == "__init__.py":
            continue
        yield ".".join(p.relative_to(REPO).with_suffix("").parts)


@pytest.mark.parametrize("mod", [
    pytest.param(mod, marks=pytest.mark.released_data if mod in RELEASED_DATA_AT_IMPORT else ())
    for mod in _all_modules()
])
def test_every_module_imports(mod):
    """Every shipped module imports by its package-qualified name, with no sys.path preparation."""
    try:
        importlib.import_module(mod)
    except ImportError as e:
        if mod in OPTIONAL_AT_IMPORT:
            pytest.skip(f"optional dependency missing: {e}")
        raise


def test_no_developer_home_paths_in_shipped_code():
    """A default pointing at one developer's home directory silently misconfigures every other install."""
    sources = list(_python_sources())
    assert REPO / "predictor" / "__init__.py" in sources, "source-policy scan omitted package __init__.py"
    offenders = []
    for p in sources:
        for i, line in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if "/home/matia" in line or "C:/Users/matia" in line or r"C:\Users\matia" in line:
                offenders.append(f"{p.relative_to(REPO)}:{i}: {line.strip()[:90]}")
    assert not offenders, "hardcoded developer paths:\n" + "\n".join(offenders)


def test_no_committed_credentials():
    """.env must never be tracked; .env.example must carry only placeholders."""
    tracked = subprocess.run(["git", "ls-files", ".env"], cwd=REPO,
                             capture_output=True, text=True).stdout.strip()
    assert not tracked, ".env is tracked by git -- credentials must never be committed"
    example = (REPO / ".env.example").read_text(encoding="utf-8")
    for line in example.splitlines():
        if "=" in line and not line.strip().startswith("#"):
            _key, _, val = line.partition("=")
            assert ("your-" in val or "example" in val or not val.strip()), \
                f".env.example may contain a real secret: {line[:60]}"


class TestProjectMetadata:
    @staticmethod
    def _pyproject():
        return tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))

    def test_console_entry_points_resolve(self):
        """Every declared console-script target must exist, or the installed command is dead."""
        for name, target in self._pyproject()["project"]["scripts"].items():
            mod_name, _, func = target.partition(":")
            mod = importlib.import_module(mod_name)
            assert callable(getattr(mod, func, None)), f"entry point {name} -> {target} is not callable"

    def test_declares_homepage(self):
        urls = self._pyproject()["project"].get("urls") or {}
        assert urls.get("Homepage"), "no Homepage URL -- users cannot find the source from the package"

    def test_requires_python_matches_syntax_used(self):
        # the codebase uses `X | None` unions at runtime, which need 3.10+
        assert self._pyproject()["project"]["requires-python"] == ">=3.10"

    def test_python_310_dev_extra_includes_a_toml_reader(self):
        dev = self._pyproject()["project"]["optional-dependencies"]["dev"]
        assert "tomli>=2.0; python_version < '3.11'" in dev


@pytest.mark.released_data
class TestVendoredData:
    """The reference data the pipeline cannot run without. Missing data is the most common way an
    install looks fine and then fails on first use."""

    def test_pfam_hmm_library_present(self):
        hmm = resources.REFS_DIR / "tf_pfam.hmm"
        assert hmm.is_file() and hmm.stat().st_size > 1_000_000, "Pfam HMM library missing/truncated"

    def test_pfam_family_map_covers_every_ssn_family(self):
        """Every SSN family must be reachable from a Pfam accession, or its members never get classified."""
        import json
        from predictor.annotate import ssn_clusters as ssn
        mapped = set(json.loads((resources.REFS_DIR / "pfam_to_family.json")
                                .read_text(encoding="utf-8")).values())
        missing = [f for f in ssn.FAMILIES if f not in mapped]
        assert not missing, f"SSN families with no Pfam route: {missing}"

    def test_every_ssn_family_has_a_member_db(self):
        from predictor.annotate import ssn_clusters as ssn
        missing = [f for f in ssn.FAMILIES
                   if not (resources.SSN_DATABASE / f"ssn_members_{ssn.family_tag(f)}.fasta").is_file()]
        assert not missing, f"SSN families with no vendored member DB: {missing}"

    def test_data_manifest_covers_packaged_databases(self):
        import json

        manifest = json.loads((resources.DATA_DIR / "MANIFEST.json").read_text(encoding="utf-8"))
        paths = {row["path"] for row in manifest["files"]}
        # Database assets remain at their original release when only software packaging changes.
        from tools.build_data_manifest import DATABASE_RELEASE
        assert manifest["database_release"] == DATABASE_RELEASE
        assert "refs/tf_pfam.hmm" in paths
        assert "ssn/database/accession_to_cluster.json" in paths
        assert len([p for p in paths if p.startswith("ssn/clusters/") and p.endswith(".fasta")]) == 71

    def test_cluster_table_and_families_agree(self):
        from predictor.annotate import ssn_clusters as ssn
        assert len(ssn.FAMILIES) == 12, f"expected the 12-family SSN, got {len(ssn.FAMILIES)}"
        for fam in ssn.FAMILIES:
            assert ssn.load_clusters(fam), f"{fam} has no clusters"


# --------------------------------------------------------------------------- release coherence
class TestReleaseCoherence:
    """One release, one version, and a distribution that contains only what it should.

    These are the facts `tools/release_check.py` verifies end to end against a built wheel. They are
    repeated here as cheap unit checks so a mismatch fails in CI in seconds rather than in a six
    minute release run.
    """

    @staticmethod
    def _versions():
        import predictor
        from predictor.provenance import runtime_provenance
        pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
        return {
            "pyproject.version": pyproject["project"]["version"],
            "predictor.__version__": predictor.__version__,
            "provenance.package_version": str(runtime_provenance()["package_version"]),
        }

    @pytest.mark.released_data
    def test_every_declared_version_agrees(self):
        versions = self._versions()
        assert len(set(versions.values())) == 1, (
            "the release version disagrees across sources:\n  "
            + "\n  ".join(f"{k:<32} {v}" for k, v in versions.items())
        )

    @pytest.mark.released_data
    def test_provenance_reports_the_shipped_data_release(self):
        from predictor.provenance import runtime_provenance
        manifest = json.loads((REPO / "predictor" / "data" / "MANIFEST.json")
                              .read_text(encoding="utf-8"))
        assert str(runtime_provenance()["data_release"]) == str(manifest["database_release"])

    @pytest.mark.released_data
    def test_package_data_patterns_cover_every_shipped_file(self):
        """A data file the patterns miss is absent from the wheel and present in every dev checkout.

        That asymmetry is the whole difficulty of packaging data: it works locally by construction.
        """
        import fnmatch
        pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
        patterns = pyproject["tool"]["setuptools"]["package-data"]["predictor"]
        manifest = json.loads((REPO / "predictor" / "data" / "MANIFEST.json")
                              .read_text(encoding="utf-8"))
        uncovered = []
        for entry in manifest["files"]:
            rel = f"data/{entry['path']}"
            if rel.endswith(".py"):
                continue                      # shipped as package source, not as package data
            if not any(fnmatch.fnmatch(rel, pat) for pat in patterns):
                uncovered.append(rel)
        assert not uncovered, (
            f"{len(uncovered)} manifest file(s) match no package-data pattern, so a wheel omits "
            f"them:\n  " + "\n  ".join(uncovered[:15])
        )

    @pytest.mark.parametrize("name", ["research_archive", "tests", "docs", ".github"])
    def test_non_runtime_trees_are_outside_the_package(self, name):
        """Only `predictor*` is packaged, so anything outside it cannot reach a wheel."""
        assert not (REPO / "predictor" / name).exists(), \
            f"{name} is inside the package directory and would be distributed"

    def test_research_archive_is_never_imported(self):
        """Historical material must be unreachable from runtime, not merely unused."""
        offenders = []
        for path in _python_sources():
            if "research_archive" in path.read_text(encoding="utf-8", errors="replace"):
                offenders.append(str(path.relative_to(REPO)))
        assert not offenders, f"runtime code references research_archive: {offenders}"

    def test_packages_find_includes_only_the_predictor_package(self):
        pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
        assert pyproject["tool"]["setuptools"]["packages"]["find"]["include"] == ["predictor*"]

    def test_the_build_frontend_is_a_declared_dev_dependency(self):
        """docs/DISTRIBUTION.md and tools/release_check.py both run `python -m build`."""
        pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
        dev = pyproject["project"]["optional-dependencies"]["dev"]
        assert any(spec.startswith("build") for spec in dev), \
            "`build` is missing from the dev extra, so the documented release command fails"
