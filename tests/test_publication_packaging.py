"""Publication deployment contracts: adapters must survive a normal wheel install."""
from pathlib import Path

from predictor import resources
from predictor.discovery import bitacora
from predictor.structure import metalnet


def test_engine_runners_are_inside_the_installable_package():
    for path in (bitacora.RUNNER, metalnet.RUNNER, bitacora.PATCH_DIR):
        assert path.exists(), path
        assert resources.PACKAGE_DIR in path.parents
    assert (bitacora.PATCH_DIR / "get_blastp_parsed_newv2.pl").is_file()


def test_folding_adapter_is_available_without_a_source_checkout(monkeypatch):
    monkeypatch.delenv("PREDICTOR_ESMFOLD_LIB", raising=False)
    lib = resources.esmfold_lib_dir()
    assert lib is not None
    assert resources.PACKAGE_DIR in Path(lib).parents
    from predictor._engines.esmfold_lib._lib import esmfold2_inputs
    assert callable(esmfold2_inputs.build_apo_homomer_input)
