"""Execute the optional GUI in Streamlit's actual rerun/widget test harness.

No real prediction or retrieval is invoked. CI installs the GUI extra on both OSes.
"""
import csv
import json
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parents[1] / "predictor" / "gui_app.py"


@pytest.fixture
def out_root(monkeypatch, tmp_path):
    from predictor import resources
    root = tmp_path / "results"
    root.mkdir()
    # The GUI reads only manifest-based resource presence; keep test results isolated.
    monkeypatch.setattr(resources, "OUTPUT_DIR", root)
    monkeypatch.setattr(resources, "audit", lambda: {"database": {"release": "test"}})
    # Production discovery is probed through WSL/BLAST/HMMER; tests must not depend on the host having them.
    from predictor.discovery import bitacora
    monkeypatch.setattr(bitacora, "engine_available", lambda: (True, ""))
    import streamlit as st
    st.cache_data.clear()
    return root


@pytest.fixture
def app(out_root):
    return AppTest.from_file(str(APP), default_timeout=30).run()


def _button(app, label):
    return next(widget for widget in app.button if widget.label == label)


def _markdown(app) -> str:
    return "\n".join(element.value for element in app.markdown)


def _make_run(root: Path, name="saved_run", *, predicted=True) -> Path:
    """A finished two-candidate run: one metal call and one saved final null (unresolved)."""
    run = root / name
    (run / "jobs").mkdir(parents=True)
    rows = [{"run_name": "G__Fur__NP_1", "family": "Fur", "gene": "fur", "locus_tag": "b0683", "aa_len": "148"},
            {"run_name": "G__MerR__NP_2", "family": "MerR", "gene": "", "locus_tag": "b0001", "aa_len": "130"}]
    with (run / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if predicted:
        (run / "jobs" / "G__Fur__NP_1").mkdir()
        (run / "jobs" / "G__Fur__NP_1" / "dossier.json").write_text("{}", encoding="utf-8")
        (run / "summary.tsv").write_text(
            "run_name\tfamily\tinducer_top\tinducer_class\tinducer_class_provenance\tagreement\tprimary_operator\t"
            "n_regulon_operons\tstatus\n"
            "G__Fur__NP_1\tFur\tFe2+\tmetal\texplicit\t1.0\tTTACTTTGTTATAAGACGTGA\t25\tok\n"
            "G__MerR__NP_2\tMerR\t\t\texplicit\t0.0\t\t25\tok\n", encoding="utf-8")
    return run


def test_initial_page_and_all_entry_modes_rerun_without_exceptions(app):
    assert not app.exception
    assert [widget.label for widget in app.radio] == ["Starting point", "Discovery engine"]
    assert app.radio[0].value == "Genome accession"
    assert app.radio[1].value.startswith("Production")
    assert app.text_input[0].label == "Genome accession (required)"
    assert app.get("button_group")[0].value == "Genome scan"
    app.get("button_group")[0].set_value("Protein list").run()
    assert not app.exception
    assert not app.radio
    assert app.text_input[0].label == "Genome accession (optional context)"
    assert [widget.label for widget in app.text_area] == ["Or paste accessions / FASTA"]
    app.get("button_group")[0].set_value("Single TF").run()
    assert not app.exception
    assert "Protein accession" in [widget.label for widget in app.text_input]
    app.get("button_group")[1].set_value("Paste sequence").run()
    assert not app.exception
    assert [widget.label for widget in app.text_area] == ["Protein sequence"]
    app.get("button_group")[0].set_value("Genome scan").run()
    assert not app.exception


def test_missing_discovery_engine_falls_back_to_the_quick_census_with_an_explanation(out_root, monkeypatch):
    from predictor.discovery import bitacora
    import streamlit as st
    monkeypatch.setattr(bitacora, "engine_available", lambda: (False, "BITACORA is not installed"))
    st.cache_data.clear()
    app = AppTest.from_file(str(APP), default_timeout=30).run()
    assert not app.exception
    assert app.radio[1].value.startswith("Quick census")
    assert any("BITACORA" in warning.value for warning in app.warning)


def test_local_genome_options_and_invalid_submission_are_actionable(app, out_root):
    for mode in ("Local genome FASTA", "GenBank", "Proteome"):
        # Starting point updates guidance immediately, before form submission.
        app.radio[0].set_value(mode).run()
        assert app.text_input[0].label == "Genome accession (optional context)"
        _button(app, "Find candidate TFs").click().run()
        assert not app.exception
        assert app.radio[0].value == mode
        assert app.text_input[0].label == "Genome accession (optional context)"
        assert [error.value for error in app.error] == ["Choose a local input file"]
    assert list(out_root.iterdir()) == [], "a rejected form must not leave an empty run folder behind"


def test_an_empty_accession_is_rejected_before_any_folder_is_made(app, out_root):
    _button(app, "Find candidate TFs").click().run()
    assert [error.value for error in app.error] == ["Enter a genome accession, for example GCF_000005845.2."]
    assert list(out_root.iterdir()) == []


def test_pasted_protein_list_becomes_an_accession_table(app, out_root, monkeypatch):
    from predictor import api
    calls = []
    monkeypatch.setattr(api, "predict_many", lambda **kwargs: calls.append(kwargs))
    app.get("button_group")[0].set_value("Protein list").run()
    app.text_area[0].set_value("WP_000240047.1\nNP_417751.1  # ZntR")
    app.text_input[0].set_value("GCF_000005845.2")
    _button(app, "Read protein list").click().run()
    # the mocked runner writes nothing; the page must still report cleanly rather than crash
    assert calls and calls[0]["discover_only"] is True
    assert calls[0]["organism"] == "GCF_000005845.2"
    assert calls[0]["source"].read_text(encoding="utf-8") == "accession\nWP_000240047.1\nNP_417751.1\n"


def test_genome_submission_uses_shared_api_and_keeps_selected_names(app, monkeypatch):
    from predictor import api, genome_run
    calls = []

    def discover(**kwargs):
        calls.append(kwargs)
        out = kwargs["out_dir"]
        with (out / "manifest.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, ["run_name", "family", "protein_id"])
            writer.writeheader()
            writer.writerows([{"run_name": "one", "family": "MerR", "protein_id": "one"},
                              {"run_name": "two", "family": "Fur", "protein_id": "two"}])

    def predict(rows, out, **kwargs):
        calls.append({"rows": rows, "out": out, **kwargs})
        summary = out / "summary.tsv"
        summary.write_text("run_name\tstatus\ntwo\tok\n", encoding="utf-8")
        return summary

    monkeypatch.setattr(api, "run_genome", discover)
    monkeypatch.setattr(genome_run, "run_candidates", predict)
    app.text_input[0].set_value("GCF_000005845.2")
    _button(app, "Find candidate TFs").click().run()
    assert not app.exception
    assert calls[0]["accession"] == "GCF_000005845.2"
    assert calls[0]["discover_only"] is True
    assert calls[0]["census_only"] is False
    # the family filter comes first; the per-regulator selection sits behind it
    assert app.multiselect[0].label == "Families" and app.multiselect[0].value == ["Fur", "MerR"]
    assert app.multiselect[1].value == ["one", "two"]
    app.multiselect[1].set_value(["two"])
    _button(app, "Build selected reports").click().run()
    assert not app.exception
    assert calls[-1]["selected_names"] == ["two"]
    assert [row["run_name"] for row in calls[-1]["rows"]] == ["one", "two"]
    # options that used to be dropped on this path now reach the shared runner
    assert calls[-1]["metalnet"] is True and calls[-1]["af3"] is True and calls[-1]["fold"] is False


def test_family_filter_narrows_the_selection(app, out_root):
    app.session_state["last_shared_genome"] = str(_make_run(out_root, predicted=False))
    app.run()
    assert not app.exception
    assert app.multiselect[1].value == ["G__Fur__NP_1", "G__MerR__NP_2"]
    app.multiselect[0].set_value(["Fur"]).run()
    assert app.multiselect[1].value == ["G__Fur__NP_1"]
    assert "1 of 2" in _markdown(app)


def test_results_show_every_inducer_class_and_never_call_a_null_organic(app, out_root):
    run = _make_run(out_root)
    app.session_state["last_shared_genome"] = str(run)
    app.run()
    assert not app.exception
    page = _markdown(app)
    assert 'class="k k-metal"' in page and 'class="k k-unresolved"' in page
    assert 'class="k k-organic"' not in page, "a saved final null must never be drawn as organic"
    assert "TTACTTTGTTATAAGACGTGA" in page
    assert "Read one report" in [h.value for h in app.subheader]


def test_saved_runs_can_be_reopened_by_path_only_inside_allowed_folders(app, out_root, tmp_path):
    inside = _make_run(out_root, "earlier")
    outside = _make_run(tmp_path / "elsewhere", "foreign")
    key = "last_shared_genome"
    app.text_input(key=f"saved_path_{key}").set_value(str(outside)).run()  # the button appears once text exists
    app.button(key=f"saved_open_path_{key}").click().run()
    assert any("only folders inside" in error.value for error in app.error)
    assert key not in app.session_state or not app.session_state[key]
    app.text_input(key=f"saved_path_{key}").set_value(f'"{inside}"').run()  # Windows "Copy as path" adds quotes
    app.button(key=f"saved_open_path_{key}").click().run()
    assert not app.exception
    assert Path(app.session_state[key]).resolve() == inside.resolve()


def test_saved_run_archive_is_prepared_only_on_explicit_action(app, monkeypatch, tmp_path):
    from predictor import gui_support
    root = tmp_path / "saved_run"
    root.mkdir()
    (root / "manifest.csv").write_text("run_name,family,protein_id\nTF,MerR,TF\n")
    (root / "RUN_REPORT.html").write_text('<a href="jobs/TF/REPORT.html">TF</a>')
    calls = []
    monkeypatch.setattr(gui_support, "context_archive", lambda path: (calls.append(path) or b"ZIP", []))
    app.session_state["last_shared_genome"] = str(root)
    app.run()
    assert not app.exception and not calls
    _button(app, "Prepare full run download").click().run()
    assert not app.exception and calls == [root]
    app.run()
    assert calls == [root]
    assert any(widget.label == "Download full run · ZIP" for widget in app.get("download_button"))
    next(widget for widget in app.toggle if widget.label == "Preview run report").set_value(True).run()
    assert not app.exception
    assert any("file links are disabled" in widget.value for widget in app.caption)


def test_single_tf_result_reads_the_stored_class_and_agreement(app, out_root):
    from predictor import api, gui_support
    bundle = out_root / "jobs" / "zntR"
    bundle.mkdir(parents=True)
    (bundle / "dossier.json").write_text(json.dumps({
        "family": "MerR", "family_support": "family_validated", "genome_accession": "NC_000913.3",
        "inducers": {"top": "Zn2+", "agreement": 1.0, "sensor_class": "metal", "finalized": True},
        "operators": {"primary": {"sequence": "TGGTCGCAGATTCAAATACT", "provenance": {"locality": "proximal"}}},
        "regulon": [], "rescan": {"hits": []}}), encoding="utf-8")
    app.get("button_group")[0].set_value("Single TF").run()
    app.session_state["last_result"] = api.Result.from_run(gui_support.result_dict(bundle))
    app.run()
    assert not app.exception
    page = _markdown(app)
    assert "Zn2+" in page and 'class="k k-metal"' in page and "exact agreement 1.00" in page
    assert any("REPORT.html is missing" in error.value for error in app.error)


def test_a_run_folder_shows_the_terminal_command_that_continues_it(app, out_root):
    run = _make_run(out_root, "ecoli")
    (run / "run_input.json").write_text(json.dumps({"mode": "genome_accession", "accession": "GCF_000005845.2"}))
    app.session_state["last_shared_genome"] = str(run)
    app.run()
    assert not app.exception
    commands = [code.value for code in app.code]
    assert any(c.startswith("tfop genome GCF_000005845.2 --out ") for c in commands)
