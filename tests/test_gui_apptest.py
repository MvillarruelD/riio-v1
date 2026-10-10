"""Execute the optional GUI in Streamlit's actual rerun/widget test harness.

No real prediction or retrieval is invoked. CI installs the GUI extra on both OSes.
"""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest


@pytest.fixture
def app(monkeypatch, tmp_path):
    from predictor import resources
    # The GUI reads only manifest-based resource presence; keep test results isolated.
    monkeypatch.setattr(resources, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(resources, "audit", lambda: {"database": {"release": "test"}})
    return AppTest.from_file(str(Path(__file__).resolve().parents[1] / "predictor" / "gui_app.py"),
                             default_timeout=20).run()


def test_initial_page_and_all_entry_modes_rerun_without_exceptions(app):
    assert not app.exception
    assert [widget.label for widget in app.radio] == ["Starting point", "Discovery engine"]
    assert app.radio[0].value == "Genome accession"
    assert app.text_input[0].label == "Genome accession (required)"
    assert app.get("button_group")[0].value == "Genome scan"
    app.get("button_group")[0].set_value("Protein list").run()
    assert not app.exception
    assert "Read protein list" in [widget.label for widget in app.button]
    assert not app.radio
    assert app.text_input[0].label == "Genome accession (optional context)"
    app.get("button_group")[0].set_value("Single TF").run()
    assert not app.exception
    assert "Protein accession" in [widget.label for widget in app.text_input]
    app.get("button_group")[1].set_value("Paste sequence").run()
    assert not app.exception
    assert [widget.label for widget in app.text_area] == ["Protein sequence"]
    app.get("button_group")[0].set_value("Genome scan").run()
    assert not app.exception


def test_local_genome_options_and_invalid_submission_are_actionable(app):
    for mode in ("Local genome FASTA", "GenBank", "Proteome"):
        # Starting point updates guidance immediately, before form submission.
        app.radio[0].set_value(mode).run()
        assert app.text_input[0].label == "Genome accession (optional context)"
        next(widget for widget in app.button if widget.label == "Find candidate TFs").click().run()
        assert not app.exception
        assert app.radio[0].value == mode
        assert app.text_input[0].label == "Genome accession (optional context)"
        assert [error.value for error in app.error] == ["Choose a local input file"]


def test_genome_submission_uses_shared_api_and_keeps_selected_names(app, monkeypatch):
    import csv
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
        summary.write_text("status\nok\n", encoding="utf-8")
        return summary

    monkeypatch.setattr(api, "run_genome", discover)
    monkeypatch.setattr(genome_run, "run_candidates", predict)
    app.text_input[0].set_value("GCF_000005845.2")
    next(widget for widget in app.button if widget.label == "Find candidate TFs").click().run()
    assert not app.exception
    assert calls[0]["accession"] == "GCF_000005845.2"
    assert calls[0]["discover_only"] is True
    assert calls[0]["census_only"] is False
    assert app.multiselect[0].value == ["one", "two"]
    app.multiselect[0].set_value(["two"])
    next(widget for widget in app.button if widget.label == "Build selected reports").click().run()
    assert not app.exception
    assert calls[-1]["selected_names"] == ["two"]
    assert [row["run_name"] for row in calls[-1]["rows"]] == ["one", "two"]


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
    next(widget for widget in app.button if widget.label == "Prepare full run download").click().run()
    assert not app.exception and calls == [root]
    app.run()
    assert calls == [root]
    assert any(widget.label == "Download full run · ZIP" for widget in app.get("download_button"))
    next(widget for widget in app.toggle if widget.label == "Preview run report").set_value(True).run()
    assert not app.exception
    assert any("file links are disabled" in widget.value for widget in app.caption)
