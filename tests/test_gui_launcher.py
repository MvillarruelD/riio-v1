"""UX contracts for the optional one-page interface.

These tests intentionally exercise the lightweight launcher rather than importing ``gui_app``: importing
a Streamlit script executes its page, which is not a stable unit-test boundary.
"""
from __future__ import annotations

import sys

from predictor import gui


def test_launcher_applies_the_product_theme(monkeypatch):
    """An installed command must not depend on a checkout-local Streamlit config file."""
    command = []

    def record(argv):
        command.extend(argv)
        return 0

    monkeypatch.setattr(gui.subprocess, "call", record)

    assert gui.main() == 0
    assert command[:4] == [sys.executable, "-m", "streamlit", "run"]
    assert "--server.headless=true" in command
    assert "--theme.base=light" in command
    assert "--theme.primaryColor=#126652" in command
    assert "--theme.backgroundColor=#F5F7F4" in command
    assert "--theme.secondaryBackgroundColor=#FFFFFF" in command
    assert "--theme.textColor=#17231F" in command


def test_source_checkout_uses_the_same_theme():
    """Direct ``streamlit run`` development sessions should match the installed launcher."""
    config = gui.Path(__file__).resolve().parents[1] / ".streamlit" / "config.toml"
    text = config.read_text(encoding="utf-8")
    for argument in gui.THEME_ARGS:
        key, value = argument.removeprefix("--theme.").split("=", 1)
        expected = f'{key} = "{value}"'
        assert expected in text
