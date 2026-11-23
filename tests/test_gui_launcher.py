"""UX contracts for the optional one-page interface.

These tests exercise the lightweight launcher rather than importing ``gui_app``: importing a Streamlit
script executes its page, which is not a stable unit-test boundary. They run unchanged on Windows and Linux.
"""
from __future__ import annotations

import socket
import sys

import pytest

from predictor import gui


@pytest.fixture
def launched(monkeypatch):
    """Record the Streamlit command instead of starting a server, and never open a real browser."""
    commands = []
    monkeypatch.setattr(gui.subprocess, "call", lambda argv: commands.append(argv) or 0)
    monkeypatch.setattr(gui, "_open_when_ready", lambda *a, **k: None)
    monkeypatch.setattr(gui, "free_port", lambda *a, **k: gui.DEFAULT_PORT)
    return commands


def test_launcher_applies_the_product_theme(launched):
    """An installed command must not depend on a checkout-local Streamlit config file."""
    assert gui.main([]) == 0
    command = launched[0]
    assert command[:4] == [sys.executable, "-m", "streamlit", "run"]
    assert "--server.headless=true" in command
    for argument in ("--theme.base=light", "--theme.primaryColor=#B4431F", "--theme.backgroundColor=#F4F1EA",
                     "--theme.secondaryBackgroundColor=#EAE6DB", "--theme.textColor=#1D1C19"):
        assert argument in command


def test_launcher_serves_this_machine_only_and_forwards_flags(launched):
    """No login protects the GUI, so it must not listen on the network unless explicitly asked."""
    assert gui.main([]) == 0
    assert "--server.address=localhost" in launched[0]
    assert "--browser.gatherUsageStats=false" in launched[0]
    assert gui.main(["--server.address=0.0.0.0", "--server.port=8600"]) == 0
    assert launched[1][-2:] == ["--server.address=0.0.0.0", "--server.port=8600"], \
        "user flags come last so they override the defaults"
    assert "--server.port=8501" not in launched[1], "an explicit port is not duplicated"


def test_port_alias_and_validation(launched, capsys):
    assert gui.main(["--port", "8600"]) == 0
    assert launched[0][-1] == "--server.port=8600"
    assert gui.main(["--port"]) == 2
    assert gui.main(["--server.port=abc"]) == 2
    assert len(launched) == 1, "bad ports never start a server"
    assert "port" in capsys.readouterr().out.lower()


def test_free_port_skips_a_busy_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as busy:
        busy.bind(("127.0.0.1", 0))
        taken = busy.getsockname()[1]
        found = gui.free_port(taken, span=20)
    assert found != taken and taken < found < taken + 20


def test_default_port_is_chosen_when_none_is_given(monkeypatch):
    commands = []
    monkeypatch.setattr(gui.subprocess, "call", lambda argv: commands.append(argv) or 0)
    monkeypatch.setattr(gui, "_open_when_ready", lambda *a, **k: None)
    monkeypatch.setattr(gui, "free_port", lambda *a, **k: 8517)
    assert gui.main([]) == 0
    assert "--server.port=8517" in commands[0]


def test_browser_is_opened_unless_suppressed_or_remote(monkeypatch):
    started = []

    class FakeThread:
        def __init__(self, target, args, daemon):
            started.append(args)

        def start(self):
            pass

    monkeypatch.setattr(gui.threading, "Thread", FakeThread)
    monkeypatch.setattr(gui.subprocess, "call", lambda argv: 0)
    monkeypatch.setattr(gui, "free_port", lambda *a, **k: 8501)
    monkeypatch.setattr(gui, "_can_open_browser", lambda: True)
    gui.main([])
    assert started and started[0][0] == "http://localhost:8501"
    started.clear()
    gui.main(["--no-browser"])
    assert not started
    monkeypatch.setattr(gui, "_can_open_browser", lambda: False)
    gui.main([])
    assert not started


def test_ssh_sessions_never_open_a_browser_and_get_forwarding_help(monkeypatch, capsys):
    monkeypatch.setenv("SSH_CONNECTION", "10.0.0.2 5000 10.0.0.9 22")
    assert gui._remote_session() and not gui._can_open_browser()
    monkeypatch.setattr(gui.subprocess, "call", lambda argv: 0)
    monkeypatch.setattr(gui, "free_port", lambda *a, **k: 8501)
    gui.main([])
    assert "ssh -L 8501:127.0.0.1:8501" in capsys.readouterr().out


def test_ctrl_c_stops_quietly(monkeypatch, capsys):
    def interrupted(argv):
        raise KeyboardInterrupt

    monkeypatch.setattr(gui.subprocess, "call", interrupted)
    monkeypatch.setattr(gui, "_open_when_ready", lambda *a, **k: None)
    monkeypatch.setattr(gui, "free_port", lambda *a, **k: 8501)
    assert gui.main([]) == 0
    assert "Stopped" in capsys.readouterr().out


def test_help_does_not_start_a_server(launched, capsys):
    assert gui.main(["--help"]) == 0
    assert not launched and "--no-browser" in capsys.readouterr().out


def test_tfop_gui_subcommand_forwards_its_arguments(monkeypatch):
    """`tfop gui` and `tfop-gui` are the same launcher on every platform."""
    from predictor import run_cli
    seen = []
    monkeypatch.setattr(gui, "main", lambda argv=None: seen.append(argv) or 0)
    assert run_cli.main(["gui", "--port", "8600", "--no-browser"]) == 0
    assert seen == [["--port", "8600", "--no-browser"]]


def test_source_checkout_uses_the_same_theme():
    """Direct ``streamlit run`` development sessions should match the installed launcher."""
    config = gui.Path(__file__).resolve().parents[1] / ".streamlit" / "config.toml"
    text = config.read_text(encoding="utf-8")
    for argument in gui.THEME_ARGS:
        key, value = argument.removeprefix("--theme.").split("=", 1)
        expected = f'{key} = "{value}"'
        assert expected in text
