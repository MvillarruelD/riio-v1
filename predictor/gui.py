"""Launch the optional one-page Streamlit interface installed as ``tfop-gui`` (also ``tfop gui``).

One launcher for Windows, Linux (including WSL) and macOS:

* serves this machine only (the GUI has no login and writes local files);
* picks the first free port from 8501, unless ``--server.port`` / ``--port`` is given;
* opens the default browser once the server answers, except over SSH or with no display;
* exits quietly on Ctrl+C.
"""
from __future__ import annotations

import os
import socket
import subprocess
import sys
import threading
import time
import webbrowser
from pathlib import Path


# Keep these in step with ``.streamlit/config.toml`` (a test compares them) so a development
# ``streamlit run`` looks identical to the installed command.
THEME_ARGS = (
    "--theme.base=light",
    "--theme.primaryColor=#B4431F",
    "--theme.backgroundColor=#F4F1EA",
    "--theme.secondaryBackgroundColor=#EAE6DB",
    "--theme.textColor=#1D1C19",
    "--theme.font=sans-serif",
)

#: The GUI runs predictions and reads/writes local files without any login, so it serves this machine
#: only (WSL2 still forwards it to the Windows browser). Usage statistics are not sent to Streamlit.
SERVER_ARGS = (
    "--server.headless=true",
    "--server.address=localhost",
    "--browser.gatherUsageStats=false",
)

DEFAULT_PORT = 8501
USAGE = """usage: tfop-gui [--port N] [--no-browser] [streamlit flags]

  --port N        listen on this port (default: first free port from 8501)
  --no-browser    do not open a browser window
  any other flag is passed to `streamlit run` and overrides the defaults,
  e.g. --server.address=0.0.0.0 to serve other computers (no login protects the GUI)."""


def _port_free(port: int, host: str = "127.0.0.1") -> bool:
    """Nothing answers on either loopback and the port can be bound.

    Binding alone is not enough on Windows, where 127.0.0.1 can be bound under a listener on 0.0.0.0 or ::1.
    """
    for family, loopback in ((socket.AF_INET, host), (socket.AF_INET6, "::1")):
        try:
            with socket.socket(family, socket.SOCK_STREAM) as client:
                client.settimeout(.25)
                if client.connect_ex((loopback, port)) == 0:
                    return False
        except OSError:                     # no IPv6 on this machine
            pass
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        if os.name == "nt":
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        else:                               # a just-closed server in TIME_WAIT is not "busy"
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            probe.bind((host, port))
        except OSError:
            return False
    return True


def free_port(start: int = DEFAULT_PORT, span: int = 40) -> int:
    """First port from ``start`` that nothing is listening on; ``start`` itself if none is found."""
    return next((p for p in range(start, start + span) if _port_free(p)), start)


def _flag_value(args: list[str], name: str) -> str | None:
    for index, arg in enumerate(args):
        if arg.startswith(name + "="):
            return arg.split("=", 1)[1]
        if arg == name and index + 1 < len(args):
            return args[index + 1]
    return None


def _remote_session() -> bool:
    """True over SSH, where a browser opened on the server would be invisible."""
    return bool(os.environ.get("SSH_CONNECTION") or os.environ.get("SSH_TTY"))


def _can_open_browser() -> bool:
    if _remote_session():
        return False
    if sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")
                                                 or os.environ.get("WSL_DISTRO_NAME")):
        return False
    return True


def _open_when_ready(url: str, host: str, port: int, timeout: float = 90.0) -> None:
    """Open ``url`` as soon as the server accepts connections (daemon thread; failure is silent)."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with socket.create_connection((host, port), timeout=1):
                break
        except OSError:
            time.sleep(.4)
    else:
        return
    try:
        opened = webbrowser.open(url)
        if not opened and os.environ.get("WSL_DISTRO_NAME"):
            subprocess.run(["explorer.exe", url], check=False, timeout=10)
    except Exception:  # a missing browser must never take the server down
        pass


def main(argv=None) -> int:
    """Start the GUI; extra Streamlit flags (e.g. ``--server.port=8600``) are passed through and win."""
    extra = sys.argv[1:] if argv is None else list(argv)
    if any(arg in ("-h", "--help") for arg in extra):
        print(USAGE)
        return 0
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print('The web interface is optional. Install it with: python -m pip install "tf-operator-predictor[gui]"')
        return 1

    open_browser = "--no-browser" not in extra
    extra = [arg for arg in extra if arg != "--no-browser"]
    extra = [f"--server.port={arg[len('--port='):]}" if arg.startswith("--port=") else arg for arg in extra]
    if "--port" in extra:  # friendlier alias for --server.port
        index = extra.index("--port")
        if index + 1 >= len(extra):
            print("--port needs a number, e.g. --port 8600")
            return 2
        extra[index:index + 2] = [f"--server.port={extra[index + 1]}"]
    explicit = _flag_value(extra, "--server.port")
    try:
        port = int(explicit) if explicit is not None else free_port()
        if not 0 < port < 65536:
            raise ValueError(port)
    except ValueError:
        print(f"Not a port number: {explicit!r}")
        return 2
    host = _flag_value(extra, "--server.address") or "localhost"
    shown = "localhost" if host in ("0.0.0.0", "::", "") else host
    url = f"http://{shown}:{port}"

    app = Path(__file__).with_name("gui_app.py")
    print(f"RIIO is starting at {url}", flush=True)
    if explicit is None and port != DEFAULT_PORT:
        print(f"(port {DEFAULT_PORT} is busy, so {port} is used)")
    if _remote_session():
        print(f"Over SSH, forward the port first:  ssh -L {port}:127.0.0.1:{port} <user>@<host>  "
              f"then open {url} on your own computer.")
    print("Press Ctrl+C in this window to stop it.", flush=True)
    if open_browser and _can_open_browser():
        threading.Thread(target=_open_when_ready,
                         args=(url, "127.0.0.1" if shown == "localhost" else shown, port), daemon=True).start()
    command = [sys.executable, "-m", "streamlit", "run", str(app), *SERVER_ARGS, *THEME_ARGS]
    if not explicit:
        command.append(f"--server.port={port}")
    command.extend(extra)
    try:
        return subprocess.call(command)
    except KeyboardInterrupt:
        print("\nStopped.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
