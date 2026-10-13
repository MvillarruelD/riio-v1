"""Launch the optional one-page Streamlit interface installed as ``tfop-gui``."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


THEME_ARGS = (
    "--theme.base=light",
    "--theme.primaryColor=#126652",
    "--theme.backgroundColor=#F5F7F4",
    "--theme.secondaryBackgroundColor=#FFFFFF",
    "--theme.textColor=#17231F",
    "--theme.font=sans-serif",
)


#: The GUI runs predictions and reads/writes local files without any login, so it serves this machine
#: only (WSL2 still forwards it to the Windows browser). Usage statistics are not sent to Streamlit.
SERVER_ARGS = (
    "--server.headless=true",
    "--server.address=localhost",
    "--browser.gatherUsageStats=false",
)


def main(argv=None) -> int:
    """Start the GUI; extra Streamlit flags (e.g. ``--server.port=8600``) are passed through and win."""
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print('The web interface is optional. Install it with: pip install "tf-operator-predictor[gui]"')
        return 1
    extra = sys.argv[1:] if argv is None else list(argv)
    app = Path(__file__).with_name("gui_app.py")
    print("Starting the RIIO GUI. Open the URL printed below in your web browser "
          "(normally http://localhost:8501); press Ctrl+C here to stop it.", flush=True)
    return subprocess.call([
        sys.executable, "-m", "streamlit", "run", str(app),
        *SERVER_ARGS, *THEME_ARGS, *extra,
    ])


if __name__ == "__main__":
    raise SystemExit(main())
