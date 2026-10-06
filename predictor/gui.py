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


def main() -> int:
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print('The web interface is optional. Install it with: pip install "tf-operator-predictor[gui]"')
        return 1
    app = Path(__file__).with_name("gui_app.py")
    return subprocess.call([
        sys.executable, "-m", "streamlit", "run", str(app),
        "--server.headless=true", *THEME_ARGS,
    ])


if __name__ == "__main__":
    raise SystemExit(main())
