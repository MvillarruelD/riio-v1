"""Convert a saved DOCX with LibreOffice, or Word/docx2pdf on Windows/macOS."""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def export_pdf(docx: Path, *, required: bool = False) -> Path | None:
    docx = docx.resolve()
    target = docx.with_suffix('.pdf')
    if not docx.is_file():
        raise FileNotFoundError(docx)
    soffice = shutil.which('soffice') or shutil.which('libreoffice')
    # Convert into a fresh directory: a stale sibling PDF cannot impersonate success.
    with tempfile.TemporaryDirectory(prefix='riio-pdf-') as temporary:
        work = Path(temporary)
        converted = work / target.name
        if soffice:
            profile = work / 'profile'
            subprocess.run([soffice, '-env:UserInstallation=' + profile.as_uri(),
                            '--headless', '--convert-to', 'pdf', '--outdir', str(work),
                            str(docx)], check=True, timeout=180, capture_output=True, text=True)
        elif sys.platform in ('win32', 'darwin'):
            try:
                from docx2pdf import convert
            except ImportError:
                if required:
                    raise RuntimeError('PDF export requires LibreOffice or docx2pdf with Microsoft Word') from None
                return None
            convert(str(docx), str(converted))
        elif required:
            raise RuntimeError('PDF export requires LibreOffice (soffice on PATH)')
        else:
            return None
        if not converted.is_file() or not converted.read_bytes().startswith(b'%PDF-'):
            raise RuntimeError(f'PDF converter did not produce a valid PDF: {target}')
        # Move only after conversion validation; retain the old PDF on any earlier failure.
        shutil.copyfile(converted, target)

    return target


def main(argv=None) -> int:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('docx', type=Path)
    args = parser.parse_args(argv)
    print(export_pdf(args.docx, required=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
