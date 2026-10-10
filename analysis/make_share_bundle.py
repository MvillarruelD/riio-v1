#!/usr/bin/env python
"""Publish a share directory from staging; never merge old and new deliverables.

Scientific interpretation belongs in the reports, not this packaging script.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import html
import os
from urllib.parse import quote, unquote, urlsplit
from html.parser import HTMLParser
import shutil
import sys
import tempfile
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
import uuid
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import manifest as MF
if os.environ.get('TFOP_RUN_TAG') or 'run_paths' in sys.modules:
    import run_paths as RP
else:
    RP = None
from project_config import PREDICTOR_ROOT

if str(PREDICTOR_ROOT) not in sys.path:
    sys.path.insert(0, str(PREDICTOR_ROOT))
from predictor.report.bundle_contract import contained_file, inventory, validate_files

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
JOBS = RP.JOBS if RP is not None else None
OUT = PROJ / f"SHARE_FOR_AUTHORS_{RP.TAG}" if RP is not None else None


def select_run() -> None:
    """Resolve required single-run context only for single-run publication."""
    global RP, JOBS, OUT
    if RP is None:
        import run_paths
        RP = run_paths
        JOBS = RP.JOBS
        OUT = PROJ / f'SHARE_FOR_AUTHORS_{RP.TAG}'


def _existing(*paths: Path) -> Path:
    for path in paths:
        if path.is_file():
            return path
    raise FileNotFoundError("missing run export: " + ", ".join(str(path) for path in paths))


def _manifest_path() -> Path:
    return RP.MANIFEST


def _rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("candidate manifest is empty")
    names = set()
    for row in rows:
        name = row.get("run_name", "")
        issue = MF.unsafe_run_name(name)
        if issue or name.casefold() in names:
            raise ValueError(f"invalid or duplicate run_name {name!r}: {issue or 'duplicate'}")
        if not row.get("fasta"):
            raise ValueError(f"missing candidate input path: {name}")
        names.add(name.casefold())
    return rows


def _coverage_note(n: int) -> str:
    canonical = RP.ANALYSIS / f"run_manifest_{RP.TAG}.csv"
    if not canonical.is_file() or canonical.resolve() == _manifest_path().resolve():
        return ""
    canonical_names = {row["run_name"].casefold() for row in _rows(canonical)}
    selected = {row["run_name"].casefold() for row in _rows(_manifest_path())}
    missing, extra = canonical_names - selected, selected - canonical_names
    if not missing and not extra:
        return ""
    return (f"This package differs from the canonical manifest: {n} selected; "
            f"{len(missing)} canonical candidates absent; {len(extra)} additional candidates. "
            "It must not be described as the complete canonical run.\n")


def build_readme(*, no_bundles: bool = False) -> str:
    rows = _rows(_manifest_path())
    predictions = list(csv.DictReader(
        _existing(RP.COMPLEMENT_TSV, RP.RESULTS / "per_regulator_predictions.tsv")
        .open(encoding="utf-8"), delimiter="\t"))
    n_metal = sum(row.get("inducer_class") == "metal" for row in predictions)
    n_protein = sum(
        row.get("inducer_class") == "metal"
        and bool((row.get("ssn_ligand") or "").strip()
                 or (row.get("gate") or "").strip() == "True")
        for row in predictions
    )
    scope = ("Individual bundles were omitted by --no-bundles. This is a summary export."
             if no_bundles else "05_per_regulator/ contains the selected individual report directories.")
    return (
        f"# Metalloregulator prediction package\n\nPrepared {dt.date.today().isoformat()} for run {RP.TAG}. "
        "Private working material for the manuscript authors; do not redistribute without agreement.\n\n"
        f"The selected manifest contains {len(rows)} candidates. The published survey contains 150; "
        "the cohorts are close but not identical. "
        f"The comparable protein-level subset calls {n_protein}/{len(predictions)} "
        f"({100*n_protein/len(predictions):.0f}%) metal-responsive. The broader fused result is "
        f"{n_metal}/{len(predictions)} ({100*n_metal/len(predictions):.0f}%) and reuses MetalNet and "
        "genomic-context evidence, so it is not an independent replication of the survey result. "
        f"{_coverage_note(len(rows))}\n\n"
        "## Files\n\n"
        "- 01_reports/: main Word report and concise Markdown summary.\n"
        "- 02_figures/: available SVG, PDF, PNG and TIFF exports.\n"
        "- 03_data/: result tables, candidate manifest and locally copied candidate inputs.\n"
        "- 04_manuscript_integration/: available composition exports.\n"
        f"- {scope}\n"
        "- 06_benchmark_context/: historical benchmark panels, when present, kept separately.\n\n"
        "## Reuse\n\n"
        "Keep this directory intact. Individual REPORT.html files open offline. The manifest in "
        "03_data/ uses paths relative to 03_data/; the original manifest is retained alongside it. "
        "share_manifest.json records file sizes and SHA-256 digests, excluding its own digest.\n\n"
        "Predicted operators and regulons are ranked hypotheses for experimental prioritisation, "
        "not validated sites. These are recorded deliverables, not a complete offline replay environment. Read the "
        "individual bundle guides for input and provenance limitations. Packaging does not "
        "validate predictions or establish computational reproducibility.\n\n"
        "Benchmark results are maintained separately under analysis/benchmarking/. No benchmark "
        "accuracy figures are inserted into this guide by the packaging code. Historical reports "
        "may still contain benchmark discussion.\n"
    )


def build_summary() -> str:
    """A short, cautious index for readers who do not start with the Word report."""
    rows = list(csv.DictReader(
        _existing(RP.COMPLEMENT_TSV, RP.RESULTS / "per_regulator_predictions.tsv")
        .open(encoding="utf-8"), delimiter="\t"))
    n = len(rows)
    metal = sum(row.get("inducer_class") == "metal" for row in rows)
    protein = sum(
        row.get("inducer_class") == "metal"
        and bool((row.get("ssn_ligand") or "").strip()
                 or (row.get("gate") or "").strip() == "True")
        for row in rows
    )
    return (
        "# Report summary\n\n"
        f"Run **{RP.TAG}** contains **{n}** candidate regulators from four bacterial genomes and "
        f"{n}/{n} recorded per-regulator bundles. The published SSN survey contains 150 candidates; "
        "the two cohorts are not identical.\n\n"
        f"The pipeline calls **{protein}/{n} ({100*protein/n:.0f}%)** candidates metal-responsive "
        "when restricted to calls supported by the SSN clade or coordination gate. This is the "
        "closest comparison with the survey's **41/150 (27%)** MetalNet-site result, but the differing "
        "cohorts make it convergence rather than validation. The full fusion calls "
        f"**{metal}/{n} ({100*metal/n:.0f}%)** metal-responsive; it also uses MetalNet and genomic "
        "context and is not an independent comparison.\n\n"
        "Each regulator bundle records the fused inducer hypothesis, ranked operator candidates, a "
        "candidate regulon and source-level evidence. Operator and regulon results are prioritisation "
        "tools, not experimentally validated calls. Read the Word report for methods, evidence tiers "
        "and limitations before reusing a result.\n"
    )


def copy(src: Path, dst: Path, label: str) -> int:
    """Copy one required deliverable; a missing input fails publication."""
    if not src.is_file():
        raise FileNotFoundError(f"missing {label}: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return 1


def assemble(destination: Path, *, no_bundles: bool = False) -> None:
    """Populate a new directory from the selected run."""
    select_run()
    rows = _rows(_manifest_path())
    copy(RP.REPORT_DOCX, destination / "01_reports" / RP.REPORT_DOCX.name, "report")
    pdf = RP.REPORT_DOCX.with_suffix('.pdf')
    if pdf.is_file():
        copy(pdf, destination / '01_reports' / pdf.name, 'PDF report')
    (destination / "01_reports" / "report_summary.md").write_text(build_summary(), encoding="utf-8")
    for source, name in ((_existing(RP.COMPLEMENT_TSV, RP.RESULTS / "per_regulator_predictions.tsv"),
                          "per_regulator_predictions.tsv"),
                         (_existing(RP.REGULON_SUMMARY, RP.RESULTS / "regulon_summary.tsv"),
                          "regulon_summary.tsv")):
        copy(source, destination / "03_data" / name, name)
    ledger = RP.RESULTS / "evidence_ledger.tsv"
    if ledger.is_file():
        copy(ledger, destination / "03_data" / ledger.name, "evidence ledger")
    copy(_manifest_path(), destination / "03_data" / "candidate_manifest.original.csv", "manifest")
    copied_rows = []
    for row in rows:
        source = MF.resolve_fasta(row["fasta"], PROJ)
        relative = f"inputs/{row['run_name']}.fasta"
        copy(source, destination / "03_data" / relative, "candidate input")
        copied_rows.append({**row, "fasta": relative})
        if not no_bundles:
            bundle = contained_file(JOBS, row["run_name"])
            issues = validate_files(bundle)
            if issues:
                raise ValueError(f"cannot share {row['run_name']}: {'; '.join(issues)}")
            target = destination / "05_per_regulator" / row["run_name"]
            shutil.copytree(bundle, target)
            issues = validate_files(target)
            if issues:
                raise ValueError(f"copied bundle failed verification: {'; '.join(issues)}")
    with (destination / "03_data" / "candidate_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(copied_rows)
    figures = sorted(p for p in RP.FIGURES.iterdir() if p.suffix.lower() in {".png", ".svg", ".pdf", ".tiff"})
    if not figures:
        raise ValueError("no figure exports are available")
    for figure in figures:
        folder = "06_benchmark_context" if figure.name.startswith("F21_") else "02_figures"
        copy(figure, destination / folder / figure.name, "figure")
        if figure.name.startswith(("F23_", "F24_")):
            copy(figure, destination / "04_manuscript_integration" / figure.name, "composition")
    (destination / "README.md").write_text(build_readme(no_bundles=no_bundles), encoding="utf-8")
    write_index(destination, rows, no_bundles=no_bundles)
    (destination / "share_manifest.json").write_text(json.dumps({
        "schema_version": 1, "run_tag": RP.TAG, "candidate_count": len(rows),
        "individual_bundles_included": not no_bundles, "artifacts": inventory(destination),
    }, indent=2), encoding="utf-8")


def write_index(destination: Path, rows: list[dict], *, no_bundles: bool = False) -> None:
    """Generate escaped, relative links only to published artifacts."""
    links = []
    for folder in ('01_reports', '03_data'):
        for path in sorted((destination / folder).glob('*')):
            if path.is_file():
                rel = path.relative_to(destination).as_posix()
                links.append(f'<li><a href="{quote(rel)}">{html.escape(path.name)}</a></li>')
    body = []
    if not no_bundles:
        for row in rows:
            name = row['run_name']
            rel = f'05_per_regulator/{name}/REPORT.html'
            if not (destination / rel).is_file():
                raise FileNotFoundError(f'missing published report: {rel}')
            body.append('<tr><td><a href="' + quote(rel) + '">' + html.escape(name) +
                        '</a></td><td>' + html.escape(row.get('family', '')) + '</td></tr>')
    index = ('<!doctype html><html lang="en"><meta charset="utf-8"><title>RIIO run reports</title>'
             '<style>body{font:16px Arial,sans-serif;max-width:1100px;margin:40px auto;padding:0 20px}'
             'td,th{padding:8px;text-align:left;border-bottom:1px solid #ddd}table{border-collapse:collapse}</style>'
             '<h1>RIIO run reports</h1><p>Recorded computational predictions for experimental prioritisation.</p>'
             '<ul>' + ''.join(links) + '</ul><table><thead><tr><th>Regulator report</th><th>Family</th></tr></thead>'
             '<tbody>' + ''.join(body) + '</tbody></table></html>')
    (destination / 'index.html').write_text(index, encoding='utf-8')


def publish(*, no_bundles: bool = False) -> Path:
    """Publish after assembly succeeds; preserve any prior release as a sibling backup."""
    select_run()
    final = OUT.resolve()
    if final.parent != PROJ.resolve() or not final.name.startswith("SHARE_FOR_AUTHORS_"):
        raise ValueError("share output must be a direct child of the project root")
    staging = Path(tempfile.mkdtemp(prefix=".share-", dir=final.parent))
    backup = None
    try:
        assemble(staging, no_bundles=no_bundles)
        if final.exists():
            backup = final.with_name(f".{final.name}.previous-{uuid.uuid4().hex[:10]}")
            final.replace(backup)
        try:
            staging.replace(final)
        except BaseException:
            if backup is not None and backup.exists() and not final.exists():
                backup.replace(final)
            raise
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(f"Shared deliverables: {final}")
    if backup is not None:
        print(f"Previous release preserved: {backup}")
    return final


def relocate_run_links(run_dir: Path, key: str) -> None:
    """Update local report links after bundles move out of a run's jobs directory."""
    import openpyxl
    for report in run_dir.glob("*.html"):
        report.write_text(report.read_text(encoding="utf-8").replace(
            'jobs/', '../../bundles/' + key + '/'), encoding="utf-8")
    for path in run_dir.glob("*.xlsx"):
        workbook = openpyxl.load_workbook(path)
        for sheet in workbook:
            for row in sheet:
                for cell in row:
                    link = cell.hyperlink
                    if link and link.target and link.target.startswith("jobs/"):
                        link.target = "../../bundles/" + key + "/" + link.target[5:]
        workbook.save(path)
        workbook.close()


def relocate_discovery_links(discovery_dir: Path, key: str) -> None:
    """The discovery page links predicted candidates as ../jobs/<TF>/; in the share they sit in bundles/."""
    report = discovery_dir / "DISCOVERY_REPORT.html"
    if report.is_file():
        report.write_text(report.read_text(encoding="utf-8").replace(
            '../jobs/', '../../../bundles/' + key + '/'), encoding="utf-8")


def validate_share_links(root: Path) -> list[str]:
    """Every HTML href/src and workbook hyperlink must resolve inside this package."""
    root = root.resolve()
    issues = []
    def check(path, target):
        parsed = urlsplit(target)
        if parsed.scheme in {"http", "https", "mailto", "data"} or not parsed.path:
            return
        resolved = (path.parent / unquote(parsed.path)).resolve()
        if parsed.scheme or parsed.netloc or not resolved.is_relative_to(root):
            issues.append(f"{path.relative_to(root)}: external local link {target}")
        elif not resolved.is_file():
            issues.append(f"{path.relative_to(root)}: missing link {target}")
    class Links(HTMLParser):
        def handle_starttag(self, tag, attrs):
            for name, value in attrs:
                if name in {"href", "src"} and value:
                    check(self.path, value)
    for path in sorted(root.rglob("*.html")):
        parser = Links()
        parser.path = path
        parser.feed(path.read_text(encoding="utf-8"))
    for path in sorted(root.rglob("*.xlsx")):
        # XLSX external hyperlinks live in worksheet relationships. Reading those
        # directly avoids instantiating every cell in 238 potentially large tables.
        # OOXML relationship targets remain relative to the workbook's directory.
        with zipfile.ZipFile(path) as book:
            for name in book.namelist():
                if name.startswith('xl/worksheets/_rels/') and name.endswith('.rels'):
                    for relationship in ET.fromstring(book.read(name)):
                        if relationship.get('Type', '').endswith('/hyperlink'):
                            check(path, relationship.get('Target', ''))
    return issues


def publish_canonical(destination: Path) -> Path:
    """Rebuild all three runs without changing any source bundle."""
    destination = destination.expanduser().resolve()
    if destination.exists():
        raise FileExistsError(f'canonical share requires a fresh destination: {destination}')
    cfg = tomllib.loads((HERE / 'canonical_runs.toml').read_text(encoding='utf-8'))
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.canonical-share-', dir=destination.parent))
    try:
        all_rows, run_links = [], []
        for key, run in cfg['runs'].items():
            root = Path(run['root'])
            rows = _rows(HERE / run['manifest'])
            if len(rows) != run['n_candidates']:
                raise ValueError(f'{key}: expected {run["n_candidates"]} candidates, found {len(rows)}')
            for row in rows:
                source = contained_file(root / 'jobs', row['run_name'])
                issues = validate_files(source)
                if issues:
                    raise ValueError(f'{row["run_name"]}: {issues}')
                target = staging / 'bundles' / key / row['run_name']
                shutil.copytree(source, target)
                all_rows.append((key, row))
            for name in ('RUN_REPORT.html', 'Supplementary_Table_S1.xlsx', 'Supplementary_Table_S1.tsv'):
                copy(root / name, staging / 'runs' / key / name, name)
            # Preserve relative jobs links used by reports without copying a second bundle tree.
            if (staging / 'runs' / key / 'RUN_REPORT.html').is_file():
                report = staging / 'runs' / key / 'RUN_REPORT.html'
                relocate_run_links(report.parent, key)
                run_links.append(f'<li><a href="runs/{key}/RUN_REPORT.html">{html.escape(key)} run report</a></li>')
            for folder in ('discovery', 'af3'):
                if (root / folder).is_dir():
                    shutil.copytree(root / folder, staging / 'runs' / key / folder)
                    if folder == 'discovery':
                        relocate_discovery_links(staging / 'runs' / key / folder, key)
                    if folder == 'discovery' and (root / folder / 'DISCOVERY_REPORT.html').is_file():
                        run_links.append(f'<li><a href="runs/{key}/discovery/DISCOVERY_REPORT.html">{html.escape(key)} discovery report</a></li>')
        overall = PROJ / 'reports'
        for path in sorted(overall.glob('run_report*.docx')) + sorted(overall.glob('run_report*.pdf')):
            copy(path, staging / 'overall_report' / path.name, 'overall report')
            run_links.append(f'<li><a href="overall_report/{quote(path.name)}">{html.escape(path.name)}</a></li>')
        rows_html = ''.join('<tr><td>' + html.escape(key) + '</td><td><a href="bundles/' + key + '/' +
                            quote(row['run_name']) + '/REPORT.html">' + html.escape(row['run_name']) +
                            '</a></td><td>' + html.escape(row.get('family', '')) + '</td></tr>' for key, row in all_rows)
        (staging / 'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
            '<title>RIIO canonical reports</title><h1>RIIO canonical reports</h1><p>Predictor ' +
            html.escape(cfg['provenance']['predictor_commit']) + '; ' + str(len(all_rows)) +
            ' recorded candidates. Computational predictions.</p><ul>' + ''.join(run_links) +
            '</ul><table><tr><th>Run</th><th>Regulator report</th><th>Family</th></tr>' + rows_html + '</table></html>', encoding='utf-8')
        link_issues = validate_share_links(staging)
        if link_issues:
            raise ValueError('Broken share links: ' + '; '.join(link_issues[:30]))
        (staging / 'share_manifest.json').write_text(json.dumps({'schema_version': 1,
            'candidate_count': len(all_rows), 'predictor_commit': cfg['provenance']['predictor_commit'],
            'artifacts': inventory(staging)}, indent=2), encoding='utf-8')
        staging.replace(destination)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return destination


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-bundles", action="store_true", help="create a summary without individual bundles")
    parser.add_argument('--canonical', type=Path, help='fresh three-run canonical share destination')
    args = parser.parse_args(argv)
    try:
        if args.canonical:
            publish_canonical(args.canonical)
            return 0
        publish(no_bundles=args.no_bundles)
    except (OSError, ValueError) as exc:
        print(f"Share publication failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
