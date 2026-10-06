"""
genome_scan.py -- the genome-first entry point: a Pfam-HMM census of candidate TFs.

Same strategy as the SSN study (Rondon et al.): scan a proteome against the family Pfam HMM library
(`predictor/data/refs/tf_pfam.hmm`, 89 models spanning 22 TF families) at the curated **gathering** thresholds, and
report every protein that hits a TF family. It is ONE batched `pyhmmer` hmmsearch over the whole proteome --
seconds for a bacterial genome, no per-protein reload -- so it is cheap enough to run on any genome.

This opens the pipeline's second avenue:

  * avenue A (this module, cheap): genome/proteome in  ->  list of candidate TFs + family + Pfam hit.
  * avenue B (`--predict`): run the existing per-TF operator pipeline on each candidate to add its
    inducer and operators. Heavier (homolog fetch + genome rescan per TF), so it is opt-in and defaults to
    the families the operator layer actually supports.

A raw nucleotide genome is gene-called in-process with `pyrodigal` first; a protein FASTA (proteome) is
scanned directly.

    tfop scan proteome.faa
    tfop scan genome.fna --families ArsR/SmtB,MerR,Fur --predict --organism GCF_000195955.2
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from dataclasses import asdict, dataclass
from html import escape
from pathlib import Path
from urllib.parse import quote, unquote

try:
    from .family_db import _PFAM_HMM, _PFAM_MAP
    from .genome_resolver import is_dna
    from . import ssn_clusters as _ssn
except ImportError:                                            # run as a script
    from predictor.annotate.family_db import _PFAM_HMM, _PFAM_MAP
    from predictor.annotate.genome_resolver import is_dna
    from predictor.annotate import ssn_clusters as _ssn

#: Families the downstream pipeline can source SSN-cluster homologs for -- i.e. every family with an SSN.
#: `scan` reports EVERY family it detects; `--predict` runs the operator pipeline on the candidates in
#: these. This was once a hand-picked 8 of the 12, chosen by "primary function is metal sensing", which
#: left MarR/SlyA and LysR-type needing an extra opt-in flag to be predicted at all. Every family now goes
#: through the same route, so the list is just `ssn.FAMILIES` and there is nothing to opt into.
PIPELINE_FAMILIES = tuple(_ssn.FAMILIES)
# Families in the census whose best-known members commonly include metal-responsive regulators.
# This is a biological grouping for overview graphics, not a confidence tier: other family members
# can sense non-metal ligands, and a dark bar does not establish metal responsiveness for a protein.
POTENTIAL_METALLOREGULATOR_FAMILIES = (
    "ArsR/SmtB", "CopY", "CsoR/FrmR", "DtxR/MntR", "Fur", "MerR", "NikR", "Rrf2",
)
_TRAIN_MIN_BP = 100_000            # Prodigal's own single-genome training floor; below this -> meta mode

_s = lambda x: x.decode() if isinstance(x, (bytes, bytearray)) else (x or "")


@dataclass
class Candidate:
    protein_id: str
    family: str | None
    pfam_acc: str | None
    pfam_name: str | None
    evalue: float
    bits: float
    length: int
    sequence: str = ""
    contig: str | None = None
    start: int | None = None        # GFF convention: 1-based, inclusive
    end: int | None = None          # GFF convention: 1-based, inclusive
    strand: str | None = None
    locus_tag: str | None = None
    gene: str | None = None
    product: str | None = None
    coordinate_source: str | None = None

    @property
    def locus(self) -> str:
        if self.contig and self.start is not None and self.end is not None:
            return f"{self.contig}:{self.start:,}-{self.end:,} ({self.strand or '?'})"
        return "Not available"


# --------------------------------------------------------------------------- input -> proteins
def read_fasta(path):
    recs, rid, buf = [], None, []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if rid is not None:
                recs.append((rid, "".join(buf)))
            rid, buf = line[1:].split()[0] if line[1:].strip() else "seq", []
        elif line.strip():
            buf.append(line.strip())
    if rid is not None:
        recs.append((rid, "".join(buf)))
    return recs


def proteins_from_gff(fna_path, gff_path, *, return_metadata=False):
    """[(protein_id, translation)] from a genome's OWN annotation: GFF CDS x genome FASTA.

    Segments sharing a protein_id are joined in genome order before translation, so a CDS split
    across several GFF lines yields the whole protein rather than its first exon.
    """
    from collections import defaultdict

    from Bio.Seq import Seq

    contigs = dict(read_fasta(fna_path))
    segments = defaultdict(list)
    metadata = defaultdict(dict)
    for line in Path(gff_path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "CDS":
            continue
        attrs = {k: unquote(v) for k, v in
                 (kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)}
        pid = attrs.get("protein_id") or attrs.get("Name") or attrs.get("ID")
        if pid and c[0] in contigs:
            segments[pid].append((c[0], int(c[3]), int(c[4]), c[6]))
            current = metadata[pid]
            current.update({k: v for k, v in {
                "locus_tag": attrs.get("locus_tag"),
                "gene": attrs.get("gene"),
                "product": attrs.get("product"),
            }.items() if v and not current.get(k)})
    prots = []
    for pid, segs in segments.items():
        segs.sort(key=lambda s: s[1])
        nt = "".join(contigs[s[0]][s[1] - 1:s[2]] for s in segs)
        if segs[0][3] == "-":
            nt = str(Seq(nt).reverse_complement())
        aa = str(Seq(nt[:len(nt) - len(nt) % 3]).translate(table=11, to_stop=True))
        if aa:
            prots.append((pid, aa))
            metadata[pid].update({
                "contig": segs[0][0], "start": min(s[1] for s in segs),
                "end": max(s[2] for s in segs), "strand": segs[0][3],
                "coordinate_source": "source GFF3",
            })
    return (prots, dict(metadata)) if return_metadata else prots


def proteins_from_input(path, *, gff=None, verbose=True, return_metadata=False):
    """(list[(id, protein)], kind).  kind in {'proteome','annotated_genome','genome'}.

    A nucleotide genome is used with its OWN annotation when a sibling `.gff` exists, and only
    gene-called when it does not. That order is load-bearing rather than a nicety: pyrodigal invents
    its own protein ids and picks its own start codons, so a re-called genome cannot be joined to
    anything keyed on accession. Measured on the survey's four genomes
    (`analysis/benchmarking/discovery/`): the accession match rate against the published 150 candidates is
    **149/150 from the annotation and 0/150 from a re-call**, and only 62/150 of the re-called gene
    models even reproduce the published sequences. Everything downstream -- regulon members,
    neighbourhood flags, any join to the collaborator's tables -- is keyed on accession, so discarding
    the annotation silently destroys the join key while still looking like a successful scan.

    Falling back to a gene call for genuinely unannotated input is engine-style degradation, not a run
    mode: the route taken is returned in `kind` and recorded in the scan manifest.
    """
    recs = read_fasta(path)
    if not recs:
        result = ([], "empty", {}) if return_metadata else ([], "empty")
        return result
    if not is_dna("".join(s for _, s in recs[:3])[:2000]):
        proteins = [(i, s.replace("*", "").strip()) for i, s in recs]
        return (proteins, "proteome", {}) if return_metadata else (proteins, "proteome")
    if gff:
        annotation = Path(gff)
    else:
        annotation = next((p for p in (Path(path).with_suffix(".gff"),
                                       Path(path).with_suffix(".gff3")) if p.exists()),
                          Path(path).with_suffix(".gff"))
    if annotation.exists():
        try:
            prots, metadata = proteins_from_gff(path, annotation, return_metadata=True)
            if prots:
                if verbose:
                    print(f"read {len(prots)} annotated protein(s) from {annotation.name} "
                          f"(accessions preserved; no gene call needed)")
                return ((prots, "annotated_genome", metadata) if return_metadata
                        else (prots, "annotated_genome"))
            if verbose:
                print(f"  {annotation.name} declares no usable CDS -- falling back to a gene call")
        except Exception as e:                                 # a malformed GFF must not sink the scan
            if verbose:
                print(f"  annotation unusable ({type(e).__name__}: {e}) -- falling back to a gene call")
    # unannotated nucleotide genome -> gene-call with pyrodigal (in-process, no external binary)
    try:
        import pyrodigal
    except ImportError as e:
        raise RuntimeError("input looks like a nucleotide genome; `pip install pyrodigal` to gene-call it, "
                           "or pass a protein FASTA (proteome).") from e
    total = sum(len(s) for _, s in recs)
    meta = total < _TRAIN_MIN_BP
    gf = pyrodigal.GeneFinder(meta=True) if meta else pyrodigal.GeneFinder(meta=False)
    if not meta:
        try:
            gf.train("TTAATTAATTAA".join(s for _, s in recs))
        except Exception:
            gf = pyrodigal.GeneFinder(meta=True)               # fall back to metagenomic mode
    prots, metadata = [], {}
    for cid, seq in recs:
        for k, gene in enumerate(gf.find_genes(seq), 1):
            pid = f"{cid}_{k}"
            prots.append((pid, gene.translate().rstrip("*")))
            metadata[pid] = {
                "contig": cid, "start": int(gene.begin), "end": int(gene.end),
                "strand": "+" if int(gene.strand) == 1 else "-",
                "locus_tag": pid, "product": "hypothetical protein",
                "coordinate_source": "Pyrodigal gene call",
            }
    if verbose:
        print(f"gene-called {len(prots)} ORFs from {len(recs)} contig(s) "
              f"({'single-genome' if not meta else 'metagenomic'} mode)")
    return (prots, "genome", metadata) if return_metadata else (prots, "genome")


# --------------------------------------------------------------------------- HMM census
def scan_proteins(records, *, families=None, cpus=0, metadata=None):
    """One hmmsearch of all proteins vs the Pfam HMM library at gathering cutoffs -> [Candidate].

    `records` is [(id, protein_seq)]. `families` (optional) keeps only those family labels. Returns one
    Candidate per protein that passes the gathering threshold, best (lowest-e) family hit, sorted by
    family then significance.
    """
    if not _PFAM_HMM.exists() or not _PFAM_MAP.exists():
        raise RuntimeError(f"Pfam HMM library missing ({_PFAM_HMM}); build it with "
                           "`python -m predictor.annotate.family_db build`.")
    import pyhmmer
    pmap = json.loads(_PFAM_MAP.read_text(encoding="utf-8"))
    metadata = metadata or {}
    alpha = pyhmmer.easel.Alphabet.amino()
    seqmap = {i: s for i, s in records}
    digital = [pyhmmer.easel.TextSequence(name=i.encode(), sequence=s).digitize(alpha)
               for i, s in records if s]
    with pyhmmer.plan7.HMMFile(str(_PFAM_HMM)) as hf:
        hmms = list(hf)
    best = {}                                                   # protein_id -> (evalue, bits, acc, name, fam)
    for hits in pyhmmer.hmmer.hmmsearch(hmms, digital, bit_cutoffs="gathering", cpus=cpus):
        name = _s(hits.query.name)
        acc = _s(hits.query.accession).split(".")[0] or name
        fam = pmap.get(acc) or pmap.get(name)
        for h in hits:
            if not h.included:
                continue
            pid = _s(h.name)
            if pid not in best or h.evalue < best[pid][0]:
                best[pid] = (float(h.evalue), float(h.score), acc, name, fam)
    out = []
    for pid, (ev, bits, acc, name, fam) in best.items():
        if families and fam not in families:
            continue
        seq = seqmap.get(pid, "")
        out.append(Candidate(pid, fam, acc, name, ev, bits, len(seq), seq,
                             **{k: v for k, v in metadata.get(pid, {}).items()
                                if k in Candidate.__dataclass_fields__}))
    fam_order = {f: i for i, f in enumerate(PIPELINE_FAMILIES)}
    out.sort(key=lambda c: (fam_order.get(c.family, 99), c.family or "~", c.evalue))
    return out


def scan(path, *, gff=None, families=None, cpus=0, verbose=True):
    """Top level: proteome/genome path -> (candidates, kind, n_proteins_scanned)."""
    prots, kind, metadata = proteins_from_input(path, gff=gff, verbose=verbose, return_metadata=True)
    t = time.time()
    cands = scan_proteins(prots, families=families, cpus=cpus, metadata=metadata)
    if verbose:
        print(f"scanned {len(prots)} proteins vs Pfam HMMs in {time.time()-t:.1f}s -> "
              f"{len(cands)} candidate TF(s)")
    return cands, kind, len(prots)


# --------------------------------------------------------------------------- outputs
def _write_locus_tracks(cands, out_dir: Path) -> dict:
    """Write interoperable GFF3 and BED overlays when genome coordinates are available."""
    located = [c for c in cands if c.contig and c.start is not None and c.end is not None]
    if not located:
        return {}
    gff = out_dir / "candidate_tfs.gff3"
    with gff.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("##gff-version 3\n")
        for c in located:
            attrs = {
                "ID": f"tfop:{c.protein_id}", "Name": c.locus_tag or c.protein_id,
                "protein_id": c.protein_id, "tf_family": c.family or "unresolved",
                "pfam_accession": c.pfam_acc or "", "hmm_evalue": f"{c.evalue:.3e}",
                "bitscore": f"{c.bits:.1f}", "prediction_status": "PREDICTED_UNVERIFIED",
                "coordinate_source": c.coordinate_source or "unknown",
            }
            if c.gene:
                attrs["gene"] = c.gene
            if c.product:
                attrs["product"] = c.product
            encoded = ";".join(f"{k}={quote(str(v), safe='._:-')}" for k, v in attrs.items())
            fh.write(f"{c.contig}\ttfop\tgene\t{c.start}\t{c.end}\t{c.bits:.1f}\t"
                     f"{c.strand or '.'}\t.\t{encoded}\n")
    bed = out_dir / "candidate_tfs.bed"
    with bed.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write('track name="TF candidates" description="Predicted TF gene loci" '
                 'visibility=pack itemRgb="On"\n')
        for c in located:
            color = ("14,79,64" if c.family in POTENTIAL_METALLOREGULATOR_FAMILIES
                     else "87,131,115")
            name = (c.locus_tag or c.protein_id).replace("\t", "_").replace(" ", "_")
            fh.write(f"{c.contig}\t{max(0, c.start - 1)}\t{c.end}\t{name}\t0\t"
                     f"{c.strand or '.'}\t{max(0, c.start - 1)}\t{c.end}\t{color}\n")
    return {"candidate_gff3": gff.name, "candidate_bed": bed.name}


def _copy_source_inputs(out_dir: Path, *, kind: str, source_path=None, annotation_path=None,
                        keep_existing: bool = False) -> dict:
    """Keep source sequence/annotation beside the overlay so a scan folder is self-contained."""
    files = {}
    source_path = Path(source_path) if source_path and Path(source_path).is_file() else None
    annotation_path = (Path(annotation_path)
                       if annotation_path and Path(annotation_path).is_file() else None)
    target_name = "source_proteome.faa" if kind == "proteome" else "source_genome.fna"
    target = out_dir / target_name
    if source_path:
        if source_path.resolve() != target.resolve():
            shutil.copyfile(source_path, target)
        files["source_sequence"] = target.name
        if not annotation_path and kind == "annotated_genome":
            annotation_path = next((p for p in (source_path.with_suffix(".gff"),
                                                 source_path.with_suffix(".gff3")) if p.is_file()), None)
    elif keep_existing and target.is_file():
        # The GUI writes the scan twice: once after discovery and again after selected predictions.
        # The second call has no temporary upload path, but the copied source still belongs in the
        # regenerated manifest.
        files["source_sequence"] = target.name
    annotation_target = out_dir / "source_annotation.gff3"
    if annotation_path:
        if annotation_path.resolve() != annotation_target.resolve():
            shutil.copyfile(annotation_path, annotation_target)
        files["source_annotation"] = annotation_target.name
    elif keep_existing and annotation_target.is_file():
        files["source_annotation"] = annotation_target.name
    return files


def write_outputs(cands, out_dir, *, kind="proteome", n_proteins=0, source="", predictions=None,
                  source_path=None, annotation_path=None):
    """Write the machine-readable data files + a self-contained HTML report for a scan.

    predictions (optional, from avenue B): {protein_id: {"inducer":.., "primary_operator":.., "bundle":..}}
    is merged into candidate_tfs.json and shown in the report.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    predictions = predictions or {}
    if source_path is None and source and Path(str(source)).is_file():
        source_path = source
    previous_manifest = {}
    manifest_path = out_dir / "scan_manifest.json"
    if manifest_path.is_file():
        try:
            previous_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            pass
    keep_existing = (previous_manifest.get("input") == source
                     and previous_manifest.get("input_kind") == kind)
    source_files = _copy_source_inputs(out_dir, kind=kind, source_path=source_path,
                                       annotation_path=annotation_path, keep_existing=keep_existing)
    track_files = _write_locus_tracks(cands, out_dir)
    records = []
    for c in cands:
        r = asdict(c)
        if c.protein_id in predictions:
            r["prediction"] = predictions[c.protein_id]
        records.append(r)
    # manifest so a downstream tool / agent can discover the data files without guessing
    manifest = {
        "tool": "tf-operator-predictor genome scan",
        "avenue": "B (candidates + operators)" if predictions else "A (candidate TF census)",
        "input": source, "input_kind": kind, "proteins_scanned": n_proteins,
        "n_candidates": len(cands),
        "families": dict(family_counts(cands)),
        "method": "pyhmmer hmmsearch vs packaged tf_pfam.hmm at Pfam gathering thresholds",
        # Whether protein_id is a real accession decides if these rows can be joined to anything else.
        "protein_ids": {
            "annotated_genome": "the genome's own GFF protein_id -- real accessions, joinable",
            "proteome": "the input FASTA's own ids -- joinable if the input carried accessions",
            "genome": "pyrodigal-invented <contig>_<n> -- NOT accessions and NOT joinable; the input "
                      "had no sibling .gff to annotate from",
        }.get(kind, "unknown"),
        "files": {
            "candidates": "candidate_tfs.json  (list; each: protein_id, family, pfam_acc/name, evalue, "
                          "bits, length, sequence, genomic locus[, prediction])",
            "candidates_tsv": "candidate_tfs.tsv", "candidates_faa": "candidate_tfs.faa",
            "report": "SCAN_REPORT.html",
            **track_files,
            **source_files,
            **({"predictions_tsv": "candidate_predictions.tsv"} if predictions else {}),
            **({"per_tf_bundles": "results/jobs/<protein_id>/  (dossier.json + REPORT.html each)"}
               if predictions else {}),
        },
    }
    (out_dir / "candidate_tfs.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    (out_dir / "scan_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    with (out_dir / "candidate_tfs.tsv").open("w", encoding="utf-8") as fh:
        fh.write("protein_id\tlocus_tag\tgene\tcontig\tstart\tend\tstrand\tfamily\tpfam_acc\t"
                 "pfam_name\tevalue\tbits\tlength\tcoordinate_source\n")
        for c in cands:
            values = (c.protein_id, c.locus_tag, c.gene, c.contig, c.start, c.end, c.strand,
                      c.family, c.pfam_acc, c.pfam_name, f"{c.evalue:.3e}", f"{c.bits:.1f}",
                      c.length, c.coordinate_source)
            fh.write("\t".join("" if value is None else str(value).replace("\t", " ")
                               for value in values) + "\n")
    with (out_dir / "candidate_tfs.faa").open("w", encoding="ascii", errors="replace") as fh:
        for c in cands:
            fh.write(f">{c.protein_id} {c.family}|{c.pfam_name}\n{c.sequence}\n")
    if predictions:
        with (out_dir / "candidate_predictions.tsv").open("w", encoding="utf-8") as fh:
            fh.write("protein_id\tfamily\tinducer\tprimary_operator\tstatus\terror\tbundle\n")
            for c in cands:
                if c.protein_id not in predictions:
                    continue
                p = predictions[c.protein_id]
                values = (c.protein_id, c.family, p.get("inducer"), p.get("primary_operator"),
                          p.get("status"), p.get("error"), p.get("bundle"))
                fh.write("\t".join("" if value is None else str(value).replace("\t", " ")
                                   for value in values) + "\n")
    write_scan_report(cands, out_dir / "SCAN_REPORT.html", kind=kind, n_proteins=n_proteins,
                      source=source, predictions=predictions)
    return out_dir


def family_counts(cands):
    from collections import Counter
    return Counter(c.family for c in cands)


def _esc(s):
    return escape(str(s), quote=True)


def write_scan_report(cands, html_path, *, kind="proteome", n_proteins=0, source="", predictions=None):
    """A self-contained, organized HTML report of the candidate-TF census."""
    predictions = predictions or {}
    counts = family_counts(cands)
    supported = sum(counts.get(f, 0) for f in PIPELINE_FAMILIES)
    top = max(counts.values()) if counts else 1
    bars = "".join(
        f'<tr><td class="fam">{_esc(f)}</td><td class="barcell"><span class="bar'
        f'{" metal" if f in POTENTIAL_METALLOREGULATOR_FAMILIES else " supported" if f in PIPELINE_FAMILIES else ""}'
        f'" style="width:{max(2, 100*k//top)}%"></span>'
        f'<span class="n">{k}</span></td></tr>'
        for f, k in counts.most_common())
    has_pred = bool(predictions)
    inducer_counts = {}
    for c in cands:
        if c.protein_id in predictions:
            inducer = predictions[c.protein_id].get("inducer") or "Unresolved"
            key = (c.family or "Unresolved family", inducer)
            inducer_counts[key] = inducer_counts.get(key, 0) + 1
    inducer_top = max(inducer_counts.values()) if inducer_counts else 1
    inducer_bars = "".join(
        f'<tr><td class="fam">{_esc(fam)}</td><td>{_esc(inducer)}</td>'
        f'<td class="barcell"><span class="bar inducer" style="width:{max(2, 100*n//inducer_top)}%"></span>'
        f'<span class="n">{n}</span></td></tr>'
        for (fam, inducer), n in sorted(inducer_counts.items(), key=lambda item: (item[0][0], -item[1], item[0][1])))
    head = ("<th>protein</th><th>locus tag</th><th>genomic locus</th><th>family</th><th>Pfam</th>"
            "<th><abbr title='Expected number of equally strong HMM matches by chance. Smaller values mean a more significant family match.'>e-value</abbr></th>"
            "<th><abbr title='HMM alignment score in bits. Larger values mean stronger support for this family match.'>bits</abbr></th><th>len</th>"
            + ("<th>inducer</th><th>primary operator</th><th>report</th>" if has_pred else ""))
    rows = []
    for c in cands:
        cls = ' class="supported"' if c.family in PIPELINE_FAMILIES else ""
        cells = (f"<td class='mono'>{_esc(c.protein_id)}</td>"
                 f"<td class='mono'>{_esc(c.locus_tag or '—')}</td>"
                 f"<td class='mono' title='{_esc(c.coordinate_source or 'No genomic coordinates supplied')}'>{_esc(c.locus)}</td>"
                 f"<td{cls}>{_esc(c.family)}</td>"
                 f"<td>{_esc(c.pfam_name)} <span class='dim'>{_esc(c.pfam_acc)}</span></td>"
                 f"<td class='mono'>{c.evalue:.2e}</td><td>{c.bits:.1f}</td><td>{c.length}</td>")
        if has_pred:
            p = predictions.get(c.protein_id, {})
            op = p.get("primary_operator") or ""
            rep = p.get("bundle")
            link = (f"<a href='{_esc((Path(rep).resolve() / 'REPORT.html').as_uri())}'>open</a>"
                    if rep else "<span class='dim'>-</span>")
            cells += (f"<td>{_esc(p.get('inducer') or '')}</td>"
                      f"<td class='mono'>{_esc(op)}</td><td>{link}</td>")
        rows.append(f"<tr>{cells}</tr>")
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="description" content="Candidate bacterial transcription-factor census">
<title>Candidate TF census - {_esc(source or kind)}</title>
<style>
 :root{{--ink:#18201c;--muted:#657069;--paper:#f3f5f1;--surface:#fbfcf9;--line:#d7ddd7;
 --accent:#176b58;--accent-dark:#0e4f40;--soft:#e1eee8}}
 *{{box-sizing:border-box}} html{{background:var(--paper)}}
 body{{font:14px/1.58 Aptos,"Segoe UI",Arial,sans-serif;color:var(--ink);max-width:1160px;margin:0 auto;
 padding:36px clamp(18px,4vw,52px) 72px;background-size:32px 32px,32px 32px,auto,auto;background-image:
 linear-gradient(rgba(24,32,28,.022) 1px,transparent 1px),linear-gradient(90deg,rgba(24,32,28,.022) 1px,transparent 1px),
 radial-gradient(circle at 90% 0,rgba(23,107,88,.10),transparent 28rem)}}
 .report-top{{display:flex;align-items:center;justify-content:space-between;gap:1rem;border-bottom:1px solid var(--line);padding-bottom:12px}}
 .eyebrow{{font-size:11px;letter-spacing:.16em;text-transform:uppercase;color:var(--accent);font-weight:750}}
 .repo-link{{color:var(--accent-dark);font-size:12px;font-weight:650;text-decoration:none;border-bottom:1px solid #96b5aa}}
 h1{{font-size:clamp(38px,6vw,64px);line-height:.98;letter-spacing:-.048em;margin:26px 0 10px;max-width:15ch}}
 h2{{font-size:23px;letter-spacing:-.025em;border-top:1px solid var(--line);padding-top:22px;margin-top:48px}}
 .sub{{color:var(--muted);margin:0;max-width:75ch}}
 .guide{{border-left:3px solid var(--accent);padding:13px 16px;background:var(--soft);color:#31423a;font-size:13px;margin:22px 0}}
 .scope{{display:block;margin-top:6px;color:var(--accent-dark);font-weight:650}}
 .cards{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin:24px 0}}
 .card{{border-top:2px solid var(--ink);padding:11px 2px 0;min-width:0}}
 .card .k{{font-size:28px;font-weight:720;font-variant-numeric:tabular-nums;letter-spacing:-.03em}}
 .card .l{{color:var(--muted);font-size:11px;letter-spacing:.06em;text-transform:uppercase}}
 .table-wrap{{overflow-x:auto;margin:.7rem 0 1.5rem;border:1px solid var(--line);border-radius:9px;background:var(--surface)}}
 table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
 th,td{{text-align:left;padding:8px 10px;border-bottom:1px solid var(--line);font-size:12px;vertical-align:top}}
 th{{background:#edf2ed;color:#3e4d45;font-size:10px;letter-spacing:.045em;text-transform:uppercase}}
 tr:last-child td{{border-bottom:0}} tr:nth-child(even){{background:#faf9f4}}
 td.supported,.fam{{font-weight:650}} td.supported{{color:var(--accent-dark)}}
 .mono{{font-family:"Cascadia Mono","SFMono-Regular",Consolas,monospace;overflow-wrap:anywhere}}
 .dim{{color:var(--muted);font-size:.82em}} a{{color:var(--accent-dark);font-weight:650}}
 .barcell{{position:relative;min-width:170px}} .bar{{display:inline-block;height:10px;background:#aeb9b3;border-radius:2px;vertical-align:middle}}
 .bar.supported{{background:#5f9b89}} .bar.metal{{background:var(--accent-dark)}}
 .bar.inducer{{background:#b8874b}} .n{{margin-left:.5rem;color:var(--muted);font-size:11px}}
 .legend{{display:flex;gap:1.2rem;flex-wrap:wrap;color:var(--muted);font-size:12px;margin:.2rem 0 .7rem}}
 .swatch{{display:inline-block;width:18px;height:8px;border-radius:2px;margin-right:.35rem;background:#aeb9b3}}
 .swatch.supported{{background:#5f9b89}}.swatch.metal{{background:var(--accent-dark)}}
 abbr{{text-decoration:underline dotted;text-underline-offset:3px;cursor:help}}
 .foot{{color:var(--muted);font-size:12px;margin-top:2rem;border-top:1px solid var(--line);padding-top:.9rem}}
 code{{background:var(--soft);padding:.08rem .3rem;border-radius:3px}}
 @media(max-width:640px){{body{{padding:24px 15px 48px}}.cards{{grid-template-columns:1fr}}h1{{font-size:42px}}}}
 @media print{{html,body{{background:white}}body{{max-width:none;padding:12mm}}.table-wrap{{overflow:visible}}}}
</style></head><body>
<header><div class="report-top"><span class="eyebrow">Bacterial regulatory genomics</span>
<a class="repo-link" href="https://github.com/MvillarruelD/riio-v1">GitHub ↗</a></div>
<h1>Genome TF census</h1>
<p class="sub">{_esc(source or '(input)')} · {_esc(kind)} · {n_proteins} proteins screened with the bundled
Pfam HMM library.</p></header>
<div class="cards">
  <div class="card"><div class="k">{len(cands)}</div><div class="l">candidate TFs</div></div>
  <div class="card"><div class="k">{len(counts)}</div><div class="l">families</div></div>
  <div class="card"><div class="k">{supported}</div><div class="l">ready for operator prediction</div></div>
</div>
<div class="guide"><b>Scope and interpretation.</b> This is a census of likely TF proteins, not experimental
validation. The census screens 22 families; the complete operator workflow currently supports 12.
<span class="scope">{' · '.join(_esc(f) for f in PIPELINE_FAMILIES)}</span></div>
<h2>Family overview</h2>
<div class="legend"><span><i class="swatch metal"></i>family can include metalloregulators</span>
<span><i class="swatch supported"></i>other full-report family</span><span><i class="swatch"></i>census only</span></div>
<p class="sub">Dark color indicates a family with known metal-responsive members. It is not evidence that every
candidate in that family senses a metal.</p><div class="table-wrap"><table><tbody>{bars}</tbody></table></div>
{f'<h2>Possible inducer distribution</h2><p class="sub">Counts reflect completed computational calls and retain unresolved results.</p><div class="table-wrap"><table><thead><tr><th>family</th><th>possible inducer</th><th>candidates</th></tr></thead><tbody>{inducer_bars}</tbody></table></div>' if has_pred else ''}
<h2>All candidates</h2><div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>
<p class="foot">Generated by <code>tfop scan</code> (tf-operator-predictor). Data files alongside this report:
<code>candidate_tfs.json</code> / <code>.tsv</code> / <code>.faa</code>, <code>scan_manifest.json</code>
{", <code>candidate_predictions.tsv</code>" if has_pred else ""}. When coordinates are available,
<code>candidate_tfs.gff3</code> and <code>candidate_tfs.bed</code> provide genome annotation/browser overlays.</p>
</body></html>"""
    Path(html_path).write_text(html, encoding="utf-8")
    return html_path


# --------------------------------------------------------------------------- self-test / CLI
def _main(argv):
    if not argv:
        print("usage: python -m predictor.annotate.genome_scan <proteome.faa|genome.fna> [families]")
        return 2
    fams = set(argv[1].split(",")) if len(argv) > 1 else None
    cands, kind, n = scan(argv[0], families=fams)
    print(f"\n[{kind}] {n} proteins -> {len(cands)} candidate TFs")
    for fam, k in family_counts(cands).most_common():
        print(f"  {str(fam):16s} {k}")
    for c in cands[:12]:
        print(f"    {c.protein_id:24s} {str(c.family):14s} {c.pfam_name:12s} e={c.evalue:.1g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
