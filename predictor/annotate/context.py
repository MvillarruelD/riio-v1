"""
context.py -- build a GenomeContext and extract the inter-operon region from a RESOLVED genome.

Takes the output of `genome_resolver` (accession + TF locus + a fetched genomic record) and produces:
  * a `GenomeContext` (genes) consumed by `signals/motif_rescan.py`, and
  * the **inter-operon region 5' of the TF** -- Snowprint's autoregulatory substrate, also where our
    motif finder looks first -- with `widen()` for the closed loop's re-extraction step.

GenBank parsing uses Biopython (CDS features -> genes). The same `Gene` / `GenomeContext` types as
`motif_rescan` are reused so the two modules connect end to end.

Run `python context.py` for a self-test (builds a synthetic GenBank, round-trips it, extracts region).
"""
from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path

try:                                    # share the genome model with motif_rescan
    from ..signals.motif_rescan import (Gene, GenomeContext,  # noqa: F401 - re-exported
                                        intergenic_intervals)
except ImportError:
    from predictor.signals.motif_rescan import Gene, GenomeContext

#: GFF feature types admitted into the gene model. Non-coding genes matter for two reasons even though
#: they are never protein targets: (1) an operon walk that cannot see an intervening tRNA/rRNA measures
#: the gap straight ACROSS it and overshoots the same-operon cutoff, silently truncating the operon;
#: (2) a motif hit inside an ncRNA gene body was scored as "intergenic" because the gene was invisible,
#: so real coding space was treated as promoter space. Structural RNA operons (rrn, tRNA arrays) are also
#: genuine regulon members for several metalloregulators.
_GENE_FEATURES = {"CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "antisense_RNA", "RNase_P_RNA", "SRP_RNA"}


#: Run of N inserted between replicons when a multi-replicon assembly is concatenated into one
#: coordinate space. **This number is load-bearing, not cosmetic.** It is what keeps the
#: concatenation from inventing relationships across a boundary, and it works because it is larger
#: than every span the pipeline is willing to reach across:
#:
#:   * the widest motif searched          `operator_logo.OPERATOR_WIDTHS`      max 41 bp
#:   * operator -> regulated gene         `_regulated_gene` / `regulated_first_genes` max_dist  400 bp
#:   * same-operon intergenic gap         `regulon.DEFAULT_OPERON_MAX_GAP`         150 bp
#:
#: A window containing N scores -inf (`motifs.pwm_scan`), so no site can be called inside or across
#: a pad; and because 500 exceeds the two distance limits, no site can be assigned a gene -- or walk
#: an operon -- on the far side of one. Verified on the SL1344 run: 0 of 2,681 sites overlapped a
#: pad and 0 were assigned a gene on another replicon.
#:
#: Raising any of those three limits past this value silently re-opens boundary leakage.
#: `tests/test_production_hardening.py` asserts the relation so that cannot happen quietly.
REPLICON_PAD = 500


@dataclass
class Region:
    accession: str
    start: int          # 0-based half-open, + strand coords
    end: int
    strand: str         # strand of the TF this region is upstream of
    sequence: str
    anchor: str = "tf_5prime"


# --------------------------------------------------------------------------- parsing
def from_genbank(text_or_path) -> GenomeContext:
    """Parse a GenBank record (text, file path, or open handle) into a GenomeContext.
    Coordinates are 0-based half-open; coding and annotated non-coding gene features become Genes
    (named by gene/locus_tag)."""
    from Bio import SeqIO

    if hasattr(text_or_path, "read"):
        handle, close = text_or_path, False
    elif "\n" in str(text_or_path) or str(text_or_path).lstrip().startswith("LOCUS"):
        handle, close = StringIO(str(text_or_path)), True
    else:
        handle, close = open(text_or_path, encoding="ascii"), True
    try:
        rec = SeqIO.read(handle, "genbank")
    finally:
        if close:
            handle.close()

    genes = []
    wrapped = []
    for f in rec.features:
        if f.type not in _GENE_FEATURES:
            continue
        q = f.qualifiers
        name = (q.get("gene") or q.get("locus_tag") or q.get("protein_id") or [""])[0]
        parts = tuple(getattr(f.location, "parts", ()) or ())
        if len(parts) > 1 and any(int(part.start) == 0 for part in parts) and any(
            int(part.end) == len(rec.seq) for part in parts
        ):
            wrapped.append(name or "unnamed_feature")
            continue
        genes.append(Gene(int(f.location.start), int(f.location.end),
                          "+" if f.location.strand != -1 else "-", name,
                          protein_id=(q.get("protein_id") or [""])[0],
                          product=(q.get("product") or [""])[0]))
    genes.sort(key=lambda g: g.start)
    circular = str(rec.annotations.get("topology", "")).lower() == "circular"
    return GenomeContext(rec.id or rec.name or "unknown", str(rec.seq), genes,
                         circular=circular, origin_wrapped_features=tuple(wrapped))


# --------------------------------------------------------------------------- GFF3 (offline mirror)
def _gff_attrs(field: str) -> dict:
    """Parse a GFF3 attribute column (key=value;key=value) into a dict (URL-unescaped values)."""
    from urllib.parse import unquote
    out = {}
    for kv in field.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def from_gff(gff_text_or_path, fasta_text_or_path, *, accession: str | None = None) -> GenomeContext:
    """Build a GenomeContext from a GFF3 annotation + its FASTA -- the OFFLINE path for mirrored
    genomes (no EFetch). CDS features become Genes (1-based inclusive GFF coords -> 0-based half-open);
    name prefers gene -> locus_tag, and protein_id/product are captured for effector inference.

    Multi-line (joined) CDS are merged by their feature ID into one min-start..max-end span.
    """
    # ---- sequence (single-contig assumption for the resolved region; first/largest record kept) ----
    def _read(src):
        if hasattr(src, "read"):
            return src.read()
        s = str(src)
        if "\n" in s or s.lstrip().startswith((">", "##gff", "#")):
            return s
        return Path(src).read_text(encoding="ascii", errors="replace")

    fa = _read(fasta_text_or_path)
    seqs, cur, acc0 = {}, None, None
    for line in fa.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            acc0 = acc0 or cur
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    acc = accession or acc0 or "unknown"
    # For multi-contig assemblies: concatenate all contigs (N-padded) and build an offset map
    # so GFF coordinates (per-contig) can be lifted to the concatenated coordinate space.
    # For single-contig or when the requested accession is an exact key, this is a no-op.
    #
    # Asking for ONE record by name deliberately gives you that record alone. But `acc` falls back
    # to the first record when no accession is given, and that fallback used to select the same
    # branch -- so a multi-replicon FASTA passed without an accession silently kept only its first
    # record and dropped the rest, with nothing in the output to say so. Take the single-record
    # branch only when the caller actually named a record, or when there is only one.
    if (accession and acc in seqs) or (acc in seqs and len(seqs) == 1):
        seq = seqs[acc]
        contig_offsets: dict[str, int] = {acc: 0}
    else:
        # sort largest-first so the primary chromosome dominates the early coordinate space
        ordered = sorted(seqs.items(), key=lambda kv: len(kv[1]), reverse=True)
        PAD = REPLICON_PAD
        parts, contig_offsets = [], {}
        offset = 0
        for cid, s in ordered:
            contig_offsets[cid] = offset
            parts.append(s)
            offset += len(s) + PAD
        seq = ("N" * PAD).join(parts) if parts else ""
        acc = acc0 or acc

    gff = _read(gff_text_or_path)
    merged = {}                              # feature-id -> [start1, stop1, strand, name, prot, product, off]
    circular = False
    for line in gff.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        if cols[2] == "region":
            circular = circular or _gff_attrs(cols[8]).get("Is_circular", "").lower() == "true"
            continue
        if cols[2] not in _GENE_FEATURES:
            continue
        seqid = cols[0]
        if seqid not in contig_offsets:
            continue
        off = contig_offsets[seqid]
        a = _gff_attrs(cols[8])
        fid = a.get("ID") or a.get("protein_id") or f"{seqid}:{cols[3]}-{cols[4]}"
        start1, stop1, strand = int(cols[3]) + off, int(cols[4]) + off, cols[6]
        name = a.get("gene") or a.get("locus_tag") or a.get("Name") or fid
        ltag = a.get("locus_tag") or ""
        prot = a.get("protein_id") or (a.get("Name") if a.get("Name", "").startswith(("WP_", "NP_", "YP_")) else "")
        prod = a.get("product", "")
        if fid in merged:
            m = merged[fid]
            m[0], m[1] = min(m[0], start1), max(m[1], stop1)
        else:
            merged[fid] = [start1, stop1, strand, name, prot, prod, ltag]

    seq_len = len(seq)
    # Replicon lengths, for the wrap guard below. The guard's threshold has to be the length of the
    # replicon the feature sits ON, not of the concatenation: measured against a 5 Mb concatenation,
    # a feature wrapping the origin of an 8.7 kb plasmid spans 0.2 % and sails through a test meant
    # to catch exactly that. Sorted offsets let a feature's replicon be found from its start.
    _rep_len = {}
    _ordered_offs = sorted(contig_offsets.items(), key=lambda kv: kv[1])
    for _i, (_cid, _off) in enumerate(_ordered_offs):
        _nxt = _ordered_offs[_i + 1][1] - REPLICON_PAD if _i + 1 < len(_ordered_offs) else seq_len
        _rep_len[_cid] = max(0, _nxt - _off)

    def _replicon_len(pos: int) -> int:
        best = seq_len
        for _cid, _off in _ordered_offs:
            if _off <= pos:
                best = _rep_len.get(_cid, seq_len)
        return best or seq_len

    genes, wrapped = [], []
    for fid, (s1, e1, st, nm, pr, pd, lt) in merged.items():
        # Drop circular-wrap artefacts: two CDS parts with the same GFF ID (one near the
        # origin, one near the end) merge to min(start)..max(end) which can span almost the
        # entire genome.  A real gene cannot occupy >80% of a circular replicon.
        if seq_len > 0 and (e1 - s1 + 1) > 0.8 * _replicon_len(s1 - 1):
            wrapped.append(nm or fid)
            continue
        genes.append(Gene(s1 - 1, e1, "+" if st != "-" else "-", nm, protein_id=pr, product=pd,
                          locus_tag=lt))
    genes.sort(key=lambda g: g.start)
    return GenomeContext(acc, seq, genes, circular=circular,
                         origin_wrapped_features=tuple(wrapped),
                         contig_offsets=dict(contig_offsets))


def genome_context_from_mirror(accession: str, *, mirror_dir=None) -> GenomeContext | None:
    """Load a mirrored genome (`<acc>.fna` + `<acc>.gff`) into a GenomeContext, fully offline.
    Returns None if the genome (or its GFF) is not in the mirror."""
    from pathlib import Path as _P
    if mirror_dir is None:
        from predictor.annotate.genome_db import MIRROR_DIR as mirror_dir
    mirror_dir = _P(mirror_dir)
    fna, gff = mirror_dir / f"{accession}.fna", mirror_dir / f"{accession}.gff"
    if not fna.exists() or not gff.exists():
        return None
    return from_gff(gff, fna, accession=accession)


def find_tf_gene(ctx: GenomeContext, tf_start: int, tf_end: int):
    """Match the resolved TF locus to a parsed gene by best span overlap; None if no overlap."""
    best, best_ov = None, 0
    for g in ctx.genes:
        ov = max(0, min(g.end, tf_end) - max(g.start, tf_start))
        if ov > best_ov:
            best, best_ov = g, ov
    return best


# --------------------------------------------------------------------------- inter-operon region
def inter_operon_region(ctx: GenomeContext, tf, *, upstream: int = 300, downstream: int = 20,
                        gene_aware: bool = True) -> Region:
    """Extract the region immediately 5' of the TF (its autoregulatory promoter).

    `tf` is a Gene or a (start, end, strand) tuple. With `gene_aware`, the upstream edge is clipped at
    the nearest neighbouring CDS so the window is the true intergenic gap, not a slice of an upstream
    gene. Strand-aware: for a '+' TF the region lies at lower coordinates; for '-', higher.
    """
    if isinstance(tf, Gene):
        s, e, strand = tf.start, tf.end, tf.strand
    else:
        s, e, strand = tf
    L = len(ctx.sequence)

    if strand == "+":
        win_start = max(0, s - upstream)
        win_end = min(L, s + downstream)
        if gene_aware:                                  # clip to nearest upstream gene end
            ups = [g.end for g in ctx.genes if g.end <= s and not (g.start == s and g.end == e)]
            if ups:
                win_start = max(win_start, max(ups))
    else:
        win_end = min(L, e + upstream)
        win_start = max(0, e - downstream)
        if gene_aware:                                  # clip to nearest downstream gene start
            downs = [g.start for g in ctx.genes if g.start >= e and not (g.start == s and g.end == e)]
            if downs:
                win_end = min(win_end, min(downs))

    return Region(ctx.accession, win_start, win_end, strand,
                  ctx.sequence[win_start:win_end], anchor="tf_5prime")


def widen(ctx: GenomeContext, region: Region, by: int = 200) -> Region:
    """Enlarge the upstream side of a region (the closed loop's re-extraction step)."""
    L = len(ctx.sequence)
    if region.strand == "+":
        ns = max(0, region.start - by)
        return Region(region.accession, ns, region.end, region.strand,
                      ctx.sequence[ns:region.end], region.anchor)
    ne = min(L, region.end + by)
    return Region(region.accession, region.start, ne, region.strand,
                  ctx.sequence[region.start:ne], region.anchor)


def locate_operator_window(ctx: GenomeContext, tf_start: int, tf_end: int, tf_strand: str):
    """The autoregulatory operator region = the intergenic gap immediately upstream of the TF gene
    (divergent promoter). Returns (center, start, end) in genome coordinates.

    (Relocated from the benchmark's live_refine_real._locate_operator so the production pipeline no longer
    depends on the validation layer; logic is unchanged.)"""
    genes = sorted(ctx.genes, key=lambda g: g.start)
    if tf_strand == "-":
        # upstream of a minus-strand gene = higher genomic coords (above tf_end)
        nxt = next((g for g in genes if g.start >= tf_end
                    and not (g.start == tf_start and g.end == tf_end)), None)
        gap_lo = tf_end
        gap_hi = nxt.start if nxt else min(len(ctx.sequence), tf_end + 300)
    else:
        # upstream of a plus-strand gene = lower genomic coords (below tf_start)
        prev = None
        for g in genes:
            if g.start >= tf_end:
                break
            prev = g
        gap_lo = (prev.end if prev else max(0, tf_start - 200))
        gap_hi = tf_start
        if gap_hi - gap_lo < 20:                          # fall back to a fixed upstream window
            gap_lo, gap_hi = tf_start - 160, tf_start
    # For autoregulatory TFs the operator overlaps the promoter, within ~75 bp of the 5' end.
    # A naive midpoint of a long gap (e.g. IS26 composite transposon, 684 bp) pushes the
    # initial fold center far from the promoter; cap at 75 bp from the 5' anchor.
    gap_len = gap_hi - gap_lo
    if tf_strand == "-":
        center = gap_lo + min(gap_len // 2, 75)    # operator near gap_lo (5' of - gene)
    else:
        center = gap_hi - min(gap_len // 2, 75)    # operator near gap_hi (5' of + gene)
    return center, gap_lo, gap_hi


# --------------------------------------------------------------------------- self-test
def _build_genbank() -> str:
    """Synthesize a small GenBank record: an upstream divergent gene, the TF, and a downstream gene,
    with a planted 'operator' tag in the TF's 5' intergenic gap."""
    import random
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    random.seed(5)
    op = "TTGACCTAGGTCAA"
    seq = list("".join(random.choice("ACGT") for _ in range(1200)))
    seq[330:330 + len(op)] = list(op)                   # operator in the 250..400 intergenic gap
    rec = SeqRecord(Seq("".join(seq)), id="SYN_0001.1", name="SYN0001",
                    description="synthetic test genome")
    rec.annotations["molecule_type"] = "DNA"
    feats = [
        (100, 250, -1, "g0_div"),       # divergent upstream gene (- strand)
        (400, 700, +1, "tfR"),          # the TF (+ strand) -> intergenic gap 250..400 is its promoter
        (760, 1050, +1, "g2"),
    ]
    for s, e, strand, name in feats:
        rec.features.append(SeqFeature(FeatureLocation(s, e, strand=strand),
                                       type="CDS", qualifiers={"gene": [name]}))
    out = StringIO()
    from Bio import SeqIO
    SeqIO.write(rec, out, "genbank")
    return out.getvalue()


def _demo() -> None:
    gb = _build_genbank()
    ctx = from_genbank(gb)
    print(f"parsed {ctx.accession}: {len(ctx.sequence)} bp, {len(ctx.genes)} genes")
    for g in ctx.genes:
        print(f"  {g.name:<8} {g.start}-{g.end}({g.strand})")

    tf = find_tf_gene(ctx, 400, 700)
    assert tf is not None and tf.name == "tfR", "TF gene not matched from locus"

    reg = inter_operon_region(ctx, tf, upstream=300, downstream=20)
    print(f"inter-operon region: {reg.start}-{reg.end}({reg.strand}) len={len(reg.sequence)} "
          f"anchor={reg.anchor}")
    # the operator at 330 must fall inside the extracted promoter region
    assert reg.start <= 330 < reg.end, "operator not inside extracted region"
    # gene-aware clip: upstream edge should not cross into g0_div (ends at 250)
    assert reg.start >= 250, "region should be clipped at the upstream gene end (250)"
    assert "TTGACCTAGGTCAA" in reg.sequence, "operator sequence missing from region"

    wide = widen(ctx, reg, by=200)
    assert wide.start < reg.start and len(wide.sequence) > len(reg.sequence), "widen failed"
    print(f"widened region: {wide.start}-{wide.end} len={len(wide.sequence)}")

    # ---- GFF3 + FASTA offline path (mirror genomes) must reproduce the same genes ----
    fasta = f">SYN_0001.1 synthetic\n{ctx.sequence}\n"
    gff_lines = ["##gff-version 3"]
    for s, e, strand, name in [(100, 250, "-", "g0_div"), (400, 700, "+", "tfR"), (760, 1050, "+", "g2")]:
        attrs = f"ID=cds-{name};gene={name};locus_tag={name};protein_id=WP_{name};product={name} protein"
        gff_lines.append(f"SYN_0001.1\tRefSeq\tCDS\t{s + 1}\t{e}\t.\t{strand}\t0\t{attrs}")
    gctx = from_gff("\n".join(gff_lines), fasta, accession="SYN_0001.1")
    print(f"from_gff parsed {gctx.accession}: {len(gctx.sequence)} bp, {len(gctx.genes)} CDS genes")
    assert gctx.sequence == ctx.sequence, "GFF-path sequence mismatch"
    assert [(g.start, g.end, g.strand, g.name) for g in gctx.genes] == \
           [(g.start, g.end, g.strand, g.name) for g in ctx.genes], "GFF genes != GenBank genes"
    tfg = find_tf_gene(gctx, 400, 700)
    assert tfg is not None and tfg.name == "tfR" and tfg.protein_id == "WP_tfR", "GFF protein_id lost"
    rg = inter_operon_region(gctx, tfg, upstream=300, downstream=20)
    assert "TTGACCTAGGTCAA" in rg.sequence, "operator missing from GFF-path region"
    print("OK: GFF3+FASTA offline path reproduces genes/coords/region and carries protein_id/product.")
    print("OK: genome parsed, TF located, inter-operon region extracted (operator inside), widen works.")


if __name__ == "__main__":
    _demo()
