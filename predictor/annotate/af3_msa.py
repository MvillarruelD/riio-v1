"""
af3_msa.py -- mine the AF3 server MSAs to deepen the family homology finder (improvement, 2026-06-15).

The conservation-first ArsR architecture (the AF3 verdict: structure can't discriminate ArsR, conservation
must) rests on the homolog finder -- which today uses a 92-sequence LOCAL BLAST DB. But every AF3 job ships
a DEEP MSA (`msas/*_unpaired_msa_*.a3m`, ~10^4 homologs with UniRef IDs + taxonomy) computed for free by
the server. This module reuses those MSAs to:

  * `parse_a3m`           -- stream (header, aligned_seq, uniref_id, taxid, genus, identity-to-query);
  * `homolog_stats`       -- depth, taxonomic diversity (taxids/genera), and the identity distribution;
  * `conservation_confidence` -- a calibrated [0,1] homolog-AVAILABILITY score (the axis we've lacked:
                            "is the evolutionary signal usable for this TF, or must something else carry
                            discovery?"). Favours MANY homologs across DIVERSE genera at READABLE identity
                            (close homologs most likely share the operator; diversity de-redundifies);
  * `write_homolog_fasta` -- dump the homolog protein sequences so they can be added to the tblastn DB /
                            used to pull more homolog promoter regions (the genome-neighborhood expansion).

Pure / offline (the .a3m is already on disk). Run `python af3_msa.py --benchmark` to tabulate the
homolog-availability axis across the AF3-folded TFs.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from glob import glob
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
_AF3 = resources.cache_path("af3_decision")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

_TAX = re.compile(r"Tax=([^=]+?)\s+TaxID=(\d+)")
_UNIREF = re.compile(r">(\w+?_\w+)")


def _identity_to_query(seq, query):
    """Fraction of query (uppercase, non-gap) columns matched. a3m: uppercase = aligned to query columns,
    lowercase = insertions (ignored). seq is already in query-column register for uppercase letters."""
    up = [c for c in seq if not c.islower()]               # aligned columns (incl. '-')
    n = m = 0
    for q, c in zip(query, up):
        if q == "-":
            continue
        n += 1
        if c == q:
            m += 1
    return m / n if n else 0.0


def parse_a3m(path, *, max_seqs=None):
    """Stream homolog records from an a3m. Yields dict(header, uniref, taxid, genus, identity). The first
    record is the query (identity 1.0)."""
    path = Path(path)
    query = None
    header = None
    buf = []
    n = 0

    def emit(h, s):
        uniref = (_UNIREF.match(h).group(1) if _UNIREF.match(h) else "")
        mt = _TAX.search(h)
        taxname, taxid = (mt.group(1), mt.group(2)) if mt else ("", "")
        genus = taxname.split()[0] if taxname else ""
        ident = 1.0 if query is None else _identity_to_query(s, query)
        return {"header": h, "uniref": uniref, "taxid": taxid, "genus": genus,
                "taxname": taxname, "identity": ident, "seq": s}

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    s = "".join(buf)
                    if query is None:
                        query = "".join(c for c in s if not c.islower())
                        yield emit(header, s)
                    else:
                        yield emit(header, s)
                        n += 1
                        if max_seqs and n >= max_seqs:
                            return
                header, buf = line, []
            else:
                buf.append(line)
        if header is not None:
            yield emit(header, "".join(buf))


def homolog_stats(a3m_path, *, max_seqs=20000):
    recs = list(parse_a3m(a3m_path, max_seqs=max_seqs))
    homs = recs[1:] if recs else []
    idents = [r["identity"] for r in homs]
    taxids = {r["taxid"] for r in homs if r["taxid"]}
    genera = {r["genus"] for r in homs if r["genus"]}
    def frac(lo, hi):
        return sum(lo <= i < hi for i in idents) / len(idents) if idents else 0.0
    n_close = frac(0.5, 0.95) * len(homs)                  # ABSOLUTE close-ortholog count (the key metric)
    return {"depth": len(homs), "n_taxids": len(taxids), "n_genera": len(genera),
            "frac_close": frac(0.5, 0.95), "n_close": round(n_close),
            "frac_twilight": frac(0.25, 0.5), "frac_near_dup": frac(0.95, 1.01),
            "mean_identity": (sum(idents) / len(idents) if idents else 0.0)}


def conservation_confidence(stats):
    """Calibrated [0,1] homolog-availability score, driven by the ABSOLUTE number of CLOSE orthologs
    (>=50% id). Operator conservation tracks close orthology, NOT distant family membership: the ArsR/SmtB
    AF3 MSAs carry ~17k homologs but ~1800 genera of DISTANT paralogs with different operators/effectors,
    so raw depth saturates and is uninformative. What separates a usable conservation signal is how many
    CLOSE orthologs (likely to retain the operator) exist; genus diversity de-redundifies them."""
    n_close = stats.get("n_close", stats["frac_close"] * stats["depth"])
    close = min(1.0, math.log10(max(1, n_close)) / math.log10(200))   # 1.0 at ~200 close orthologs
    div = 1.0 / (1.0 + math.exp(-(math.log10(max(1, stats["n_genera"])) - 1.0) * 2.0))
    return round(0.75 * close + 0.25 * div, 3)


def write_homolog_fasta(a3m_path, out_path, *, min_identity=0.25, max_identity=0.95, max_seqs=20000):
    """Dump de-gapped homolog protein sequences in the usable identity band (drop near-dups and the
    unalignable tail) -- ready to add to the tblastn DB for genome-neighborhood homolog expansion."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out_path, "w", encoding="utf-8") as out:
        for r in parse_a3m(a3m_path, max_seqs=max_seqs):
            if r["identity"] >= 1.0:                        # query
                continue
            if not (min_identity <= r["identity"] <= max_identity):
                continue
            seq = "".join(c for c in r["seq"] if c.isalpha()).upper()
            if len(seq) < 40:
                continue
            out.write(f">{r['uniref'] or 'hom'}|{r['taxname']}|id{r['identity']:.2f}\n{seq}\n")
            n += 1
    return n


def _op_msa(tf_id):
    hits = glob(str(_AF3 / f"*{tf_id.lower()}_op_af3decision" / "msas" / "*unpaired*.a3m"))
    return hits[0] if hits else None


def benchmark():
    tfs = {"CueR_Ecoli": "MerR", "ZiaR_Synechocystis": "ArsR/SmtB", "NmtR_Mtb": "ArsR/SmtB",
           "SmtB_Synechococcus": "ArsR/SmtB", "CzrA_Saureus": "ArsR/SmtB"}
    print(f"{'TF':<22}{'family':<11}{'depth':>7}{'genera':>8}{'n_close':>8}"
          f"{'fclose':>7}{'twil':>7}{'CONF':>7}")
    rows = {}
    for tf, fam in tfs.items():
        p = _op_msa(tf)
        if not p:
            print(f"{tf:<22}{fam:<11}  (no AF3 MSA)"); continue
        s = homolog_stats(p)
        conf = conservation_confidence(s)
        rows[tf] = {**s, "family": fam, "conservation_confidence": conf}
        print(f"{tf:<22}{fam:<11}{s['depth']:>7}{s['n_genera']:>8}{s['n_close']:>8}"
              f"{s['frac_close']:>7.2f}{s['frac_twilight']:>7.2f}{conf:>7.2f}")
    import json
    resources.cache_path("homolog_availability.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8")
    print("\nsaved -> results/homolog_availability.json")
    print("CONF = homolog-availability/conservation-confidence axis (high => conservation can carry "
          "discovery; low => need structure/other evidence).")
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", action="store_true")
    ap.add_argument("--a3m", default=None)
    ap.add_argument("--fasta-out", default=None)
    a = ap.parse_args(argv)
    if a.benchmark:
        benchmark()
    elif a.a3m:
        s = homolog_stats(a.a3m)
        print("stats:", s, "conf:", conservation_confidence(s))
        if a.fasta_out:
            n = write_homolog_fasta(a.a3m, a.fasta_out)
            print(f"wrote {n} homolog sequences -> {a.fasta_out}")


if __name__ == "__main__":
    main()
