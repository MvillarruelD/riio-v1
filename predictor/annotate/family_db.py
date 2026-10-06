"""
family_db.py -- holistic TF-family classification (Phase 0.2).

Classifies a TF into ANY family, not just the two archetypes. Two complementary, cross-platform paths:

  HMM (primary)  -- pyhmmer hmmscan against a Pfam HMM library built from the 89 Pfam domains that
                    appear across the 531-TF operators_db. Domain-based, family-level, principled.
  blastp (fallback / corroboration) -- nearest neighbour against a broad reference of 307 TFs
                    (operators_db with sequences + the curated ArsR/SmtB & MerR set), spanning 22
                    families. Always available; links divergent families the MMseqs prefilter misses.

Both references/libraries are built once from the existing curated data (no manual curation). The Pfam
HMMs are fetched from the InterPro API. Family labels come from operators_db's own `protein_family`.

  build       python family_db.py build      # broad blastp reference + Pfam HMM library
  classify    python family_db.py classify <fasta|seq>
"""
from __future__ import annotations

import collections
import glob
import gzip
import json
import shutil
import sys
import tempfile
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from predictor import resources

from predictor.annotate.genome_resolver import is_dna, _blast, _run, _write_fasta

_REPO = Path(__file__).resolve().parents[2]

#: E-value at or below which a Pfam HMM hit, or the nearest reference by blastp, supports a family
#: call without a low-confidence flag.
FAMILY_EVALUE = 1e-3
_OPDB_TFS = _REPO.parent / "5.8 Promoters" / "operators_db" / "transcription_factors"
_CURATED_TFS = _REPO.parent / "5.8 Promoters" / "arsr_merr_families" / "transcription_factors"
_REFDIR = resources.REFS_DIR
_ALL_FAA = _REFDIR / "all_tfs.faa"
_ALL_JSON = _REFDIR / "all_tfs.json"

#: Shortest sequence allowed into the blastp nearest-neighbour reference. The smallest genuine
#: one-component regulators (ArsR/SmtB, CsoR) run ~90-110 aa, so anything below this is a fragment.
MIN_REFERENCE_LEN = 90
_PFAM_HMM = _REFDIR / "tf_pfam.hmm"
_PFAM_MAP = _REFDIR / "pfam_to_family.json"

_INTERPRO_HMM = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}?annotation=hmm"


@dataclass
class FamilyCall:
    family: str | None
    pfam_id: str | None = None
    pfam_name: str | None = None
    method: str = "none"            # 'hmm' | 'blastp' | 'structure' | 'none'
    confidence: float = 0.0         # [0,1]
    nearest_tf: str | None = None
    evidence: str = ""
    flags: list = field(default_factory=list)
    struct_family: str | None = None   # third axis (annotate.family_struct); None if no structure/foldseek
    struct_tm: float = 0.0             # query-normalised TM of the winning reference fold
    nearest_identity: float | None = None  # actual blastp identity; never a family-confidence surrogate
    nearest_evalue: float | None = None


# --------------------------------------------------------------------------- reference build
def _iter_tf_yamls():
    import yaml
    for src, root in (("operators_db", _OPDB_TFS), ("curated", _CURATED_TFS)):
        for f in sorted(glob.glob(str(root / "**" / "*.yaml"), recursive=True)):
            try:
                d = yaml.safe_load(open(f, encoding="utf-8"))
            except Exception:
                continue
            if isinstance(d, dict) and d.get("id"):
                yield src, d


def build_reference(out_faa: Path = _ALL_FAA, out_json: Path = _ALL_JSON) -> int:
    """Broad blastp reference (FASTA + metadata) from operators_db + curated TFs, spanning families."""
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    meta, records, seen = {}, [], set()
    skipped = []
    for src, d in _iter_tf_yamls():
        seq = ((d.get("sequence") or {}).get("protein") or "").replace("\n", "").strip()
        if not seq:
            continue
        tid = d["id"]
        if tid in seen:
            continue
        # Fragments cannot serve as a nearest-neighbour reference: they are too short to produce a
        # meaningful blastp alignment, and any hit they do win is spurious -- which then becomes a
        # confident-looking family call. Five curated records shared one 24-aa fragment; excluding them
        # here keeps the shipped reference honest even while the upstream yamls are incomplete.
        # `env/verify_tf_reference.py` reports every record filtered out this way.
        if len(seq) < MIN_REFERENCE_LEN:
            skipped.append((tid, len(seq)))
            continue
        seen.add(tid)
        cls = d.get("classification") or {}
        meta[tid] = {
            "family": cls.get("protein_family"),
            "family_group": cls.get("family_group"),
            "pfam": [(p.get("id") if isinstance(p, dict) else p)
                     for p in (cls.get("pfam") or []) if p],
            "uniprot": (d.get("identifiers") or {}).get("uniprot"),
            "organism": (d.get("organism") or {}).get("species"),
            "source": src,
        }
        records.append((tid, seq))
    # newline="\n" everywhere a vendored reference is written -- see tf_record.build_reference for why
    with out_faa.open("w", encoding="ascii", newline="\n") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")
    out_json.write_text(json.dumps(meta), encoding="utf-8", newline="\n")
    if skipped:
        print(f"  build_reference: skipped {len(skipped)} fragment record(s) "
              f"(<{MIN_REFERENCE_LEN} aa): {', '.join(f'{t}({n}aa)' for t, n in skipped)}")
    return len(records)


# --------------------------------------------------------------------------- Pfam HMM library
def _pfam_family_votes():
    """Collect {pfam_acc -> Counter(family)} and {pfam_acc -> name} from operators_db."""
    votes = collections.defaultdict(collections.Counter)
    names = {}
    for _src, d in _iter_tf_yamls():
        cls = d.get("classification") or {}
        fam = cls.get("protein_family")
        for p in (cls.get("pfam") or []):
            acc = p.get("id") if isinstance(p, dict) else p
            nm = p.get("name") if isinstance(p, dict) else None
            if acc and fam:
                votes[acc][fam] += 1
                if nm:
                    names[acc] = nm
    return votes, names


def build_pfam_hmm_library(out_hmm: Path = _PFAM_HMM, out_map: Path = _PFAM_MAP,
                           rebuild: bool = False) -> int:
    """Fetch the operators_db Pfam HMMs from InterPro into one HMM file; write pfam->family map."""
    out_hmm.parent.mkdir(parents=True, exist_ok=True)
    votes, names = _pfam_family_votes()
    fam_of = {acc: c.most_common(1)[0][0] for acc, c in votes.items()}
    if out_hmm.exists() and out_map.exists() and not rebuild:
        return sum(1 for _ in out_hmm.read_text(encoding="ascii", errors="replace").splitlines()
                   if _.startswith("NAME"))
    pmap, n_ok = {}, 0
    with out_hmm.open("w", encoding="ascii", errors="replace", newline="\n") as out:
        for acc in sorted(fam_of):
            try:
                raw = urllib.request.urlopen(_INTERPRO_HMM.format(acc=acc), timeout=30).read()
                txt = gzip.decompress(raw).decode("ascii", "replace")
            except Exception:
                continue
            out.write(txt if txt.endswith("\n") else txt + "\n")
            fam = fam_of[acc]
            pmap[acc] = fam
            if acc in names:
                pmap[names[acc]] = fam            # hmmscan reports the HMM NAME, not the accession
            n_ok += 1
    out_map.write_text(json.dumps(pmap), encoding="utf-8", newline="\n")
    return n_ok


# --------------------------------------------------------------------------- classifiers
def _to_protein(seq: str) -> str:
    if not is_dna(seq):
        return seq.replace("\n", "").strip()
    from Bio.Seq import Seq
    return str(Seq(seq.replace("\n", "").strip()).translate(to_stop=True))


def hmm_classify(prot: str):
    """pyhmmer hmmscan of one protein vs the Pfam HMM library -> best (acc/name, family, e, bits)."""
    if not _PFAM_HMM.exists() or not _PFAM_MAP.exists():
        return None
    try:
        import pyhmmer
    except ImportError:
        return None
    pmap = json.loads(_PFAM_MAP.read_text(encoding="utf-8"))
    alpha = pyhmmer.easel.Alphabet.amino()
    q = pyhmmer.easel.TextSequence(name=b"query", sequence=prot).digitize(alpha)
    with pyhmmer.plan7.HMMFile(str(_PFAM_HMM)) as hf:
        hmms = list(hf)
    def _s(x):
        return x.decode() if isinstance(x, (bytes, bytearray)) else (x or "")

    best = None
    for top in pyhmmer.hmmer.hmmscan([q], hmms):
        for hit in top:
            name = _s(hit.name)
            acc = _s(hit.accession).split(".")[0] or name
            fam = pmap.get(acc) or pmap.get(name)
            if fam and (best is None or hit.evalue < best[3]):
                best = (acc, name, fam, hit.evalue, hit.score)
    return best


def blastp_classify(prot: str, *, exclude_id=None, reference=None, evalue: float = 100.0):
    """Nearest curated/operators_db TF by blastp -> (nearest_id, family, pident, evalue)."""
    faa = reference or (str(_ALL_FAA) if _ALL_FAA.exists() else None)
    if faa is None:
        build_reference()
        faa = str(_ALL_FAA)
    meta = json.loads(_ALL_JSON.read_text(encoding="utf-8"))
    tmp = Path(tempfile.mkdtemp(prefix="famclass_"))
    try:
        shutil.copy(faa, tmp / "db.faa")
        _run([_blast("makeblastdb"), "-in", "db.faa", "-dbtype", "prot", "-out", "db"], cwd=tmp)
        _write_fasta(tmp / "q.faa", [("query", prot)])
        _run([_blast("blastp"), "-query", "q.faa", "-db", "db", "-evalue", evalue,
              "-outfmt", "6 sseqid pident evalue bitscore", "-max_target_seqs", 50,
              "-out", "h.tsv"], cwd=tmp)
        best = {}
        for line in (tmp / "h.tsv").read_text().splitlines():
            if not line.strip():
                continue
            sid, pid, ev, bits = line.split("\t")
            if sid == exclude_id:
                continue
            bits = float(bits)
            if sid not in best or bits > best[sid][2]:
                best[sid] = (float(pid) / 100.0, float(ev), bits)
        if not best:
            return None
        top = sorted(best.items(), key=lambda kv: -kv[1][2])[0]
        sid, (pid, ev, _b) = top
        return sid, meta.get(sid, {}).get("family"), pid, ev
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _structure_axis(cif, uniprot_acc):
    """Third evidence axis: fold similarity to the family-representative structures.

    `annotate.family_struct` is NOT part of this release -- it needs Foldseek and a curated
    reference-fold set, neither of which ships here -- so this axis is inert by design and returns
    None. The hook is kept because `FamilyCall` carries the struct fields and the two sequence axes
    are written to combine with it; installing the module enables the axis with no other change.

    Never raises: an absent module, an absent structure and an absent Foldseek are all normal.
    """
    if not (cif or uniprot_acc):
        return None
    try:
        from predictor.annotate import family_struct           # type: ignore[attr-defined]
    except Exception:
        return None
    try:
        sc = family_struct.structure_classify(cif, uniprot_acc=uniprot_acc)
        return sc if sc.available else None
    except Exception:
        return None


def classify_family(seq: str, *, exclude_id=None, prefer: str = "hmm",
                    cif: str | None = None, uniprot_acc: str | None = None) -> FamilyCall:
    """Holistic family call from three independent axes, in decreasing order of reliability:

      1. **Pfam HMM**   -- domain-level, principled, spans every family in `pfam_to_family.json`.
      2. **blastp**     -- nearest curated TF; catches families the HMM library misses.
      3. **structure**  -- fold similarity to a family-representative set (Foldseek). Requires both a
         structure (`cif`, or a `uniprot_acc` to pull from AFDB) and the optional
         `annotate.family_struct` module, which is NOT shipped in this release. Absent it, the call
         rests on axes 1 and 2 and every `struct_*` field is None.

    The call is **agnostic of the SSN cluster DB** -- none of the three axes is limited to the families
    that happen to have clusters (that is `ssn_clusters.FAMILIES`, a different, smaller set).

    Structure **corroborates but never overrides** a confident HMM call, for one measured reason: its
    reference set is CLOSED (12 folds), so a family outside it is silently forced onto the nearest fold
    (e.g. SoxS/AraC-XylS -> TetR/AcrR at qTM 0.35), whereas the HMM library knows AraC/XylS and calls it
    correctly. Structure's real job is the twilight zone -- when HMM and blastp both come up empty, a
    confident fold match (qTM >= 0.5) is the only remaining evidence, which is exactly the case Rondon
    et al. argue for.
    """
    prot = _to_protein(seq)
    hmm = hmm_classify(prot) if prefer == "hmm" else None
    # BLAST is corroborating/fallback evidence. A confident HMM call must remain usable when BLAST+
    # is not installed (or its subprocess fails); previously the unconditional call made the primary,
    # self-contained pyhmmer route fail before its result could be returned.
    bp = None
    bp_error = None
    try:
        bp = blastp_classify(prot, exclude_id=exclude_id)
    except Exception as exc:
        bp_error = f"{type(exc).__name__}: {exc}"
    sc = _structure_axis(cif, uniprot_acc)
    s_fam = sc.family if sc else None
    s_tm = sc.tm if sc else 0.0

    if hmm and hmm[3] < FAMILY_EVALUE:
        acc, name, fam, ev, bits = hmm
        conf = max(0.0, min(1.0, -__import__("math").log10(max(ev, 1e-300)) / 50.0))
        flags = [] if conf > 0.2 else ["weak HMM match"]
        ev_txt = f"Pfam {acc}/{name} e={ev:.1g}"
        if s_fam:
            if s_fam == fam:
                conf = min(1.0, conf + 0.1)                    # independent axes agree
                ev_txt += f"; structure agrees ({sc.target} qTM={s_tm:.2f})"
            else:
                flags.append(f"structure disagrees: fold nearest {s_fam} ({sc.target} qTM={s_tm:.2f}) "
                             f"-- HMM kept (structural reference set is closed)")
        if bp_error:
            flags.append("blastp corroboration unavailable; family call rests on Pfam HMM")
        return FamilyCall(fam, acc, name, "hmm", round(conf, 3),
                          bp[0] if bp else None, ev_txt, flags, s_fam, s_tm,
                          bp[2] if bp else None, bp[3] if bp else None)
    if bp:
        sid, fam, pid, ev = bp
        flags = [] if ev <= FAMILY_EVALUE else [f"low-confidence family (nearest {sid} at e={ev:.1g})"]
        ev_txt = f"nearest {sid} id={pid:.0%} e={ev:.1g}"
        conf = pid
        # a weak blastp call that the fold corroborates is worth more than either alone
        if s_fam and s_fam == fam:
            conf = min(1.0, conf + 0.15)
            ev_txt += f"; structure agrees ({sc.target} qTM={s_tm:.2f})"
        elif s_fam and sc.confident and ev > 1e-3:
            # blastp is weak AND structure is confident and says otherwise -> prefer the fold
            return FamilyCall(s_fam, None, None, "structure", round(s_tm, 3), sid,
                              f"fold {sc.target} qTM={s_tm:.2f} (blastp weak: {sid} e={ev:.1g})",
                              ["structure-led call: sequence evidence weak; reference set is closed",
                               f"blastp suggested {fam}"], s_fam, s_tm, pid, ev)
        return FamilyCall(fam, None, None, "blastp", round(conf, 3), sid, ev_txt, flags, s_fam, s_tm,
                          pid, ev)
    if s_fam and sc.confident:
        # the twilight-zone case: no HMM, no blastp -- the fold is the only evidence left
        return FamilyCall(s_fam, None, None, "structure", round(s_tm, 3), None,
                          f"fold {sc.target} qTM={s_tm:.2f} (no HMM/blastp evidence)",
                          ["structure-only call; the reference fold set is closed, so a family "
                           "outside it cannot be named"], s_fam, s_tm)
    flags = ["no family evidence (HMM + blastp both empty)"]
    if bp_error:
        flags.append("blastp unavailable; install BLAST+ to enable nearest-neighbour fallback")
    if sc and s_fam:
        flags.append(f"structure inconclusive: nearest {s_fam} qTM={s_tm:.2f} < {family_struct_tm_floor()}")
    return FamilyCall(None, flags=flags, struct_family=s_fam, struct_tm=s_tm)


def family_struct_tm_floor() -> float:
    """The qTM below which a fold match is inconclusive.

    0.5 is the fallback used whenever the optional `family_struct` module is absent, which is the
    normal case in this release. It only ever reaches a message string, never a decision.
    """
    try:
        from predictor.annotate import family_struct           # type: ignore[attr-defined]
        return family_struct.TM_CONFIDENT
    except Exception:
        return 0.5


# --------------------------------------------------------------------------- CLI / self-test
def _main(argv):
    if argv and argv[0] == "build":
        n = build_reference()
        print(f"broad reference: {n} TFs -> {_ALL_FAA.name}")
        m = build_pfam_hmm_library()
        print(f"Pfam HMM library: {m} HMMs -> {_PFAM_HMM.name}")
    elif argv and argv[0] == "classify":
        seq = argv[1]
        if Path(seq).exists():
            seq = "".join(l.strip() for l in Path(seq).read_text().splitlines()
                          if not l.startswith(">"))
        print(classify_family(seq))
    else:
        _demo()


def _demo() -> None:
    import yaml
    if not _ALL_FAA.exists():
        build_reference()
    have_hmm = _PFAM_HMM.exists()
    print(f"broad reference present; Pfam HMM library present: {have_hmm}")

    # holistic: a TetR-family TF must classify TetR/AcrR, NOT be forced into the two archetypes
    tetr = next((d for _s, d in _iter_tf_yamls()
                 if (d.get("classification") or {}).get("protein_family") == "TetR/AcrR"
                 and (d.get("sequence") or {}).get("protein")), None)
    if tetr is None:
        print("  note: no TetR/AcrR yaml with a sequence found in operators_db/curated -- skipping")
    else:
        fc = classify_family(tetr["sequence"]["protein"], exclude_id=tetr["id"])
        print(f"TetR TF {tetr['id']} -> family={fc.family} via {fc.method} ({fc.evidence})")
        assert fc.family and "TetR" in fc.family, f"expected TetR/AcrR, got {fc.family}"

    # archetypes still correct
    adh = yaml.safe_load(open(_CURATED_TFS / "merr" / "AdhR_Bsubtilis.yaml", encoding="utf-8"))
    fa = classify_family(adh["sequence"]["protein"], exclude_id="AdhR_Bsubtilis")
    print(f"AdhR -> family={fa.family} via {fa.method} ({fa.evidence})")
    assert fa.family == "MerR", f"expected MerR, got {fa.family}"

    lys = next((d for _s, d in _iter_tf_yamls()
                if (d.get("classification") or {}).get("protein_family") == "LysR-type (LTTR)"
                and (d.get("sequence") or {}).get("protein")), None)
    if lys is None:
        print("  note: no LysR-type (LTTR) yaml with a sequence found in operators_db/curated -- skipping")
    else:
        fl = classify_family(lys["sequence"]["protein"], exclude_id=lys["id"])
        print(f"LysR TF {lys['id']} -> family={fl.family} via {fl.method}")
        assert fl.family and "LysR" in fl.family, f"expected LysR, got {fl.family}"
    print("OK: holistic family classification across families (TetR, MerR, LysR).")


if __name__ == "__main__":
    _main(sys.argv[1:])
