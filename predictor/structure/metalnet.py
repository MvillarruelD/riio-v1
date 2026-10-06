"""
metalnet.py -- MetalNet2 as a FAMILY-AGNOSTIC metal-binding-site gate (the fifth inducer source).

Why this exists
---------------
`structure.metal_site.coordination_gate` encodes real metal-site chemistry, but only for the families we
wrote chemistry for: MerR, ArsR/SmtB, NikR, CsoR/FrmR, Rrf2 and Fur. For every other family it returns
`has_site=None` -- silent. Measured on the 150-candidate run that is **82 candidates (55 %)** with no
structural opinion at all: GntR (32), TetR (24), MarR (17), DtxR (4), CopY (3), LysR (2). Their inducer
call rests entirely on SSN cluster membership, genomic context, Ligify and regulon composition, with
nothing structural available to confirm or veto it.

MetalNet2 (Cheng et al.; https://github.com/wangchulab/MetalNet2) has no such family scope. It finds
metal sites from sequence + MSA alone: co-evolving Cys/His/Glu/Asp ("CHED") pairs from the MSA
Transformer's contact head, scored by an AutoGluon classifier over ESM-2 embeddings, then filtered to
connected constellations. That is evidence about *this protein*, independent of both our chemistry and
the SSN label -- which is exactly what the 82 are missing.

What this module is
-------------------
The pipeline-side adapter. The engine itself runs in its own conda env (`metalnet`) because MetalNet2's
classifiers are AutoGluon 0.8.0 pickles requiring Python 3.9 / numpy 1.25 / pandas 1.5.3 / sklearn 1.2.2,
which cannot be loaded by this interpreter. Same shape as the Folddisco bridge in
`structure.metal_site_search`: locate the engine, shell out, parse, and degrade gracefully when it is
absent (`available=False`, `has_site=None` -- never a fabricated verdict, and never an exception that
reaches the run).

Resolution of the engine python: `$METALNET_PYTHON` -> `<conda root>/envs/metalnet/python[.exe]` -> None.
Set up with `docs/ENGINES.md` (§ MetalNet2).

MSA
---
MetalNet needs an alignment, and the coordinate frame of every position it reports is the a3m's first
record -- so the wrong a3m is worse than none. Sources, in order, with the one used recorded in the
result so a reader can tell how the verdict was obtained:

  1. `results/metalnet/msa/<key>.a3m`  -- already fetched for this exact sequence;
  2. the AF3 server a3m on disk, when this TF was folded (`annotate.af3_msa._op_msa`) -- free, offline;
  3. ColabFold's UniRef a3m (`annotate.msa_homologs.colabfold_a3m(..., prefer="uniref")`) -- the same
     MSA MetalNet2 was trained and benchmarked on, so this is the faithful route; needs the network;
  4. a LOCAL a3m built by MMseqs2 against the vendored SSN family member DB. Fully offline, but it is a
     single-family alignment rather than a UniRef search, so it is shallower and less diverse than what
     the model saw in training. Recorded as `msa_source="local_family"`; treat those calls as weaker.

Results are cached per sequence under `results/metalnet/`, keyed by sequence + engine parameters, so a
re-run of the batch costs nothing.

  python -m predictor.structure.metalnet --status
  python -m predictor.structure.metalnet --seq <fasta|sequence> --family MerR
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]

RUNNER = resources.PACKAGE_DIR / "_engines" / "metalnet_runner.py"
ENGINE_DIR = Path(os.environ.get("METALNET2_DIR", _REPO / "env" / "tools" / "metalnet2"))
CACHE = resources.cache_path("metalnet")
MSA_CACHE = CACHE / "msa"

#: Bumped whenever the engine parameters or the summary schema change, so cached verdicts computed under
#: the old settings are recomputed instead of silently reused.
CACHE_VERSION = "mn2-v1"

#: Minimum alignment depth we are willing to run on. Below this the MSA Transformer's contact head has
#: nothing to co-evolve over and the CHED pairs are noise; upstream's own metalloproteome analysis
#: (`app/metalloproteome/analysis/analyze_msa_depth.ipynb`) reports the same sensitivity. We abstain
#: rather than emit a low-confidence verdict, because an abstention is handled correctly downstream and
#: a wrong metal call is not.
MIN_MSA_DEPTH = 16


# --------------------------------------------------------------------------- engine resolution
def engine_python() -> str | None:
    """The interpreter of the `metalnet` env: `$METALNET_PYTHON` -> conda envs/metalnet -> None."""
    env = os.environ.get("METALNET_PYTHON")
    if env and Path(env).exists():
        return env
    roots = [Path(sys.prefix).parent.parent, Path(sys.prefix).parent, Path.home() / "miniconda3",
             Path.home() / "anaconda3", Path.home() / "miniforge3"]
    for root in roots:
        for exe in (root / "envs" / "metalnet" / "python.exe", root / "envs" / "metalnet" / "bin" / "python"):
            if exe.exists():
                return str(exe)
    return None


@lru_cache(maxsize=1)
def status() -> dict:
    """Engine readiness. Pure inspection, never raises -- `tfop setup` and the pipeline both use it."""
    py = engine_python()
    out = {"engine": "MetalNet2", "python": py, "runner": str(RUNNER), "repo": str(ENGINE_DIR),
           "runner_ok": RUNNER.exists(), "repo_ok": ENGINE_DIR.is_dir(), "ready": False, "detail": None}
    if not py:
        out["detail"] = "no `metalnet` conda env (set $METALNET_PYTHON) -- see docs/ENGINES.md"
        return out
    if not RUNNER.exists():
        out["detail"] = f"runner missing: {RUNNER}"
        return out
    try:
        r = subprocess.run([py, str(RUNNER), "--status"], capture_output=True, text=True, timeout=300,
                           env={**os.environ, "KMP_DUPLICATE_LIB_OK": "TRUE"})
        info = json.loads(r.stdout or "{}")
        out.update({k: info[k] for k in ("model_1_ok", "model_2_ok", "packages") if k in info})
        out["ready"] = bool(info.get("ready"))
        if not out["ready"]:
            missing = [k for k, v in (info.get("packages") or {}).items() if str(v).startswith("MISSING")]
            out["detail"] = ("model weights not installed" if not info.get("model_1_ok")
                             else f"missing packages: {', '.join(missing)}" if missing else "not ready")
    except Exception as e:
        out["detail"] = f"{type(e).__name__}: {e}"
    return out


# --------------------------------------------------------------------------- MSA sourcing
def _key(seq: str) -> str:
    return hashlib.sha1(f"{CACHE_VERSION}|{(seq or '').upper()}".encode()).hexdigest()[:16]


def _clean_a3m(text: str) -> str:
    """ColabFold's a3m can carry NUL bytes (upstream strips them with `sed s/\\x00//g`). They are never
    meaningful, and left in place they survive into the aligned rows and shift every column."""
    return text.replace("\x00", "")


def _local_family_a3m(seq: str, family: str | None, out_path: Path, *, max_seqs: int = 512,
                      verbose: bool = False) -> Path | None:
    """Build an a3m from the VENDORED SSN family member DB with MMseqs2 -- the fully offline route.

    Only query columns are emitted (uppercase or '-'); insertions are dropped, which is lossless for this
    consumer because MetalNet's `parse_a3m_file` deletes lowercase columns anyway. Rows are de-duplicated
    and ordered best-e-value first; the engine then picks 64 maximally-diverse rows from them."""
    try:
        from predictor.annotate import ssn_clusters as _ssn
        from predictor.annotate.homologs import find_mmseqs, _run
    except Exception:
        return None
    tag = None
    try:
        tag = _ssn.family_tag(family) if family else None
    except Exception:
        tag = None
    db = resources.SSN_DATABASE / f"ssn_members_{tag}.fasta" if tag else None
    if db is None or not db.exists():
        if verbose:
            print(f"  metalnet: no vendored member DB for family {family!r} -> no local MSA")
        return None
    import tempfile
    seq = seq.upper()
    with tempfile.TemporaryDirectory(prefix="mn2_msa_") as td:
        tmp = Path(td)
        (tmp / "q.fasta").write_text(f">query\n{seq}\n", encoding="ascii")
        m8 = tmp / "hits.m8"
        try:
            _run(find_mmseqs(), ["easy-search", tmp / "q.fasta", db, m8, tmp / "t",
                                 "-s", 7.5, "--max-seqs", max_seqs * 2, "-e", 1e-3, "--threads", 4,
                                 "--format-output", "target,qstart,qaln,taln,evalue", "-v", 1], cwd=tmp)
        except Exception as e:
            if verbose:
                print(f"  metalnet: local MSA search failed ({type(e).__name__}: {e})")
            return None
        rows, seen = [], {seq}
        for line in m8.read_text(encoding="ascii", errors="replace").splitlines():
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            target, qstart, qaln, taln, evalue = parts[0], int(parts[1]), parts[2], parts[3], float(parts[4])
            row = ["-"] * len(seq)
            qi = qstart - 1                                  # mmseqs qstart is 1-based
            for qc, tc in zip(qaln, taln):
                if qc == "-":                                # insertion relative to the query -> dropped
                    continue
                if 0 <= qi < len(row):
                    row[qi] = tc.upper()
                qi += 1
            aligned = "".join(row)
            if aligned in seen:
                continue
            seen.add(aligned)
            rows.append((evalue, target, aligned))
            if len(rows) >= max_seqs:
                break
    if not rows:
        # Not hypothetical: the vendored TetR member DB holds 332 sequences (the other eleven families
        # hold 7k-75k) and E. coli AcrR hits none of them, so TetR candidates get NO offline alignment.
        # Say so rather than returning a bare None that reads downstream as "MetalNet has no opinion".
        if verbose:
            print(f"  metalnet: no MMseqs2 hits for this query in the vendored {tag} member DB "
                  f"({db.name}) -> no local MSA; this family needs the ColabFold route")
        return None
    rows.sort(key=lambda r: r[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write(f">query\n{seq}\n")
        for _e, target, aligned in rows:
            fh.write(f">{target}\n{aligned}\n")
    if verbose:
        print(f"  metalnet: local family MSA -> {out_path} ({len(rows)} rows, {tag})")
    return out_path


def a3m_for(seq: str, *, tf_id: str | None = None, family: str | None = None,
            allow_colabfold: bool = False, allow_local: bool = True,
            verbose: bool = False) -> tuple[Path | None, str | None]:
    """Resolve an alignment for `seq`. Returns (path, source) or (None, None)."""
    MSA_CACHE.mkdir(parents=True, exist_ok=True)
    cached = MSA_CACHE / f"{_key(seq)}.a3m"
    src_file = MSA_CACHE / f"{_key(seq)}.source"
    if cached.exists() and cached.stat().st_size > 0:
        try:
            return cached, json.loads(src_file.read_text(encoding="utf-8"))
        except Exception:
            return cached, "cache"

    def _remember(path: Path, source: str) -> tuple[Path, str]:
        try:
            cached.write_text(_clean_a3m(Path(path).read_text(encoding="utf-8", errors="replace")),
                              encoding="utf-8", newline="\n")
            src_file.write_text(json.dumps(source), encoding="utf-8")
            return cached, source
        except Exception:
            return Path(path), source

    if tf_id:
        try:
            from predictor.annotate.af3_msa import _op_msa
            hit = _op_msa(tf_id)
            if hit:
                return _remember(Path(hit), "af3")
        except Exception:
            pass
    if allow_colabfold:
        try:
            from predictor.annotate import msa_homologs as _mh
            got = _mh.colabfold_a3m(seq, MSA_CACHE / f"{_key(seq)}.colabfold.a3m", prefer="uniref",
                                    verbose=verbose)
            if got:
                return _remember(Path(got), "colabfold")
        except Exception as e:
            if verbose:
                print(f"  metalnet: ColabFold MSA failed ({type(e).__name__}: {e})")
    if allow_local:
        got = _local_family_a3m(seq, family, MSA_CACHE / f"{_key(seq)}.local.a3m", verbose=verbose)
        if got:
            return _remember(Path(got), "local_family")
    return None, None


# --------------------------------------------------------------------------- prediction
def _unavailable(reason: str) -> dict:
    """The one shape a caller must handle when the engine cannot speak. `has_site=None` means 'no
    opinion' and is treated exactly like a silent coordination gate downstream -- never as a negative."""
    return {"available": False, "has_site": None, "reason": reason, "engine": "MetalNet2"}


def predict_site(seq: str, *, tf_id: str | None = None, family: str | None = None,
                 allow_colabfold: bool = False, allow_local: bool = True, use_cache: bool = True,
                 timeout: int = 3600, verbose: bool = False) -> dict:
    """MetalNet2's verdict for one sequence. Never raises: any failure returns `has_site=None` with a
    reason, so a run without the engine behaves exactly as it did before this source existed."""
    seq = (seq or "").upper()
    if not seq:
        return _unavailable("empty sequence")
    CACHE.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE / f"{_key(seq)}.json"
    if use_cache and cache_file.exists():
        try:
            return json.loads(cache_file.read_text(encoding="utf-8"))
        except Exception:
            pass

    st = status()
    if not st.get("ready"):
        return _unavailable(f"engine unavailable: {st.get('detail') or 'not ready'}")

    a3m, source = a3m_for(seq, tf_id=tf_id, family=family, allow_colabfold=allow_colabfold,
                          allow_local=allow_local, verbose=verbose)
    if a3m is None:
        return _unavailable("no MSA available (needs a cached/AF3 a3m, ColabFold, or a vendored "
                            "family member DB)")

    import tempfile
    with tempfile.TemporaryDirectory(prefix="mn2_run_") as td:
        work = Path(td)
        sid = f"mn_{_key(seq)}"
        (work / "jobs.json").write_text(
            json.dumps([{"seq_id": sid, "seq": seq, "a3m": str(a3m)}]), encoding="utf-8")
        out = work / "out.json"
        try:
            r = subprocess.run([st["python"], str(RUNNER), "--jobs", str(work / "jobs.json"),
                                "--out", str(out), "--work", str(work / "scratch"), "--quiet"],
                               capture_output=True, text=True, timeout=timeout,
                               env={**os.environ, "KMP_DUPLICATE_LIB_OK": "TRUE"})
        except Exception as e:
            return _unavailable(f"engine call failed ({type(e).__name__}: {e})")
        if not out.exists():
            tail = (r.stderr or r.stdout or "")[-800:]
            return _unavailable(f"engine produced no output (rc={r.returncode}): {tail}")
        payload = json.loads(out.read_text(encoding="utf-8"))

    res = payload["results"].get(sid, {})
    got = {"available": True, "engine": "MetalNet2", "msa_source": source,
           "params": payload.get("params"), **res}
    # Shallow alignments are not evidence in either direction -- abstain rather than assert.
    if res.get("ok") and int(res.get("n_msa_used") or 0) < MIN_MSA_DEPTH:
        got["has_site"] = None
        got["reason"] = (f"MSA too shallow ({res.get('n_msa_used')} rows < {MIN_MSA_DEPTH}) -> "
                         f"co-evolution unusable; abstaining")
    elif not res.get("ok"):
        got["has_site"] = None
        got["reason"] = res.get("reason") or "engine could not score this sequence"
    if use_cache:
        try:
            cache_file.write_text(json.dumps(got, indent=2), encoding="utf-8", newline="\n")
        except Exception:
            pass
    return got


def predict_sites(items, *, allow_colabfold: bool = False, allow_local: bool = True,
                  use_cache: bool = True, timeout: int | None = None, verbose: bool = False) -> list:
    """Batch form of `predict_site` -- ONE engine process for many sequences.

    Worth using for anything above a couple of TFs: a single call spends ~a minute loading ESM-2 650M
    and the AutoGluon bag before it scores anything, and that cost is per PROCESS, not per sequence. The
    150-candidate run is the reason this exists.

    `items` are dicts with `seq` and optional `tf_id` / `family` (a bare string is accepted as `seq`).
    Returns a list of results aligned with `items`. Cached sequences are served from disk and never sent
    to the engine, so re-running a partially-completed batch only does the remaining work."""
    norm = [{"seq": (it if isinstance(it, str) else it.get("seq") or "").upper(),
             "tf_id": None if isinstance(it, str) else it.get("tf_id"),
             "family": None if isinstance(it, str) else it.get("family")} for it in items]
    out: list = [None] * len(norm)
    CACHE.mkdir(parents=True, exist_ok=True)

    todo = []
    for i, it in enumerate(norm):
        if not it["seq"]:
            out[i] = _unavailable("empty sequence")
            continue
        cf = CACHE / f"{_key(it['seq'])}.json"
        if use_cache and cf.exists():
            try:
                out[i] = json.loads(cf.read_text(encoding="utf-8"))
                continue
            except Exception:
                pass
        todo.append(i)
    if not todo:
        return out

    st = status()
    if not st.get("ready"):
        for i in todo:
            out[i] = _unavailable(f"engine unavailable: {st.get('detail') or 'not ready'}")
        return out

    jobs, sources, index = [], {}, {}
    for i in todo:
        it = norm[i]
        a3m, source = a3m_for(it["seq"], tf_id=it["tf_id"], family=it["family"],
                              allow_colabfold=allow_colabfold, allow_local=allow_local, verbose=verbose)
        if a3m is None:
            out[i] = _unavailable("no MSA available (needs a cached/AF3 a3m, ColabFold, or a vendored "
                                  "family member DB)")
            continue
        sid = f"mn_{_key(it['seq'])}"
        index[sid], sources[sid] = i, source
        jobs.append({"seq_id": sid, "seq": it["seq"], "a3m": str(a3m)})
    if not jobs:
        return out

    import tempfile
    with tempfile.TemporaryDirectory(prefix="mn2_batch_") as td:
        work = Path(td)
        (work / "jobs.json").write_text(json.dumps(jobs), encoding="utf-8")
        res_path = work / "out.json"
        cmd = [st["python"], str(RUNNER), "--jobs", str(work / "jobs.json"), "--out", str(res_path),
               "--work", str(work / "scratch")] + ([] if verbose else ["--quiet"])
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=timeout or (600 + 120 * len(jobs)),
                               env={**os.environ, "KMP_DUPLICATE_LIB_OK": "TRUE"})
        except Exception as e:
            for sid, i in index.items():
                out[i] = _unavailable(f"engine call failed ({type(e).__name__}: {e})")
            return out
        if not res_path.exists():
            tail = (r.stderr or r.stdout or "")[-800:]
            for sid, i in index.items():
                out[i] = _unavailable(f"engine produced no output (rc={r.returncode}): {tail}")
            return out
        payload = json.loads(res_path.read_text(encoding="utf-8"))

    for sid, i in index.items():
        res = payload["results"].get(sid, {})
        got = {"available": True, "engine": "MetalNet2", "msa_source": sources[sid],
               "params": payload.get("params"), **res}
        if res.get("ok") and int(res.get("n_msa_used") or 0) < MIN_MSA_DEPTH:
            got["has_site"] = None
            got["reason"] = (f"MSA too shallow ({res.get('n_msa_used')} rows < {MIN_MSA_DEPTH}) -> "
                             f"co-evolution unusable; abstaining")
        elif not res.get("ok"):
            got["has_site"] = None
            got["reason"] = res.get("reason") or "engine could not score this sequence"
        out[i] = got
        if use_cache:
            try:
                (CACHE / f"{_key(norm[i]['seq'])}.json").write_text(
                    json.dumps(got, indent=2), encoding="utf-8", newline="\n")
            except Exception:
                pass
    return out


def summary(res: dict | None) -> str:
    """One-line human-readable form for logs and the dossier."""
    if not res or not res.get("available"):
        return f"MetalNet: unavailable ({(res or {}).get('reason', 'not run')})"
    if res.get("has_site") is None:
        return f"MetalNet: no opinion ({res.get('reason')})"
    if res.get("has_site"):
        return (f"MetalNet: metal site -- {res.get('n_site_residues')} residues "
                f"({', '.join(res.get('site_residues') or [])}), {res.get('n_site_pairs')} pairs, "
                f"type suggestion {res.get('metal_type_top') or '-'} [advisory], msa={res.get('msa_source')}")
    return (f"MetalNet: no metal site ({res.get('n_coevo_pairs')} CHED pairs tested, "
            f"max prob {res.get('max_prob')}), msa={res.get('msa_source')}")


def _self_test() -> None:
    """Offline + engine-agnostic. Pins the contract every caller depends on: the adapter never raises,
    an absent engine abstains rather than denies, and the cache key is sequence-scoped."""
    st = status()
    print(f"engine: ready={st['ready']} python={st.get('python')} ({st.get('detail') or 'ok'})")

    # the keys any caller may rely on when the engine cannot speak
    assert {"available", "has_site", "reason", "engine"} <= set(_unavailable("x")), \
        "the unavailable shape is part of the contract"
    assert _unavailable("x")["has_site"] is None, "unavailable must ABSTAIN, never deny"
    assert predict_site("")["has_site"] is None, "an empty sequence must not raise"

    assert _key("MKV") == _key("mkv"), "the key is case-insensitive"
    assert _key("MKV") != _key("MKVA"), "the key must separate different sequences"
    assert _clean_a3m(">q\nMK\x00V\n") == ">q\nMKV\n", "NUL bytes shift every column -- must be stripped"

    from predictor.structure import metal_site as _ms
    assert _ms.metalnet_verdict({"available": True, "has_site": True}) is True
    assert _ms.metalnet_verdict({"available": False, "has_site": False}) is None, \
        "an unavailable engine is silence, not a negative"

    for state, want in ((_unavailable("no engine"), "unavailable"),
                        ({"available": True, "has_site": None, "reason": "shallow"}, "no opinion"),
                        ({"available": True, "has_site": True, "n_site_residues": 4,
                          "site_residues": ["C1"], "n_site_pairs": 2}, "metal site"),
                        ({"available": True, "has_site": False, "n_coevo_pairs": 3,
                          "max_prob": 0.1}, "no metal site")):
        assert want in summary(state), f"summary({state}) should mention {want!r}"
    print("OK: MetalNet adapter abstains cleanly when the engine is absent and never raises.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="MetalNet2 metal-binding-site gate")
    ap.add_argument("--seq", help="protein sequence or a fasta path")
    ap.add_argument("--family", default=None)
    ap.add_argument("--tf-id", default=None)
    ap.add_argument("--colabfold", action="store_true", help="allow the ColabFold MSA fetch (network)")
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--batch", help="fasta of TFs -> warm the cache in ONE engine process")
    ap.add_argument("--family-of", help="json {tf_id: family} for --batch (optional; the family only "
                                        "selects the vendored member DB for the offline MSA route)")
    a = ap.parse_args(argv)
    if a.self_test:
        _self_test()
        return 0
    if a.batch:
        fams = json.loads(Path(a.family_of).read_text(encoding="utf-8")) if a.family_of else {}
        items, name, buf = [], None, []
        for line in Path(a.batch).read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                if name:
                    items.append({"seq": "".join(buf), "tf_id": name, "family": fams.get(name)})
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if name:
            items.append({"seq": "".join(buf), "tf_id": name, "family": fams.get(name)})
        res = predict_sites(items, allow_colabfold=a.colabfold, use_cache=not a.no_cache, verbose=True)
        for it, r in zip(items, res):
            print(f"{it['tf_id']:<28}{summary(r)}")
        n = sum(1 for r in res if r.get("has_site"))
        print(f"\n{n}/{len(res)} with a predicted metal site "
              f"({sum(1 for r in res if r.get('has_site') is None)} abstentions)")
        return 0
    if a.status or not a.seq:
        print(json.dumps(status(), indent=2))
        return 0
    s = a.seq
    if Path(s).exists():
        s = "".join(l.strip() for l in Path(s).read_text().splitlines() if not l.startswith(">"))
    res = predict_site(s, tf_id=a.tf_id, family=a.family, allow_colabfold=a.colabfold,
                       use_cache=not a.no_cache, verbose=True)
    print(json.dumps(res, indent=2))
    print(summary(res))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
