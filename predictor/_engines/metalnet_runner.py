"""
metalnet_runner.py -- the ENGINE side of the MetalNet2 metal-binding-site gate.

This script runs INSIDE the dedicated `metalnet` conda env, never in the pipeline's own interpreter.
MetalNet2 ships its two classifiers as AutoGluon 0.8.0 stacked ensembles pickled under Python 3.9 with
numpy 1.25 / pandas 1.5.3 / scikit-learn 1.2.2; those pickles cannot be loaded by the pipeline env
(3.11 / numpy 2.x), so the engine is bridged as a subprocess exactly like Folddisco is bridged through
WSL. `predictor/structure/metalnet.py` is the pipeline-side adapter; this file is what it execs.

The science is UPSTREAM's, not ours: every step imports the corresponding function from the vendored
`env/tools/metalnet2/` checkout rather than reimplementing it, so a MetalNet2 update changes our answer
too. What this file adds is only orchestration:

  * ONE process for a whole batch. Upstream's `run_prediction_workflow.py` shells out four times per
    invocation and reloads the MSA Transformer (100M) and ESM-2 650M each time; folding the batch into a
    single process makes a 150-TF run a single pair of model loads.
  * no `os.system` / `rm -rf` -- upstream's workflow is POSIX-only and this project runs on Windows.
  * a strict query/MSA agreement check. Every downstream position is an index into the a3m's first
    record, so an a3m built for a DIFFERENT protein yields confident predictions at meaningless
    positions. `--jobs` carries the query sequence and the runner refuses any job whose a3m query row
    disagrees with it.
  * a JSON summary per sequence instead of four intermediate TSVs.

Pipeline: a3m -> hamming-filtered 64 seqs -> MSA-Transformer coevolution -> CHED pairs above the
coevolution threshold -> ESM-2 650M per-residue embeddings -> AutoGluon model 1 (metal-binding pair
probability) -> graph filter -> [AutoGluon model 2: metal type]. Defaults are upstream's tuned values
(`model/train/cmd.sh`): num_seq 64, coevo 0.1, prob 0.8, node-num 4, esm2/avg.

  python metalnet_runner.py --status
  python metalnet_runner.py --jobs jobs.json --out out.json --work <scratch dir>

`jobs.json` is [{"seq_id": str, "seq": str, "a3m": path}, ...].
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from argparse import Namespace
from collections import Counter
from pathlib import Path

_HERE = Path(__file__).resolve().parent
#: The vendored upstream checkout. Overridable so a user can point at their own clone.
MN2 = Path(__file__).resolve().parents[2] / "env" / "tools" / "metalnet2"

# upstream's tuned production settings -- model/train/cmd.sh + predict_pairs.py defaults
DEFAULTS = {"num_seq": 64, "msa_filter_type": "hamming", "coevo_threshold": 0.1,
            "prob_threshold": 0.8, "node_num_threshold": 4,
            "encode_method": "esm2", "encode_strategy": "avg"}
MAX_CHAIN_LENGTH = 1023           # extract_coevo_pairs.max_chain_length (MSA Transformer column limit)


def _mn2_dir() -> Path:
    import os
    return Path(os.environ.get("METALNET2_DIR") or MN2)


def _add_upstream_to_path() -> Path:
    """Put `model/scripts/` on sys.path so upstream's `utils` package resolves, and return the dir."""
    scripts = _mn2_dir() / "model" / "scripts"
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    return scripts


class _posix_path_shim:
    """Let a Linux-trained fastai learner unpickle on Windows.

    MetalNet2's model 1 is an AutoGluon bag that includes `NeuralNetFastAI_BAG_L1`, and a fastai
    `Learner` pickles its `path` as a `pathlib.PosixPath`. Unpickling one on Windows raises
    `NotImplementedError: cannot instantiate 'PosixPath' on your system` before any of the model weights
    are read, so the whole ensemble fails to load. Mapping PosixPath onto WindowsPath for the duration of
    the load is the standard fastai remedy; the path is a training-time artefact that is never used at
    inference, and everything else in the pickle is numeric. Scoped to the load and to Windows so nothing
    else in the process sees a mangled `pathlib`.
    """

    def __enter__(self):
        import pathlib
        self._saved = None
        if sys.platform.startswith("win"):
            self._saved = pathlib.PosixPath
            pathlib.PosixPath = pathlib.WindowsPath
        return self

    def __exit__(self, *exc):
        import pathlib
        if self._saved is not None:
            pathlib.PosixPath = self._saved
        return False


def model_paths() -> dict:
    d = _mn2_dir()
    return {"repo": d,
            "model_1": d / "model" / "train" / "MetalNet_AutogluonModels",
            "model_2": d / "extra" / "train" / "MetalNet_metal_type_AutogluonModels"}


# --------------------------------------------------------------------------- status
def status() -> dict:
    """What this engine can do here. Never raises -- the adapter reports whatever comes back."""
    p = model_paths()
    out = {"python": sys.version.split()[0], "repo": str(p["repo"]), "repo_ok": p["repo"].is_dir(),
           "model_1_ok": (p["model_1"] / "predictor.pkl").exists(),
           "model_2_ok": (p["model_2"] / "predictor.pkl").exists(),
           "packages": {}, "ready": False}
    for mod in ("numpy", "pandas", "sklearn", "torch", "esm", "networkx", "autogluon.tabular"):
        try:
            out["packages"][mod] = getattr(__import__(mod, fromlist=["__version__"]),
                                           "__version__", "installed")
        except Exception as e:
            out["packages"][mod] = f"MISSING ({type(e).__name__})"
    out["ready"] = bool(out["repo_ok"] and out["model_1_ok"]
                        and not any(str(v).startswith("MISSING") for v in out["packages"].values()))
    return out


# --------------------------------------------------------------------------- stage 1: coevolution
def _coevo_pairs(jobs, params, *, cuda=None, verbose=True):
    """Upstream `extract_coevo_pairs`: filter each a3m to `num_seq` maximally-diverse rows, run the MSA
    Transformer's contact head, keep CHED pairs above `coevo_threshold`. Returns (DataFrame, per-job info).

    Failures are per-job: a job whose MSA is missing, disagrees with its query, or is too long is dropped
    with a reason and the rest of the batch still runs."""
    import pandas as pd
    _add_upstream_to_path()
    import extract_coevo_pairs as ecp                                    # upstream, unmodified

    transformer, batch_converter = ecp.load_msa_transformer(cuda)
    use_gpu = cuda is not None
    frames, info = [], {}
    for job in jobs:
        sid, seq, a3m = job["seq_id"], job["seq"].upper(), job.get("a3m")
        rec = {"ok": False, "reason": None, "n_msa_total": 0, "n_msa_used": 0}
        info[sid] = rec
        try:
            if not a3m or not Path(a3m).exists():
                rec["reason"] = "no a3m"
                continue
            full = ecp.parse_a3m_file(a3m)
            rec["n_msa_total"] = len(full)
            if not full:
                rec["reason"] = "empty a3m"
                continue
            # The a3m's FIRST record defines the coordinate frame for every position reported below.
            # An a3m built for another protein would still produce confident-looking pairs, at positions
            # that mean nothing -- so disagreement is a hard error, not a warning.
            query_row = full[0].replace("-", "").upper()
            if query_row != seq:
                rec["reason"] = (f"a3m query row does not match the submitted sequence "
                                 f"(a3m {len(query_row)} aa, query {len(seq)} aa)")
                continue
            if len(seq) > MAX_CHAIN_LENGTH:
                rec["reason"] = f"sequence too long for the MSA Transformer ({len(seq)} > {MAX_CHAIN_LENGTH})"
                continue
            aln = ecp.filter_msa(a3m, params["num_seq"], params["msa_filter_type"])
            rec["n_msa_used"] = len(aln)
            mtx = ecp.get_coevo_matrix(aln, transformer, batch_converter, use_gpu)
            df = ecp.get_coevo_pairs(aln[0], mtx, params["coevo_threshold"])
            df["seq_id"] = sid
            frames.append(df)
            rec["ok"] = True
            if verbose:
                print(f"  [coevo] {sid}: {rec['n_msa_used']}/{rec['n_msa_total']} msa rows -> "
                      f"{len(df)} CHED pairs", flush=True)
        except Exception as e:                                            # one bad job must not kill the batch
            rec["reason"] = f"{type(e).__name__}: {e}"
            if verbose:
                traceback.print_exc()
    if frames:
        out = pd.concat(frames, ignore_index=True)
        out["resi_seq_posi_1"] = out["resi_seq_posi_1"].astype(int)
        out["resi_seq_posi_2"] = out["resi_seq_posi_2"].astype(int)
    else:
        out = pd.DataFrame(columns=["resi_1", "resi_seq_posi_1", "resi_2", "resi_seq_posi_2",
                                    "coevo_value", "seq_id"])
    return out, info


# --------------------------------------------------------------------------- stage 2: ESM-2 embeddings
def _esm2_encodings(jobs, work: Path, *, cuda=None, esm2_model=None) -> dict:
    """Upstream `utils/esm2.py` (fair-esm's extract.py) -> one `<seq_id>.pt` of per-residue layer-33
    representations per sequence. Returns {seq_id: path}. Called once for the whole batch.

    Loaded by path rather than imported: `utils/esm2.py` is a script, and importing it as `utils.esm2`
    would shadow the `esm` package name inside upstream's own `utils` package on some layouts."""
    import importlib.util
    scripts = _add_upstream_to_path()
    spec = importlib.util.spec_from_file_location("_mn2_esm2", scripts / "utils" / "esm2.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    fasta = work / "batch.fasta"
    fasta.write_text("".join(f">{j['seq_id']}\n{j['seq'].upper()}\n" for j in jobs), encoding="utf-8")
    out_dir = work / "esm2"
    mod.run(Namespace(model_location=esm2_model or _maybe_download_esm2(),
                      fasta_file=fasta, output_dir=out_dir,
                      toks_per_batch=4096, repr_layers=[-1], include=["per_tok"],
                      truncation_seq_length=2700, nogpu=cuda is None))
    return {j["seq_id"]: str((out_dir / f"{j['seq_id']}.pt").resolve())
            for j in jobs if (out_dir / f"{j['seq_id']}.pt").exists()}


def _maybe_download_esm2(model_name: str = "esm2_t33_650M_UR50D") -> str:
    """Upstream `run_prediction_workflow.maybe_download_esm2`, verbatim in effect: the fair-esm torch-hub
    checkpoint, fetched once into the torch cache and reused offline afterwards."""
    import torch
    from esm import pretrained
    path = Path(torch.hub.get_dir()) / "checkpoints" / f"{model_name}.pt"
    if not path.exists():
        pretrained._download_model_and_regression_data(model_name)
    return str(path)


# --------------------------------------------------------------------------- stage 3 + 4: the classifiers
def _predict_pairs(df_pairs, encodings: dict, params: dict):
    """Upstream `predict_pairs.predict`: AutoGluon model 1 gives each CHED pair a metal-binding
    probability, then the graph filter keeps only pairs that sit in a connected constellation (leaf
    residues removed; a rescue path when the surviving graph is smaller than `node_num_threshold`)."""
    _add_upstream_to_path()
    import predict_pairs as pp                                            # upstream, unmodified
    from autogluon.tabular import TabularPredictor
    from utils.encode_utils import encode                                 # upstream, unmodified

    enc = encode(df=df_pairs, dict_file=encodings, encode_method=params["encode_method"],
                 strategy=params["encode_strategy"])
    if enc.empty:
        return enc
    with _posix_path_shim():
        model = TabularPredictor.load(str(model_paths()["model_1"]))
        return pp.predict(predictor=model, df_encoded_pairs=enc,
                          prob_threshold=params["prob_threshold"],
                          node_num_threshold=params["node_num_threshold"])


def _predict_metal_type(df_pos, encodings: dict, params: dict):
    """Upstream `extra/scripts/predict_metal_type.predict` on the pairs the graph filter kept. The label
    is a PDB chemical-component code (ZN, CU, FE, SF4, ...) and can be multi-valued ("ZN;MN").

    Deliberately advisory downstream: the pipeline uses MetalNet for WHETHER a site exists, not WHICH
    metal it binds -- and CA/MG in particular are usually structural rather than sensed."""
    import numpy as np
    import pandas as pd
    _add_upstream_to_path()
    from autogluon.tabular import TabularPredictor
    from utils.encode_utils import encode

    enc = encode(df=df_pos, dict_file=encodings, encode_method=params["encode_method"],
                 strategy=params["encode_strategy"])
    if enc.empty:
        return pd.DataFrame(columns=["seq_id", "resi_seq_posi_1", "resi_seq_posi_2", "pred"])
    with _posix_path_shim():
        model = TabularPredictor.load(str(model_paths()["model_2"]))
        data_x = pd.DataFrame(np.stack(enc["x"]))
        out = enc.drop(columns=["x"]).copy().reset_index()
        out["pred"] = model.predict(pd.DataFrame(data_x), as_pandas=False)
    return out


# --------------------------------------------------------------------------- summary
def _summarise(sid: str, df_seq, coevo_info: dict, type_rows) -> dict:
    """One sequence's verdict. `has_site` is upstream's own definition of a predicted metal-binding
    protein (`app/metalloproteome/cmd.sh` keeps the rows with `filter_by_graph == 1`)."""
    pos = df_seq[df_seq["filter_by_graph"] == 1] if len(df_seq) else df_seq
    residues, pairs = {}, []
    for _, r in (pos.iterrows() if len(pos) else []):
        for aa, p in ((r["resi_1"], int(r["resi_seq_posi_1"])), (r["resi_2"], int(r["resi_seq_posi_2"]))):
            residues[p] = aa
        pairs.append({"resi_1": r["resi_1"], "posi_1": int(r["resi_seq_posi_1"]),
                      "resi_2": r["resi_2"], "posi_2": int(r["resi_seq_posi_2"]),
                      "coevo": round(float(r["coevo_value"]), 4), "prob": round(float(r["prob"]), 4)})
    votes = Counter()
    for lab in (type_rows or []):
        for part in str(lab).split(";"):
            if part and part != "nan":
                votes[part] += 1
    return {
        "ok": bool(coevo_info.get("ok")),
        "reason": coevo_info.get("reason"),
        "n_msa_total": coevo_info.get("n_msa_total", 0),
        "n_msa_used": coevo_info.get("n_msa_used", 0),
        "n_coevo_pairs": int(len(df_seq)),
        "n_site_pairs": int(len(pos)),
        "has_site": bool(len(pos) > 0),
        "max_prob": round(float(df_seq["prob"].max()), 4) if len(df_seq) else 0.0,
        "site_residues": [f"{residues[p]}{p + 1}" for p in sorted(residues)],   # 1-based, as upstream plots
        "n_site_residues": len(residues),
        "site_pairs": pairs,
        "metal_type_votes": dict(votes.most_common()),
        "metal_type_top": (votes.most_common(1)[0][0] if votes else None),
    }


def predict(jobs, work: Path, *, params=None, cuda=None, metal_type=True, verbose=True) -> dict:
    import pandas as pd
    params = {**DEFAULTS, **(params or {})}
    work.mkdir(parents=True, exist_ok=True)

    empty = pd.DataFrame(columns=["resi_1", "resi_seq_posi_1", "resi_2", "resi_seq_posi_2",
                                  "coevo_value", "seq_id", "prob", "filter_by_graph"])
    df_pairs, coevo_info = _coevo_pairs(jobs, params, cuda=cuda, verbose=verbose)
    runnable = [j for j in jobs if coevo_info.get(j["seq_id"], {}).get("ok")]
    results = {}
    if not runnable or df_pairs.empty:
        for j in jobs:
            results[j["seq_id"]] = _summarise(j["seq_id"], empty, coevo_info.get(j["seq_id"], {}), None)
        return {"params": params, "results": results}

    encodings = _esm2_encodings(runnable, work, cuda=cuda)
    df_pairs = df_pairs[df_pairs["seq_id"].isin(encodings)]
    df_out = _predict_pairs(df_pairs, encodings, params)
    if not len(df_out):
        df_out = empty

    type_by_seq = {}
    if metal_type and len(df_out) and model_paths()["model_2"].exists():
        pos = df_out[df_out["filter_by_graph"] == 1]
        if len(pos):
            try:
                t = _predict_metal_type(pos, encodings, params)
                for sid, g in t.groupby("seq_id"):
                    type_by_seq[sid] = list(g["pred"])
            except Exception:
                if verbose:
                    traceback.print_exc()                   # advisory only -- never fails the gate

    for j in jobs:
        sid = j["seq_id"]
        results[sid] = _summarise(sid, df_out[df_out["seq_id"] == sid],
                                  coevo_info.get(sid, {}), type_by_seq.get(sid))
    return {"params": params, "results": results}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs", help="json: [{seq_id, seq, a3m}, ...]")
    ap.add_argument("--out", help="json output path")
    ap.add_argument("--work", help="scratch dir for the fasta + esm2 tensors")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--cuda", type=int, default=None)
    ap.add_argument("--no-metal-type", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_', '-')}", type=type(v), default=v)
    a = ap.parse_args(argv)

    if a.status or not a.jobs:
        print(json.dumps(status(), indent=2))
        return 0
    jobs = json.loads(Path(a.jobs).read_text(encoding="utf-8"))
    work = Path(a.work) if a.work else Path(a.out).resolve().parent / "_metalnet_work"
    res = predict(jobs, work, params={k: getattr(a, k) for k in DEFAULTS},
                  cuda=a.cuda, metal_type=not a.no_metal_type, verbose=not a.quiet)
    Path(a.out).write_text(json.dumps(res, indent=2), encoding="utf-8")
    if not a.quiet:
        n_ok = sum(1 for r in res["results"].values() if r["ok"])
        n_site = sum(1 for r in res["results"].values() if r["has_site"])
        print(f"metalnet: {n_ok}/{len(jobs)} scored, {n_site} with a predicted metal site -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
