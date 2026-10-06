"""Integrity of the vendored TF reference — the blastp nearest-neighbour database.

`predictor/data/refs/all_tfs.faa` classifies a novel TF whenever the Pfam HMM has no confident hit. A wrong
sequence there fails silently: the tool still answers, it just answers with the family of whatever the
corrupt record happens to resemble.

The shipped reference had eight such records. Five distinct TFs shared one 24-aa fragment; SoxR_Ecoli and
ZntR_Ecoli shared a single 107-aa sequence that was neither protein; TipA_Slividans was a 73-aa fragment
of a 253-aa protein. These checks are offline and cheap; `env/verify_tf_reference.py` adds the
authoritative UniProt comparison, which needs network and so is not run here.
"""
import json
from collections import defaultdict
from pathlib import Path

import pytest

from predictor.annotate import family_db as FD

REPO = Path(__file__).resolve().parents[1]


def _read_fasta(p: Path) -> dict[str, str]:
    seqs, cur = {}, None
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        ln = ln.strip()
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += ln
    return seqs


@pytest.fixture(scope="module")
def reference():
    faa = FD._ALL_FAA
    if not faa.is_file():
        pytest.skip("broad reference not built")
    meta = json.loads(FD._ALL_JSON.read_text(encoding="utf-8")) if FD._ALL_JSON.is_file() else {}
    return _read_fasta(faa), meta


@pytest.mark.released_data
def test_reference_is_populated(reference):
    seqs, _ = reference
    assert len(seqs) > 250, f"reference unexpectedly small ({len(seqs)} sequences)"


@pytest.mark.released_data
def test_operator_index_is_shipped_and_populated():
    """The vendored operator index must ship with content, not as an empty stub."""
    from predictor.annotate import known_operators as KO
    assert KO._INDEX.is_file(), "operators_by_family.json is not shipped"
    idx = json.loads(KO._INDEX.read_text(encoding="utf-8"))
    assert idx.get("by_family"), "the vendored operator index is empty"


def test_rebuilding_the_index_cannot_empty_it(tmp_path):
    """`build_index()` must never replace a good index with an empty rebuild.

    Its source YAMLs (`_TF_DIRS`) live outside the repo, so on any plain install they are absent and the
    rebuild yields {}. Writing that destroyed the 212 KB vendored index -- measured on a fresh
    `pip install` from GitHub, where merely running the self-test shrank it to 30 bytes. That made
    self-verification corrupt the install it was verifying.
    """
    from predictor.annotate import known_operators as KO

    good = {"by_family": {"MerR": ["ACGTACGTACGT"]}, "by_tf": {"X": ["ACGTACGTACGT"]}}
    target = tmp_path / "operators_by_family.json"
    target.write_text(json.dumps(good), encoding="utf-8")

    # both source lists must be blanked: _TF_DIRS supplies the TF->family map, _OP_DIRS the operators
    original_tf, original_op = KO._TF_DIRS, KO._OP_DIRS
    try:
        KO._TF_DIRS = [tmp_path / "definitely-absent"]      # force an empty rebuild
        KO._OP_DIRS = [tmp_path / "definitely-absent"]
        KO.build_index(out=target)
    finally:
        KO._TF_DIRS, KO._OP_DIRS = original_tf, original_op

    assert json.loads(target.read_text(encoding="utf-8")) == good, \
        "an empty rebuild overwrote a populated operator index"


def test_rebuilding_the_curated_reference_cannot_empty_it(tmp_path):
    """Same guarantee for `tf_record.build_reference()`, which regenerates `curated_tfs.faa` from the
    curated YAML tree -- also outside the repo, also absent on a plain install. Measured on a fresh
    `pip install`: 53 sequences -> 0, from running the self-test alone."""
    from predictor.annotate import tf_record as TR

    faa = tmp_path / "curated_tfs.faa"
    jsn = tmp_path / "curated_tfs.json"
    faa.write_text(">TF_One\nMKVKCLLAVFLIVLIAAEHCQ\n", encoding="ascii")
    jsn.write_text(json.dumps({"TF_One": {"family": "MerR"}}), encoding="utf-8")

    original = TR._CURATED_TFS
    try:
        TR._CURATED_TFS = tmp_path / "definitely-absent"
        TR.build_reference(out_faa=faa, out_json=jsn)
    finally:
        TR._CURATED_TFS = original

    assert faa.read_text(encoding="utf-8").count(">") == 1, \
        "an empty rebuild emptied a populated curated reference"


@pytest.mark.released_data
def test_no_fragments(reference):
    """A sequence below MIN_REFERENCE_LEN cannot produce a meaningful blastp alignment; any hit it wins
    is spurious and becomes a confident-looking family call."""
    seqs, _ = reference
    short = {k: len(v) for k, v in seqs.items() if len(v) < FD.MIN_REFERENCE_LEN}
    assert not short, f"fragment records in the reference: {short}"


@pytest.mark.released_data
def test_no_shared_sequences_between_distinct_proteins(reference):
    """Two ids may share a sequence only when they are the same protein under different names -- i.e.
    they carry the same UniProt accession. Anything else means at least one record is wrong."""
    seqs, meta = reference
    by_seq = defaultdict(list)
    for k, v in seqs.items():
        by_seq[v].append(k)
    bad = []
    for ids in by_seq.values():
        if len(ids) > 1:
            accs = {(meta.get(i) or {}).get("uniprot") for i in ids}
            if len(accs) != 1 or None in accs:
                bad.append(sorted(ids))
    assert not bad, f"distinct TFs sharing one sequence: {bad}"


@pytest.mark.released_data
def test_sequences_are_protein_not_dna(reference):
    """A record holding a nucleotide sequence would blast as nonsense against a protein database."""
    dna = set("ACGTN")
    offenders = [k for k, v in reference[0].items()
                 if len(v) > 30 and set(v.upper()) <= dna]
    assert not offenders, f"records that look like DNA, not protein: {offenders}"


@pytest.mark.released_data
def test_curated_records_all_have_a_family(reference):
    """Every CURATED record must carry a family -- that is what curation means.

    The secondary `operators_db` source has a small unlabelled tail (11/307 at the time of writing) whose
    sequences the Pfam HMM also cannot place confidently (18-40 bits). That is a coverage gap, not a
    correctness bug: an unlabelled record simply contributes no family call, which degrades gracefully
    rather than producing a wrong one. The threshold guards against that tail growing.
    """
    seqs, meta = reference
    curated_missing = [k for k in seqs
                       if (meta.get(k) or {}).get("source") == "curated"
                       and not (meta.get(k) or {}).get("family")]
    assert not curated_missing, f"curated records with no family label: {curated_missing}"

    unlabelled = [k for k in seqs if not (meta.get(k) or {}).get("family")]
    assert len(unlabelled) <= 15, (
        f"{len(unlabelled)}/{len(seqs)} reference records carry no family label "
        f"(was 11): {unlabelled[:12]}")


@pytest.mark.released_data
def test_family_labels_agree_with_the_pfam_hmm(reference):
    """Where the HMM is confident, the curated family label must agree with it.

    A mislabelled record poisons every query whose nearest neighbour it is. `SmtB_Msmegmatis` was filed as
    ArsR/SmtB while hitting Pfam FUR at 122 bits (vs 16 for the ArsR-mapped domain) and sharing only 22%
    identity with the genuine mycobacterial ArsR Rv2358; it has been re-filed as Fur.

    Only strong hits are checked (>=50 bits): a weak hit is not evidence the curator was wrong, and some
    families legitimately share a DNA-binding domain (ModE routes to the LysR-type HTH_1, for instance).
    """
    seqs, meta = reference
    disagreements = []
    for tid, seq in sorted(seqs.items()):
        filed = (meta.get(tid) or {}).get("family")
        if not filed:
            continue
        hit = FD.hmm_classify(seq)
        if not hit or hit[4] < 50.0:
            continue
        if hit[2] != filed:
            disagreements.append(f"{tid}: filed {filed!r}, HMM {hit[2]!r} ({hit[0]} {hit[4]:.0f} bits)")
    # Shared-DBD co-classifications are legitimate; a hard failure only for a large, unexplained set.
    assert len(disagreements) <= 12, (
        f"{len(disagreements)} curated records disagree with a strong Pfam HMM hit:\n  "
        + "\n  ".join(disagreements[:20]))


def test_regenerated_references_are_written_with_lf(tmp_path):
    """Vendored references must be written LF-only, on every platform.

    `.gitattributes` checks these files out as LF, but Python text mode translates to CRLF on Windows.
    A rebuild therefore left `curated_tfs.faa`/`.json` permanently "modified" in git with ZERO content
    difference -- and `tfop selftest` triggers a rebuild on any machine that has the curated YAML tree.
    That noise is not cosmetic: these are exactly the files whose accidental truncation the two tests
    above exist to prevent, so a standing spurious diff invites someone to "clean it up" and commit a
    mangled reference.
    """
    from predictor.annotate import tf_record as TR

    faa, jsn = tmp_path / "curated_tfs.faa", tmp_path / "curated_tfs.json"
    src = tmp_path / "curated"
    (src / "merr").mkdir(parents=True)
    (src / "merr" / "TF_One.yaml").write_text(
        "id: TF_One\n"
        "classification:\n  protein_family: MerR\n"
        "sequence:\n  protein: MKVKCLLAVFLIVLIAAEHCQ\n",
        encoding="utf-8")

    original = TR._CURATED_TFS
    try:
        TR._CURATED_TFS = src
        n = TR.build_reference(out_faa=faa, out_json=jsn)
    finally:
        TR._CURATED_TFS = original

    assert n == 1, f"fixture should yield one record, got {n}"
    for p in (faa, jsn):
        raw = p.read_bytes()
        assert b"\r\n" not in raw, f"{p.name} was written with CRLF; pass newline='\n'"


@pytest.mark.released_data
def test_documented_reference_size_matches_the_shipped_file():
    """The reference size is stated in prose in two modules; it must equal what actually ships.

    Both statements were wrong before 2026-09-02 -- `tf_record.classify` said 312 TFs and
    `family_db`'s module docstring said 317, while the file holds 307. Neither number was ever
    checked, so the docstrings drifted independently of the data and of each other. A reader
    calibrating how much to trust a blastp family call was reading a made-up denominator.
    """
    faa = FD._ALL_FAA
    jsn = FD._ALL_JSON
    n_records = sum(1 for ln in faa.read_text(encoding="utf-8").splitlines() if ln.startswith(">"))
    n_meta = len(json.loads(jsn.read_text(encoding="utf-8")))
    assert n_records == n_meta, f"FASTA has {n_records} records, metadata has {n_meta}"

    from predictor.annotate import tf_record
    stated = f"{n_records}-TF broad reference"
    assert stated in (tf_record.classify.__doc__ or ""), (
        f"tf_record.classify docstring must say {stated!r}; the shipped reference holds {n_records}")
    assert f"broad reference of {n_records} TFs" in (FD.__doc__ or ""), (
        f"family_db docstring must say 'broad reference of {n_records} TFs'")
