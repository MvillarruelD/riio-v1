"""bitacora.py -- run the survey's own discovery tool, rather than a lookalike of it.

BITACORA (Vizueta, Sanchez-Gracia & Rozas) is what Rondon et al. used to find the 150 candidate
metalloregulators: per genome, per seed profile, a combined BLASTP / TBLASTN / HMMER search with
gene-model curation. Our `annotate.genome_scan` census reaches the same candidate set
(149/150 by accession; the 2026-08-19 discovery decision, recorded in the analysis repo before
commit 5815559), but reaching the same answer is not the same as reproducing the published method,
and reproducing it is the point of this module.

Same bridge pattern as `structure.metalnet`: the engine is an external checkout under `env/tools/`
(gitignored, user-installed per `docs/ENGINES.md`), `env/tools/bitacora_runner.py` is the thin runner we
own, and this adapter is what the pipeline calls. When the engine is absent this module ABSTAINS with
a stated reason -- it never silently substitutes the census, because the two carry different
provenance.

Three things the port has to get right, all of them learned the hard way:

1. **Emit full-length proteins, keyed by accession.** BITACORA's headline output is
   `*_proteins_trimmed.fasta`, cut down to the matched domain -- which is why the survey's manifest
   lists IdeR (Rv2711) at 154 aa against its real 230, with `start`/`end` describing the trimmed span.
   Trimmed sequences cannot be joined to anything, so the trim is carried as an ANNOTATION on a
   full-length record, never as the record itself.
2. **The seed profiles are already on disk.** All 15 are 100 % contained in the vendored SSN member
   DBs (`predictor/data/ssn/database/ssn_members_<FAMILY>.fasta`); measured, 0 sequences of theirs that we lack. So
   nothing is downloaded to run this -- the profiles are built from what we ship.
3. **Record which engine produced a candidate list.** `discovery_engine` rides on every result.

    python -m predictor.discovery.bitacora --self-test          # offline; passes without the engine
    python -m predictor.discovery.bitacora --check              # is the engine installed?
"""
from __future__ import annotations

import os
import shutil
import sys
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]

#: The upstream checkout. Gitignored and user-installed -- see `docs/ENGINES.md` section BITACORA.
ENGINE_DIR = Path(os.environ.get("BITACORA_DIR", _REPO / "env" / "tools" / "bitacora"))
#: The runner we own (committed), which normalises BITACORA's outputs into full-length records.
RUNNER = resources.PACKAGE_DIR / "_engines" / "bitacora_runner.py"
#: Seed profiles, built from the vendored SSN member DBs by `build_profiles()`.
PROFILE_DIR = Path(os.environ.get("BITACORA_PROFILES", resources.cache_path("bitacora_profiles")))
CACHE_DIR = resources.cache_path("bitacora")

#: BITACORA keeps a BLASTP HSP covering `BITACORA_QCOV` of the query OR `BITACORA_SCOV` of the
#: subject. Upstream is 2/3 and 0.8, and that rule -- not the E-value -- is what determines the
#: candidate set: of the 16 candidates the survey lists and our 140-candidate run did not recover,
#: ELEVEN were matched at our own E-value and removed afterwards by this filter.
#:
#: Measured over all four genomes with the E-value held at 1e-5, relaxing only this rule to
#: 0.50/0.50 recovers 17 more of the survey's candidates for 6 additional proteins the survey does
#: not list. Loosening the E-value to 1e-3 instead costs 6.9 unlisted proteins per recovery -- about
#: twenty times worse. See `env/tools/bitacora_patches/README.md`.
COVERAGE_QCOV = "0.50"
COVERAGE_SCOV = "0.50"

#: The engine checkout is gitignored and user-installed, so a patch applied there is lost on
#: reinstall and invisible to anyone reproducing the run. The patched files are tracked here instead
#: and copied over the install before every search.
PATCH_DIR = resources.PACKAGE_DIR / "_engines" / "bitacora_patches"


def ensure_patched() -> list[str]:
    """Copy tracked patches over the user-installed engine. Returns the files it refreshed.

    Silent when the engine is absent -- that is the normal offline state, not an error.
    """
    import filecmp
    import shutil
    dest = ENGINE_DIR / "Scripts"
    if not dest.is_dir() or not PATCH_DIR.is_dir():
        return []
    done = []
    for src in sorted(PATCH_DIR.glob("*.pl")):
        tgt = dest / src.name
        if not tgt.exists() or not filecmp.cmp(src, tgt, shallow=False):
            shutil.copy2(src, tgt)
            done.append(src.name)
    return done

#: E-value for BOTH searches -- BITACORA takes a single `-e` and applies it to BLASTP and HMMER
#: alike, so this is one knob, not two.
#:
#: Promoted from 1e-3 on 2026-08-21 after the seed profiles were fixed (SEED_CLADES + clade-built
#: HMMs). Measured on H37Rv against the survey's published 37, and confirmed on M. avium:
#:
#:     1e-3 -> 67 candidates, all-recall 100 %, site-recall 10/10, precision 55 %
#:     1e-4 -> 45 candidates, all-recall  95 %, site-recall 10/10, precision 78 %
#:     1e-5 -> 36 candidates, all-recall  89 %, site-recall 10/10, precision 92 %
#:
#: All-recall falling is not damage here: EVERY candidate dropped at every threshold has
#: `has_metalnet_site = no` -- three `TetR_SczA` and one `PbrR;ecZntR` at 1e-5, none with a predicted
#: metal site, all of them rows the inducer stage abstains on or calls non-metal anyway. Site-bearing
#: recall is 100 % at all three. Tightening further is NOT safe by the same evidence: it has not been
#: measured, and 1e-5 already sits at the survey's own candidate count.
SEARCH_EVALUE = "1e-5"

#: The survey's 15 seed profiles -> the family whose vendored member DB supplies their sequences.
#: From `analysis/collaborator_mining/README.md`; NikR_NikR returned zero candidates in the survey and
#: is kept so the profile set matches the published one exactly.
SEED_PROFILES = {
    "BsCzrA": "ArsR", "mtNmtR": "ArsR",
    "PbrR": "MerR", "ecZntR": "MerR",
    "ecFur": "Fur", "ecZur": "Fur",
    "CsoR_bsCsoR": "CsoR", "DtxR_bsMntR": "DtxR", "Rrf2_ecIscR": "Rrf2",
    "MarR_AdcR": "MarR", "GntR_LldR": "GntR", "TetR_SczA": "TetR",
    "CopY_saMecI": "CopY", "LysR_ModE": "LysR", "NikR_NikR": "NikR",
}


@dataclass
class BitacoraCandidate:
    """One candidate, always full-length; the trim is an annotation on it."""
    protein_id: str
    sequence: str                      # FULL-LENGTH -- never the trimmed domain
    profile: str = ""
    family: str = ""
    contig: str = ""
    start: int | None = None
    end: int | None = None
    strand: str = ""
    #: what BITACORA matched, as an offset into `sequence`. None when it reported no trim.
    trim_start: int | None = None
    trim_end: int | None = None
    evidence: dict = field(default_factory=dict)

    @property
    def trimmed_sequence(self) -> str:
        """The domain BITACORA matched. Provided for parity with their tables; never an identity."""
        if self.trim_start is None or self.trim_end is None:
            return ""
        return self.sequence[self.trim_start:self.trim_end]


@dataclass
class BitacoraResult:
    status: str                        # 'ok' | 'engine_absent' | 'no_profiles' | 'failed'
    candidates: tuple = ()
    engine: str = "bitacora"
    reason: str = ""
    genome: str = ""
    #: How the searched proteins were obtained: 'proteome' | 'annotated_genome' | 'genome'. The last
    #: means genes were CALLED rather than read, so the identifiers are ours and cannot be joined to
    #: an external table by accession -- the trap that made the survey's manifest unjoinable.
    protein_source: str = ""

    @property
    def available(self) -> bool:
        return self.status == "ok"


#: BITACORA is a Perl/shell pipeline that shells out to BLAST+ and HMMER binaries. HMMER has no
#: Windows build, so on Windows the engine runs through WSL -- the same bridge pattern Folddisco uses.
#: `BITACORA_BIN` is the WSL directory holding hmmsearch/blastp (a conda env's bin/ is the no-sudo way
#: to get them, and is what `docs/ENGINES.md` documents).
_DISTRO = os.environ.get("BITACORA_WSL_DISTRO", "Ubuntu")
_WSL_BIN = os.environ.get("BITACORA_BIN", "$HOME/miniforge3/envs/bitacora/bin")


def _win2wsl(p) -> str:
    """Windows path -> its /mnt/<drive> WSL form (absolute paths only)."""
    p = str(Path(p).resolve())
    if len(p) > 1 and p[1] == ":":
        return "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")
    return p.replace("\\", "/")


def _wsl_exe() -> str | None:
    for cand in (os.environ.get("SystemRoot", r"C:\Windows") + r"\System32\wsl.exe", "wsl.exe", "wsl"):
        if Path(cand).exists() or cand in ("wsl.exe", "wsl"):
            return cand
    return None


def _wsl_has(tool: str) -> bool:
    """Is `tool` runnable inside WSL with BITACORA_BIN on PATH?"""
    exe = _wsl_exe()
    if exe is None:
        return False
    try:
        r = subprocess.run([exe, "-d", _DISTRO, "-e", "bash", "-lc",
                            f'export PATH={_WSL_BIN}:$PATH; command -v {tool}'],
                           capture_output=True, text=True, timeout=60)
        return r.returncode == 0 and bool(r.stdout.strip())
    except Exception:
        return False


def engine_available() -> tuple[bool, str]:
    """(installed?, why not). Checked the same way `tfop setup` checks every other engine.

    The binaries are looked for on the Windows PATH first and then inside WSL, because that is where
    HMMER can actually live on this platform. Reporting "missing" while a working WSL install sits
    behind the bridge would be a false negative, and the whole point of an engine check is that its
    answer can be trusted.
    """
    if not ENGINE_DIR.is_dir():
        return False, (f"BITACORA is not installed at {ENGINE_DIR} "
                       f"(set $BITACORA_DIR, or see docs/ENGINES.md section BITACORA)")
    script = next(ENGINE_DIR.glob("runBITACORA*.sh"), None)
    if script is None:
        return False, f"{ENGINE_DIR} exists but holds no runBITACORA*.sh -- incomplete checkout"
    if not RUNNER.exists():
        return False, f"the runner {RUNNER} is missing from this checkout"
    missing = [t for t in ("blastp", "hmmsearch")
               if not shutil.which(t) and not _wsl_has(t)]
    if missing:
        return False, (f"BITACORA needs {', '.join('`' + m + '`' for m in missing)} (BLAST+ and "
                       f"HMMER). Neither the Windows PATH nor WSL ({_DISTRO}, PATH={_WSL_BIN}) has "
                       f"them -- see docs/ENGINES.md section BITACORA")
    return True, ""


#: Which Pfam accessions define each seed profile's family. BITACORA requires a `<name>_db.hmm`
#: alongside every `<name>_db.fasta`, and this is where that HMM comes from -- the SAME bundled Pfam
#: library `annotate.genome_scan` censuses with (packaged `tf_pfam.hmm`), so the two discovery routes
#: search on the same profiles rather than on two different notions of the family.
#:
#: Override per profile with `BITACORA_HMM_<PROFILE>` pointing at your own .hmm file: the point of
#: routing through BITACORA is that the profile set is yours to change, and hard-coding ours would
#: take that away.
def _pfam_accessions(family_tag: str) -> list[str]:
    """Pfam accessions for a census family TAG (`ArsR`), via the full family name (`ArsR/SmtB`).

    The two vocabularies differ and silently returning nothing for the mismatch is how seven of the
    fifteen profiles ended up with no HMM on the first attempt. `ssn_clusters.family_tag` is the one
    place that knows the mapping, so it is inverted here rather than re-guessed.
    """
    import json

    from predictor.annotate import ssn_clusters as _ssn
    from predictor.annotate.family_db import _PFAM_MAP

    full = next((f for f in _ssn.FAMILIES if _ssn.family_tag(f) == family_tag), family_tag)
    try:
        m = json.loads(Path(_PFAM_MAP).read_text(encoding="utf-8"))
    except Exception:
        return []
    out = []
    for acc, v in m.items():
        fam = v.get("family") if isinstance(v, dict) else v
        if fam in (full, family_tag) and acc.upper().startswith("PF"):
            out.append(acc)
    return sorted(set(out))


#: The ONE Pfam model that defines each seed profile: its DNA-binding domain.
#:
#: The first implementation wrote every Pfam accession mapped to the profile's FAMILY into the
#: profile HMM, and hmmsearch reports a hit when a protein matches ANY model in the file. So
#: `TetR_SczA` carried 32 models -- PF00440 (TetR_N, the DBD of essentially every TetR) plus 31
#: accessory ligand-binding domains that are shared with unrelated regulators -- and returned 83 hits
#: on H37Rv where the survey published 6. `GntR_LldR` carried 9, including an aminotransferase
#: domain (PF00155). This was the HMMER half of the over-call; the BLASTP half was the family-pool
#: FASTA fixed in SEED_CLADES.
#:
#: The DBD is the right choice because it is what makes a protein a member of the family at all,
#: while the accessory domains are what it happens to bind -- and those are precisely the domains
#: shared with proteins that are not regulators of this family.
#:
#: NOTE what this does NOT do: Pfam models are family-level, so paired profiles necessarily share a
#: DBD (`BsCzrA`/`mtNmtR` are both PF01022; `ecFur`/`ecZur` both PF01475). Clade specificity comes
#: from the BLASTP half via SEED_CLADES, never from here. Do not expect this table to separate them.
PROFILE_DBD = {
    "BsCzrA": "PF01022",        # HTH_5, the ArsR/SmtB DBD
    "mtNmtR": "PF01022",
    "PbrR": "PF00376",          # MerR HTH
    "ecZntR": "PF00376",
    "ecFur": "PF01475",         # FUR
    "ecZur": "PF01475",
    "CsoR_bsCsoR": "PF02583",   # CsoR / transcriptional repressor of metal efflux
    "DtxR_bsMntR": "PF01325",   # Fe_dep_repress, the DtxR/MntR DBD
    "Rrf2_ecIscR": "PF02082",   # Rrf2
    "MarR_AdcR": "PF01047",     # MarR
    "GntR_LldR": "PF00392",     # GntR HTH -- not PF00155/PF07702, which are ligand-binding
    "TetR_SczA": "PF00440",     # TetR_N -- not the 31 accessory domains
    "CopY_saMecI": "PF03965",   # Penicillinase_R
    "LysR_ModE": "PF00126",     # HTH_1, the LysR-type DBD
    "NikR_NikR": "PF08753",     # NikR_C
}


def _write_profile_hmm(profile: str, family: str, dest: Path) -> str:
    """Write `<profile>_db.hmm`, extracting this family's models from the bundled Pfam library.

    Returns a short status string. An explicit `BITACORA_HMM_<PROFILE>` env var wins, so a caller can
    substitute a hand-built or differently-thresholded profile without touching this code.
    """
    import os as _os
    import shutil as _sh

    target = dest / f"{profile}_db.hmm"
    override = _os.environ.get(f"BITACORA_HMM_{profile.upper()}")
    if override and Path(override).exists():
        _sh.copyfile(override, target)
        return f"from ${{BITACORA_HMM_{profile.upper()}}}"

    from predictor.annotate.family_db import _PFAM_HMM
    wanted = {PROFILE_DBD[profile]} if profile in PROFILE_DBD else set(_pfam_accessions(family))
    if not wanted or not Path(_PFAM_HMM).exists():
        return "no Pfam model for this family"
    # HMMER3 flatfile: records run from "HMMER3/" to "//". Keep whole records whose ACC matches.
    text = Path(_PFAM_HMM).read_text(encoding="utf-8", errors="replace")
    kept, cur, acc = [], [], None
    for line in text.splitlines(keepends=True):
        if line.startswith("HMMER3"):
            cur, acc = [line], None
            continue
        cur.append(line)
        if line.startswith("ACC "):
            acc = line.split()[1].split(".")[0].strip()
        elif line.startswith("//"):
            if acc in wanted:
                kept.append("".join(cur))
            cur, acc = [], None
    if not kept:
        return "no Pfam model matched"
    target.write_text("".join(kept), encoding="utf-8", newline="\n")

    # A Pfam DBD identifies the FAMILY, which is not what a seed profile is. PF00440 matches all 52
    # TetR-family proteins in H37Rv where the survey published 6 -- the model cannot separate the
    # SczA clade from the rest of the family, because it was never built to. So use the Pfam model
    # only as the alignment frame, and build the profile's real model from the CLADE's own sequences:
    # hmmalign the clade against the DBD, then hmmbuild. Measured on that same TetR case, 52 -> 14.
    # Needs nothing but HMMER, which the engine already requires.
    built = _build_clade_hmm(profile, target, dest)
    return built or f"{len(kept)} Pfam model(s) (family-level; clade model not built)"


def _build_clade_hmm(profile: str, pfam_hmm: Path, dest: Path) -> str:
    """Rebuild `<profile>_db.hmm` from the clade FASTA, in the Pfam DBD's alignment frame.

    Returns a status string, or "" if the model could not be built -- in which case the caller keeps
    the Pfam model, which is over-inclusive but correct as far as it goes.
    """
    clade_fa = dest / f"{profile}_db.fasta"
    if not clade_fa.exists():
        return ""
    on_windows = os.name == "nt"
    sto, tmp_hmm = "/tmp/_tfop_clade.sto", "/tmp/_tfop_clade.hmm"
    q_hmm = _win2wsl(pfam_hmm) if on_windows else str(pfam_hmm)
    q_fa = _win2wsl(clade_fa) if on_windows else str(clade_fa)
    q_out = _win2wsl(dest / f"{profile}_db.hmm") if on_windows else str(dest / f"{profile}_db.hmm")
    inner = (f'export PATH={_WSL_BIN}:$PATH; '
             f'hmmalign --amino --trim -o {sto} "{q_hmm}" "{q_fa}" && '
             f'hmmbuild --amino -n {profile}_clade {tmp_hmm} {sto} >/dev/null && '
             f'cp {tmp_hmm} "{q_out}"')
    cmd = ([os.environ.get("SystemRoot", r"C:\Windows") + r"\System32\wsl.exe",
            "-d", _DISTRO, "-e", "bash", "-lc", inner] if on_windows
           else ["bash", "-lc", inner])
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except Exception:
        return ""
    if proc.returncode != 0:
        return ""
    n = sum(1 for ln in clade_fa.read_text(encoding="utf-8", errors="replace").splitlines()
            if ln.startswith(">"))
    return f"clade model from {n} seq(s)"


#: Each seed profile is the SSN CLADE its name refers to -- not the family pool.
#:
#: The family pool was the first implementation and it is measurably wrong. Measured on H37Rv
#: (`analysis/benchmarking/discovery/`): recall stayed perfect at 37/37, but precision was 37/181 = 20 %,
#: because a 74,698-sequence BLASTP database at e-value 1e-3 matches nearly any HTH protein. Worse,
#: it made the PAIRED profiles identical -- `BsCzrA` and `mtNmtR` were byte-for-byte the same file,
#: as were `ecFur`/`ecZur` and `PbrR`/`ecZntR` -- so 15 profiles collapsed to 9 and the metal
#: specificity that makes them different seeds (Zn vs Ni, Fe vs Zn, Cd/Pb vs Zn) was erased. The
#: survey's whole design is that these are distinct anchors.
#:
#: Every mapping below was read off `ssn_clusters.cluster_table()` by hand rather than matched by
#: name: `CsoR_bsCsoR` happens to equal its cluster id, `MarR_AdcR` does not, and a fuzzy matcher
#: silently resolved four profiles onto an unrelated unnamed clade.
SEED_CLADES = {
    "BsCzrA": "ArsR_c2",        # CzrA/SmtB/ZiaR, Zn2+
    "mtNmtR": "ArsR_c7",        # NmtR, Ni2+ -- distinct from CzrA, which is the point
    "PbrR": "MerR_2",           # CadR/PbrR691, Cd2+/Pb2+
    "ecZntR": "MerR_9",         # ecZntR, Zn2+
    "ecFur": "Fur_ecFur",       # Fe2+
    "ecZur": "Fur_ecZur",       # Zn2+
    "CsoR_bsCsoR": "CsoR_bsCsoR",   # CsoR/RicR, Cu+
    "DtxR_bsMntR": "DtxR_bsMntR",   # MntR, Mn2+
    "Rrf2_ecIscR": "Rrf2_ecIscR",   # IscR, redox[2Fe-2S]
    "MarR_AdcR": "MarR_c1",     # AdcR, Zn2+
    "GntR_LldR": "GntR_c3",     # PdhR/LldR, 2-oxoacid
    "TetR_SczA": "TetR_c1_tetR",    # KstR2/TetR -- the only TetR clade we hold
    "CopY_saMecI": "CopY_c4",   # BlaI/MecI
    "LysR_ModE": "LysR_ModE",   # molybdate
    "NikR_NikR": "NikR_c1",     # Ni2+ -- returned zero candidates in the survey; kept for parity
}


def build_profiles(dest: Path | None = None, *, force: bool = False) -> dict:
    """Write the 15 seed-profile FASTAs, one SSN clade each (see `SEED_CLADES`).

    Nothing is downloaded: the clade FASTAs are the ones `ssn_clusters` already points at. This is
    still not a byte-identical replay of the survey's profiles -- ours are SSN clades, theirs were
    their own curated seed sets -- and the returned dict says so per profile.
    """
    from predictor.annotate import ssn_clusters as _sc

    dest = Path(dest or PROFILE_DIR)
    dest.mkdir(parents=True, exist_ok=True)
    by_id = {c["cluster_id"]: c for c in _sc.cluster_table()}
    out = {}
    for profile, family in sorted(SEED_PROFILES.items()):
        cid = SEED_CLADES.get(profile, "")
        info = by_id.get(cid)
        target = dest / f"{profile}_db.fasta"
        if info is None:
            out[profile] = {"status": "missing_clade", "family": family, "clade": cid}
            continue
        src = Path(info["fasta"])
        if not src.is_absolute():
            src = _REPO / src
        if not src.exists():
            out[profile] = {"status": "missing_clade_fasta", "family": family, "clade": cid,
                            "path": str(src)}
            continue
        rec = {"status": "cached", "family": family, "clade": cid,
               "representative": info.get("representative", ""), "inducer": info.get("inducer"),
               "n_seqs": info.get("n_seqs"), "path": str(target),
               "note": "one SSN clade; NOT a byte-identical copy of the survey's own seed set"}
        if not target.exists() or force:
            shutil.copyfile(src, target)
            rec["status"] = "written"
        rec["hmm"] = ("present" if (dest / f"{profile}_db.hmm").exists() and not force
                      else _write_profile_hmm(profile, family, dest))
        out[profile] = rec
    return out


def _proteome_for(genome_fasta, gff, out_dir: Path, *, verbose: bool = False):
    """Write the proteins BITACORA should search, and say where they came from.

    Reuses `annotate.genome_scan.proteins_from_input`, which prefers a genome's own annotation and
    only calls genes when there is none. Cached next to the run so a repeat costs nothing.
    """
    from predictor.annotate import genome_scan

    cache = Path(out_dir) / "proteins.faa"
    kind_file = Path(out_dir) / "protein_source.txt"
    if cache.exists() and cache.stat().st_size > 0 and kind_file.exists():
        return cache, kind_file.read_text(encoding="utf-8").strip()
    try:
        recs, kind = genome_scan.proteins_from_input(str(genome_fasta), gff=gff, verbose=verbose)
    except TypeError:                       # older signature without `gff`
        recs, kind = genome_scan.proteins_from_input(str(genome_fasta), verbose=verbose)
    except Exception as exc:
        if verbose:
            print(f"  bitacora: could not read proteins: {exc}")
        return None, ""
    if not recs:
        return None, ""
    # Unix newlines explicitly: this FASTA is read by BLAST+ and HMMER inside WSL, and CRLF makes
    # the last residue of every line part of the sequence alphabet as far as they are concerned.
    with cache.open("w", encoding="utf-8", newline="\n") as fh:
        n = 0
        for acc, seq in recs:                # genome_scan yields (protein_id, translation) pairs
            acc, seq = str(acc).split()[0] if acc else "", (seq or "").strip()
            if not acc or not seq:
                continue
            fh.write(f">{acc}\n{seq}\n")
            n += 1
    # Never hand the engine an empty FASTA and let it report the failure: BITACORA's own message for
    # that case says nothing about WHY, and the first time this happened the cause was here.
    if not n:
        if verbose:
            print(f"  bitacora: read {len(recs)} record(s) but none had both an id and a sequence")
        return None, ""
    kind_file.write_text(kind, encoding="utf-8")
    return cache, kind


def find_candidates(genome_fasta, gff=None, *, profiles=None, out_dir=None,
                    timeout: int = 7200, verbose: bool = False) -> BitacoraResult:
    """Run BITACORA over one genome. Abstains -- never raises -- when the engine is absent."""
    ok, why = engine_available()
    if not ok:
        if verbose:
            print(f"  bitacora: {why}")
        return BitacoraResult(status="engine_absent", reason=why, genome=str(genome_fasta))

    prof_dir = Path(profiles or PROFILE_DIR)
    if not prof_dir.is_dir() or not any(prof_dir.glob("*_db.fasta")):
        # Profiles are deterministic derivatives of packaged clade FASTAs. A clean wheel should be
        # usable without knowing about this cache-construction step, so materialize it on first use.
        # The directory remains an explicit override when supplied by a caller.
        built = build_profiles(prof_dir)
        if not any(prof_dir.glob("*_db.fasta")):
            missing = [name for name, row in built.items() if row.get("status", "").startswith("missing")]
            return BitacoraResult(
                status="no_profiles", genome=str(genome_fasta),
                reason=f"could not build seed profiles in {prof_dir}; missing: {', '.join(missing) or 'unknown'}")

    out_dir = Path(out_dir or (CACHE_DIR / Path(genome_fasta).stem))
    out_dir.mkdir(parents=True, exist_ok=True)

    # BITACORA runs in PROTEIN mode (`bitacora_runner.run_bitacora` says why), so what it needs is a
    # proteome, not the genome -- passing `--genome` alone made this function fail 100 % of the time
    # with the engine PRESENT, which the offline self-test could not see. Deriving the proteome here,
    # with the same annotation-preferring reader the census uses, is also what keeps the two discovery
    # routes comparable: both see the same proteins under the same identifiers, so a difference in
    # their output is a difference in METHOD rather than in what was searched.
    prot_fasta, protein_source = _proteome_for(genome_fasta, gff, out_dir, verbose=verbose)
    if prot_fasta is None:
        return BitacoraResult(status="failed", genome=str(genome_fasta),
                              reason=f"no proteins could be read from {genome_fasta}")

    # sys.executable, not "python": this repo is normally driven by an interpreter that is not on
    # PATH under that name, and the failure it produced was a Windows Store shim, not a clean error.
    # Keep the user-installed engine in step with the tracked patches, and put the coverage rule in
    # the environment the runner will carry into WSL.
    ensure_patched()
    os.environ.setdefault("BITACORA_QCOV", COVERAGE_QCOV)
    os.environ.setdefault("BITACORA_SCOV", COVERAGE_SCOV)

    cmd = [os.environ.get("BITACORA_PY", sys.executable), str(RUNNER),
           "--proteins", str(prot_fasta), "--profiles", str(prof_dir),
           "--engine", str(ENGINE_DIR), "--out", str(out_dir),
           "--evalue", os.environ.get("BITACORA_EVALUE", SEARCH_EVALUE),
           "--name", Path(genome_fasta).stem[:24] or "tfop"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return BitacoraResult(status="failed", genome=str(genome_fasta),
                              reason=f"BITACORA exceeded {timeout}s")
    if proc.returncode != 0:
        return BitacoraResult(status="failed", genome=str(genome_fasta),
                              reason=(proc.stderr or proc.stdout or "")[-500:])
    res = _parse_runner_output(out_dir, genome=str(genome_fasta))
    res.protein_source = protein_source
    return res


def _family_of_profiles(profile_field: str) -> str:
    """Family for a `profile` cell, which the runner writes as `A|B` when both profiles hit.

    A bare `SEED_PROFILES.get()` returns "" for those, silently dropping the family of exactly the
    proteins the search was most confident about -- two seed profiles agreeing is the strongest
    evidence BITACORA produces, not a reason to discard it. Where the profiles disagree on family the
    first is taken and both remain visible in `profile`, because the family here is a routing hint
    that `annotate.ssn_clusters` re-decides downstream, not a claim.
    """
    fams = [SEED_PROFILES.get(p.strip(), "") for p in profile_field.split("|") if p.strip()]
    return next((f for f in fams if f), "")


def _parse_runner_output(out_dir: Path, *, genome: str) -> BitacoraResult:
    """Read the runner's normalised TSV. Kept separate so it is testable without the engine."""
    import csv

    table = Path(out_dir) / "candidates.tsv"
    if not table.exists():
        return BitacoraResult(status="failed", genome=genome,
                              reason=f"the runner wrote no candidates.tsv in {out_dir}")
    cands = []
    with table.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r.get("protein_id") or not r.get("sequence"):
                continue
            def _int(k):
                v = (r.get(k) or "").strip()
                return int(v) if v.lstrip("-").isdigit() else None
            cands.append(BitacoraCandidate(
                protein_id=r["protein_id"].strip(), sequence=r["sequence"].strip().upper(),
                profile=(r.get("profile") or "").strip(),
                family=_family_of_profiles(r.get("profile") or ""),
                contig=(r.get("contig") or "").strip(),
                start=_int("start"), end=_int("end"), strand=(r.get("strand") or "").strip(),
                trim_start=_int("trim_start"), trim_end=_int("trim_end"),
                evidence={k: r[k] for k in ("evalue", "bits", "search") if r.get(k)}))
    return BitacoraResult(status="ok", candidates=tuple(cands), genome=genome)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    import tempfile

    ok, why = engine_available()
    print(f"engine installed: {ok}" + (f"  ({why})" if not ok else ""))

    # abstention must be explicit, and must NOT fall back to the census: different provenance
    if not ok:
        r = find_candidates("nonexistent.fna")
        assert r.status == "engine_absent" and not r.available
        assert r.candidates == (), "an absent engine must yield no candidates, not borrowed ones"
        assert r.reason, "abstention must say why"

    assert len(SEED_PROFILES) == 15, "the survey used 15 seed profiles"
    assert SEED_PROFILES["NikR_NikR"] == "NikR", "kept although it found 0 candidates upstream"

    # full-length is the record; the trim is an annotation on it
    c = BitacoraCandidate(protein_id="NP_217227.1", sequence="M" + "A" * 229,
                          trim_start=10, trim_end=164)
    assert len(c.sequence) == 230, "IdeR is 230 aa; the manifest's 154 is BITACORA's trim"
    assert len(c.trimmed_sequence) == 154
    assert c.protein_id == "NP_217227.1", "the accession is the identity, not the sequence"

    # the parser is exercised without the engine
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "candidates.tsv").write_text(
            "protein_id\tsequence\tprofile\tcontig\tstart\tend\tstrand\ttrim_start\ttrim_end\tevalue\n"
            "NP_1.1\tMAKV\tecFur\tNC_1\t100\t400\t+\t1\t3\t1e-40\n"
            "\tMISSINGID\tecFur\tNC_1\t1\t2\t+\t\t\t\n",
            encoding="utf-8")
        res = _parse_runner_output(d, genome="test.fna")
    assert res.available and len(res.candidates) == 1, "rows with no id must be dropped, not guessed"
    got = res.candidates[0]
    assert got.protein_id == "NP_1.1" and got.family == "Fur", "profile -> family mapping"
    assert got.trimmed_sequence == "AK"

    missing = _parse_runner_output(Path(tempfile.gettempdir()) / "definitely-not-here", genome="x")
    assert missing.status == "failed" and missing.reason

    print("OK: abstains without the engine, keeps full-length records keyed by accession, "
          "maps profiles to families, and parses/validates the runner's table.")


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        _demo()
    elif "--check" in sys.argv:
        ok, why = engine_available()
        print(f"BITACORA engine: {'OK' if ok else 'MISSING'}")
        if not ok:
            print(f"  {why}")
        sys.exit(0 if ok else 1)
    else:
        print(__doc__)
        sys.exit(2)
