"""Concurrent cache readers must never observe incomplete genome publication."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import threading
import zipfile

import pytest

from predictor.annotate import genome_mirror as mirror


ACCESSION = "GCF_TEST_CACHE"
OLD_FASTA = ">record\nACGTACGT\n"
NEW_FASTA = ">record\nACGTACGTACGT\n"
GFF = "##gff-version 3\nrecord\ttest\tCDS\t1\t6\t.\t+\t0\tID=cds\n"


def archive(root, name="archive", fasta=NEW_FASTA, gff=GFF):
    path = root / (name + ".zip")
    with zipfile.ZipFile(path, "w") as stream:
        prefix = f"ncbi_dataset/data/{ACCESSION}/"
        stream.writestr(prefix + "genomic.fna", fasta)
        if gff is not None:
            stream.writestr(prefix + "genomic.gff", gff)
    return path


@pytest.mark.parametrize("existing", [False, True])
def test_reader_sees_old_complete_file_or_cache_miss_during_copy(tmp_path, monkeypatch, existing):
    target = tmp_path / "mirror"
    target.mkdir()
    final = target / (ACCESSION + ".fna")
    if existing:
        final.write_text(OLD_FASTA)
        (target / (ACCESSION + ".gff")).write_text(GFF)
    source = archive(tmp_path)
    entered, proceed = threading.Event(), threading.Event()
    original = mirror.shutil.copyfile

    def paused_copy(source, destination):
        if Path(source).suffix == ".fna":
            Path(destination).write_bytes(b"")
            entered.set()
            assert proceed.wait(5), "test did not release writer"
        return original(source, destination)

    monkeypatch.setattr(mirror.shutil, "copyfile", paused_copy)
    with ThreadPoolExecutor(max_workers=1) as pool:
        writing = pool.submit(mirror.ingest_zip, source, mirror_dir=target)
        try:
            assert entered.wait(5)
            fna, _ = mirror.ensure_genome(ACCESSION, allow_download=False, mirror_dir=target)
            if existing:
                assert fna == final
                assert final.read_text() == OLD_FASTA
            else:
                assert fna is None
                assert not final.exists()
        finally:
            proceed.set()
        assert writing.result(timeout=5) == [ACCESSION]
    assert final.read_text() == NEW_FASTA
    assert not list(target.glob("*.tmp"))


def test_two_writers_publish_only_complete_files(tmp_path, monkeypatch):
    target = tmp_path / "mirror"
    target.mkdir()
    final = target / (ACCESSION + ".fna")
    final.write_text(OLD_FASTA)
    first = archive(tmp_path, "first")
    second_text = ">record\nTTTTCCCCAAAAGGGG\n"
    second = archive(tmp_path, "second", fasta=second_text)
    barrier = threading.Barrier(3)
    proceed = threading.Event()
    original = mirror.shutil.copyfile

    def paused_copy(source, destination):
        if Path(source).suffix == ".fna":
            Path(destination).write_bytes(b"")
            barrier.wait(timeout=5)
            assert proceed.wait(5)
        return original(source, destination)

    monkeypatch.setattr(mirror.shutil, "copyfile", paused_copy)
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(mirror.ingest_zip, source, mirror_dir=target) for source in (first, second)]
        try:
            barrier.wait(timeout=5)
            assert final.read_text() == OLD_FASTA
            assert mirror.ensure_genome(ACCESSION, allow_download=False, mirror_dir=target)[0] == final
        finally:
            proceed.set()
        assert all(future.result(timeout=5) == [ACCESSION] for future in futures)
    assert final.read_text() in (NEW_FASTA, second_text)
    assert not list(target.glob("*.tmp"))


@pytest.mark.parametrize("invalid", ["", ">record\n", "service error", ">\nACGT\n", ">empty\n>full\nACGT\n"])
def test_corrupt_archive_cannot_replace_known_good_files(tmp_path, invalid):
    target = tmp_path / "mirror"
    target.mkdir()
    final = target / (ACCESSION + ".fna")
    final.write_text(OLD_FASTA)
    annotation = target / (ACCESSION + ".gff")
    annotation.write_text(GFF)
    with pytest.raises(ValueError, match="empty or malformed FASTA"):
        mirror.ingest_zip(archive(tmp_path, fasta=invalid), mirror_dir=target)
    assert final.read_text() == OLD_FASTA
    assert annotation.read_text() == GFF
    assert not list(target.glob("*.tmp"))


def test_failed_copy_cleans_temporary_without_truncating_old_file(tmp_path, monkeypatch):
    target = tmp_path / "mirror"
    target.mkdir()
    final = target / (ACCESSION + ".fna")
    final.write_text(OLD_FASTA)

    def fail_copy(source, destination):
        Path(destination).write_bytes(b"partial")
        raise OSError("disk failure")

    monkeypatch.setattr(mirror.shutil, "copyfile", fail_copy)
    with pytest.raises(OSError, match="disk failure"):
        mirror.ingest_zip(archive(tmp_path, gff=None), mirror_dir=target)
    assert final.read_text() == OLD_FASTA
    assert not list(target.glob("*.tmp"))


def test_existing_empty_cache_is_not_reused_and_online_lookup_can_repair(tmp_path, monkeypatch):
    target = tmp_path / "mirror"
    target.mkdir()
    final = target / (ACCESSION + ".fna")
    final.write_text("")
    assert mirror.ensure_genome(ACCESSION, allow_download=False, mirror_dir=target) == (None, None)
    source = archive(tmp_path)
    called = []

    def fetch(accession, *, mirror_dir):
        called.append(accession)
        return {"ingested": mirror.ingest_zip(source, mirror_dir=mirror_dir)}

    monkeypatch.setattr(mirror, "fetch_accession", fetch)
    fna, gff = mirror.ensure_genome(ACCESSION, mirror_dir=target)
    assert called == [ACCESSION]
    assert fna.read_text() == NEW_FASTA
    assert gff.read_text() == GFF


def test_empty_annotation_is_rejected_before_publication(tmp_path):
    target = tmp_path / "mirror"
    with pytest.raises(ValueError, match="empty annotation"):
        mirror.ingest_zip(archive(tmp_path, gff=""), mirror_dir=target)
    assert not list(target.iterdir())
