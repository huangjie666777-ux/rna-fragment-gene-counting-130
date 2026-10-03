from __future__ import annotations

import io
import os
from pathlib import Path

import pysam
import pytest

from app.bam import coverage_intervals, iter_fragments
from app.counter import count_fragments
from app.gtf import parse_gtf
from app.models import CountError, StrandMode
from app.storage import JobRepository


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session", autouse=True)
def examples():
    out = ROOT / "examples"
    if not (out / "reads.bam").exists():
        from scripts.make_examples import main
        main(str(out))
    return out


@pytest.fixture
def gtf_genes(examples):
    return parse_gtf((examples / "annotation.gtf").read_bytes())


@pytest.fixture
def bam_bytes(examples):
    return (examples / "reads.bam").read_bytes()


def test_gtf_coordinate_conversion_and_merge():
    data = b'1\t.\texon\t1\t10\t.\t+\t.\tgene_id "g";\n1\t.\texon\t5\t20\t.\t+\t.\tgene_id "g";\n'
    genes = parse_gtf(data)
    assert [(e.start, e.end) for e in genes["g"].exons] == [(0, 20)]


def test_gtf_requires_gene_id():
    with pytest.raises(CountError):
        parse_gtf(b'1\t.\texon\t1\t10\t.\t+\t.\tgene_name "x";\n')


def test_gtf_consistent_strand_and_contig():
    with pytest.raises(CountError):
        parse_gtf(
            b'1\t.\texon\t1\t10\t.\t+\t.\tgene_id "g";\n'
            b'1\t.\texon\t21\t30\t.\t-\t.\tgene_id "g";\n'
        )


def test_cigar_coverage():
    seg = pysam.AlignedSegment()
    seg.cigartuples = [(0, 10), (3, 100), (7, 5), (2, 4), (8, 6), (1, 3), (4, 2)]
    seg.reference_start = 0
    assert coverage_intervals(seg) == [(0, 10), (110, 115), (119, 125)]


def test_counts_forward(gtf_genes, bam_bytes):
    result = count_fragments(iter_fragments(io.BytesIO(bam_bytes)), gtf_genes, StrandMode.FORWARD, 10)
    assert result.total == 12
    assert result.filtered == 6
    assert result.ambiguous == 1
    assert result.unassigned == 3
    assert result.assigned == 2
    assert result.counts["GENE_A"] == 2
    assert result.counts["GENE_B"] == 0
    assert result.counts["GENE_C"] == 0
    assert sum(result.counts.values()) == result.assigned


def test_counts_reverse(gtf_genes, bam_bytes):
    result = count_fragments(iter_fragments(io.BytesIO(bam_bytes)), gtf_genes, StrandMode.REVERSE, 10)
    assert result.counts["GENE_A"] == 0
    assert result.counts["GENE_B"] == 1
    assert result.counts["GENE_C"] == 1
    assert result.unassigned == 4


def test_counts_unstranded(gtf_genes, bam_bytes):
    result = count_fragments(iter_fragments(io.BytesIO(bam_bytes)), gtf_genes, StrandMode.UNSTRANDED, 10)
    assert result.counts["GENE_A"] == 2
    assert result.counts["GENE_B"] == 1
    assert result.counts["GENE_C"] == 1
    assert result.ambiguous == 1


def test_minmapq_boundary(gtf_genes, bam_bytes):
    result = count_fragments(iter_fragments(io.BytesIO(bam_bytes)), gtf_genes, StrandMode.FORWARD, 9)
    assert result.filtered == 5
    # frag_lowmapq now passes filters but its reads don't match a same-strand gene here


def _make_bam(records):
    header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as tmp:
        path = tmp.name
    try:
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            for seg in records:
                bam.write(seg)
        return Path(path).read_bytes()
    finally:
        os.unlink(path)


def _seg(name, flag, start=10):
    seg = pysam.AlignedSegment()
    seg.query_name = name
    seg.query_sequence = "A" * 20
    seg.flag = flag
    seg.reference_id = 0
    seg.reference_start = start
    seg.mapping_quality = 60
    seg.cigartuples = [(0, 20)]
    seg.set_tag("NH", 1)
    return seg


def test_reappearing_qname_rejected():
    data = _make_bam([
        _seg("a", 0x1 | 0x40),
        _seg("b", 0x1 | 0x40),
        _seg("a", 0x1 | 0x80),
    ])
    with pytest.raises(CountError):
        list(iter_fragments(io.BytesIO(data)))


def test_truncated_bam_rejected():
    data = _make_bam([_seg("a", 0x1 | 0x40)])
    with pytest.raises(CountError):
        list(iter_fragments(io.BytesIO(data[: len(data) // 2])))


def test_repository_persistence_and_reload(tmp_path, gtf_genes, bam_bytes):
    repo = JobRepository(tmp_path)
    job_id = repo.create(10, StrandMode.FORWARD)
    result = count_fragments(iter_fragments(io.BytesIO(bam_bytes)), gtf_genes, StrandMode.FORWARD, 10)
    repo.save_result(job_id, gtf_genes, result)
    reloaded = JobRepository(tmp_path)
    meta = reloaded.get(job_id)
    assert meta["status"] == "complete"
    assert meta["gene_counts"]["GENE_A"] == 2
    assert reloaded.counts_path(job_id).read_text().splitlines()[1].split("\t") == ["GENE_A", "2"]


def test_http_roundtrip_and_download(tmp_path, examples, monkeypatch):
    monkeypatch.setenv("RNACOUNT_DATA_DIR", str(tmp_path))
    import importlib
    import app.config as config
    importlib.reload(config)
    import app.storage as storage
    importlib.reload(storage)
    import app.main as main
    importlib.reload(main)
    from fastapi.testclient import TestClient

    with TestClient(main.app) as client:
        with (examples / "annotation.gtf").open("rb") as gtf, (examples / "reads.bam").open("rb") as bam:
            response = client.post(
                "/jobs",
                files={"gtf": ("a.gtf", gtf, "text/plain"), "bam": ("a.bam", bam, "application/octet-stream")},
                data={"min_mapq": "10", "strand": "reverse"},
            )
        assert response.status_code == 201, response.text
        job_id = response.json()["job_id"]

        got = client.get(f"/jobs/{job_id}")
        assert got.status_code == 200
        tsv = client.get(f"/jobs/{job_id}/counts.tsv")
        assert tsv.status_code == 200
        body = tsv.text
        assert body.startswith("gene_id\tcount\n")
        assert "GENE_B\t1" in body

        # invalid MAPQ rejected by validation
        with (examples / "annotation.gtf").open("rb") as gtf, (examples / "reads.bam").open("rb") as bam:
            bad = client.post(
                "/jobs",
                files={"gtf": ("a.gtf", gtf, "text/plain"), "bam": ("a.bam", bam, "application/octet-stream")},
                data={"min_mapq": "300", "strand": "forward"},
            )
        assert bad.status_code == 422


def test_concurrent_jobs_do_not_overwrite(tmp_path, examples, monkeypatch):
    import importlib
    import app.config as config
    monkeypatch.setenv("RNACOUNT_DATA_DIR", str(tmp_path))
    importlib.reload(config)
    import app.storage as storage
    importlib.reload(storage)
    import app.main as main
    importlib.reload(main)
    from fastapi.testclient import TestClient

    import concurrent.futures

    def submit(strand):
        with TestClient(main.app) as client:
            with (examples / "annotation.gtf").open("rb") as gtf, (examples / "reads.bam").open("rb") as bam:
                resp = client.post(
                    "/jobs",
                    files={"gtf": ("a.gtf", gtf), "bam": ("a.bam", bam)},
                    data={"min_mapq": "10", "strand": strand},
                )
            return resp

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        responses = list(pool.map(submit, ["forward", "reverse", "unstranded"]))
    payloads = [r.json() for r in responses]
    assert all(r.status_code == 201 for r in responses)
    ids = {p["job_id"] for p in payloads}
    assert len(ids) == 3
    # Each job keeps its own strand-specific result
    by_strand = {p["strand"]: p for p in payloads}
    assert by_strand["forward"]["gene_counts"]["GENE_A"] == 2
    assert by_strand["reverse"]["gene_counts"]["GENE_B"] == 1
    assert by_strand["unstranded"]["gene_counts"]["GENE_B"] == 1
