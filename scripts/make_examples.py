"""Generate a small GTF + QNAME-grouped paired-end BAM for demo and tests.

Run:  .venv/bin/python -m scripts.make_examples examples
"""
from __future__ import annotations

import sys
from pathlib import Path

import pysam

GTF_TEXT = """##description demo annotation
chr1\tdemo\texon\t11\t40\t.\t+\t.\tgene_id "GENE_A";
chr1\tdemo\texon\t31\t60\t.\t+\t.\tgene_id "GENE_A";
chr1\tdemo\texon\t101\t130\t.\t+\t.\tgene_id "GENE_A";
chr1\tdemo\texon\t201\t230\t.\t-\t.\tgene_id "GENE_B";
chr1\tdemo\texon\t301\t330\t.\t+\t.\tgene_id "GENE_C";
"""

# Flags
R1 = 0x1 | 0x40
R2 = 0x1 | 0x80
R1_REV = R1 | 0x10
R2_REV = R2 | 0x10
R2_DUP = R2_REV | 0x400
R1_QC = R1 | 0x200


def _seg(name, flag, contig, start, cigar, mapq=60, nh=1):
    seg = pysam.AlignedSegment()
    seg.query_name = name
    seg.query_sequence = "A" * 30
    seg.flag = flag
    seg.reference_id = 0
    seg.reference_start = start
    seg.mapping_quality = mapq
    seg.cigartuples = cigar
    seg.set_tag("NH", nh)
    return seg


M = (0, 30)
SPLICED = [(0, 20), (3, 60), (0, 10)]  # 20M60N10M


RECORDS = [
    # frag_ok: R1 forward on GENE_A exon1, R2 forward exon2 -> forward+ assigned GENE_A
    _seg("frag_ok", R1, "chr1", 10, [M]),
    _seg("frag_ok", R2_REV, "chr1", 35, [M]),
    # frag_multi: R1 covers exon1/exon2 span incl intron, R2 on GENE_C -> ambiguous (A,C)
    _seg("frag_multi", R1, "chr1", 10, SPLICED),
    _seg("frag_multi", R2_REV, "chr1", 305, [M]),
    # frag_b: R1 reverse-complement flag means read strand '-' ; reverse protocol -> matches '+' GENE_C
    _seg("frag_revc", R1_REV, "chr1", 305, [M]),
    _seg("frag_revc", R2, "chr1", 305, [M]),
    # frag_bminus: R1 '+' on reverse-strand GENE_B; reverse protocol -> assigned GENE_B
    _seg("frag_bminus", R1, "chr1", 205, [M]),
    _seg("frag_bminus", R2_REV, "chr1", 205, [M]),
    # frag_dup: filtered (duplicate)
    _seg("frag_dup", R1, "chr1", 10, [M]),
    _seg("frag_dup", R2_DUP, "chr1", 10, [M]),
    # frag_qc: filtered (QC fail)
    _seg("frag_qc", R1_QC, "chr1", 10, [M]),
    _seg("frag_qc", R2_REV, "chr1", 10, [M]),
    # frag_lowmapq: filtered (MAPQ 9 < default 10)
    _seg("frag_lowmapq", R1, "chr1", 10, [M], mapq=9),
    _seg("frag_lowmapq", R2_REV, "chr1", 10, [M], mapq=60),
    # frag_multi_hit: filtered (NH=2 on R2)
    _seg("frag_multi_hit", R1, "chr1", 10, [M]),
    _seg("frag_multi_hit", R2_REV, "chr1", 10, [M], nh=2),
    # frag_unmapped: filtered (R2 unmapped)
    _seg("frag_unmapped", R1, "chr1", 10, [M]),
    ("frag_unmapped", R2 | 0x4),
    # frag_missing: filtered (R1 only)
    _seg("frag_missing", R1, "chr1", 305, [M]),
    # frag_nowhere: passes filters but hits no exon -> unassigned
    _seg("frag_nowhere", R1, "chr1", 500, [M]),
    _seg("frag_nowhere", R2_REV, "chr1", 500, [M]),
    # frag_overlap: mates overlap each other on GENE_A -> counted once, not twice
    _seg("frag_overlap", R1, "chr1", 10, [M]),
    _seg("frag_overlap", R2_REV, "chr1", 12, [M]),
]


def main(out_dir: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    gtf_path = out / "annotation.gtf"
    gtf_path.write_text(GTF_TEXT.lstrip(), encoding="utf-8")

    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"LN": 1000, "SN": "chr1"}],
    }
    bam_path = out / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for item in RECORDS:
            if isinstance(item, tuple):
                name, flag = item
                seg = pysam.AlignedSegment()
                seg.query_name = name
                seg.query_sequence = "A" * 30
                seg.flag = flag
                seg.reference_id = -1
                seg.reference_start = -1
                seg.mapping_quality = 0
                seg.cigartuples = []
            else:
                seg = item
            bam.write(seg)
    print(f"wrote {gtf_path} and {bam_path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "examples")
