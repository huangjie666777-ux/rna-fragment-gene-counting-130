"""Fragment-level gene counting."""
from __future__ import annotations

from bisect import bisect_right
from typing import Dict, List, Optional, Set, Tuple

from .models import FragmentResult, Gene, ReadView, StrandMode


def _interval_genes(
    index: Dict[str, Dict[str, Tuple[List[Tuple[int, int]], List[List[str]]]]],
    contig: str,
    strand: str,
    intervals: List[Tuple[int, int]],
) -> Set[str]:
    strand_index = index.get(contig, {}).get(strand)
    if strand_index is None:
        return set()
    segments, members = strand_index
    # segments sorted by start; remember end of previous segment for overlap test
    hits: Set[str] = set()
    starts = [seg[0] for seg in segments]
    for qstart, qend in intervals:
        if qend <= qstart:
            continue
        # candidate segments whose start < qend; scan back while they overlap
        j = bisect_right(starts, qend - 1) - 1
        while j >= 0:
            seg_start, seg_end = segments[j]
            if seg_end <= qstart:
                break
            # actual overlap guaranteed by boundary alignment
            hits.update(members[j])
            j -= 1
    return hits


def _fragment_strand(r1: ReadView, mode: StrandMode) -> Optional[str]:
    if mode is StrandMode.UNSTRANDED:
        return None
    if mode is StrandMode.FORWARD:
        return r1.strand
    return "-" if r1.strand == "+" else "+"


def count_fragments(
    fragments,
    genes: Dict[str, Gene],
    mode: StrandMode,
    min_mapq: int,
) -> FragmentResult:
    from .gtf import build_interval_index

    index = build_interval_index(genes)
    counts: Dict[str, int] = {gid: 0 for gid in sorted(genes)}
    result = FragmentResult(counts=counts)

    for group in fragments:
        result.total += 1

        r1 = r2 = None
        for read in group:
            if read.read_number == 1 and r1 is None:
                r1 = read
            elif read.read_number == 2 and r2 is None:
                r2 = read
            else:
                r1 = None  # duplicate/extra end -> malformed fragment -> filter
                break
        # filtering happens before any gene assignment
        if (
            r1 is None
            or r2 is None
            or len(group) != 2
            or not r1.mapped
            or not r2.mapped
            or r1.duplicate
            or r2.duplicate
            or r1.qc_fail
            or r2.qc_fail
            or r1.mapq < min_mapq
            or r2.mapq < min_mapq
            or r1.nh > 1
            or r2.nh > 1
        ):
            result.filtered += 1
            continue

        # Assignment happens only after all filters above have passed.
        genes_hit: Set[str] = set()
        if mode is StrandMode.UNSTRANDED:
            for read in (r1, r2):
                for st in ("+", "-"):
                    genes_hit.update(_interval_genes(index, read.contig, st, read.intervals))
        else:
            strand = _fragment_strand(r1, mode)
            for read in (r1, r2):
                genes_hit.update(_interval_genes(index, read.contig, strand, read.intervals))

        if len(genes_hit) == 1:
            counts[next(iter(genes_hit))] += 1
            result.assigned += 1
        elif len(genes_hit) > 1:
            result.ambiguous += 1
        else:
            result.unassigned += 1

    assert (
        result.assigned + result.ambiguous + result.unassigned + result.filtered
        == result.total
    )
    return result
