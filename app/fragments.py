"""Paired-end fragment extraction from a QNAME-grouped BAM.

Only primary alignments (no secondary/supplementary flags) are examined.
Each QNAME must form exactly one R1 primary and one R2 primary; every other
condition (missing mate, unmapped, duplicate, QC fail, low MAPQ, NH > 1)
filters the whole fragment.  Coverage only ever comes from CIGAR M/=/X
operations; D and N merely advance the reference cursor.
"""

from __future__ import annotations

from dataclasses import dataclass

import pysam

from .errors import JobRejected
from .limits import MAX_FRAGMENTS

# Pysam flag constants.  AlignedSegment has is_secondary/is_supplementary,
# but the flags are listed explicitly to document the contract.
PAIRED = 0x1
PROPER_PAIR = 0x2
READ_UNMAPPED = 0x4
REVERSE = 0x10
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0x100
QCFAIL = 0x200
DUP = 0x400
SUPPLEMENTARY = 0x800

FILTER_FLAGS = READ_UNMAPPED | QCFAIL | DUP


@dataclass(frozen=True)
class FragmentCoverage:
    """Everything assignment needs for one fragment."""

    contig: str
    r1_is_reverse: bool
    intervals: tuple[tuple[int, int], ...]  # M/=/X reference coverage


def _coverage_intervals(aln: pysam.AlignedSegment) -> list[tuple[int, int]]:
    """Return reference intervals covered by CIGAR M/=/X operations."""
    if not aln.cigartuples:
        raise JobRejected(
            f"mapped alignment for QNAME {aln.query_name!r} has no CIGAR"
        )
    ref_pos = aln.reference_start
    intervals: list[tuple[int, int]] = []
    for operation, length in aln.cigartuples:
        # M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8
        if operation in (0, 7, 8):  # M, =, X produce coverage
            start = ref_pos
            ref_pos += length
            if intervals and intervals[-1][1] == start:
                intervals[-1] = (intervals[-1][0], ref_pos)
            else:
                intervals.append((start, ref_pos))
        elif operation in (2, 3):  # D, N advance reference coordinates only
            ref_pos += length
        elif operation in (1, 4, 5, 6):  # I, S, H, P touch no reference base
            continue
        else:  # pragma: no cover - pysam only defines 0..8
            raise JobRejected(
                f"unsupported CIGAR operation {operation} for QNAME "
                f"{aln.query_name!r}"
            )
    return intervals


def _alignment_is_passing(aln: pysam.AlignedSegment, min_mapq: int) -> bool:
    if aln.flag & FILTER_FLAGS:
        return False
    if aln.mapping_quality < min_mapq:
        return False
    nh = aln.get_tag("NH") if aln.has_tag("NH") else 1
    if nh != 1:
        return False
    return True


def iter_fragments(
    bam_path: str, min_mapq: int
) -> list[FragmentCoverage]:
    """Read every fragment, ordered by first QNAME appearance.

    Returns one entry per QNAME group.  Groups that fail any filter are
    returned as ``None`` entries so callers can account them as filtered.
    """
    try:
        bamfile = pysam.AlignmentFile(bam_path, mode="rb", check_sq=False)
    except (OSError, ValueError, pysam.utils.SamtoolsError) as exc:
        raise JobRejected(f"cannot open BAM: {exc}") from exc

    with bamfile:
        order: list[str] = []
        groups: dict[str, dict[int, pysam.AlignedSegment]] = {}
        try:
            for aln in bamfile.fetch(until_eof=True):
                if aln.flag & (SECONDARY | SUPPLEMENTARY):
                    continue
                if not aln.query_name:
                    raise JobRejected("BAM contains a record without QNAME")
                if aln.flag & READ1:
                    slot = 1
                elif aln.flag & READ2:
                    slot = 2
                else:
                    # Primary alignment that is neither R1 nor R2.
                    slot = 0
                bucket = groups.get(aln.query_name)
                if bucket is None:
                    bucket = {}
                    groups[aln.query_name] = bucket
                    order.append(aln.query_name)
                if slot in bucket:
                    raise JobRejected(
                        f"BAM QNAME {aln.query_name!r} repeats primary "
                        + ("R1" if slot == 1 else "R2" if slot == 2 else "alignment")
                    )
                bucket[slot] = aln
                if len(order) > MAX_FRAGMENTS:
                    raise JobRejected(
                        f"BAM exceeds the limit of {MAX_FRAGMENTS} fragments"
                    )
        except (OSError, ValueError, pysam.utils.SamtoolsError, StopIteration) as exc:
            if isinstance(exc, JobRejected):
                raise
            raise JobRejected(f"truncated or invalid BAM: {exc}") from exc

    fragments: list[FragmentCoverage | None] = []
    for qname in order:
        bucket = groups[qname]
        r1 = bucket.get(1)
        r2 = bucket.get(2)
        # Exactly one primary R1 and one primary R2; no stray primaries,
        # and both mates must be paired reads.
        if (
            r1 is None
            or r2 is None
            or set(bucket) != {1, 2}
            or not (r1.flag & PAIRED and r2.flag & PAIRED)
        ):
            fragments.append(None)
            continue
        if not (_alignment_is_passing(r1, min_mapq) and _alignment_is_passing(r2, min_mapq)):
            fragments.append(None)
            continue
        try:
            intervals = _coverage_intervals(r1)
            intervals.extend(_coverage_intervals(r2))
        except JobRejected:
            raise
        fragments.append(
            FragmentCoverage(
                contig=r1.reference_name,
                r1_is_reverse=bool(r1.flag & REVERSE),
                intervals=tuple(intervals),
            )
        )
    return fragments
