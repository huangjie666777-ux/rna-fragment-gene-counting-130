"""Paired-end fragment grouping and read-level filtering."""
from __future__ import annotations

import os
import tempfile
from typing import BinaryIO, Dict, Iterator, List, Optional, Union

import pysam

from .config import MAX_FRAGMENTS
from .models import CountError, ReadView

# pysam flags
F_PAIRED = 0x1
F_PROPER_PAIR = 0x2
F_UNMAP = 0x4
F_MUNMAP = 0x8
F_REVERSE = 0x10
F_MREVERSE = 0x20
F_READ1 = 0x40
F_READ2 = 0x80
F_SECONDARY = 0x100
F_QCFAIL = 0x200
F_DUP = 0x400
F_SUPPLEMENTARY = 0x800

COVERS_REF = frozenset((0, 7, 8))  # M, =, X
ADVANCES_REF = frozenset((0, 2, 3, 7, 8))  # M, D, N, =, X
_CODES = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"}


def coverage_intervals(alignment: pysam.AlignedSegment) -> List[tuple]:
    """0-based half-open reference intervals covered by M/=/X only.

    D and N advance the reference position without producing coverage.
    """
    ref_pos = alignment.reference_start
    intervals: List[tuple] = []
    for op, length in alignment.cigartuples or []:
        if op in COVERS_REF:
            intervals.append((ref_pos, ref_pos + length))
        if op in ADVANCES_REF:
            ref_pos += length
        elif op not in (1, 4, 5, 6):  # I/S/H/P don't advance reference
            raise CountError(f"unsupported CIGAR op {_CODES.get(op, op)} for {alignment.query_name}")
    return intervals


def _nh(alignment: pysam.AlignedSegment) -> int:
    try:
        value = alignment.get_tag("NH")
    except KeyError:
        return 1
    if not isinstance(value, int) or value < 1:
        raise CountError(f"invalid NH tag for {alignment.query_name}: {value!r}")
    return value


def _as_read_view(alignment: pysam.AlignedSegment) -> ReadView:
    mapped = alignment.reference_id >= 0 and not (alignment.flag & F_UNMAP)
    return ReadView(
        qname=alignment.query_name,
        read_number=1 if alignment.is_read1 else 2,
        contig=alignment.reference_name if mapped else None,
        strand="-" if (mapped and alignment.is_reverse) else "+",
        intervals=coverage_intervals(alignment) if mapped else [],
        mapq=alignment.mapping_quality,
        nh=_nh(alignment) if mapped else 1,
        duplicate=bool(alignment.flag & F_DUP),
        qc_fail=bool(alignment.flag & F_QCFAIL),
        mapped=mapped,
    )


def iter_fragments(bam_data: Union[bytes, BinaryIO]) -> Iterator[List[ReadView]]:
    """Yield primary-alignment groups keyed by QNAME (contiguous).

    Secondary and supplementary alignments are ignored. Non-contiguous
    QNAME groups (a group reappearing after a different group was seen)
    reject the whole file. Non-consecutive duplicates of identical qname
    (re-groups) also reject.
    """
    raw = bam_data if isinstance(bam_data, (bytes, bytearray)) else bam_data.read()
    tmp = tempfile.NamedTemporaryFile(prefix="upload-", suffix=".bam", delete=False)
    tmp_path = tmp.name
    try:
        with tmp:
            tmp.write(bytes(raw))
        with pysam.AlignmentFile(tmp_path, "rb") as bam:
            seen_groups: Dict[str, object] = {}
            current_name: Optional[str] = None
            current: List[ReadView] = []
            group_count = 0

            for alignment in bam.fetch(until_eof=True):
                flag = alignment.flag
                if flag & (F_SECONDARY | F_SUPPLEMENTARY):
                    continue
                name = alignment.query_name
                if not name:
                    raise CountError("alignment record is missing QNAME")
                if name != current_name:
                    if current_name is not None:
                        yield current
                        current = []
                    if name in seen_groups:
                        raise CountError(
                            f"QNAME group {name!r} reappears non-contiguously; whole file rejected"
                        )
                    seen_groups[name] = None
                    current_name = name
                    group_count += 1
                    if group_count > MAX_FRAGMENTS:
                        raise CountError(f"BAM exceeds the limit of {MAX_FRAGMENTS} fragments")
                current.append(_as_read_view(alignment))
            if current_name is not None:
                yield current
    except CountError:
        raise
    except (OSError, ValueError, StopIteration) as exc:
        # pysam raises ValueError / OSError for truncated or invalid BAM
        raise CountError(f"invalid or truncated BAM: {exc}") from exc
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
