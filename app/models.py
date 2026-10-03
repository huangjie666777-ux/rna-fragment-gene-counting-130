"""Shared value objects."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


class StrandMode(str, Enum):
    FORWARD = "forward"
    REVERSE = "reverse"
    UNSTRANDED = "unstranded"


@dataclass(frozen=True)
class Exon:
    start: int  # 0-based, half-open
    end: int


@dataclass
class Gene:
    gene_id: str
    contig: str
    strand: str  # '+' or '-'
    exons: List[Exon] = field(default_factory=list)


class CountError(Exception):
    """Raised for any malformed/over-limit input; the whole job is rejected."""


@dataclass
class ReadView:
    qname: str
    read_number: int          # 1 or 2
    contig: Optional[str]
    strand: str               # '+' / '-'
    intervals: List[tuple]    # reference coverage [(start, end), ...]
    mapq: int
    nh: int
    duplicate: bool
    qc_fail: bool
    mapped: bool


@dataclass
class FragmentResult:
    total: int = 0
    filtered: int = 0
    ambiguous: int = 0
    unassigned: int = 0
    assigned: int = 0
    counts: Optional[dict] = None
