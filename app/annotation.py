"""GTF parsing, gene validation and exon interval indexing.

GTF uses 1-based closed coordinates; everything here is converted to
0-based half-open intervals, and overlapping exons of the same gene are
merged.  Exons of different genes are sliced into an elementary partition
so that any covered base maps to exactly one (possibly multi-gene) set.
"""

from __future__ import annotations

import re
from bisect import bisect_right
from dataclasses import dataclass, field
from pathlib import Path

from .errors import JobRejected
from .limits import MAX_EXONS

_GENE_ID_RE = re.compile(r'gene_id\s+"([^"]+)"')


@dataclass(frozen=True)
class GeneMeta:
    gene_id: str
    contig: str
    strand: str  # '+' or '-'


@dataclass
class ContigIndex:
    starts: list[int] = field(default_factory=list)
    ends: list[int] = field(default_factory=list)
    gene_sets: list[frozenset[str]] = field(default_factory=list)

    def genes_overlapping(self, start: int, end: int) -> frozenset[str]:
        """Return the set of genes whose exons overlap [start, end)."""
        hits: set[str] = set()
        i = bisect_right(self.ends, start)
        while i < len(self.starts) and self.starts[i] < end:
            hits.update(self.gene_sets[i])
            i += 1
        return frozenset(hits)


@dataclass
class AnnotationIndex:
    genes: dict[str, GeneMeta]
    contigs: dict[str, ContigIndex]
    exon_count: int

    def all_gene_ids(self) -> list[str]:
        return sorted(self.genes)


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    intervals.sort()
    merged: list[tuple[int, int]] = []
    for start, end in intervals:
        if merged and start <= merged[-1][1]:  # touching intervals merge too
            if end > merged[-1][1]:
                merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def _build_contig_index(
    exons_by_contig: dict[str, list[tuple[str, int, int]]]
) -> dict[str, ContigIndex]:
    index: dict[str, ContigIndex] = {}
    for contig, raw_exons in exons_by_contig.items():
        per_gene: dict[str, list[tuple[int, int]]] = {}
        for gene_id, start, end in raw_exons:
            per_gene.setdefault(gene_id, []).append((start, end))
        # Sweep merged per-gene exons into elementary partitions.
        events: list[tuple[int, int, str]] = []
        for gene_id, intervals in per_gene.items():
            for start, end in _merge_intervals(intervals):
                events.append((start, 1, gene_id))
                events.append((end, -1, gene_id))
        # Ends (-1) before starts (+1) at the same coordinate, since half-open
        # intervals do not actually touch at that base.
        events.sort(key=lambda ev: (ev[0], ev[1]))
        cidx = ContigIndex()
        active: set[str] = set()
        prev: int | None = None
        for pos, delta, gene_id in events:
            if prev is not None and pos > prev and active:
                cidx.starts.append(prev)
                cidx.ends.append(pos)
                cidx.gene_sets.append(frozenset(active))
            if delta == 1:
                active.add(gene_id)
            else:
                active.discard(gene_id)
            prev = pos
        index[contig] = cidx
    return index


def parse_gtf(path: str | Path) -> AnnotationIndex:
    path = Path(path)
    try:
        raw = path.read_bytes()
    except OSError as exc:  # pragma: no cover - defensive
        raise JobRejected(f"cannot read GTF file: {exc}") from exc
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise JobRejected("GTF must be UTF-8 encoded") from exc

    genes: dict[str, GeneMeta] = {}
    exons_by_contig: dict[str, list[tuple[str, int, int]]] = {}
    exon_count = 0

    for lineno, line in enumerate(text.splitlines(), start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise JobRejected(f"GTF line {lineno}: expected 9 tab-separated fields")
        seqname, _source, feature, start_s, end_s, _score, strand, _frame, attrs = fields
        if feature != "exon":
            continue
        if not seqname:
            raise JobRejected(f"GTF line {lineno}: empty seqname")
        if strand not in ("+", "-"):
            raise JobRejected(
                f"GTF line {lineno}: exon strand must be '+' or '-', got {strand!r}"
            )
        try:
            start_1based = int(start_s)
            end_1based = int(end_s)
        except ValueError as exc:
            raise JobRejected(f"GTF line {lineno}: non-integer coordinates") from exc
        if start_1based < 1 or end_1based < start_1based:
            raise JobRejected(
                f"GTF line {lineno}: invalid 1-based interval {start_s}-{end_s}"
            )
        match = _GENE_ID_RE.search(attrs)
        if not match or not match.group(1):
            raise JobRejected(f"GTF line {lineno}: exon is missing required gene_id")
        gene_id = match.group(1)

        existing = genes.get(gene_id)
        if existing is None:
            genes[gene_id] = GeneMeta(gene_id, seqname, strand)
        elif existing.contig != seqname or existing.strand != strand:
            raise JobRejected(
                f"GTF line {lineno}: gene {gene_id!r} was previously seen on "
                f"{existing.contig}:{existing.strand}, now {seqname}:{strand}"
            )

        # Convert 1-based closed [start, end] to 0-based half-open [s, e).
        exons_by_contig.setdefault(seqname, []).append(
            (gene_id, start_1based - 1, end_1based)
        )
        exon_count += 1
        if exon_count > MAX_EXONS:
            raise JobRejected(f"GTF exceeds the limit of {MAX_EXONS} exons")

    if exon_count == 0:
        raise JobRejected("GTF contains no exon features")

    return AnnotationIndex(
        genes=genes,
        contigs=_build_contig_index(exons_by_contig),
        exon_count=exon_count,
    )
