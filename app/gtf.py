"""GTF parsing: exon-only, 1-based closed -> 0-based half-open, merged per gene."""
from __future__ import annotations

from typing import Dict, List, Tuple

from .config import MAX_EXONS
from .models import CountError, Exon, Gene

VALID_STRANDS = set("+-")


def _parse_attributes(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    # GTF attributes: key "value"; key "value";
    i = 0
    n = len(text)
    while i < n:
        while i < n and text[i].isspace():
            i += 1
        key_start = i
        while i < n and (text[i].isalnum() or text[i] in "_-"):
            i += 1
        key = text[key_start:i]
        while i < n and text[i].isspace():
            i += 1
        if not key or i >= n or text[i] != '"':
            # skip malformed remainder token
            while i < n and text[i] != ";":
                i += 1
            if i < n:
                i += 1
            continue
        i += 1  # opening quote
        val_start = i
        buf: List[str] = []
        escaped = False
        while i < n:
            ch = text[i]
            if escaped:
                buf.append(ch)
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == '"':
                break
            else:
                buf.append(ch)
            i += 1
        if i >= n:
            raise CountError("GTF attribute value missing closing quote")
        _ = val_start
        i += 1  # closing quote
        while i < n and text[i] != ";":
            i += 1
        if i < n:
            i += 1
        attrs[key] = "".join(buf)
    return attrs


def parse_gtf(data: bytes) -> Dict[str, Gene]:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CountError(f"GTF is not valid UTF-8: {exc}") from exc

    genes: Dict[str, Gene] = {}
    exon_count = 0

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise CountError(
                f"GTF line {lineno}: expected 9 tab-separated fields, got {len(fields)}"
            )
        seqname, source, feature, start_s, end_s, score, strand, frame, attrs_s = fields
        if feature != "exon":
            continue
        if strand not in VALID_STRANDS:
            raise CountError(f"GTF line {lineno}: exon strand must be '+' or '-', got {strand!r}")
        try:
            start_1 = int(start_s)
            end_1 = int(end_s)
        except ValueError as exc:
            raise CountError(f"GTF line {lineno}: non-integer exon coordinates") from exc
        if start_1 < 1 or end_1 < start_1:
            raise CountError(f"GTF line {lineno}: invalid 1-based exon coordinates {start_1}-{end_1}")

        attrs = _parse_attributes(attrs_s)
        gene_id = attrs.get("gene_id")
        if not gene_id:
            raise CountError(f"GTF line {lineno}: exon is missing required gene_id attribute")

        gene = genes.get(gene_id)
        if gene is None:
            gene = Gene(gene_id=gene_id, contig=seqname, strand=strand)
            genes[gene_id] = gene
        else:
            if gene.contig != seqname:
                raise CountError(
                    f"GTF line {lineno}: gene {gene_id!r} exons appear on multiple sequences "
                    f"({gene.contig!r} vs {seqname!r})"
                )
            if gene.strand != strand:
                raise CountError(
                    f"GTF line {lineno}: gene {gene_id!r} exons appear on both strands "
                    f"({gene.strand!r} vs {strand!r})"
                )

        gene.exons.append(Exon(start=start_1 - 1, end=end_1))  # inclusive 1-based -> half-open
        exon_count += 1
        if exon_count > MAX_EXONS:
            raise CountError(f"GTF exceeds the limit of {MAX_EXONS} exon features")

    if not genes:
        raise CountError("GTF contains no exon features")

    for gene in genes.values():
        gene.exons = _merge_exons(gene.exons)
    return genes


def _merge_exons(exons: List[Exon]) -> List[Exon]:
    ordered = sorted(exons, key=lambda e: (e.start, e.end))
    merged: List[Exon] = []
    for exon in ordered:
        if merged and exon.start < merged[-1].end:
            if exon.end > merged[-1].end:
                merged[-1] = Exon(merged[-1].start, exon.end)
        elif merged and exon.start == merged[-1].end and exon.end == merged[-1].end:
            # zero-length duplicate, nothing to merge
            continue
        else:
            merged.append(exon)
    return merged


def build_interval_index(
    genes: Dict[str, Gene]
) -> Dict[str, Dict[str, Tuple[List[Tuple[int, int]], List[List[str]]]]]:
    """Per-contig non-overlapping exon index for fast gene lookup.

    Returns contig -> (sorted (start, end) segments, gene_id lists per segment).
    Segments are split at every exon boundary so overlapping genes of the
    same strand share a segment and both gene ids are reported.
    """
    # contig -> strand -> list of (start, end, gene_id)
    by_contig: Dict[str, Dict[str, List[Tuple[int, int, str]]]] = {}
    for gene in genes.values():
        by_contig.setdefault(gene.contig, {}).setdefault(gene.strand, [])
        for exon in gene.exons:
            by_contig[gene.contig][gene.strand].append((exon.start, exon.end, gene.gene_id))

    # contig -> strand -> (segments, gene_id lists)
    index: Dict[str, Dict[str, Tuple[List[Tuple[int, int]], List[List[str]]]]] = {}
    for contig, strand_map in by_contig.items():
        index.setdefault(contig, {})
        for strand, exons_on_strand in strand_map.items():
            boundaries = sorted({b for start, end, _ in exons_on_strand for b in (start, end)})
            segments: List[Tuple[int, int]] = []
            members: List[List[str]] = []
            for left, right in zip(boundaries, boundaries[1:]):
                if left >= right:
                    continue
                hits = {
                    gid
                    for start, end, gid in exons_on_strand
                    if start <= left and right <= end
                }
                if hits:
                    segments.append((left, right))
                    members.append(sorted(hits))
            index[contig][strand] = (segments, members)
    return index
