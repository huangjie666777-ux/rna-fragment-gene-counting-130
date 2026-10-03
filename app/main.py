"""HTTP layer for RNA paired-end gene counting."""
from __future__ import annotations

import io
from typing import Any, Dict

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse

from .bam import iter_fragments
from .config import MAX_FILE_BYTES
from .counter import count_fragments
from .gtf import parse_gtf
from .models import CountError, StrandMode
from .storage import repository

app = FastAPI(
    title="RNA paired-end gene counter",
    version="1.0.0",
    description="Submit a GTF annotation and a QNAME-grouped paired-end BAM to count fragments per gene.",
)


async def _read_limited(file: UploadFile, label: str) -> bytes:
    chunks = []
    size = 0
    while True:
        chunk = await file.read(1024 * 1024)
        if not chunk:
            break
        size += len(chunk)
        if size > MAX_FILE_BYTES:
            raise HTTPException(
                status_code=413,
                detail=f"{label} exceeds the {MAX_FILE_BYTES // (1024 * 1024)} MiB upload limit",
            )
        chunks.append(chunk)
    data = b"".join(chunks)
    if not data:
        raise HTTPException(status_code=400, detail=f"{label} must not be empty")
    return data


def _public_meta(meta: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "job_id": meta["job_id"],
        "status": meta["status"],
        "min_mapq": meta["min_mapq"],
        "strand": meta["strand"],
        "created_at": meta.get("created_at"),
        "finished_at": meta.get("finished_at"),
        "num_genes": meta.get("num_genes"),
        "total_fragments": meta.get("total_fragments"),
        "filtered": meta.get("filtered"),
        "ambiguous": meta.get("ambiguous"),
        "unassigned": meta.get("unassigned"),
        "assigned": meta.get("assigned"),
        "gene_counts": meta.get("gene_counts"),
    }


@app.post("/jobs", status_code=201)
async def submit_job(
    gtf: UploadFile = File(..., description="GTF annotation (exon features with gene_id)"),
    bam: UploadFile = File(..., description="QNAME-grouped paired-end BAM"),
    min_mapq: int = Form(..., ge=0, le=255),
    strand: StrandMode = Form(..., description="forward | reverse | unstranded"),
) -> Dict[str, Any]:
    gtf_data = await _read_limited(gtf, "GTF")
    bam_data = await _read_limited(bam, "BAM")

    job_id = repository.create(min_mapq=min_mapq, strand=strand)
    try:
        genes = parse_gtf(gtf_data)
        result = count_fragments(iter_fragments(io.BytesIO(bam_data)), genes, strand, min_mapq)
    except CountError as exc:
        repository.fail(job_id, str(exc))
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except HTTPException:
        repository.fail(job_id, "upload rejected")
        raise
    except Exception as exc:  # never publish partial results
        repository.fail(job_id, repr(exc))
        raise HTTPException(status_code=400, detail=f"input rejected: {exc}") from exc

    repository.save_result(job_id, genes, result)
    meta = repository.get(job_id)
    return _public_meta(meta)


@app.get("/jobs")
def list_jobs() -> Dict[str, Any]:
    jobs = [_public_meta(meta) for meta in repository.list_jobs() if meta.get("status") == "complete"]
    return {"jobs": jobs}


@app.get("/jobs/{job_id}")
def get_job(job_id: str) -> Dict[str, Any]:
    meta = repository.get(job_id)
    if meta is None or meta.get("status") != "complete":
        raise HTTPException(status_code=404, detail="job not found")
    return _public_meta(meta)


@app.get("/jobs/{job_id}/counts.tsv")
def download_counts(job_id: str) -> FileResponse:
    meta = repository.get(job_id)
    if meta is None or meta.get("status") != "complete":
        raise HTTPException(status_code=404, detail="job not found")
    path = repository.counts_path(job_id)
    if path is None:
        raise HTTPException(status_code=404, detail="counts not available")
    return FileResponse(
        path,
        media_type="text/tab-separated-values",
        filename=f"{job_id}_counts.tsv",
    )


@app.get("/health")
def health() -> Dict[str, str]:
    return {"status": "ok"}
