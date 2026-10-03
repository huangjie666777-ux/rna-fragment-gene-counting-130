"""On-disk job repository with atomic writes and restart recovery."""
from __future__ import annotations

import json
import os
import tempfile
import threading
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import uuid4

from .config import ensure_data_dir
from .models import FragmentResult, Gene, StrandMode

META_NAME = "meta.json"
COUNTS_NAME = "counts.tsv"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


class JobRepository:
    def __init__(self, root: Optional[Path] = None) -> None:
        self.root = root or ensure_data_dir()
        self.root.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._jobs: Dict[str, Dict[str, Any]] = {}
        self._load_existing()

    def _load_existing(self) -> None:
        for job_dir in sorted(self.root.iterdir()):
            meta_path = job_dir / META_NAME
            if not job_dir.is_dir() or not meta_path.is_file():
                continue
            try:
                with meta_path.open("r", encoding="utf-8") as handle:
                    meta = json.load(handle)
                if isinstance(meta, dict) and "job_id" in meta:
                    self._jobs[meta["job_id"]] = meta
            except (OSError, json.JSONDecodeError):
                continue

    def create(self, min_mapq: int, strand: StrandMode) -> str:
        job_id = uuid4().hex
        meta = {
            "job_id": job_id,
            "status": "processing",
            "min_mapq": min_mapq,
            "strand": strand.value,
            "created_at": _utc_now(),
        }
        with self._lock:
            self._jobs[job_id] = meta
        return job_id

    def save_result(
        self,
        job_id: str,
        genes: Dict[str, Gene],
        result: FragmentResult,
    ) -> None:
        job_dir = self.root / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        ordered_gene_ids = sorted(genes)

        tsv_lines = ["gene_id\tcount"]
        tsv_lines.extend(f"{gid}\t{result.counts[gid]}" for gid in ordered_gene_ids)
        tsv_content = "\n".join(tsv_lines) + "\n"
        self._atomic_write(job_dir / COUNTS_NAME, tsv_content.encode("utf-8"))

        with self._lock:
            meta = dict(self._jobs[job_id])
        meta.update(
            {
                "status": "complete",
                "finished_at": _utc_now(),
                "num_genes": len(genes),
                "total_fragments": result.total,
                "filtered": result.filtered,
                "ambiguous": result.ambiguous,
                "unassigned": result.unassigned,
                "assigned": result.assigned,
                "gene_counts": {gid: result.counts[gid] for gid in ordered_gene_ids},
            }
        )
        self._atomic_write(job_dir / META_NAME, json.dumps(meta, indent=2, sort_keys=True).encode("utf-8"))
        with self._lock:
            self._jobs[job_id] = meta

    def fail(self, job_id: str, message: str) -> None:
        with self._lock:
            self._jobs.pop(job_id, None)
        job_dir = self.root / job_id
        if job_dir.exists():
            for child in job_dir.iterdir():
                child.unlink()
            job_dir.rmdir()

    @staticmethod
    def _atomic_write(path: Path, content: bytes) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(prefix=".tmp-", dir=path.parent)
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(tmp_name, path)
        except BaseException:
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
            raise

    def get(self, job_id: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            meta = self._jobs.get(job_id)
            return dict(meta) if meta else None

    def list_jobs(self) -> List[Dict[str, Any]]:
        with self._lock:
            return [dict(meta) for meta in self._jobs.values()]

    def counts_path(self, job_id: str) -> Optional[Path]:
        path = self.root / job_id / COUNTS_NAME
        return path if path.is_file() else None


repository = JobRepository()
