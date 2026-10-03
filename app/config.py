"""Process-wide limits and paths."""
from __future__ import annotations

from pathlib import Path
import os

MAX_FILE_BYTES = 20 * 1024 * 1024  # 20 MiB per upload
MAX_EXONS = 50_000
MAX_FRAGMENTS = 100_000

DATA_DIR = Path(os.environ.get("RNACOUNT_DATA_DIR", "data"))


def ensure_data_dir() -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    return DATA_DIR
