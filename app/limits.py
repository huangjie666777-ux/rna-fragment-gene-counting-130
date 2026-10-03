"""Hard limits shared by the HTTP layer and the processing pipeline."""

from __future__ import annotations

MAX_FILE_BYTES = 20 * 1024 * 1024  # 20 MiB per uploaded file
MAX_EXONS = 50_000
MAX_FRAGMENTS = 100_000

STRAND_MODES = ("forward", "reverse", "unstranded")
