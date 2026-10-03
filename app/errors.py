"""Errors that cause an entire job to be rejected atomically."""

from __future__ import annotations


class JobRejected(Exception):
    """Raised when input is illegal, truncated or exceeds a hard limit.

    No partial counts are ever persisted for such jobs.
    """

    def __init__(self, message: str, status_code: int = 400) -> None:
        super().__init__(message)
        self.message = message
        self.status_code = status_code
