"""Importable logic for the BacDive-AI robustness benchmark.

Thin Snakemake adapters live in ``benchmark/scripts/``; everything testable
lives here. Random draws use :func:`derive_seed`, a stable SHA-256 hash of the
inputs, never Python's ``hash()``.
"""

import hashlib

__all__ = ["derive_seed"]


def derive_seed(*parts):
    """Derive a stable 32-bit seed from config seed and draw coordinates."""
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode()).hexdigest()
    return int(digest[:8], 16)
