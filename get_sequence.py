#!/usr/bin/env python3
"""
get_sequence.py — Nucleotide sequence utility functions
========================================================
Author : Vipin Menon, BIG Lab, Hanyang University
Date   : September 2019 (modernized 2026)

Provides helper functions for DNA sequence manipulation used by CGD.
"""


# Precompute the complement translation table once at module level
_COMPLEMENT_TABLE = str.maketrans("ATCG", "TAGC")


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Args:
        sequence: A DNA string containing only A, T, C, G characters (uppercase).

    Returns:
        The reverse complement as an uppercase string.

    Example:
        >>> reverse_complement("ATCG")
        'CGAT'
    """
    return sequence[::-1].translate(_COMPLEMENT_TABLE)


# ---------------------------------------------------------------------------
# Backward-compatible alias (matches the original function name used in CGD.py)
# ---------------------------------------------------------------------------
reverseComp = reverse_complement
