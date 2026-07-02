"""
Regression tests for CGD.py / get_sequence.py.

Run with:  pytest tests/
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from get_sequence import reverse_complement
from CGD import parse_fasta

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_reverse_complement():
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("") == ""


def test_parse_fasta_multi_record_sequences_do_not_leak(tmp_path):
    """
    Regression test for a bug in the original (pre-2026) CGD.py where the
    sequence-accumulator list was created once per file instead of once per
    record, so every record after the first was silently prefixed with all
    prior sequences. The current parser must keep records fully independent.
    """
    fasta = tmp_path / "multi.fa"
    fasta.write_text(
        ">seq1\nACGTACGT\n"
        ">seq2\nTTTTGGGG\n"
        ">seq3\nCCCCAAAA\n"
    )
    records = dict(parse_fasta(str(fasta)))
    assert records["seq1"] == "ACGTACGT"
    assert records["seq2"] == "TTTTGGGG"
    assert records["seq3"] == "CCCCAAAA"


def test_parse_fasta_missing_file():
    import pytest

    with pytest.raises(FileNotFoundError):
        list(parse_fasta("/no/such/file.fa"))


def test_cli_output_flag_writes_to_requested_path(tmp_path):
    fasta = REPO_ROOT + "/input.fa"
    out_file = tmp_path / "custom_output.txt"
    result = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "CGD.py"), "-b", fasta, "-o", str(out_file)],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, result.stderr
    assert out_file.exists()
    content = out_file.read_text()
    assert content.startswith("ID\tStart\tEnd\tStrand\tSequence\tCGDi\n")
    assert len(content.splitlines()) > 1
