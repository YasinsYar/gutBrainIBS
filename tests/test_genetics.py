from __future__ import annotations

from src.utils.genetics import is_palindromic, strip_ensg_version


def test_palindromic_detection() -> None:
    assert is_palindromic("A", "T")
    assert is_palindromic("C", "G")
    assert not is_palindromic("A", "C")


def test_strip_ensg_version() -> None:
    assert strip_ensg_version("ENSG00000123456.7") == "ENSG00000123456"
