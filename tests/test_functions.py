import pytest
from pathlib import Path
import numpy
from offline_bases import Bases


b = Bases()
b_path = Path(b.context["local_bases_fname_pattern"].format("TEST"))
sv_path = Path(b.context["local_sv_fname_pattern"].format("TEST"))


def test_read_svd_01():
    """Test case when no cases are passed."""
    array = b.read_local_svd([], "TEST", 0.0001)
    numpy.testing.assert_allclose(array, numpy.empty([0, 0]))


def test_read_svd_02():
    """Test case when only one case is passed."""

    # Generate test data
    in_array = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    numpy.save(b_path, in_array)
    numpy.savetxt(sv_path, numpy.array([2, 1]))
    out_array = numpy.array([[0.22, 0.12], [0.42, 0.22], [0.62, 0.32], [0.82, 0.42]])

    # Actual test
    array = b.read_local_svd(["."], "TEST", 0.0001)
    numpy.testing.assert_allclose(array, out_array)

    # Remove test data
    b_path.unlink()
    sv_path.unlink()


def test_read_svd_03():
    """Test case when more than one case are passed."""

    # Generate test data
    in_array = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    numpy.save(b_path, in_array)
    numpy.savetxt(sv_path, numpy.array([2, 1]))
    out_array = numpy.array(
        [
            [0.22, 0.12, 0.22, 0.12],
            [0.42, 0.22, 0.42, 0.22],
            [0.62, 0.32, 0.62, 0.32],
            [0.82, 0.42, 0.82, 0.42],
        ]
    )
    # Actual test
    array = b.read_local_svd([".", "."], "TEST", 0.0001)
    numpy.testing.assert_allclose(array, out_array)

    # Remove test data
    b_path.unlink()
    sv_path.unlink()


def test_read_svd_04():
    """Test cutoff svd filter (only one mode passes)."""

    # Generate test data
    in_array = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    numpy.save(b_path, in_array)
    numpy.savetxt(sv_path, numpy.array([2, 1]))
    out_array = numpy.array([[0.22], [0.42], [0.62], [0.82]])

    # Actual test
    array = b.read_local_svd(["."], "TEST", 1.5)
    numpy.testing.assert_allclose(array, out_array)

    # Remove test data
    b_path.unlink()
    sv_path.unlink()


def test_read_svd_05():
    """Test cutoff svd filter (no mode passes)."""

    # Generate test data
    in_array = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    numpy.save(b_path, in_array)
    numpy.savetxt(sv_path, numpy.array([2, 1]))

    # Actual test
    array = b.read_local_svd(["."], "TEST", 2.5)
    numpy.testing.assert_allclose(array, numpy.empty([4, 0]))

    # Remove test data
    b_path.unlink()
    sv_path.unlink()
