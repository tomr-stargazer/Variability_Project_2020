import pytest
import numpy as np
import astropy.table

from ..apply_new_error_correction import revise_errorbars


def test_revise_errorbars():
    # create a simple table

    expected = astropy.table.Table()
    expected_errs = np.array(
        [1.000170040995895, 1.7149899090577925, 2.9289756462650978]
    )

    uncorrected_errs = 0.025 + np.arange(3)
    test_table = astropy.table.Table()

    bands = ["J", "H", "K"]
    for b in bands:
        test_table[f"{b}APERMAG3ERR"] = uncorrected_errs

        expected[f"{b}APERMAG3ERR"] = expected_errs

    actual = revise_errorbars(test_table, 1, 2)

    assert (actual == expected).all()

    # np.testing.assert_allclose(actual, expected)
