# -*- coding: utf-8 -*-
"""
Verify that the type hinted element functions accept scalars, lists and
numpy arrays and give the same results.
"""

import sys
from pathlib import Path

import numpy as np

script_dir = Path(__file__).parent.parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

from calfem.core import spring1e, spring1s, bar1e, beam2e


def test_spring1e_input_types():
    k1 = spring1e(100)
    k2 = spring1e([100])
    k3 = spring1e(np.array([100]))

    assert k1.shape == (2, 2)
    np.testing.assert_allclose(k1, [[100, -100], [-100, 100]])
    np.testing.assert_allclose(k2, k1)
    np.testing.assert_allclose(k3, k1)


def test_spring1s_input_types():
    force1 = spring1s(100, [0.1, 0.2])
    force2 = spring1s(100, np.array([0.1, 0.2]))

    np.testing.assert_allclose(force1, 10.0)
    np.testing.assert_allclose(force2, force1)


def test_bar1e_with_and_without_load():
    ke1 = bar1e([0, 2], [210e9, 0.01])
    ke2, fe2 = bar1e([0, 2], [210e9, 0.01], [1000])

    assert ke1.shape == (2, 2)
    assert ke2.shape == (2, 2)
    assert fe2.shape == (2, 1)
    np.testing.assert_allclose(ke2, ke1)


def test_beam2e():
    ke_beam = beam2e([0, 2], [0, 0], [210e9, 0.01, 1e-4])

    assert ke_beam.shape == (6, 6)
    np.testing.assert_allclose(ke_beam, ke_beam.T)
