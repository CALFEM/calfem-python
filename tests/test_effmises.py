# -*- coding: utf-8 -*-
"""
Tests for effmises, effective von Mises stress for all analysis types.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

script_dir = Path(__file__).parent.parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

import calfem.core as cfc
import calfem.core_compat as cfc_compat


MODULES = [cfc, cfc_compat]


@pytest.mark.parametrize("module", MODULES)
@pytest.mark.parametrize("es, ptype", [
    ([[100.0, 0.0, 0.0]], 1),                   # plane stress, 3 components
    ([[100.0, 0.0, 0.0, 0.0]], 1),              # plane stress, with sigz
    ([[100.0, 0.0, 0.0, 0.0]], 2),              # plane strain
    ([[100.0, 0.0, 0.0, 0.0]], 3),              # axisymmetry
    ([[100.0, 0.0, 0.0, 0.0, 0.0, 0.0]], 4),    # three dimensional
])
def test_uniaxial_stress(module, es, ptype):
    np.testing.assert_allclose(module.effmises(np.array(es), ptype), [100.0])


@pytest.mark.parametrize("module", MODULES)
def test_pure_shear(module):
    expected = np.sqrt(3.0) * 50.0
    np.testing.assert_allclose(module.effmises(np.array([[0.0, 0.0, 50.0]]), 1), [expected])
    np.testing.assert_allclose(module.effmises(np.array([[0.0, 0.0, 0.0, 50.0]]), 2), [expected])
    for i in range(3, 6):
        es = np.zeros((1, 6))
        es[0, i] = 50.0
        np.testing.assert_allclose(module.effmises(es, 4), [expected])


@pytest.mark.parametrize("module", MODULES)
def test_hydrostatic_stress_is_zero(module):
    np.testing.assert_allclose(module.effmises(np.array([[70.0, 70.0, 70.0, 0.0]]), 2), [0.0], atol=1e-12)
    np.testing.assert_allclose(module.effmises(np.array([[70.0, 70.0, 70.0, 0.0, 0.0, 0.0]]), 4), [0.0], atol=1e-12)


@pytest.mark.parametrize("module", MODULES)
def test_one_value_per_element(module):
    es = np.array([[100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                   [0.0, 200.0, 0.0, 0.0, 0.0, 0.0]])
    np.testing.assert_allclose(module.effmises(es, 4), [100.0, 200.0])


@pytest.mark.parametrize("module", MODULES)
@pytest.mark.parametrize("es, ptype", [
    (np.zeros((1, 5)), 1),
    (np.zeros((1, 3)), 2),
    (np.zeros((1, 3)), 3),
    (np.zeros((1, 4)), 4),
    (np.zeros((1, 3)), 5),
])
def test_invalid_input_raises(module, es, ptype):
    with pytest.raises(ValueError):
        module.effmises(es, ptype)
