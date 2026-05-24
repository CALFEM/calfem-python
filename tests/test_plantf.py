# -*- coding: utf-8 -*-
"""
Regression tests for plantf stress vector handling.
"""

import sys
from pathlib import Path

import numpy as np

script_dir = Path(__file__).parent.parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

import calfem.core as cfc
import calfem.core_compat as cfc_compat


EX = [0.0, 0.1, 0.0]
EY = [0.0, 0.0, 0.1]
EP = [2, 0.1]
EXPECTED_COLUMN = np.array(
    [[-2.5e9], [-2.5e9], [2.5e9], [0.0], [0.0], [2.5e9]]
)


def plane_strain_stress():
    D = cfc.hooke(2, 200e9, 0.3)
    strain = np.array([[1.0], [1.0], [1.0], [0.0]])
    return D @ strain


def test_plantf_accepts_four_component_plane_strain_stress_shapes():
    stress = plane_strain_stress()

    variants = [
        stress,
        np.asarray(stress).T,
        np.asarray(stress).reshape(-1),
        np.asarray(stress).reshape(-1).tolist(),
    ]

    for variant in variants:
        fe = cfc.plantf(EX, EY, EP, variant)
        np.testing.assert_allclose(fe, EXPECTED_COLUMN)


def test_plantf_accepts_three_component_stress_shapes():
    stress = np.asarray(plane_strain_stress()).reshape(-1)[[0, 1, 3]]

    variants = [
        stress,
        stress.reshape(1, 3),
        stress.reshape(3, 1),
        stress.tolist(),
    ]

    for variant in variants:
        fe = cfc.plantf(EX, EY, EP, variant)
        np.testing.assert_allclose(fe, EXPECTED_COLUMN)


def test_core_compat_plantf_accepts_four_component_plane_strain_stress():
    stress = plane_strain_stress()

    fe = cfc_compat.plantf(EX, EY, EP, stress)

    np.testing.assert_allclose(fe, EXPECTED_COLUMN.reshape(-1))
