# -*- coding: utf-8 -*-
"""
Tests for the isoparametric solid elements plani4e, plani4s, soli8e and soli8s.

The element routines are checked against hand-computed values, an
independent reference implementation and basic physical requirements
(symmetry, rigid body modes, load resultants and constant strain patch tests).
"""

import sys
import itertools
from pathlib import Path

# Add src directory to system path to use local calfem package
script_dir = Path(__file__).parent.parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

import numpy as np
import pytest
from calfem.core import hooke, plani4e, plani4s, soli8e, soli8s

E = 210e9
NU = 0.3

# Distorted (non-rectangular) quadrilateral
EX4 = np.array([0.0, 2.1, 2.4, -0.2])
EY4 = np.array([0.0, 0.3, 1.8, 1.5])
AREA4 = 0.5*abs((EX4[0] - EX4[2])*(EY4[1] - EY4[3])
                - (EX4[1] - EX4[3])*(EY4[0] - EY4[2]))

# Unit cube in CALFEM node order and a distorted brick
CUBE = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                 [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
BRICK = CUBE @ np.array([[2.0, 0.3, 0.1],
                         [0.2, 1.5, -0.2],
                         [0.1, 0.25, 3.0]]).T + np.array([1.0, -2.0, 0.5])
BRICK += 0.08*np.random.default_rng(0).standard_normal((8, 3))


def gauss(ir):
    if ir == 1:
        return [0.0], [2.0]
    if ir == 2:
        return [-1/np.sqrt(3), 1/np.sqrt(3)], [1.0, 1.0]
    return [-np.sqrt(0.6), 0.0, np.sqrt(0.6)], [5/9, 8/9, 5/9]


def ref_plani4e(ex, ey, t, D, ir, q):
    """Independent implementation of the 4-node plane element."""
    xi_n = np.array([-1, 1, 1, -1])
    eta_n = np.array([-1, -1, 1, 1])
    X = np.column_stack([ex, ey])
    K = np.zeros((8, 8))
    f = np.zeros(8)
    g, w = gauss(ir)
    for (a, wa), (b, wb) in itertools.product(zip(g, w), zip(g, w)):
        N = (1 + xi_n*a)*(1 + eta_n*b)/4
        dN = np.vstack([xi_n*(1 + eta_n*b)/4, eta_n*(1 + xi_n*a)/4])
        J = dN @ X
        dNx = np.linalg.solve(J, dN)
        B = np.zeros((3, 8))
        Nm = np.zeros((2, 8))
        B[0, 0::2] = dNx[0]
        B[1, 1::2] = dNx[1]
        B[2, 0::2] = dNx[1]
        B[2, 1::2] = dNx[0]
        Nm[0, 0::2] = N
        Nm[1, 1::2] = N
        dv = np.linalg.det(J)*wa*wb*t
        K += B.T @ D @ B*dv
        f += Nm.T @ q*dv
    return K, f


def ref_soli8e(X, D, ir, q):
    """Independent implementation of the 8-node brick element."""
    xn = np.array([-1, 1, 1, -1, -1, 1, 1, -1])
    yn = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    zn = np.array([-1, -1, -1, -1, 1, 1, 1, 1])
    K = np.zeros((24, 24))
    f = np.zeros(24)
    g, w = gauss(ir)
    for (a, wa), (b, wb), (c, wc) in itertools.product(*[list(zip(g, w))]*3):
        N = (1 + xn*a)*(1 + yn*b)*(1 + zn*c)/8
        dN = np.vstack([xn*(1 + yn*b)*(1 + zn*c),
                        yn*(1 + xn*a)*(1 + zn*c),
                        zn*(1 + xn*a)*(1 + yn*b)])/8
        J = dN @ X
        dNx = np.linalg.solve(J, dN)
        B = np.zeros((6, 24))
        Nm = np.zeros((3, 24))
        for k in range(3):
            B[k, k::3] = dNx[k]
            Nm[k, k::3] = N
        B[3, 0::3] = dNx[1]
        B[3, 1::3] = dNx[0]
        B[4, 0::3] = dNx[2]
        B[4, 2::3] = dNx[0]
        B[5, 1::3] = dNx[2]
        B[5, 2::3] = dNx[1]
        dv = np.linalg.det(J)*wa*wb*wc
        K += B.T @ D @ B*dv
        f += Nm.T @ q*dv
    return K, f


def n_zero_eigenvalues(K):
    ev = np.linalg.eigvalsh((K + K.T)/2)
    tol = 1e-9*ev.max()
    assert (ev > -tol).all(), "Stiffness matrix has negative eigenvalues"
    return int((abs(ev) < tol).sum())


# ---------------------------------------------------------------- plani4e

def test_plani4e_unit_square_hand_value():
    """Unit square, nu=0, E=t=1: K11 = 1/3 + G/3 = 1/2."""
    D = hooke(1, 1.0, 0.0)
    Ke, fe = plani4e([0, 1, 1, 0], [0, 0, 1, 1], [1, 1.0, 2], D, [2.0, -4.0])
    assert np.isclose(Ke[0, 0], 0.5)
    assert np.allclose(np.asarray(fe).ravel(), np.tile([0.5, -1.0], 4))


@pytest.mark.parametrize("ir", [1, 2, 3])
@pytest.mark.parametrize("ptype, Dsize", [(1, 3), (1, 4), (2, 3), (2, 4)])
def test_plani4e_against_reference(ptype, Dsize, ir):
    t = 0.05
    q = np.array([3.0, -7.0])
    D4 = np.asarray(hooke(2, E, NU))
    if ptype == 1:
        Dref = np.asarray(hooke(1, E, NU))
        Din = Dref if Dsize == 3 else D4
    else:
        Dref = D4[np.ix_([0, 1, 3], [0, 1, 3])]
        Din = Dref if Dsize == 3 else D4

    Ke, fe = plani4e(EX4, EY4, [ptype, t, ir], Din, q)
    Ke = np.asarray(Ke)
    fe = np.asarray(fe).ravel()
    Kr, fr = ref_plani4e(EX4, EY4, t, Dref, ir, q)

    assert Ke.shape == (8, 8)
    assert np.allclose(Ke, Kr, rtol=1e-12, atol=1e-9*abs(Kr).max())
    assert np.allclose(fe, fr)
    assert np.allclose(Ke, Ke.T, atol=1e-9*abs(Ke).max())
    # Resultant of the body force equals q*A*t
    assert np.allclose([fe[0::2].sum(), fe[1::2].sum()], q*AREA4*t)


@pytest.mark.parametrize("ir, n_zero", [(1, 5), (2, 3), (3, 3)])
def test_plani4e_rigid_body_modes(ir, n_zero):
    """3 rigid body modes, plus 2 hourglass modes for 1-point integration."""
    Ke, _ = plani4e(EX4, EY4, [1, 0.1, ir], hooke(1, E, NU))
    Ke = np.asarray(Ke)
    assert n_zero_eigenvalues(Ke) == n_zero
    for u in [np.tile([1.0, 0.0], 4), np.tile([0.0, 1.0], 4),
              np.column_stack([-EY4, EX4]).ravel()]:
        assert np.allclose(Ke @ u, 0, atol=1e-8*abs(Ke).max())


@pytest.mark.parametrize("convert", [list, lambda a: a[None, :], np.matrix])
def test_plani4e_input_types(convert):
    D = hooke(1, E, NU)
    K_ref, _ = plani4e(EX4, EY4, [1, 0.1, 2], D)
    K, _ = plani4e(convert(EX4), convert(EY4), [1, 0.1, 2], D)
    assert np.allclose(K, K_ref)


def test_plani4e_invalid_integration_rule():
    with pytest.raises(ValueError):
        plani4e(EX4, EY4, [1, 0.1, 4], hooke(1, E, NU))


# ---------------------------------------------------------------- plani4s

EPS0_2D = np.array([1e-3, -4e-4, 6e-4])  # epsx epsy gamxy


def linear_displacement_2d(ex, ey, eps):
    """Nodal displacements for a constant in-plane strain field."""
    ux = eps[0]*ex + eps[2]/2*ey
    uy = eps[2]/2*ex + eps[1]*ey
    return np.column_stack([ux, uy]).ravel()


@pytest.mark.parametrize("ir", [1, 2, 3])
@pytest.mark.parametrize("ptype", [1, 2])
def test_plani4s_constant_strain_patch_3x3(ptype, ir):
    D = np.asarray(hooke(1, E, NU))
    ed = linear_displacement_2d(EX4, EY4, EPS0_2D)
    es, et = plani4s(EX4, EY4, [ptype, 0.1, ir], D, ed)
    assert es.shape == (ir*ir, 3)
    assert et.shape == (ir*ir, 3)
    assert np.allclose(et, np.tile(EPS0_2D, (ir*ir, 1)))
    assert np.allclose(es, np.tile(D @ EPS0_2D, (ir*ir, 1)))


@pytest.mark.parametrize("Dsize", [4, 6])
def test_plani4s_plane_stress_full_D(Dsize):
    """With a 4x4/6x6 D, sigz = 0 and epsz = -nu/E*(sigx + sigy)."""
    D = np.asarray(hooke(2, E, NU) if Dsize == 4 else hooke(4, E, NU))
    ed = linear_displacement_2d(EX4, EY4, EPS0_2D)
    es, et = plani4s(EX4, EY4, [1, 0.1, 2], D, ed)
    sig = np.asarray(hooke(1, E, NU)) @ EPS0_2D
    assert es.shape == (4, Dsize)
    assert np.allclose(es[:, [0, 1, 3]], sig)
    assert np.allclose(es[:, 2], 0, atol=1e-12*E)
    assert np.allclose(et[:, [0, 1, 3]], EPS0_2D)
    assert np.allclose(et[:, 2], -NU/E*(sig[0] + sig[1]))
    assert np.allclose(es[:, 4:], 0) and np.allclose(et[:, 4:], 0)


@pytest.mark.parametrize("Dsize", [4, 6])
def test_plani4s_plane_strain_full_D(Dsize):
    """With a 4x4/6x6 D, epsz = 0 and sigz = nu*(sigx + sigy)."""
    D = np.asarray(hooke(2, E, NU) if Dsize == 4 else hooke(4, E, NU))
    ed = linear_displacement_2d(EX4, EY4, EPS0_2D)
    es, et = plani4s(EX4, EY4, [2, 0.1, 2], D, ed)
    e_full = np.zeros(Dsize)
    e_full[[0, 1, 3]] = EPS0_2D
    assert es.shape == (4, Dsize)
    assert np.allclose(et, e_full)
    assert np.allclose(es, D @ e_full)
    assert np.allclose(es[:, 2], NU*(es[:, 0] + es[:, 1]))


def test_plani4s_integration_point_order():
    """u = v = x*y on the unit square gives epsx = y, epsy = x at the
    integration points, ordered as in plani4e."""
    x = np.array([0.0, 1.0, 1.0, 0.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    ed = np.column_stack([x*y, x*y]).ravel()
    _, et = plani4s(x, y, [1, 1.0, 2], hooke(1, E, NU), ed)
    a = 0.5/np.sqrt(3)
    assert np.allclose(et[:, 0], 0.5 + a*np.array([-1, -1, 1, 1]))
    assert np.allclose(et[:, 1], 0.5 + a*np.array([-1, 1, -1, 1]))


@pytest.mark.parametrize("ptype", [1, 2])
def test_plani4s_consistent_with_plani4e(ptype):
    """Strain energy ed^T*Ke*ed equals t*sum(es.et)*detJ*w
    (unit square: detJ = 1/4, ir=2: w = 1)."""
    t = 0.2
    D = np.asarray(hooke(1, E, NU))
    x = np.array([0.0, 1.0, 1.0, 0.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    ed = np.random.default_rng(2).standard_normal(8)*1e-3
    Ke, _ = plani4e(x, y, [ptype, t, 2], D)
    es, et = plani4s(x, y, [ptype, t, 2], D, ed)
    assert np.isclose(ed @ np.asarray(Ke) @ ed, t*np.sum(es*et)/4, rtol=1e-12)


def test_plani4s_rigid_body_gives_zero_stress():
    ed = (np.column_stack([-EY4 + 1.0, EX4 + 2.0])*1e-3).ravel()
    es, et = plani4s(EX4, EY4, [1, 0.1, 3], hooke(1, E, NU), ed)
    assert np.allclose(et, 0, atol=1e-15)
    assert np.allclose(es, 0, atol=1e-15*E)


@pytest.mark.parametrize("convert", [list, lambda a: a[None, :], np.matrix])
def test_plani4s_input_types(convert):
    D = hooke(1, E, NU)
    ed = linear_displacement_2d(EX4, EY4, EPS0_2D)
    es_ref, et_ref = plani4s(EX4, EY4, [1, 0.1, 2], D, ed)
    es, et = plani4s(convert(EX4), convert(EY4), [1, 0.1, 2], D, convert(ed))
    assert np.allclose(es, es_ref)
    assert np.allclose(et, et_ref)


def test_plani4s_invalid_arguments():
    D = hooke(1, E, NU)
    with pytest.raises(ValueError):
        plani4s(EX4, EY4, [1, 0.1, 4], D, np.zeros(8))
    with pytest.raises(ValueError):
        plani4s(EX4, EY4, [3, 0.1, 2], D, np.zeros(8))


# ----------------------------------------------------------------- soli8e

def test_soli8e_unit_cube_hand_value():
    """Unit cube, nu=0, E=1: K11 = 1/9 + 2*G/9 = 2/9."""
    D = hooke(4, 1.0, 0.0)
    Ke, fe = soli8e(*CUBE.T, [2], D, [1.0, 2.0, -8.0])
    assert np.isclose(Ke[0, 0], 2/9)
    assert np.allclose(np.asarray(fe).ravel(), np.tile([1.0, 2.0, -8.0], 8)/8)


@pytest.mark.parametrize("ir", [1, 2, 3])
def test_soli8e_against_reference(ir):
    D = np.asarray(hooke(4, E, NU))
    q = np.array([1.0, -2.0, 4.0])
    Ke, fe = soli8e(*BRICK.T, [ir], D, q)
    Ke = np.asarray(Ke)
    fe = np.asarray(fe).ravel()
    Kr, fr = ref_soli8e(BRICK, D, ir, q)

    assert Ke.shape == (24, 24)
    assert np.allclose(Ke, Kr, rtol=1e-12, atol=1e-9*abs(Kr).max())
    assert np.allclose(fe, fr)
    assert np.allclose(Ke, Ke.T, atol=1e-9*abs(Ke).max())


def test_soli8e_body_force_resultant():
    """Resultant of the body force equals q*V (exact volume with ir=3)."""
    q = np.array([1.0, -2.0, 4.0])
    _, f_vol = ref_soli8e(BRICK, np.eye(6), 3, np.array([1.0, 0.0, 0.0]))
    V = f_vol[0::3].sum()
    for ir in (2, 3):
        _, fe = soli8e(*BRICK.T, [ir], hooke(4, E, NU), q)
        fe = np.asarray(fe).ravel()
        assert np.allclose([fe[k::3].sum() for k in range(3)], q*V,
                           rtol=1e-10)


@pytest.mark.parametrize("ir, n_zero", [(1, 18), (2, 6), (3, 6)])
def test_soli8e_rigid_body_modes(ir, n_zero):
    """6 rigid body modes, plus 12 hourglass modes for 1-point integration."""
    Ke = np.asarray(soli8e(*BRICK.T, [ir], hooke(4, E, NU)))
    assert n_zero_eigenvalues(Ke) == n_zero
    x, y, z = BRICK.T
    zero = 0*x
    modes = [np.tile([1.0, 0.0, 0.0], 8), np.tile([0.0, 1.0, 0.0], 8),
             np.tile([0.0, 0.0, 1.0], 8),
             np.column_stack([-y, x, zero]).ravel(),
             np.column_stack([zero, -z, y]).ravel(),
             np.column_stack([z, zero, -x]).ravel()]
    for u in modes:
        assert np.allclose(Ke @ u, 0, atol=1e-8*abs(Ke).max())


def test_soli8e_return_values():
    D = hooke(4, E, NU)
    assert not isinstance(soli8e(*CUBE.T, [2], D), tuple)
    assert len(soli8e(*CUBE.T, [2], D, [0, 0, 1])) == 2


@pytest.mark.parametrize("convert", [list, lambda a: a[None, :]])
def test_soli8e_input_types(convert):
    D = hooke(4, E, NU)
    K_ref = soli8e(*BRICK.T, [2], D)
    K = soli8e(*[convert(c) for c in BRICK.T], [2], D)
    assert np.allclose(K, K_ref)


def test_soli8e_invalid_integration_rule():
    with pytest.raises(ValueError):
        soli8e(*CUBE.T, [4], hooke(4, E, NU))


# ----------------------------------------------------------------- soli8s

EPS0 = np.array([1e-3, -2e-4, 5e-4, 3e-4, -1e-4, 2e-4])  # xx yy zz xy xz yz


def linear_displacement(X, eps):
    """Nodal displacements u = H x for a constant strain field."""
    H = np.array([[eps[0], eps[3]/2, eps[4]/2],
                  [eps[3]/2, eps[1], eps[5]/2],
                  [eps[4]/2, eps[5]/2, eps[2]]])
    return (X @ H.T).ravel()


@pytest.mark.parametrize("ir", [1, 2, 3])
def test_soli8s_constant_strain_patch(ir):
    D = np.asarray(hooke(4, E, NU))
    ed = linear_displacement(BRICK, EPS0)
    et, es, eci = soli8s(*BRICK.T, [ir], D, ed)
    ngp = ir**3
    assert np.asarray(et).shape == (ngp, 6)
    assert np.asarray(es).shape == (ngp, 6)
    assert np.allclose(et, np.tile(EPS0, (ngp, 1)))
    assert np.allclose(es, np.tile(D @ EPS0, (ngp, 1)))


def test_soli8s_rigid_body_gives_zero_stress():
    x, y, z = BRICK.T
    ed = (np.column_stack([-y + 1.0, x, 0*x + 2.0])*1e-3).ravel()
    et, es, _ = soli8s(*BRICK.T, [2], hooke(4, E, NU), ed)
    assert np.allclose(et, 0, atol=1e-15)
    assert np.allclose(es, 0, atol=1e-15*E)


def test_soli8s_integration_point_coordinates():
    g = 1/np.sqrt(3)
    _, _, eci = soli8s(*CUBE.T, [2], hooke(4, E, NU), np.zeros(24))
    expected = 0.5 + 0.5*g*np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1],
                                     [-1, 1, -1], [-1, -1, 1], [1, -1, 1],
                                     [1, 1, 1], [-1, 1, 1]])
    assert np.allclose(eci, expected)


def test_soli8s_consistent_with_soli8e():
    """Strain energy ed^T*Ke*ed equals the sum of es.et*detJ*w over the
    integration points (unit cube: detJ = 1/8, ir=2: w = 1)."""
    D = np.asarray(hooke(4, E, NU))
    ed = np.random.default_rng(1).standard_normal(24)*1e-3
    Ke = np.asarray(soli8e(*CUBE.T, [2], D))
    et, es, _ = soli8s(*CUBE.T, [2], D, ed)
    energy = np.sum(np.asarray(es)*np.asarray(et))/8
    assert np.isclose(ed @ Ke @ ed, energy, rtol=1e-12)


@pytest.mark.parametrize("convert", [list, lambda a: a[None, :], np.matrix])
def test_soli8s_displacement_input_types(convert):
    D = hooke(4, E, NU)
    ed = linear_displacement(BRICK, EPS0)
    et_ref, es_ref, _ = soli8s(*BRICK.T, [2], D, ed)
    et, es, _ = soli8s(*BRICK.T, [2], D, convert(ed))
    assert np.allclose(et, et_ref)
    assert np.allclose(es, es_ref)


def test_soli8s_invalid_integration_rule():
    with pytest.raises(ValueError):
        soli8s(*CUBE.T, [4], hooke(4, E, NU), np.zeros(24))
