# -*- coding: utf-8 -*-
"""
Test script for the red function
"""

import sys
import os
from pathlib import Path

# Add src directory to system path to use local calfem package
script_dir = Path(__file__).parent.parent
src_dir = script_dir / "src"
sys.path.insert(0, str(src_dir))

import numpy as np
from calfem.core import red

def test_red_basic():
    """Test basic matrix reduction"""
    print("Testing red() - basic case...")
    
    # Create a simple 3x3 matrix
    K = np.array([[1, 2, 3],
                  [4, 5, 6],
                  [7, 8, 9]])
    
    # Remove second DOF (1-indexed)
    bc = np.array([2])
    K_red = red(K, bc)
    
    # Expected result: 2x2 matrix with rows/cols 0 and 2
    expected = np.array([[1, 3],
                        [7, 9]])
    
    assert K_red.shape == (2, 2), f"Expected shape (2, 2), got {K_red.shape}"
    assert np.allclose(K_red, expected), f"Expected:\n{expected}\nGot:\n{K_red}"
    print("✓ Basic test passed")

def test_red_multiple_dofs():
    """Test reduction with multiple DOFs"""
    print("\nTesting red() - multiple DOFs...")
    
    # Create a 5x5 matrix
    K = np.arange(1, 26).reshape(5, 5)
    
    # Remove DOFs 2 and 4 (1-indexed)
    bc = np.array([2, 4])
    K_red = red(K, bc)
    
    # Should have rows/cols 0, 2, 4 (0-indexed)
    expected = K[np.ix_([0, 2, 4], [0, 2, 4])]
    
    assert K_red.shape == (3, 3), f"Expected shape (3, 3), got {K_red.shape}"
    assert np.allclose(K_red, expected), f"Expected:\n{expected}\nGot:\n{K_red}"
    print("✓ Multiple DOFs test passed")

def test_red_2d_bc_array():
    """Test with 2D boundary condition array (typical CALFEM format)"""
    print("\nTesting red() - 2D boundary condition array...")
    
    # Create a 4x4 stiffness matrix
    K = np.array([[10, -5,  0,  0],
                  [-5, 10, -5,  0],
                  [ 0, -5, 10, -5],
                  [ 0,  0, -5, 10]])
    
    # Boundary conditions in typical CALFEM format: [dof, value]
    bc = np.array([[1, 0],
                   [4, 0]])
    
    K_red = red(K, bc)
    
    # Should remove first and last DOFs, leaving middle 2x2
    expected = np.array([[10, -5],
                        [-5, 10]])
    
    assert K_red.shape == (2, 2), f"Expected shape (2, 2), got {K_red.shape}"
    assert np.allclose(K_red, expected), f"Expected:\n{expected}\nGot:\n{K_red}"
    print("✓ 2D boundary condition array test passed")

def test_red_symmetric_matrix():
    """Test that symmetry is preserved"""
    print("\nTesting red() - symmetric matrix preservation...")
    
    # Create a symmetric stiffness matrix
    K = np.array([[4, -2,  0, -1],
                  [-2, 6, -3,  0],
                  [0, -3,  5, -2],
                  [-1, 0, -2,  3]])
    
    # Verify input is symmetric
    assert np.allclose(K, K.T), "Input matrix should be symmetric"
    
    # Remove DOF 3
    bc = np.array([3])
    K_red = red(K, bc)
    
    # Check that result is still symmetric
    assert np.allclose(K_red, K_red.T), "Reduced matrix should remain symmetric"
    print("✓ Symmetry preservation test passed")

def test_red_realistic_fem():
    """Test with a realistic FEM scenario"""
    print("\nTesting red() - realistic FEM example...")
    
    # Simple 2-element spring system (3 DOFs total)
    # Fixed at DOF 1, free at DOF 2, fixed at DOF 3
    k1 = 1000  # N/m
    k2 = 2000  # N/m
    
    K = np.array([[k1,    -k1,      0],
                  [-k1, k1+k2,    -k2],
                  [0,     -k2,     k2]])
    
    # Fix DOFs 1 and 3
    bc = np.array([1, 3])
    K_red = red(K, bc)
    
    # Should have only the middle DOF
    expected = np.array([[k1 + k2]])
    
    assert K_red.shape == (1, 1), f"Expected shape (1, 1), got {K_red.shape}"
    assert np.allclose(K_red, expected), f"Expected:\n{expected}\nGot:\n{K_red}"
    print("✓ Realistic FEM test passed")

def test_red_input_types():
    """Test that red works with different input types"""
    print("\nTesting red() - different input types...")
    
    # Test with list inputs
    K_list = [[1, 2], [3, 4]]
    bc_list = [1]
    K_red1 = red(K_list, bc_list)
    
    # Test with numpy arrays
    K_array = np.array([[1, 2], [3, 4]])
    bc_array = np.array([1])
    K_red2 = red(K_array, bc_array)
    
    # Both should give the same result
    assert np.allclose(K_red1, K_red2), "Results should be the same regardless of input type"
    
    expected = np.array([[4]])
    assert np.allclose(K_red1, expected), f"Expected {expected}, got {K_red1}"
    print("✓ Input types test passed")

def run_all_tests():
    """Run all red function tests"""
    print("="*60)
    print("Testing red() function")
    print("="*60)
    
    try:
        test_red_basic()
        test_red_multiple_dofs()
        test_red_2d_bc_array()
        test_red_symmetric_matrix()
        test_red_realistic_fem()
        test_red_input_types()
        
        print("\n" + "="*60)
        print("All tests passed! ✓")
        print("="*60)
        return True
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return False
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
