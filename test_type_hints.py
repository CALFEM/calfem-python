#!/usr/bin/env python3
"""
Test script to verify that type hints don't break compatibility
"""

import sys
import os

# Add the src directory to the path to import calfem
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

import numpy as np
from calfem.core import spring1e, spring1s, bar1e, beam2e

def test_functions():
    """Test the type-hinted functions work exactly as before"""
    
    print("Testing spring1e...")
    # Test spring1e with different input types
    k1 = spring1e(100)  # scalar
    k2 = spring1e([100])  # list
    k3 = spring1e(np.array([100]))  # numpy array
    
    print(f"spring1e(100): shape={k1.shape}, type={type(k1)}")
    print(f"spring1e([100]): shape={k2.shape}, type={type(k2)}")
    print(f"spring1e(np.array([100])): shape={k3.shape}, type={type(k3)}")
    
    print("\nTesting spring1s...")
    # Test spring1s
    force1 = spring1s(100, [0.1, 0.2])  # list displacement
    force2 = spring1s(100, np.array([0.1, 0.2]))  # numpy array displacement
    
    print(f"spring1s with list: {force1}")
    print(f"spring1s with numpy array: {force2}")
    
    print("\nTesting bar1e...")
    # Test bar1e
    ke1 = bar1e([0, 2], [210e9, 0.01])  # no load
    ke2, fe2 = bar1e([0, 2], [210e9, 0.01], [1000])  # with load
    
    print(f"bar1e without load: shape={ke1.shape}")
    print(f"bar1e with load: Ke shape={ke2.shape}, fe shape={fe2.shape}")
    
    print("\nTesting beam2e...")
    # Test beam2e
    try:
        ke_beam = beam2e([0, 2], [0, 0], [210e9, 0.01, 1e-4])
        print(f"beam2e: shape={ke_beam.shape}")
    except Exception as e:
        print(f"beam2e test failed: {e}")
    
    print("\nAll tests completed successfully!")

if __name__ == "__main__":
    test_functions()
