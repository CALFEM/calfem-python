# -*- coding: utf-8 -*-
"""
This is very simple test of the calfem package.

Make sure all examples run without errors.
"""

import os, sys, re, subprocess
import example_outputs as eo
import numpy as np

from pprint import pprint

def extract_numeric_values(output_text):
    """Extract numeric values from output text using regex."""
    import re
    pattern = r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?'
    return [float(x) for x in re.findall(pattern, output_text)]

def compare_outputs(output1, output2, rtol=1e-5, atol=1e-8):
    """Compare outputs with tolerance."""
    import numpy as np
    values1 = extract_numeric_values(output1)
    values2 = extract_numeric_values(output2)
    
    if len(values1) != len(values2):
        return False
    
    return np.allclose(values1, values2, rtol=rtol, atol=atol)

def test_examples():

    examples_dir = "examples"

    examples = [
        "exs_bar2.py",
        "exs_bar2_la.py",
        "exs_bar2_lb.py",
        "exs_beam1.py",
        "exs_beam2.py",
        "exs_beambar2.py",
        "exs_flw_diff2.py",
        "exs_flw_temp1.py",
        "exs_flw_temp2.py",
        "exs_spring.py",
        "exm_stress_2d.py",
        "exm_stress_2d_materials.py",
        "exm_flow_model.py"
    ]

    # Set environment variable to avoid blocking of plots

    os.environ["CFV_NO_BLOCK"] = "YES"

    # Assume 0 return codes 

    for example in examples:
        print(f"Running: {example} ", end="")

        env = os.environ.copy()
        env["CFV_NO_BLOCK"] = "YES"
        
        proc = subprocess.run(
            [sys.executable, f"examples/{example}"],
            env=env,
            capture_output=True,
            text=True
        )
        
        # Check return code
        assert proc.returncode == 0, f"Example {example} failed with output: {proc.stderr}"
        
        # Compare numeric values within tolerance
        actual_values = extract_numeric_values(proc.stdout)
        expected_values = eo.examples[example]
        
        assert len(actual_values) == len(expected_values), \
            f"Expected {len(expected_values)} values, got {len(actual_values)}"
        
        assert np.allclose(actual_values, expected_values, rtol=1e-5, atol=1e-8), \
            "Numeric values differ significantly"
        
        print(f"PASSED.")

if __name__ == "__main__":
    test_examples()
