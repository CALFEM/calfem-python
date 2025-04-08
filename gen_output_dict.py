# -*- coding: utf-8 -*-
"""
This script generates a dictionary of example outputs.
"""

import os, sys, re, subprocess
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

def gen_output_examples():

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
        "exs_spring.py"
    ]

    # Set environment variable to avoid blocking of plots


    os.environ["CFV_NO_BLOCK"] = "YES"
    env = os.environ.copy()

    # Assume 0 return codes 

    return_codes = 0

    example_dict = {}

    for example in examples:
        print(f"Running: {example}")

        proc = subprocess.run(
            [sys.executable, f"examples/{example}"],
            env=env,
            capture_output=True,
            text=True
        )

        # Check return code
        assert proc.returncode == 0, f"Example {example} failed with output: {proc.stderr}"

        actual_values = extract_numeric_values(proc.stdout)

        example_dict[example] = actual_values
    

    # Save the example dictionary to a file

    with open("example_outputs.py", "w") as f:
        f.write("# Example outputs\n")
        f.write("examples = {\n")
        for example, values in example_dict.items():
            f.write(f"    '{example}': {values},\n")
        f.write("}\n")

if __name__ == "__main__":
    gen_output_examples()