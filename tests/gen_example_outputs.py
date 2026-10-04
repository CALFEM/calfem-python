# -*- coding: utf-8 -*-
"""
Generate the reference values in example_outputs.py used by test_examples.py.

Run this script after an intended change in the output of an example:

    python tests/gen_example_outputs.py
"""

import sys
from pathlib import Path

tests_dir = Path(__file__).parent
sys.path.insert(0, str(tests_dir))

from test_examples import extract_numeric_values, run_example

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
    "exm_stress_2d_materials.py",
    "exm_stress_2d.py",
    "exm_flow_model.py",
]


def gen_output_examples():
    example_dict = {}

    for example in examples:
        print(f"Running: {example}")
        proc = run_example(example)
        assert proc.returncode == 0, f"Example {example} failed with output: {proc.stderr}"
        example_dict[example] = extract_numeric_values(proc.stdout)

    with open(tests_dir / "example_outputs.py", "w") as f:
        f.write("# Example outputs\n")
        f.write("examples = {\n")
        for example, values in example_dict.items():
            f.write(f"    '{example}': {values},\n")
        f.write("}\n")


if __name__ == "__main__":
    gen_output_examples()
