# -*- coding: utf-8 -*-
"""
Run examples and compare their printed numeric output with stored reference
values in example_outputs.py.

Regenerate the reference values with tests/gen_example_outputs.py after an
intended change in the output of an example.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

tests_dir = Path(__file__).parent
root_dir = tests_dir.parent
src_dir = root_dir / "src"
sys.path.insert(0, str(src_dir))
sys.path.insert(0, str(tests_dir))

import example_outputs as eo


def extract_numeric_values(output_text):
    """Extract numeric values from output text using regex."""
    pattern = r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?'
    return [float(x) for x in re.findall(pattern, output_text)]


def run_example(example):
    """Run an example with the local calfem package and non-blocking plots."""
    env = os.environ.copy()
    env["CFV_NO_BLOCK"] = "YES"
    env["MPLBACKEND"] = "Agg"
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(src_dir), env.get("PYTHONPATH")]))

    return subprocess.run(
        [sys.executable, str(Path("examples") / example)],
        cwd=root_dir,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )


@pytest.mark.parametrize("example", sorted(eo.examples))
def test_example_output(example):
    proc = run_example(example)

    assert proc.returncode == 0, f"Example {example} failed with output: {proc.stderr}"

    actual_values = extract_numeric_values(proc.stdout)
    expected_values = eo.examples[example]

    assert len(actual_values) == len(expected_values), \
        f"Expected {len(expected_values)} values, got {len(actual_values)}"

    np.testing.assert_allclose(actual_values, expected_values, rtol=1e-5, atol=1e-8)
