# -*- coding: utf-8 -*-
"""
This is very simple test of the calfem package.

Make sure all examples run without errors.
"""

import os


def test_examples():

    examples_dir = "./examples"

    examples = [
        "exd_beam2_b.py",
        "exd_beam2_m.py",
        "exd_beam2_t.py",
        "exd_beam2_tr.py",
        "exm_circle_bsplines.py",
        "exm_flow_model.py",
        "exm_geometry.py",
        "exm_stress_2d.py",
        "exm_stress_2d_export.py",
        "exm_stress_2d_materials.py",
        "exm_stress_2d_pyvtk.py",
        "exm_structured_mesh.py",
        "exm_temp_2d_markers.py",
        "exm_temp_2d_splines_arcs.py",
        "exm_tutorial_1.py",
        "exm_tutorial_2.py",
        "exn_bar2g.py",
        "exn_bar2m.py",
        "exn_beam2.py",
        "exn_beam2_b.py",
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

    if os.path.exists("test_examples.log"):
        os.remove("test_examples.log")

    os.environ["CFV_NO_BLOCK"] = "YES"

    return_codes = 0

    for example in examples:
        print(f"Running: {example}", end="")

        echo_string = f"echo ------ {example} "
        os.system(echo_string + "-"*(40-len(example)) +
                  " >> test_examples.log 2>&1")

        example_path = os.path.join(examples_dir, example)

        return_code = os.system(
            f"python {example_path} >> test_examples.log 2>&1")

        return_codes += return_code

        print(f" return_code = {return_code}")

    assert return_codes == 0
