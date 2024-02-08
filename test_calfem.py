# -*- coding: utf-8 -*-
"""
This is very simple test of the calfem package.

Make sure all examples run without errors.
"""

import os
import example_output as eo


def test_examples():

    examples_dir = "./examples"

    examples = [
        # "exd_beam2_b.py",
        # "exd_beam2_m.py",
        # "exd_beam2_t.py",
        # "exd_beam2_tr.py",
        # "exm_circle_bsplines.py",
        # "exm_flow_model.py",
        # "exm_geometry.py",
        # "exm_stress_2d.py",
        # "exm_stress_2d_export.py",
        # "exm_stress_2d_materials.py",
        # "exm_stress_2d_pyvtk.py",
        # "exm_structured_mesh.py",
        # "exm_temp_2d_markers.py",
        # "exm_temp_2d_splines_arcs.py",
        # "exm_tutorial_1.py",
        # "exm_tutorial_2.py",
        # "exn_bar2g.py",
        # "exn_bar2m.py",
        # "exn_beam2.py",
        # "exn_beam2_b.py",
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

    # Remove old log files

    if os.path.exists("test_examples.log"):
        os.remove("test_examples.log")

    if os.path.exists("test_examples_warnings.log"):
        os.remove("test_examples_warnings.log")

    if os.path.exists("test_examples_output.log"):
        os.remove("test_examples_output.log") 

    # Set environment variable to avoid blocking of plots

    os.environ["CFV_NO_BLOCK"] = "YES"

    # Assume 0 return codes 

    return_codes = 0

    for example in examples:
        print(f"Running: {example}", end="")

        echo_string = f"echo ## EXAMPLE: {example} "
        os.system(echo_string + "-"*(40-len(example)) +
                  " >> test_examples.log 2>&1")

        example_path = os.path.join(examples_dir, example)

        return_code = os.system(
            f"python {example_path} >> test_examples.log 2>&1")

        return_codes += return_code

        print(f" return_code = {return_code}")

    # Parse for warnings and errors

    example_output = {}

    with open("test_examples.log", "r") as f:

        with open("test_examples_warnings.log", "w") as w:

            log = f.readlines()

            current_example = ""
            for line in log:
                if "## EXAMPLE:" in line:
                    current_example = line.split("## EXAMPLE: ")[1].split()[0].strip()
                    example_output[current_example] = ""
                    continue

                if "Warning:" in line:
                    line_items = line.split(":")
                    if len(line_items)>2:
                        w.write(f"{current_example}: {line_items[-2]}: {line_items[-1]} at line: {line_items[-3]}\n")

                if current_example!="":
                    example_output[current_example] += line

    # Compare output

    with open("test_examples_output.log", "w") as f:            
        for example, output in example_output.items():
            if example in eo.examples.keys():
                f.write(f"Comparing output of {example}")
                if output.strip()!=eo.examples[example].strip():
                    f.write(f"Example {example} failed!\n")
                    f.write(f"Example output:\n{output}\n")
                    f.write(f"Expected output:\n{eo.examples[example]}\n")
                    return_codes += 1
                else:
                    f.write(f" --- PASSED!\n")
    
    assert return_codes == 0

if __name__ == "__main__":
    test_examples()