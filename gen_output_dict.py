# -*- coding: utf-8 -*-
"""
This script generates a dictionary of example outputs.
"""

import os, sys
from pprint import pprint

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

    # Assume 0 return codes 

    return_codes = 0

    for example in examples:
        print(f"Running: {example}", end="")

        echo_string = f"echo ## EXAMPLE: {example} "
        os.system(echo_string + "-"*(40-len(example)) +
                  " >> run_examples.log 2>&1")

        example_path = os.path.join(examples_dir, example)
        python_executable = sys.executable
        
        return_code = os.system(
            f'"{python_executable}" {example_path} >> run_examples.log 2>&1')
        
        if return_code == 0:
            print(" --- PASSED!")
        else:
            print(" --- FAILED!")

        return_codes += return_code

    assert return_codes == 0

    example_dict = {}

    with open("run_examples.log", "r") as f:
        lines = f.readlines()

        current_example = ""

        for line in lines:
            if "## EXAMPLE:" in line:
                current_example = line.split()[2].strip()
                example_dict[current_example] = ""
            elif current_example!="":
                example_dict[current_example] = example_dict[current_example] + line

    if os.path.exists("example_output.py"):
        os.rename("example_output.py", "example_output.py.bak")

    with open("example_output.py", "w") as f:
        f.write("examples = {\n")
        for example, output in example_dict.items():
            f.write(f"    \"{example}\": '''{output}''',\n")
        f.write("}\n")


if __name__ == "__main__":
    gen_output_examples()