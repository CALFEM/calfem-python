# -*- coding: utf-8 -*-
"""
This is very simple test of the calfem package.

Make sure all examples run without errors.
"""

import os, sys
import example_output as eo
import difflib
from pprint import pprint

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
        python_executable = sys.executable
        
        return_code = os.system(
            f'"{python_executable}" {example_path} >> test_examples.log 2>&1')
        
        if return_code == 0:
            print(" --- PASSED!")
        else:
            print(" --- FAILED!")

        return_codes += return_code

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
                    example_output[current_example] += line.rstrip() + "\n"

    # Compare output

    if False:
        with open("test_examples_output.log", "w") as f:            
            for example, output in example_output.items():
                if example in eo.examples.keys():
                    print(f"Comparing output of {example}", end="")

                    f.write(f"Comparing output of {example}")
                    d = difflib.Differ()

                    lines = eo.examples[example].splitlines()
                    stripped_lines = [line.rstrip() for line in lines]

                    diff = d.compare(output.splitlines(), stripped_lines)
                    diff_list = list(diff)

                    is_identical = all(line.startswith('  ') for line in diff_list)

                    if not is_identical:

                        f.write(f"\n---------- {example} -----------\n")
                        f.write("\n".join(diff_list))
                        f.write(f"\n---------- {example} -----------\n")
                        return_codes += 1
                        print(" --- FAILED!")
                    else:
                        f.write(f" --- PASSED!\n")
                        print(" --- PASSED!")
    
    assert return_codes == 0

if __name__ == "__main__":
    try: 
        test_examples()
    except AssertionError:
        print("Test failed!")
        sys.exit(1)