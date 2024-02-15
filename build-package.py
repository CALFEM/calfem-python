
import os, sys

def update_setup(package_name, package_version, package_deps):

    with open("setup-template.py", "r") as f:
        setup_template = f.read()

    with open("setup.py", "w") as f:
        f.write(setup_template.format(package_name=package_name, package_version=package_version, package_depends=package_deps)) 

def build_package():
    os.system("python -m build --wheel")   

if __name__ == "__main__":

    package_version = "3.6.6"

    update_setup("calfem-python", package_version, "'numpy', 'visvis', 'pyvtk', 'matplotlib', 'scipy', 'gmsh', 'qtpy', 'vedo', 'tabulate'")

    build_package()

    update_setup("calfem-python-small", package_version, "'numpy', 'visvis', 'matplotlib', 'scipy', 'gmsh', 'tabulate'")

    build_package()

