#!/usr/bin/env python

import setuptools
from sphinx.setup_command import BuildDoc

name = "neighbormodels"
version = "0.0"
release = "0.0.1"

setuptools.setup(
    name=name,
    author="James K. Glasbrenner",
    author_email="jglasbr2@gmu.edu",
    license="MIT",
    version=release,
    description=(
        "A tool for building crystal lattice models with distance-dependent "
        "neighbor interactions."
    ),
    python_requires=">=3.6",
    packages=setuptools.find_packages("src"),
    package_dir={"": "src"},
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "pymatgen",
    ],
    cmdclass={"build_sphinx": BuildDoc},
    command_options={
        "build_sphinx": {
            "project": ("setup.py", name),
            "version": ("setup.py", version),
            "release": ("setup.py", release),
            "source_dir": ("setup.py", "docs"),
            "build_dir": ("setup.py", "docs/_build"),
        }
    },
)
