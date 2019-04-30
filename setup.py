#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict

import setuptools
from sphinx.setup_command import BuildDoc

with open('README.rst', 'rt', encoding='utf8') as f:
    readme = f.read()

name = "neighbormodels"
version = "0.2"
release = "0.2.0"

setuptools.setup(
    name=name,
    author="James K. Glasbrenner",
    author_email="jglasbr2@gmu.edu",
    license="MIT",
    version=release,
    url="https://github.com/jkglasbrenner/datamaterials-neighbormodels",
    project_urls=OrderedDict((
        ("Documentation", "https://neighbormodels.readthedocs.io"),
        ("Code", "https://github.com/jkglasbrenner/datamaterials-neighbormodels"),
    )),
    description=(
        "A tool for building crystal lattice models with distance-dependent "
        "neighbor interactions."
    ),
    long_description=readme,
    python_requires=">=3.6",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "pymatgen",
    ],
    extras_require={
        "docs": [
            "sphinx",
            "sphinx_rtd_theme",
        ],
    },
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
