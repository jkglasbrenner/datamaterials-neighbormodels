#!/usr/bin/env python

import setuptools

setuptools.setup(
    name="datamaterials_neighbors",
    author="James K. Glasbrenner",
    author_email="jglasbr2@gmu.edu",
    license="MIT",
    version="0.0.1",
    description=(
        "Counts the number of pairwise neighbors in a lattice within "
        "and between sub-lattices."
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
)
