#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open("README.md") as file:
    readme = file.read()

setup(
    name="molscrub",
    author="Forli Lab",
    version="0.1.1",
    license="GPL-v3",
    author_email="forli@scripps.edu",
    url="https://github.com/forlilab/molscrub",
    description="Enumerate states of small organic molecules: 3D coordinates, tautomers, pH-based adjustments",
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'scrub.py=scrubber.main:main',
        ]
    },
    package_data={
        "scrubber": ["data/*"]
    },
    include_package_data=True,
    zip_safe=True,
    install_requires=[
        "rdkit>=2022.03.1"
    ],
    python_requires=">=3.8",
)