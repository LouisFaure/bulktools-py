import os
import pathlib
from setuptools import setup, find_packages
import subprocess
from os import path

this_directory = path.abspath(path.dirname(__file__))

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

from distutils.command.build import build


class CustomBuild(build):
    def run(self):
        build.run(self)
        subprocess.call(["make", "-C", "fqcounter"])


setup(
    name="bulktools",
    description="Python package running bulk pipeline",
    version='0.1',
    setup_requires=["setuptools"],
    package_dir={"bulktools": "bulktools"},
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    cmdclass={
        "build": CustomBuild,
    },
    entry_points={
        'console_scripts': ['bt=bulktools.pipeline:main'],
    },
)