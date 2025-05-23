from setuptools import setup, find_packages
from mylotools.version import __version__

setup(
    name="mylotools",
    version=__version__,
    packages=find_packages(),
    install_requires=[
        "plotly",
        "matplotlib",
        "biopython",
    ],
    scripts=['bin/mylotools'],  # This tells setup.py to install your bin script
)

