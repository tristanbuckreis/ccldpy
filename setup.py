from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize
import numpy as np

with open("C:/Users/trist/My Jupyter Notebooks/CCLD5/README.md","r") as fh:
    long_description = fh.read()

sourcefiles = ['C:/Users/trist/My Jupyter Notebooks/CCLD5/ccldpy.pyx',
               'C:/Users/trist/My Jupyter Notebooks/CCLD5/ccldpy.c']
extensions = [Extension("ccldpy",sourcefiles)]

setup(
    name="ccldpy",
    version="0.0.1",
    author="Tristan E. Buckreis",
    author_email="tristanbuckreis@ucla.edu",
    description="Python package for simulating earthquake rupture surface representation.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=cythonize(extensions),
    include_dirs=[np.get_include()],
)