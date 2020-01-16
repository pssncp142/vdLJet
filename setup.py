from setuptools import setup, find_namespace_packages
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension("vdLJet",
            sources=["src/vdLJet.pyx"],
            include_dirs = [np.get_include()],
            extra_compile_args=['-fopenmp'],
            extra_link_args=['-fopenmp'],
            )
]

setup(name="vdLModel",  
    package_dir={'': 'src'},
    packages=find_namespace_packages(where='src'),
    ext_modules=cythonize(ext_modules)
    )
