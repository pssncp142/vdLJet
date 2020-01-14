from distutils.core import setup
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

setup(name="vdLModel", ext_modules=cythonize(ext_modules))
