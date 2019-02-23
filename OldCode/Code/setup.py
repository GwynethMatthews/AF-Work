from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("Attempt3_1.pyx",annotate=True),
)