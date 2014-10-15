from distutils.core import setup, Extension
import numpy.distutils.misc_util

pyspa = Extension("_spa", ["_spa.c", "spa.c"])

setup(
    ext_modules=[pyspa],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)
