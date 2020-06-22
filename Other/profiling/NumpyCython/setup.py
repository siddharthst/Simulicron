from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name='Recombination function in Cython',
  ext_modules=[Extension('_Cyrecombination_', ['Cyrecombination.pyx'],)],
  cmdclass={'build_ext': build_ext},
)