#!/usr/bin/env python
from setuptools import setup
# from setuptools import find_packages
# find_packages

# under windows, run
# python.exe setup.py bdist --format msi
# to create a windows installer

version_short = '0.1'
version_long = '0.1.0'

if __name__ == '__main__':
    setup(name='cc_fit',
          version=version_long,
          description='Cole-Cole fit routines',
          author='Maximilian Weigand',
          license='GPL-3',
          author_email='mweigand@geo.uni-bonn.de',
          url='http://www.geo.uni-bonn.de/~mweigand',
          # find_packages() somehow does not work under Win7 when creating a
          # msi # installer
          # packages=find_packages(),
          package_dir={'': 'lib'},
          packages=['lib_cc_fit', ],
          scripts=['src/cc_fit.py', ],
          install_requires=['numpy', 'scipy', 'matplotlib'],
          use_2to3=True,
          )
