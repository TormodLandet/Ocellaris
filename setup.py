from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
from codecs import open
import os, sys


here = os.path.abspath(os.path.dirname(__file__))


# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    

# Get the version
for line in open(os.path.join(here, 'ocellaris', '__init__.py'), encoding='utf-8'):
    if line.startswith('__version__'):
        version = line.split('=')[1].strip()[1:-1]


# Which packages we depend on
dependencies = ['PyYAML', 'h5py', 'numpy'] #'fenics-dolfin']


# No need to install dependencies on ReadTheDocs
if os.environ.get('READTHEDOCS') == 'True':
    dependencies = []


# Make the setup.py test command work
class PyTest(TestCommand):
    description = 'Run Ocellaris\' tests with pytest'
    user_options = [('run-unit-tests=', 'u', "Run unit tests"),
                    ('run-regression-tests=', 'r', "Run regression tests")]
    
    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.run_unit_tests = True
        self.run_regression_tests = True
    
    def finalize_options(self):
        TestCommand.finalize_options(self)
        assert self.run_unit_tests or self.run_regression_tests
    
    def run_tests(self):
        import pytest
        args = []
        if self.verbose:
            args.append('-vs')
        if self.run_unit_tests:
            args.append(os.path.join(here, 'tests/'))
        if self.run_regression_tests:
            args.append(os.path.join(here, 'cases/regression_tests.py'))
        
        errno = pytest.main(args)
        sys.exit(errno)


# Give setuptools/pip informattion about the Ocellaris package
setup(
    name='ocellaris',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=version,

    description='A discontinuous Galerkin FEM solver for multiphase free surface flows',
    long_description=long_description,

    # The project's main homepage.
    url='https://bitbucket.org/trlandet/ocellaris',

    # Author details
    author='Tormod Landet',
    author_email='tormod@landet.net',

    # Choose your license
    license='Apache 2.0',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: Apache Software License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: C++'
    ],

    # What does your project relate to?
    keywords='fem fenics cfd dg navier-stokes multi-phase flow',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed.
    install_requires=dependencies,

    # If there are data files included in your packages that need to be
    # installed, specify them here.
    package_data={
        'ocellaris': ['cpp/*.h', 'cpp/*/*.h'],
    },

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'ocellaris=ocellaris.__main__:run_from_console',
            'ocellaris_inspector=ocellaris_post.inspector.__main__:main',
        ],
    },
    
    # Configure the "test" command
    tests_require=['pytest'],
    cmdclass = {'test': PyTest},
)
