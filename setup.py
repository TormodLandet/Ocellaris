from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
from codecs import open
import os
import sys


here = os.path.abspath(os.path.dirname(__file__))


# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


# Get the version
for line in open(os.path.join(here, 'ocellaris', '__init__.py'), encoding='utf-8'):
    if line.startswith('__version__'):
        version = line.split('=')[1].strip()[1:-1]


# List packages we depend on
FENICS_VERSION = ">=2018.1.0.dev0,<2018.2"
dependencies = ['PyYAML',
                'h5py',
                'numpy',
                'meshio>=2.0.0',
                'raschii>=1.0.2',
                'yschema>=1.0.2']
dependencies.append('fenics-dolfin%s' % FENICS_VERSION)


# No need to install dependencies on ReadTheDocs
if os.environ.get('READTHEDOCS') == 'True':
    dependencies = []


# Make the setup.py test command work
class PyTest(TestCommand):
    description = 'Run Ocellaris\' tests with pytest'
    user_options = [('skip-unit-tests=', 'u', "Skip unit tests"),
                    ('skip-regression-tests=', 'r', "Skip regression tests"),
                    ('skip-demo-tests=', 'd', "Skip demo tests")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.skip_unit_tests = False
        self.skip_regression_tests = False
        self.skip_demo_tests = False

    def finalize_options(self):
        TestCommand.finalize_options(self)

    def run_tests(self):
        import pytest
        args = ['-v', '--durations=10']
        if self.verbose:
            args.append('-s')

        if not self.skip_unit_tests:
            args.append(os.path.join(here, 'tests/'))
        if not self.skip_regression_tests:
            args.append(os.path.join(here, 'cases/regression_tests.py'))
        if not self.skip_demo_tests:
            args.append(os.path.join(here, 'demos/'))

        if (self.skip_unit_tests and self.skip_regression_tests and self.skip_demo_tests):
            print('WARNING: You skipped all tests!')
            sys.exit(0)
        else:
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
        'Development Status :: 4 - Beta',

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
        'ocellaris': ['cpp/*.h', 'cpp/*/*.h', 'input_file_schema.yml'],
    },
    zip_safe=False,

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'ocellaris=ocellaris.__main__:run_from_console',
            'ocellaris_inspector=ocellaris_post.inspector.__main__:main',
            'ocellaris_logstats=ocellaris_post.logstats:main',
        ],
    },

    # Configure the "test" command
    tests_require=['pytest'],
    cmdclass={'test': PyTest},
)
