"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages, Extension, Command
import os
import warnings

# io.open is needed for projects that support Python 2.7
# It ensures open() defaults to text mode with universal newlines,
# and accepts an argument to specify the text encoding
# Python 3 only projects can skip this import
from io import open

# stuff to check presence of compiler:
import distutils.sysconfig
import distutils.ccompiler
from distutils.errors import CompileError, DistutilsPlatformError

import itertools

from rsciio.version import __version__

setup_path = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(setup_path, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

# Extensions. Add your extension here:
raw_extensions = [
    Extension(
        "rsciio.bruker.unbcf_fast", [os.path.join("rsciio", "bruker", "unbcf_fast.pyx")]
    ),
]

cleanup_list = []
for leftover in raw_extensions:
    path, ext = os.path.splitext(leftover.sources[0])
    if ext in (".pyx", ".py"):
        cleanup_list.append("".join([os.path.join(setup_path, path), ".c*"]))
        if os.name == "nt":
            bin_ext = ".cpython-*.pyd"
        else:
            bin_ext = ".cpython-*.so"
        cleanup_list.append("".join([os.path.join(setup_path, path), bin_ext]))


def count_c_extensions(extensions):
    c_num = 0
    for extension in extensions:
        # if first source file with extension *.c or *.cpp exists
        # it is cythonised or pure c/c++ extension:
        sfile = extension.sources[0]
        path, ext = os.path.splitext(sfile)
        if os.path.exists(path + ".c") or os.path.exists(path + ".cpp"):
            c_num += 1
    return c_num


def cythonize_extensions(extensions):
    try:
        from Cython.Build import cythonize

        return cythonize(extensions, language_level="3")
    except ImportError:
        warnings.warn(
            """WARNING: cython required to generate fast c code is not found on this system.
Only slow pure python alternative functions will be available.
To use fast implementation of some functions writen in cython either:
a) install cython and re-run the installation,
b) try alternative source distribution containing cythonized C versions of fast code,
c) use binary distribution (i.e. wheels, egg)."""
        )
        return []


def no_cythonize(extensions):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


# to cythonize, or not to cythonize... :
if len(raw_extensions) > count_c_extensions(raw_extensions):
    extensions = cythonize_extensions(raw_extensions)
else:
    extensions = no_cythonize(raw_extensions)


# to compile or not to compile... depends if compiler is present:
compiler = distutils.ccompiler.new_compiler()
assert isinstance(compiler, distutils.ccompiler.CCompiler)
distutils.sysconfig.customize_compiler(compiler)
try:
    compiler.compile(
        [os.path.join(setup_path, "rsciio", "tests", "bruker_data", "test_compilers.c")]
    )
except (CompileError, DistutilsPlatformError):
    warnings.warn(
        """WARNING: C compiler can't be found.
Only slow pure python alternative functions will be available.
To use fast implementation of some functions writen in cython/c either:
a) check that you have compiler (EXACTLY SAME as your python
distribution was compiled with) installed,
b) use binary distribution of hyperspy (i.e. wheels, egg, (only osx and win)).
Installation will continue in 5 sec..."""
    )
    extensions = []
    from time import sleep

    sleep(5)  # wait 5 secs for user to notice the message


class Recythonize(Command):

    """cythonize all extensions"""

    description = "(re-)cythonize all changed cython extensions"

    user_options = []

    def initialize_options(self):
        """init options"""
        pass

    def finalize_options(self):
        """finalize options"""
        pass

    def run(self):
        # if there is no cython it is supposed to fail:
        from Cython.Build import cythonize

        global raw_extensions
        global extensions
        cythonize(extensions, language_level="3")


install_requires = [
    "dask[array]>=2.11",
    "python-dateutil",
    "h5py>=2.3",
    "imageio",
    "numba>=0.52",
    "numpy>=1.17.1",
    "pint>=0.8",
    "python-box>=6.0,<7.0",
    "pyyaml",
    "scipy>=1.1",
    "sparse",
]

extras_require = {
    "blockfile": ["scikit-image"],
    "mrcz": ["blosc>=1.5", "mrcz>=0.3.6"],
    "scalebar_export": ["matplotlib-scalebar", "matplotlib>=3.1.3"],
    "tiff": ["tifffile>=2020.2.16", "imagecodecs>=2020.1.31"],
    "usid": ["pyUSID"],
    "zspy": ["zarr"],
    "tests": [
        "pytest>=3.6",
        "pytest-xdist",
        "pytest-rerunfailures",
        "pytest-cov",
    ],  # for testing
    "docs": [
        "pydata-sphinx-theme",
        "sphinxcontrib-towncrier",
        # TODO: Remove explicit dependency on sphinx when pydata-sphinx-theme >= 0.13
        #  is available, and 0.13 as minimial supported version of pydata-sphinx-theme
        "sphinx~=5.3",
        # pin towncrier until https://github.com/sphinx-contrib/sphinxcontrib-towncrier/issues/60 is fixed
        "towncrier<22.8",
    ],  # for building the docs
}

# Don't include "tests" and "docs" requirements since "all" is designed to be
# used for user installation.
runtime_extras_require = {
    x: extras_require[x] for x in extras_require.keys() if x not in ["tests", "docs"]
}
extras_require["all"] = list(itertools.chain(*list(runtime_extras_require.values())))

extras_require["dev"] = list(itertools.chain(*list(extras_require.values())))


# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    # This is the name of your project. The first time you publish this
    # package, this name will be registered for you. It will determine how
    # users can install this project, e.g.:
    #
    # $ pip install sampleproject
    #
    # And where it will live on PyPI: https://pypi.org/project/sampleproject/
    #
    # There are some restrictions on what makes a valid project name
    # specification here:
    # https://packaging.python.org/specifications/core-metadata/#name
    name="RosettaSciIO",  # Required
    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    #
    # For a discussion on single-sourcing the version across setup.py and the
    # project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,  # Required
    # This is a one-line description or tagline of what your project does. This
    # corresponds to the "Summary" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#summary
    description="Scientific file formats",  # Optional
    # This is an optional longer description of your project that represents
    # the body of text which users will see when they visit PyPI.
    #
    # Often, this is the same as your README, so you can just read it in from
    # that file directly (as we have already done above)
    #
    # This field corresponds to the "Description" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-optional
    long_description=long_description,  # Optional
    # Denotes that our long_description is in Markdown; valid values are
    # text/plain, text/x-rst, and text/markdown
    #
    # Optional if long_description is written in reStructuredText (rst) but
    # required for plain-text or Markdown; if unspecified, "applications should
    # attempt to render [the long_description] as text/x-rst; charset=UTF-8 and
    # fall back to text/plain if it is not valid rst" (see link below)
    #
    # This field corresponds to the "Description-Content-Type" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-content-type-optional
    long_description_content_type="text/markdown",  # Optional (see note above)
    # This should be a valid link to your project's main homepage.
    #
    # This field corresponds to the "Home-Page" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#home-page-optional
    url="https://github.com/hyperspy/rosettasciio",  # Optional
    # This should be your name or the name of the organization which owns the
    # project.
    author="The HyperSpy Develoopers",  # Optional
    # This should be a valid email address corresponding to the author listed
    # above.
    # author_email='pypa-dev@googlegroups.com',  # Optional
    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",
        # Indicate who your project is intended for
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        # Pick your license as you wish
        "License :: OSI Approved :: MIT License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # These classifiers are *not* checked by 'pip install'. See instead
        # 'python_requires' below.
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    #
    # Note that this is a string of words separated by whitespace, not a list.
    keywords="IO scientific format",  # Optional
    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    ext_modules=extensions,  # For Cython code
    cmdclass={
        "recythonize": Recythonize,
    },
    packages=find_packages(exclude=["contrib", "docs", "tests"]),  # Required
    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. If you
    # do not support Python 2, you can simplify this to '>=3.5' or similar, see
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires=">=3.6, <4",
    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=install_requires,  # Optional
    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    # Similar to `install_requires` above, these must be valid existing
    # projects.
    extras_require=extras_require,
    # If there are data files included in your packages that need to be
    # installed, specify them here.
    #
    # If using Python 2.6 or earlier, then these have to be included in
    # MANIFEST.in as well.
    package_data={
        "rsciio": [
            "*/specifications.yaml",
            "tests/blockfile_data/*.blo",
            "tests/dens_data/*.dens",
            "tests/dm_stackbuilder_plugin/test_stackbuilder_imagestack.dm3",
            "tests/dm3_1D_data/*.dm3",
            "tests/dm3_2D_data/*.dm3",
            "tests/dm3_3D_data/*.dm3",
            "tests/dm4_1D_data/*.dm4",
            "tests/dm4_2D_data/*.dm4",
            "tests/dm4_3D_data/*.dm4",
            "tests/dm3_locale/*.dm3",
            "tests/FEI_new/*.emi",
            "tests/FEI_new/*.ser",
            "tests/FEI_old/*.emi",
            "tests/FEI_old/*.ser",
            "tests/FEI_old/*.npy",
            "tests/FEI_old/*.tar.gz",
            "tests/impulse_data/*.csv",
            "tests/impulse_data/*.log",
            "tests/impulse_data/*.npy",
            "tests/msa_files/*.msa",
            "tests/hdf5_files/*.hdf5",
            "tests/hdf5_files/*.hspy",
            "tests/JEOL_files/*",
            "tests/JEOL_files/Sample/00_View000/*",
            "tests/JEOL_files/InvalidFrame/*",
            "tests/JEOL_files/InvalidFrame/Sample/00_Dummy-Data/*",
            "tests/tiff_files/*.zip",
            "tests/tiff_files/*.tif",
            "tests/tiff_files/*.tif.gz",
            "tests/tiff_files/*.dm3",
            "tests/tvips_files/*.tvips",
            "tests/npz_files/*.npz",
            "tests/unf_files/*.unf",
            "tests/bruker_data/*.bcf",
            "tests/bruker_data/*.json",
            "tests/bruker_data/*.npy",
            "tests/bruker_data/*.spx",
            "tests/ripple_files/*.rpl",
            "tests/ripple_files/*.raw",
            "tests/panta_rhei_files/*.prz",
            "tests/phenom_data/*.elid",
            "tests/sur_data/*.sur",
            "tests/sur_data/*.pro",
            "tests/emd_files/*.emd",
            "tests/emd_files/fei_emd_files.zip",
            "tests/protochips_data/*.npy",
            "tests/protochips_data/*.csv",
            "tests/nexus_files/*.nxs",
            "tests/empad_data/*.xml",
        ]
    },
    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # https://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    #
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],  # Optional
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # `pip` to create the appropriate form of executable for the target
    # platform.
    #
    # For example, the following would provide a command called `sample` which
    # executes the function `main` from this package when invoked:
    # entry_points={  # Optional
    #     'console_scripts': [
    #         'sample=sample:main',
    #     ],
    # },
    # List additional URLs that are relevant to your project as a dict.
    #
    # This field corresponds to the "Project-URL" metadata fields:
    # https://packaging.python.org/specifications/core-metadata/#project-url-multiple-use
    #
    # Examples listed include a pattern for specifying where the package tracks
    # issues, where the source is hosted, where to say thanks to the package
    # maintainers, and where to support the project financially. The key is
    # what's used to render the link text on PyPI.
    project_urls={  # Optional
        "Bug Reports": "https://github.com/hyperspy/rosettasciio/issues",
        # 'Funding': 'https://donate.pypi.org',
        # 'Say Thanks!': 'https://saythanks.io/to/example',
        "Source": "https://github.com/hyperspy/rosettasciio/",
    },
)
