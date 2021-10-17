from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages


with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='pyscal',
    version='2.10.15',
    author='Sarath Menon',
    author_email='sarath.menon@pyscal.org',
    description='Python library written in C++ for calculation of local atomic structural environment',
    long_description=readme,
    # tell setuptools to look for any packages under 'src'
    packages=find_packages('src'),
    # tell setuptools that all packages will be under the 'src' directory
    # and nowhere else
    package_dir={'':'src'},
    # add an extension module named 'python_cpp_example' to the package
    headers=["src/pyscal/atom.h", "src/pyscal/system.h", "lib/voro++/voro++.hh"],
    ext_modules=[
        Pybind11Extension(
            "pyscal.catom",
            ["src/pyscal/atom.cpp", "src/pyscal/atom_binding.cpp"],
            language='c++',
            include_dirs=['lib/voro++']
        ),
        Pybind11Extension(
            "pyscal.csystem",
            ["src/pyscal/system.cpp", "src/pyscal/system_binding.cpp", "src/pyscal/atom.cpp", "lib/voro++/voro++.cc"],
            language='c++',
            include_dirs=['lib/voro++']
        ),
    ],
    # add custom build_ext command
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    download_url = 'https://anaconda.org/conda-forge/pyscal',
    url = 'https://pyscal.org',
    install_requires=['numpy', 'ase', 'plotly', 'ipywidgets'],
    classifiers=[
        'Programming Language :: Python :: 3'
    ]
)
