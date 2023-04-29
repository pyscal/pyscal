from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='pyscal3',
    version='3.0.0',
    author='Sarath Menon',
    author_email='sarath.menon@pyscal.org',
    description='Python library written in C++ for calculation of local atomic structural environment',
    long_description=readme,
    long_description_content_type = "text/markdown",
    # tell setuptools to look for any packages under 'src'
    packages=find_packages('src'),
    # tell setuptools that all packages will be under the 'src' directory
    # and nowhere else
    package_dir={'':'src'},
    # add an extension module named 'python_cpp_example' to the package
    #headers=["src/pyscal/atom.h", "src/pyscal/system.h", "lib/voro++/voro++.hh"],
    #headers=["src/pyscal/atom.h", "src/pyscal/modsystem.h"],
    headers=["src/pyscal/modsystem.h"],
    ext_modules=[
        #Pybind11Extension(
        #    "pyscal.catom",
        #    ["src/pyscal/atom.cpp", "src/pyscal/atom_binding.cpp"],
        #    language='c++',
        #    include_dirs=['lib/voro++']
        #),
        Pybind11Extension(
            "pyscal.csystem",
            ["src/pyscal/neighbor.cpp", "src/pyscal/sh.cpp", 
            "src/pyscal/solids.cpp", "src/pyscal/voronoi.cpp", 
            "src/pyscal/modsystem_binding.cpp", "lib/voro++/voro++.cc"],
            language='c++',
            include_dirs=['lib/voro++'],
            extra_compile_args=['-O3'],
        ),

        #Pybind11Extension(
        #    "pyscal.csystem",
        #    ["src/pyscal/system.cpp", "src/pyscal/system_binding.cpp", "src/pyscal/atom.cpp", "lib/voro++/voro++.cc"],
        #    language='c++',
        #    include_dirs=['lib/voro++']
        #),
    ],
    # add custom build_ext command
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    download_url = 'https://anaconda.org/conda-forge/pyscal',
    url = 'https://pyscal.org',
    install_requires=['pybind11', 'numpy', 'ase'],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    include_package_data=True,
    #package_data={'': ['*.yaml']},
)
