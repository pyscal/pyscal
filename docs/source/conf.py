# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import sphinx_rtd_theme
import shutil

sys.path.insert(0, os.path.abspath('../../src/pyscal/'))

def skip(app, what, name, obj, would_skip, options):
    if name in ( '__init__',):
        return False
    return would_skip  
def setup(app):
    app.connect('autodoc-skip-member', skip) 

#copy ipynb here
if os.path.exists("examples"):
    shutil.rmtree("examples")
shutil.copytree("../../examples", "examples")

# -- Project information -----------------------------------------------------

project = u'pyscal'
copyright = u'2019, Sarath Menon'
author = u'Sarath Menon'

# The short X.Y version
version = u'3.0.0'
# The full version, including alpha/beta/rc tags
release = u'3.0.0'


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    "sphinx.ext.extlinks", 
    "sphinx.ext.intersphinx", 
    "sphinx.ext.todo", 
    "sphinx.ext.viewcode",
    "myst_parser", 
    "sphinx_copybutton", 
    "sphinx_design", 
    "sphinx_inline_tabs",
    "nbsphinx",
]

myst_enable_extensions = ["dollarmath", "amsmath"]
html_theme = 'furo'
html_extra_path = ['../_static' ]
html_static_path = ["../_static"]

html_theme_options = {
    "light_logo": "pyscal_logo1.png",
    "dark_logo": "pyscal_logo2.png",
}


source_suffix = ['.rst', '.md']
nbsphinx_execute = 'never'
nbsphinx_allow_errors = True

exclude_patterns = []