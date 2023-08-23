# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import subprocess

# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "pysorb"
copyright = "2023, Marcus V. Pereira"
author = "Marcus V. Pereira"

master_doc = "index"

# -- General configuration ---------------------------------------------------

# generate doxygen xml output
subprocess.run(["doxygen", "Doxyfile.in"])
# generate rst from doxygen for each class, file, namespace
# subprocess.run(['breathe-apidoc', 'build/xml', '-o', 'developer/_auto'])

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "breathe",
    "myst_parser",
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
]

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "attrs_inline",
]

autoclass_content = "class"
# autofunction_content = "function"
autosummary_generate = True

nbsphinx_prolog = """
pysorb
"""

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build"]

# breathe: generate docs from doxygen xml output
breathe_projects = {"pysorb": "build/xml"}
breathe_default_project = "pysorb"
breathe_default_members = ("members", "undoc-members")

myst_update_mathjax = False
myst_dmath_double_inline = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

source_suffix = ['.rst', '.md']