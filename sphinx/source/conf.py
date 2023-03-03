# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FEM.edu'
copyright = '2023, Peter Mackenzie-Helnwein'
author = 'Peter Mackenzie-Helnwein'
release = '0.4'

# -- add current folder to path
# https://www.jetbrains.com/pycharm/guide/tutorials/sphinx_sites/documentation/
import os
import sys
# sys.path.insert(0, os.path.abspath("../../src/elements"))
# sys.path.insert(0, os.path.abspath("../../src/materials"))
# sys.path.insert(0, os.path.abspath("../../src/solver"))
sys.path.insert(0, os.path.abspath("../../src"))
sys.path.insert(0, os.path.abspath("."))


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autosectionlabel',
    "myst_parser",
    "sphinx.ext.autodoc",
    ]

extensions.append('sphinx_git')

templates_path = ['_templates']
exclude_patterns = []

autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 6



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

