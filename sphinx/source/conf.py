# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project   = 'FEM.edu'
copyright = '2023, Peter Mackenzie-Helnwein'
author    = 'Peter Mackenzie-Helnwein'
release   = 'FEM.edu-v1.0.1'

# -- add current folder to path
# https://www.jetbrains.com/pycharm/guide/tutorials/sphinx_sites/documentation/
import os
import sys
# sys.path.insert(0, os.path.abspath("../../src/elements"))
# sys.path.insert(0, os.path.abspath("../../src/materials"))
# sys.path.insert(0, os.path.abspath("../../src/solver"))
sys.path.insert(0, os.path.abspath("../../src"))
sys.path.insert(0, os.path.abspath("."))

from sphinx_gallery.sorting import FileNameSortKey, ExplicitOrder

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions  = ['sphinx.ext.autosectionlabel']
extensions += ["myst_parser"]
extensions += ["sphinx.ext.autodoc"]
extensions += ['sphinx_git']
extensions += ['sphinx_substitution_extensions']
extensions += ['sphinx_gallery.gen_gallery']

templates_path = ['_templates']
exclude_patterns = []

autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 6

# -- Sphinx-Gallery configuration
sphinx_gallery_conf = {
    'examples_dirs': ['../../galleries/examples', '../../galleries/tutorials'],   # path to your example scripts
    'gallery_dirs': ['auto_examples','auto_tutorials'],  # path to where to save gallery generated output
    'ignore_pattern': r'__init__\.py|runall\.py|read_profile_stats.py',
    # specify that examples should be ordered according to filename
    'within_subsection_order': FileNameSortKey,
    'subsection_order': ExplicitOrder(['../../galleries/examples/beams',
                                       '../../galleries/examples/trusses',
                                       '../../galleries/examples/frames',
                                       '../../galleries/examples/plates',
                                       '../../galleries/examples/solids',
                                       '../../galleries/examples/mixed',
                                       '../../galleries/tutorials/model_creation',
                                       '../../galleries/tutorials/direct_stiffness',
                                       '../../galleries/tutorials/linear_fem',
                                       '../../galleries/tutorials/large_deformation',
                                       '../../galleries/tutorials/plasticity',
                                       ]),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# define global abbreviations
rst_prolog = """
.. |Application|  replace:: **FEM.edu**
.. |PackageName|  replace:: **FEM.edu**
.. |Author|   replace:: Peter Mackenzie-Helnwein
"""

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'
