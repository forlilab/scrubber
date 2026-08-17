# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../molscrub/'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'molscrub'
copyright = '2026, ForliLab'
author = 'ForliLab'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx_design'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

pygments_style = 'sphinx'


html_logo = "images/logo_v3.png"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

autodoc_default_options = {
    'members': True,
    'special-members': '__call__',
    'undoc-members': False,
    'inherited-members': False,
    'show-inheritance': True,
}
autoclass_content = 'both'


html_theme_options = {
    'show_toc_level': 2,
    'repository_url': 'https://github.com/forlilab/molscrub',
    'use_repository_button': True,     # add a "link to repository" button
    'navigation_with_keys': False,
}
