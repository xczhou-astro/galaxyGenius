# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
project = 'GalaxyGenius'
copyright = '2024, Your Name'
author = 'Your Name'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'myst_parser'
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------
autodoc_member_order = 'bysource'

# GitHub Pages Configuration
html_baseurl = 'https://xczhou-astro.github.io/galaxygenius/'
html_context = {
    'display_github': True,
    'github_user': 'xczhou-astro',
    'github_repo': 'galaxygenius',
    'github_version': 'main/docs/source/',
}