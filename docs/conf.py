# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ORFanage'
copyright = '2023, Ales Varabyou'
author = 'Ales Varabyou'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
#html_static_path = ['_static']
#
html_theme_options = {
    "logo": "./logo.png",
    "logo_text_align": "center",
    "description": "Logo designed by Julia Wang",
    "github_user": "alevar",
    "github_repo": "ORFanage",
    "travis_button": False,  # Circle now
    "codecov_button": False  # README badge now
}
