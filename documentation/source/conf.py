# Configuration file for the Sphinx documentation builder.
#
# -- imports -----------------------------------------------------------------
import glob
import os

import docutils

# -- Project information -----------------------------------------------------

project = 'LFRic Apps'
author = 'Simulation IT'
copyright = 'Met Office'
release = '0.1.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_sitemap',
    'sphinx_design',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
]

# Enable equation referencing and cross-referencing
mathjax3_config = { "tex": { "tags": "ams", "packages": {"[+]": ["ams"]}, } }

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'
# html_title = "LFRic Apps"

# Generate the sitemap info, this will need updating when we have versioned docs
html_baseurl = os.environ.get("SPHINX_HTML_BASE_URL", "https://metoffice.github.io/lfric_apps")
sitemap_locales = [None]
sitemap_url_scheme = "{link}"

# Hide the link which shows the rst markup
html_show_sourcelink = False

html_theme_options = {
    "announcement": "This documentation is under construction. "
                    "Thank you for your patience while we add content!",
    "navigation_with_keys": True,
    "use_edit_page_button": True,
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_align": "content",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/MetOffice/lfric_apps",
            "icon": "fa-brands fa-github"
        },
        {
            "name": "GitHub Discussions",
            "url": "https://github.com/MetOffice/simulation-systems/discussions",
            "icon": "far fa-comments",
        }
    ],
    "logo": {
        "text": "LFRic Apps",
        "image_light": "_static/MO_SQUARE_black_mono_for_light_backg_RBG.png",
        "image_dark": "_static/MO_SQUARE_for_dark_backg_RBG.png",
    },
    "secondary_sidebar_items": {
        "**/*": ["page-toc", "edit-this-page", "show-glossary"],
        "index": [],
    },
    "footer_start": ["crown-copyright"],
    "footer_center": ["show-accessibility"],
    "footer_end": ["sphinx-version", "theme-version"],
    "primary_sidebar_end": []
}

html_sidebars = {
    "index": []
}

# Provides the Edit on GitHub link in the generated docs.
html_context = {
    "display_github": True,
    "github_user": "MetOffice",
    "github_repo": "lfric_apps",
    "github_version": "main",
    "doc_path": "documentation/source"
}

# Enable numbered references to e.g. figures.
#
numfig = True

# Exclude files from Sphinx processing
exclude_patterns = ['common_links.rst']

# Create rst_epilog variable to allow concatenation of epilog parts which will
# be included at the end of every rst file.
rst_prolog = ""

# Add the contents of the common_links file to the epilog.
with open('common_links.rst') as file:
    rst_prolog += file.read()


def superscript_substitution_role(name, rawtext, text, lineno, inliner,
                                  options={}, content=[]):
    node = docutils.nodes.superscript()
    node2 = docutils.nodes.substitution_reference(refname=text)
    node += [node2]
    return [node], []


def setup(app):
    app.add_role('superscript_substitution', superscript_substitution_role)

# Options for intersphinx mapping extension
# To discover the available objects in the target project (e.g. psyclone), run:
# python -m sphinx.ext.intersphinx https://psyclone.readthedocs.io/en/stable/objects.inv
intersphinx_mapping = {
    'psyclone': ('https://psyclone.readthedocs.io/en/stable/', None),
    'simsys': ('https://metoffice.github.io/simulation-systems/', None),
    'lfric_core': ('https://metoffice.github.io/lfric_core/', None)
}

# -- Options for LaTeX/PDF output ---------------------------------------------

latex_engine = 'pdflatex'

# The Science Guide PDF branding template can be switched by setting the
# LFRIC_SCIENCE_LATEX_STYLE environment variable to the name (without
# extension) of a different .sty file placed in _templates/latex/. See
# _templates/latex/README.md for details on swapping in an official template.
science_latex_style = os.environ.get('LFRIC_SCIENCE_LATEX_STYLE', 'lfric_science')

# The Met Office logo used on the title page of the combined Science Guide PDF.
latex_logo = '_static/MO_SQUARE_black_mono_for_light_backg_RBG.png'

# Files that need to be copied into the LaTeX build directory.
latex_additional_files = [
    f"_templates/latex/{science_latex_style}.sty",
    latex_logo,
]

# Disable Sphinx's default fancy chapter styling (fncychap), since the
# Science Guide template defines its own chapter/section styling.
# Colours match those used in the web documentation theme (see
# _static/custom.css), and the title page carries a disclaimer noting the
# template is provisional (defined in the .sty file, inserted here via the
# 'maketitle' key).
latex_elements = {
    'fncychap': '',
    'preamble': r'\usepackage{%s}' % science_latex_style,
    'sphinxsetup': ('TitleColor={HTML}{0F79BE}, '
                    'InnerLinkColor={HTML}{0F79BE}, '
                    'OuterLinkColor={HTML}{0F79BE}'),
    'maketitle': (r'\newcommand\sphinxbackoftitlepage{\lfricScienceDisclaimer}'
                  r'\sphinxmaketitle'),
}

# Build one standalone PDF per Science Guide section, plus one combined
# "paper" containing all sections. New sections added under science_guide/
# are picked up automatically, with no further configuration required.
_science_guide_dir = os.path.join(os.path.dirname(__file__), 'science_guide')
_science_guide_sections = sorted(
    os.path.basename(os.path.dirname(path))
    for path in glob.glob(os.path.join(_science_guide_dir, '*', 'index.rst'))
)

latex_documents = [
    (
        f"science_guide/{section}/index",
        f"lfric_science_{section}.tex",
        f"LFRic Apps Science Guide: {section.replace('_', ' ').title()}",
        author,
        'howto',
    )
    for section in _science_guide_sections
]

# The combined "paper" is built from a separate, hand-curated toctree
# (science_guide/paper_index.rst), not the glob-based one used for web
# navigation, so that new/placeholder sections aren't pulled in automatically.
latex_documents.append(
    (
        'science_guide/paper_index',
        'lfric_science_paper.tex',
        'LFRic Apps Science Guide',
        author,
        'manual',
    )
)
