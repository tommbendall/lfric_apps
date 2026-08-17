# Configuration file for the Sphinx documentation builder.
#
# -- imports -----------------------------------------------------------------
import glob
import os
import subprocess

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
    # Converts SVG figures to PDF (via Inkscape) for the LaTeX/PDF builder,
    # since pdflatex cannot include SVGs directly.
    'sphinx.ext.imgconverter',
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
#
# Everything from here to the end of the file implements the standalone
# Science Guide PDF builds: one branded, UMDP-style PDF per section under
# 'science_guide/' (see '_science_guide_sections' below), plus one combined
# "paper" containing all sections. This is a distinct, self-contained
# feature layered on top of the normal Sphinx HTML build above it -- none
# of it is required for (or affects) 'make html'.
#
# See 'source/_templates/latex/README.md' for how the branding template
# ('lfric_science.sty') works and how to swap in a different one, and
# 'documentation/Makefile's 'science_pdf' target for how to build a single
# section's PDF on its own (rather than every section via 'make latexpdf').
#
# Key mechanism: Sphinx builds every document listed in 'latex_documents'
# (all sections + the combined paper) in a single 'sphinx-build -b latex'
# pass, sharing one 'latex_elements' preamble. Since each PDF needs
# different title-page metadata (doc number, title, owner, etc.),
# '_science_guide_titlepage_tex()' generates LaTeX that switches on
# '\jobname' (which .tex file pdflatex is currently compiling) to select
# the right values for each document -- see that function's docstring.

latex_engine = 'pdflatex'

# The Science Guide PDF branding template can be switched by setting the
# LFRIC_SCIENCE_LATEX_STYLE environment variable to the name (without
# extension) of a different .sty file placed in _templates/latex/. See
# _templates/latex/README.md for details on swapping in an official template.
science_latex_style = os.environ.get('LFRIC_SCIENCE_LATEX_STYLE', 'lfric_science')

# The Met Office logo used on the title page of the combined Science Guide PDF.
latex_logo = '_static/MO_SQUARE_black_mono_for_light_backg_RBG.png'

# The Momentum logo used alongside it on the title page.
latex_momentum_logo = '_static/momentum_logo.png'

# Files that need to be copied into the LaTeX build directory.
latex_additional_files = [
    f"_templates/latex/{science_latex_style}.sty",
    latex_logo,
    latex_momentum_logo,
]


def _git_describe():
    '''Return a version string derived from the nearest git tag, used for the
    'LFRic-Apps Version' field on each Science Guide PDF's title page.'''
    try:
        return subprocess.check_output(
            ['git', 'describe', '--tags', '--always'],
            cwd=os.path.dirname(__file__),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return 'unknown'


def _git_last_updated(path):
    '''Return the date (YYYY-MM-DD) of the last commit that touched 'path',
    used for the 'Last Updated' field on a Science Guide PDF's title page.'''
    try:
        date = subprocess.check_output(
            ['git', 'log', '-1', '--format=%cd', '--date=short', '--', path],
            cwd=os.path.dirname(__file__),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        return date or 'unknown'
    except Exception:
        return 'unknown'


def _git_last_release():
    '''Return the most recent official release tag (e.g. 'vn3.2'), used for
    the 'LFRic-Apps Version' field on each Science Guide PDF's title page.
    Unlike '_git_describe()', this omits the commit-count/hash suffix added
    for commits after the tag, so it always shows the last official release
    rather than the exact development commit.'''
    try:
        return subprocess.check_output(
            ['git', 'describe', '--tags', '--abbrev=0'],
            cwd=os.path.dirname(__file__),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return 'unknown'


_lfric_apps_version = _git_last_release()

# Build one standalone PDF per Science Guide section, plus one combined
# "paper" containing all sections. New sections added under science_guide/
# are picked up automatically, with no further configuration required.
_science_guide_dir = os.path.join(os.path.dirname(__file__), 'science_guide')
_science_guide_sections = sorted(
    os.path.basename(os.path.dirname(path))
    for path in glob.glob(os.path.join(_science_guide_dir, '*', 'index.rst'))
)

# UMDP-style document numbers shown on each section's title page. Sections
# not listed here show 'TBD' until they are formally numbered.
_science_guide_doc_numbers = {
    'nudging': '083',
}

# Title-page owner/contributors for each section. Sections not listed default
# to a placeholder owner with no additional named contributors.
_science_guide_default_owner = 'TBD'
_science_guide_owners = {
    'nudging': 'Thomas Bendall',
}
_science_guide_contributors = {
    'nudging': ['Thomas Bendall', 'Mohit Dalvi'],
}


def _science_guide_titlepage_tex():
    '''Generate LaTeX that sets the title-page macros (doc number, title,
    version, last-updated, owner, contributors) for each Science Guide PDF,
    keyed on '\\jobname' (the name of the .tex file pdflatex is currently
    compiling). This is needed because 'latex_elements' is shared across
    every PDF built in a single `sphinx-build -b latex` run, so the same
    preamble text must select different values depending on which document
    is currently being built.
    '''
    blocks = []
    for section in _science_guide_sections:
        jobname = f"lfric_science_{section}"
        title = section.replace('_', ' ').title()
        number = _science_guide_doc_numbers.get(section, 'TBD')
        last_updated = _git_last_updated(
            os.path.join('science_guide', section))
        owner = _science_guide_owners.get(section, _science_guide_default_owner)
        contributors = _science_guide_contributors.get(section, [owner])
        contributors_tex = ', '.join(contributors)
        blocks.append(
            r'\ifdefstring{\lfricJobName}{%s}{%%' % jobname + '\n'
            r'\renewcommand{\lfricDocNumber}{%s}%%' % number + '\n'
            r'\renewcommand{\lfricSectionTitle}{%s}%%' % title + '\n'
            r'\renewcommand{\lfricVersion}{%s}%%' % _lfric_apps_version + '\n'
            r'\renewcommand{\lfricLastUpdated}{%s}%%' % last_updated + '\n'
            r'\renewcommand{\lfricOwner}{%s}%%' % owner + '\n'
            r'\renewcommand{\lfricContributors}{%s}%%' % contributors_tex + '\n'
            r'}{}'
        )

    # The combined "paper" document covering every section.
    paper_last_updated = _git_last_updated('science_guide')
    paper_owner = ', '.join(sorted(set(_science_guide_owners.values()))) \
        or _science_guide_default_owner
    paper_contributors = paper_owner
    blocks.append(
        r'\ifdefstring{\lfricJobName}{lfric_science_paper}{%' + '\n'
        r'\renewcommand{\lfricDocNumber}{N/A}%' + '\n'
        r'\renewcommand{\lfricSectionTitle}{LFRic Apps Science Guide}%' + '\n'
        r'\renewcommand{\lfricVersion}{%s}%%' % _lfric_apps_version + '\n'
        r'\renewcommand{\lfricLastUpdated}{%s}%%' % paper_last_updated + '\n'
        r'\renewcommand{\lfricOwner}{%s}%%' % paper_owner + '\n'
        r'\renewcommand{\lfricContributors}{%s}%%' % paper_contributors + '\n'
        r'}{}'
    )
    return '\n'.join(blocks)


# Disable Sphinx's default fancy chapter styling (fncychap), since the
# Science Guide template defines its own chapter/section styling. Colours
# match those used in the web documentation theme (see _static/custom.css).
# The title page (defined in the .sty file) is fully custom-branded; the
# per-document metadata it displays (doc number, version, last-updated,
# owner, contributors) is set here via '_science_guide_titlepage_tex()'.
latex_elements = {
    'fncychap': '',
    # Suppress Sphinx's own default font packages (which load plain,
    # unscaled 'tgtermes'/'tgheros'). 'lfric_science.sty' loads its own
    # TeX Gyre fonts ('tgpagella' body, 'tgheros' scaled for sans/headings);
    # leaving Sphinx's default in place causes an "Option clash for
    # package tgheros" error since it would be loaded twice with
    # different options.
    'fontpkg': '',
    'preamble': (
        r'\usepackage{%s}' % science_latex_style + '\n'
        + _science_guide_titlepage_tex()
    ),
    'sphinxsetup': ('TitleColor={HTML}{0F79BE}, '
                    'InnerLinkColor={HTML}{0F79BE}, '
                    'OuterLinkColor={HTML}{0F79BE}'),
    'maketitle': r'\sphinxmaketitle',
}

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
