# Science Guide PDFs

This directory contains the LaTeX style file (`lfric_science.sty`) used to
brand the standalone PDF builds of the Science Guide: one PDF per section
under `documentation/source/science_guide/`, plus a combined "paper"
containing every section. These are built from the same `.rst` content as
the web documentation, via Sphinx's LaTeX/PDF builder.

This README is aimed at a developer who wants to build these PDFs, add a
new section to them, or change the branding template.

## Building a Science Guide PDF

No Python virtual environment is required, provided the LFRic and LaTeX
modules are loaded, for example:

```bash
cd documentation
ml use ~lfricadmin/lmod
module load lfric
module load latex
make clean
make science_pdf SECTION=nudging
```

`SECTION` must match a directory name under `source/science_guide/` (e.g.
`nudging`). The resulting PDF is written to
`documentation/build/latex/lfric_science_<SECTION>.pdf`.

`make science_pdf` builds only the requested section's PDF. To build every
section's PDF plus the combined paper in one go, use Sphinx's own target
instead:

```bash
make latexpdf
```

which writes all PDFs to `documentation/build/latex/`.

## Adding a new Science Guide section

New sections under `source/science_guide/<SECTION>/index.rst` are picked up
automatically as a standalone PDF target -- no changes are needed in
`conf.py` just to build one. However, a few optional per-section settings
in `documentation/source/conf.py` control what appears on that section's
title page:

- `_science_guide_doc_numbers`: the UMDP-style document number shown next
  to the title (e.g. `'nudging': '083'`). Sections not listed show `TBD`.
- `_science_guide_owners`: the named owner shown on the title page.
  Sections not listed show `_science_guide_default_owner` (`'TBD'`).
- `_science_guide_contributors`: the list of contributors shown on the
  title page. Sections not listed default to just the owner.

To add a new section to the combined "paper" PDF, also add it to the
hand-curated toctree in `source/science_guide/paper_index.rst` (this is
separate from the automatic per-section glob, so that new/placeholder
sections aren't pulled into the paper before they're ready).

The "LFRic-Apps Version" and "Last Updated" fields are generated
automatically (from the latest git release tag, and the last commit date
touching that section's directory, respectively) and don't need any
per-section configuration.

## Current template

`lfric_science.sty` is the template. It provides:

- brand colours matching the web documentation theme
  (`documentation/source/_static/custom.css`)
- styled chapter/section/subsection headings
- a branded header and footer on every page, including a Crown copyright
  notice and page number
- a custom title page (logos, a banner with the document title/number, the
  section title, key metadata, and the contributors list), modelled on the
  UMDP paper format

**This is a provisional template**, intended to be replaced by a final LaTeX
template once one becomes available -- see "Switching to a different template" below.

It is intentionally kept in a single, isolated `.sty` file so that it can be
replaced without needing to change any of the science content or the way
`conf.py` builds the PDF targets.

## Switching to a different template

The template in use is controlled by a single setting in
`documentation/source/conf.py`:

```python
science_latex_style = os.environ.get("LFRIC_SCIENCE_LATEX_STYLE", "lfric_science")
```

To switch to an official (or otherwise different) template:

1. Add the new `.sty` file (and any accompanying assets, e.g. a logo) to this
   `_templates/latex/` directory.
2. Either:
   - rename the new file to `lfric_science.sty`, replacing this one, or
   - give it a different name and set the `LFRIC_SCIENCE_LATEX_STYLE`
     environment variable (or edit the default in `conf.py`) to point to it.
3. Update `latex_additional_files` in `conf.py` if the new template requires
   extra files (e.g. fonts or images) to be copied into the LaTeX build
   directory.
4. Rebuild the PDFs with `make latexpdf` and check the output.

No changes to the `.rst` science content are required when switching
templates.

## Implementation notes / gotchas

These are recorded here (and alongside the relevant code in `conf.py` and
`lfric_science.sty`) since the standalone-PDF capability is a distinct,
self-contained addition to the normal Sphinx HTML build, and is easy to
break if these points are missed:

- **Per-document title-page metadata**: Sphinx builds every PDF listed in
  `latex_documents` (all sections + the combined "paper") in a single
  `sphinx-build -b latex` pass, sharing one `latex_elements['preamble']`.
  `conf.py`'s `_science_guide_titlepage_tex()` works around this by
  generating LaTeX that switches on `\jobname` (captured as `\lfricJobName`
  in `lfric_science.sty`, since `\ifdefstring` can't compare `\jobname`
  directly) to select the right doc number/title/owner/etc. per PDF.
- **`fontpkg` must be set to `''`** in `conf.py`'s `latex_elements`. Sphinx's
  own default `fontpkg` loads plain `tgheros`; `lfric_science.sty` loads it
  again with `scale=0.92`, so leaving Sphinx's default in place causes
  `Option clash for package tgheros`.
- **tex-gyre packages take `scale=`, not `scaled=`**. `scaled=` is the older
  PSNFSS/`helvet` convention and isn't a key `tgheros.sty` understands;
  using it raises `Package keyval Error: scaled undefined`.
- **SVG figures need `sphinx.ext.imgconverter`** (listed in `conf.py`'s
  `extensions`). `pdflatex` can't embed `.svg` directly; without this
  extension (which uses Inkscape) the LaTeX build fails with
  `Unknown graphics extension: .svg`.
- **`science_pdf` still runs the full `-b latex` build** before compiling
  just the requested section's PDF, since Sphinx generates every section's
  `.tex` file together in one pass -- there's no per-section Sphinx build.
  Only the final `latexmk` compile step is limited to one section.

