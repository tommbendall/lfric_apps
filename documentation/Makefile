# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	
.PHONY: clean
clean:
	@echo "***** Removing, $(BUILDDIR) directory *****"
	@rm -rf $(BUILDDIR)

# Build a single Science Guide section's standalone PDF, e.g.:
#   make science SECTION=nudging
# 'SECTION' must match a directory name under 'source/science_guide/' (and
# thus a 'lfric_science_<SECTION>.tex'/.pdf jobname, per conf.py). This
# still runs the full LaTeX build (all sections' .tex sources are generated
# together), but then compiles only the requested section's PDF via
# latexmk, rather than every section's PDF.
.PHONY: science_pdf
science_pdf:
	@if [ -z "$(SECTION)" ]; then \
		echo "Usage: make science_pdf SECTION=<section-name>  (e.g. make science_pdf SECTION=nudging)"; \
		exit 1; \
	fi
	@$(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex" $(SPHINXOPTS) $(O)
	@$(MAKE) -C "$(BUILDDIR)/latex" lfric_science_$(SECTION).pdf
	@echo "***** Science Guide PDF built: $(BUILDDIR)/latex/lfric_science_$(SECTION).pdf *****"

