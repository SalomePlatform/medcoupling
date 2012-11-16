# Copyright (C) 2007-2012  CEA/DEN, EDF R&D, OPEN CASCADE
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

# author: Guillaume Boulant (EDF/R&D, 2012)

include $(top_srcdir)/adm_local/unix/make_common_starter.am

tuidocdir = $(docdir)/tui/MED
devdocdir = $(docdir)/dev/MED
medopdocdir=$(devdocdir)

.PHONY : latex

if SPHINX_IS_OK

medopdoc_DATA=html/index.html
html/index.html:$(RSTFILES)
	make htm

endif

EXTRA_DIST += html

SPHINXOPTS      =
SOURCEDIR       = $(srcdir)
SPHINXBUILD     = sphinx-build
PAPEROPT_a4     = -D latex_paper_size=a4
ALLSPHINXOPTS   = -d doctrees $(PAPEROPT_a4) $(SPHINXOPTS) $(SOURCEDIR)
SPHINX_PYTHONPATH = $(prefix)/lib/python$(PYTHON_VERSION)/site-packages/salome:$(KERNEL_ROOT_DIR)/lib/python$(PYTHON_VERSION)/site-packages/salome:$(KERNEL_ROOT_DIR)/bin/salome:$(OMNIORB_ROOT)/lib/python$(PYTHON_VERSION)/site-packages

CONF_PATH = $(top_builddir)/src/MEDOP/doc/sphinx

SPHINX_LD_LIBRARY_PATH = $(OMNIORB_ROOT)/lib

htm:
	mkdir -p html doctrees
	PYTHONPATH=$(SPHINX_PYTHONPATH):${PYTHONPATH}; \
	LD_LIBRARY_PATH=$(SPHINX_LD_LIBRARY_PATH):${LD_LIBRARY_PATH}; \
	$(SPHINXBUILD) -c $(CONF_PATH) -b html $(ALLSPHINXOPTS) html
	@echo
	@echo "Build finished. The HTML pages are in html."

latex:
	mkdir -p latex doctrees
	PYTHONPATH=$(SPHINX_PYTHONPATH):${PYTHONPATH}; \
	LD_LIBRARY_PATH=$(SPHINX_LD_LIBRARY_PATH):${LD_LIBRARY_PATH}; \
	$(SPHINXBUILD) -c $(CONF_PATH) -b latex $(ALLSPHINXOPTS) latex
	@echo
	@echo "Build finished; the LaTeX files are in latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

html:
	mkdir -p $@

RSTFILES= \
	index.rst               \
	medop-specifications.rst \
	medop-develguide.rst     \
	medop-userguide.rst

# Files on the development of the prototype (version 2010)
RSTFILES+=\
	medop-prototype-overview.rst     \
	medop-prototype-develguide.rst   \
	medop-prototype-medmem.rst

# Files containing definitions and references
RSTFILES+=\
	medop-definitions.rst \
	medop-references.rst  \
	medop-workingnotes-2010.rst \
	medop-workingnotes-2011.rst \
	medop-workingnotes-2012.rst

EXTRA_DIST+= $(RSTFILES)

EXTRA_DIST+= \
	conf.py \
	_static \
	images

install-data-local:
	$(INSTALL) -d $(DESTDIR)$(medopdocdir)
	if test -d "html"; then b=; else b="$(srcdir)/"; fi; \
	cp -rf $$b"html"/* $(DESTDIR)$(medopdocdir) ; \
	if test -f $$b"latex"/medopdoc.pdf; then cp -f $$b"latex"/medopdoc.pdf $(DESTDIR)$(medopdocdir) ; fi;

uninstall-local:
	chmod -R +w $(DESTDIR)$(medopdocdir)
	rm -rf $(DESTDIR)$(medopdocdir)/*

clean-local:
	-chmod -R +w html latex doctrees
	rm -rf html latex doctrees
	if test -d "html"; then rm -rf html ; fi

pdf: latex
	cd latex; make all-pdf
	cp latex/*.pdf html/_static/documents/.

WIKIRELATIVEPATH="pub/userdata/gboulant/data/salome/smodules/xmed/doc"
WIKIREMOTEPATH="nepal@nepal:/nepal/$(WIKIRELATIVEPATH)"
deploy: htm pdf
	rsync -avz html/ $(WIKIREMOTEPATH)
	@echo "The documentation can be read at URL : http://nepal.der.edf.fr/$(WIKIRELATIVEPATH)/index.html"
