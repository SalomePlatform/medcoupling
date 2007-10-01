#  Copyright (C) 2003  CEA/DEN, EDF R&D
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

include $(top_srcdir)/adm_local/unix/make_common_starter.am

SUBDIRS = MED
#EXTRA_DIST+= MED

dev_docs:
	cp -fr $(srcdir)/MED ./INPUT; \
	cp -fr ./MED/doxyfile ./INPUT; \
	cp -fr ./MED/sources/static/tree.js ./INPUT/sources/static; \
	cd INPUT; \
	sed 's|../../../share/salome|$(root_srcdir)|' ./doxyfile > ./doxyfile1; \
	sed 's|../../build/salome|$(top_builddir)|' ./doxyfile1 > ./doxyfile2; \
	mv -f doxyfile2 doxyfile1; \
	echo "DOXYGEN SUPPORT PYTHON - @DOXYGEN_WITH_PYTHON@"; \
	if( test "x@DOXYGEN_WITH_PYTHON@" = "xyes"); then \
	  sed 's|python_extension_must_be_here|*.py|' ./doxyfile1 > ./doxyfile2; \
	  mv -f doxyfile2 doxyfile1; \
	  $(doxygen) -u ./doxyfile1; \
	else \
	  sed 's|python_extension_must_be_here||' ./doxyfile1 > ./doxyfile2; \
	  mv -f doxyfile2 doxyfile1; \
	fi; \
	mv -f doxyfile1 doxyfile; \
	$(doxygen) ./doxyfile; \
	cd ..; \
	cp -fr $(srcdir)/MED/sources/static/*.* ./MED/
	cp -fr $(srcdir)/MED/sources/ MED/
	cp -fr $(srcdir)/MED/HTML/ MED/
	rm -fr INPUT

#doctuidir= $(docdir)/tui/MED
#nodist_doctui_DATA= MED/doxyfile

#doctuistaticdir= $(docdir)/tui/MED/sources/static
#nodist_doctuistatic_DATA= MED/sources/static/tree.js
