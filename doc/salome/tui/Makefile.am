#  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
#
#  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
#  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
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
	sed 's|../../../share/salome|$(top_srcdir)|' ./doxyfile > ./doxyfile1; \
	sed 's|../../build/salome|$(top_builddir)|' ./doxyfile1 > ./doxyfile2; \
	mv -f doxyfile2 doxyfile1; \
	echo "DOXYGEN SUPPORT PYTHON - @DOXYGEN_WITH_PYTHON@"; \
	if( test "x@DOXYGEN_WITH_PYTHON@" = "xyes"); then \
	  sed 's|python_extension_must_be_here|*.py|' ./doxyfile1 > ./doxyfile2; \
	  mv -f doxyfile2 doxyfile1; \
	  $(DOXYGEN) -u ./doxyfile1; \
	else \
	  sed 's|python_extension_must_be_here||' ./doxyfile1 > ./doxyfile2; \
	  mv -f doxyfile2 doxyfile1; \
	fi; \
	mv -f doxyfile1 doxyfile; \
	$(DOXYGEN) ./doxyfile; \
	cd ..; \
	cp -fr $(srcdir)/MED/sources/static/*.* ./MED/
	cp -fr $(srcdir)/MED/sources/ MED/
	cp -fr $(srcdir)/MED/HTML/ MED/
	rm -fr INPUT

#doctuidir= $(docdir)/tui/MED
#nodist_doctui_DATA= MED/doxyfile

#doctuistaticdir= $(docdir)/tui/MED/sources/static
#nodist_doctuistatic_DATA= MED/sources/static/tree.js
