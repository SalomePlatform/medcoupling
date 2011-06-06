# Copyright (C) 2007-2011  CEA/DEN, EDF R&D, OPEN CASCADE
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

#  MED MEDMEM : MED files in memory
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

DOX_INPUT_FILE = Doxyfile_med_user

guidocdir = $(docdir)/gui/MED
guidoc_DATA = images/head.png

install-data-local : html-local
	@if test -d doc_ref_user; then 					\
	$(INSTALL) -d $(DESTDIR)$(docdir)/gui/MED;	     		\
	fi
	@if test -d doc_ref_user/html ; then				\
	  for filen in `find doc_ref_user/html -maxdepth 1 -type f` ; do\
	    echo "Installing $${filen}" ; 				\
	    cp -rp $${filen} $(DESTDIR)$(docdir)/gui/MED ;		\
	  done ;							\
	fi ;

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)/gui/MED

clean-local:
	rm -rf doc_ref_user log_user

EXTRA_DIST += figures \
	      main.dox \
	      Geometric2D.dox \
	      biblio.dox \
	      barycoords.dox \
	      dualmesh.dox \
	      extractor.dox \
	      field.dox \
	      grid.dox \
	      interpkernel.dox \
	      medcoupling.dox \
	      medfilebrowser.dox \
	      medloader.dox \
	      medmem.dox \
	      medsplitter.dox \
	      mesh.dox \
	      meshing.dox \
	      polygon.dox \
	      remapping.dox \
	      support.dox \
	      tools.dox \
	      static/footer.html \
	      static/doxygen.css \
	      images