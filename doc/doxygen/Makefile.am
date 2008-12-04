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
#  MED MEDMEM : MED files in memory
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

DOX_INPUT_FILE = Doxyfile_med_user

userdocdir = $(prefix)/doc/user

install-data-local : html-local
	$(INSTALL) -d $(DESTDIR)$(userdocdir)
	cp -rp doc_ref_user/html $(DESTDIR)$(userdocdir)

uninstall-local:
	rm -rf $(DESTDIR)$(userdocdir)

clean-local:
	rm -rf doc_ref_user log_user

EXTRA_DIST += figures \
	      main.dox \
	      MED_class.dox \
	      biblio.dox \
	      field.dox \
	      grid.dox \
	      medsplitter.dox \
	      mesh.dox \
	      meshing.dox \
	      polygon.dox \
	      remapping.dox \
	      support.dox \
	      tools.dox
