#  Copyright (C) 2003  CEA/DEN, EDF R&D
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

include $(top_srcdir)/adm_local/unix/make_common_starter.am

EXTRA_DIST += images static

dev_docs: doxyfile
	echo "Running doxygen in directory: "`pwd`; \
	$(DOXYGEN) $<;

clean-local:
	-rm -fr MED doxygen.bak

install-data-local:
	if test -d MED; then \
	  $(INSTALL) -d $(DESTDIR)$(docdir)/tui ; \
	  cp -rp MED $(DESTDIR)$(docdir)/tui ; \
	fi;

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)/tui/MED
