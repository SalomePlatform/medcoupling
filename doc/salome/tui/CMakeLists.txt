# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
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

SALOME_CONFIGURE_FILE(doxyfile.in doxyfile)
SALOME_CONFIGURE_FILE(static/header.html.in static/header.html)

ADD_CUSTOM_TARGET(dev_docs ${DOXYGEN_EXECUTABLE})

INSTALL(CODE "FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/MED)")
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/MED 
        DESTINATION ${SALOME_INSTALL_DOC}/tui
        PATTERN "*.md5" EXCLUDE
        PATTERN "*.map" EXCLUDE)
INSTALL(FILES images/head.png DESTINATION ${SALOME_INSTALL_DOC}/tui/MED) 

SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES MED)
