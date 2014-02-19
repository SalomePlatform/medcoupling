# Copyright (C) 2012-2014  CEA/DEN, EDF R&D, OPEN CASCADE
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

SET(input ${CMAKE_CURRENT_SOURCE_DIR}/SALOMEconfig.h.in)
SET(output ${CMAKE_CURRENT_BINARY_DIR}/SALOMEconfig.h)
MESSAGE(STATUS "Creation of ${output}")
CONFIGURE_FILE(${input} ${output})
INSTALL(FILES ${output} DESTINATION ${SALOME_INSTALL_HEADERS})

FILE(GLOB CMAKE_CONFIG_FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake_files/*.cmake")
INSTALL(FILES ${CMAKE_CONFIG_FILES} DESTINATION ${SALOME_INSTALL_CMAKE_LOCAL})
