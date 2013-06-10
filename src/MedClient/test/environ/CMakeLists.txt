# Copyright (C) 2012-2013  CEA/DEN, EDF R&D
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

ADD_SUBDIRECTORY(csh)

SET(MedClientTestsConfig 
#  runTestMedCorba 
  runContainer
)

FOREACH(MedClientTestConfig ${MedClientTestsConfig})
  SET(input ${CMAKE_CURRENT_SOURCE_DIR}/${MedClientTestConfig}.in)
  SET(output ${CMAKE_CURRENT_BINARY_DIR}/${MedClientTestConfig})
  MESSAGE(STATUS "Creation of ${output}")
  CONFIGURE_FILE(${input} ${output})
  INSTALL(FILES ${output} PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ DESTINATION Tests/environ/csh)
ENDFOREACH(MedClientTestConfig ${MedClientTestsConfig})
