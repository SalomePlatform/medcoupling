# Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config_appli_model.xml.in ${CMAKE_CURRENT_BINARY_DIR}/config_appli_model.xml @ONLY)

SET(MED_RESOURCES_FILES
  ${CMAKE_CURRENT_BINARY_DIR}/config_appli_model.xml
  appli-splashscreen.jpg
  SalomeApp.xml
  )

INSTALL(FILES ${MED_RESOURCES_FILES} DESTINATION ${MED_salomebin_BINS}/med/appligen)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/appligen.sh.in ${CMAKE_CURRENT_BINARY_DIR}/appligen.sh @ONLY)

SET(MED_RESOURCES_FILES2
  ${CMAKE_CURRENT_BINARY_DIR}/appligen.sh
  )

INSTALL(FILES ${MED_RESOURCES_FILES2} DESTINATION ${MED_salomebin_BINS}/med/appligen PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
