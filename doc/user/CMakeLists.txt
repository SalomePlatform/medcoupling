# Copyright (C) 2012-2019  CEA/DEN, EDF R&D
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

INCLUDE(UseSphinx)

SALOME_CONFIGURE_FILE(conf.py.in conf.py)

IF(MEDCOUPLING_BUILD_FRENCH_DOC)
  SET(BUILD_LANGUAGE_SWITCHER "true")
  SALOME_CONFIGURE_FILE(static/switchers.js.in static/switchers.js)
  ADD_MULTI_LANG_DOCUMENTATION(TARGET_NAME usr_doc MODULE html LANGUAGES fr)
ELSE()
  SET(BUILD_LANGUAGE_SWITCHER "false")
  SALOME_CONFIGURE_FILE(static/switchers.js.in static/switchers.js)
  ADD_MULTI_LANG_DOCUMENTATION(TARGET_NAME usr_doc MODULE html)
ENDIF()

#ADD_DEPENDENCIES(usr_doc dev_docs)

INSTALL(CODE "EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" --build ${PROJECT_BINARY_DIR} --target usr_doc)")
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${MEDCOUPLING_INSTALL_DOC}/user)
IF(MEDCOUPLING_BUILD_FRENCH_DOC)
  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html_fr DESTINATION ${MEDCOUPLING_INSTALL_DOC}/user)
ENDIF()

SET(make_clean_files html html_fr doctrees locale potfiles static)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${make_clean_files}")
