# Copyright (C) 2012-2025  CEA, EDF
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this library; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email :
# webmaster.salome@opencascade.com
#

include(UseSphinx)

salome_configure_file(conf.py.in conf.py)

if(MEDCOUPLING_BUILD_FRENCH_DOC)
  set(BUILD_LANGUAGE_SWITCHER "true")
  salome_configure_file(static/switchers.js.in static/switchers.js)
  add_multi_lang_documentation(TARGET_NAME usr_doc MODULE html LANGUAGES fr)
else()
  set(BUILD_LANGUAGE_SWITCHER "false")
  salome_configure_file(static/switchers.js.in static/switchers.js)
  add_multi_lang_documentation(TARGET_NAME usr_doc MODULE html)
endif()

install(
  CODE "EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" --build ${PROJECT_BINARY_DIR} --target usr_doc)"
)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
        DESTINATION ${MEDCOUPLING_INSTALL_DOC}/user)
if(MEDCOUPLING_BUILD_FRENCH_DOC)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html_fr
          DESTINATION ${MEDCOUPLING_INSTALL_DOC}/user)
endif()

set(make_clean_files html html_fr doctrees locale potfiles static)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                    "${make_clean_files}")
