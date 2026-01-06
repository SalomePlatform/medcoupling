# Copyright (C) 2012-2026  CEA, EDF
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

salome_configure_file(conf.py.in conf.py)

set(_cmd_options
    -c
    ${CMAKE_CURRENT_BINARY_DIR}
    -b
    html
    -d
    doctrees
    -D
    latex_paper_size=a4
    ${CMAKE_CURRENT_SOURCE_DIR}
    html)
salome_generate_environment_script(_cmd env_script "${SPHINX_EXECUTABLE}"
                                   "${_cmd_options}")

add_custom_target(tutorial_doc COMMAND ${_cmd})
add_dependencies(doxygen_documentation tutorial_doc)

# INSTALL(CODE "EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" --build
# ${PROJECT_BINARY_DIR} --target tutorial_doc)")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
        DESTINATION ${MEDCOUPLING_INSTALL_DOC}/tutorial)
install(DIRECTORY images/
        DESTINATION ${MEDCOUPLING_INSTALL_DOC}/tutorial/images)
install(DIRECTORY data/ DESTINATION ${MEDCOUPLING_INSTALL_DOC}/tutorial/data)

set(make_clean_files html doctrees)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                    "${make_clean_files}")
