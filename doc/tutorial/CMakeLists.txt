# Copyright (C) 2012-2016  CEA/DEN, EDF R&D
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

SALOME_CONFIGURE_FILE(conf.py.in conf.py)

SET(_cmd_options -c ${CMAKE_CURRENT_BINARY_DIR} -b html -d doctrees -D latex_paper_size=a4 ${CMAKE_CURRENT_SOURCE_DIR} html)
SALOME_GENERATE_ENVIRONMENT_SCRIPT(_cmd env_script "${SPHINX_EXECUTABLE}" "${_cmd_options}")      

ADD_CUSTOM_TARGET(tutorial_doc COMMAND ${_cmd})
ADD_DEPENDENCIES(usr_docs tutorial_doc)

#INSTALL(CODE "EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" --build ${PROJECT_BINARY_DIR} --target tutorial_doc)")
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ DESTINATION ${MEDCOUPLING_INSTALL_DOC}/tutorial)
INSTALL(DIRECTORY images/ DESTINATION ${MEDCOUPLING_INSTALL_DOC}/tutorial/images)
INSTALL(DIRECTORY data/ DESTINATION ${MEDCOUPLING_INSTALL_DOC}/tutorial/data)

SET(make_clean_files html doctrees)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${make_clean_files}")
