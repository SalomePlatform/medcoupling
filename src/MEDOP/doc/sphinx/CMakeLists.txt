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

IF(SPHINX_STATUS)
  SET(srcdir ${CMAKE_CURRENT_SOURCE_DIR})
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)
  #FILE(GLOB MEDOP_DOC_RST "${CMAKE_CURRENT_SOURCE_DIR}/*.rst")
  
  #FOREACH(SRCFNAME ${MEDOP_DOC_RST})
  #  GET_FILENAME_COMPONENT(FNAME ${SRCFNAME} NAME)
  #  CONFIGURE_FILE(${SRCFNAME} ${CMAKE_CURRENT_BINARY_DIR}/${FNAME} COPY_ONLY)
  #ENDFOREACH(FNAME ${MEDOP_DOC_RST})
  
  #FILE(MAKE_DIRECTORY _static)
  #CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/_static/xmed.css ${CMAKE_CURRENT_BINARY_DIR}/_static/xmed.css COPY_ONLY)

  #FILE(GLOB MEDOP_DOC_PNG_RST "${CMAKE_CURRENT_SOURCE_DIR}/images/*.png")
  #FOREACH(SRCFNAME ${MEDOP_DOC_PNG_RST})
  #  GET_FILENAME_COMPONENT(FNAME ${SRCFNAME} NAME)
  #  CONFIGURE_FILE(${SRCFNAME} ${CMAKE_CURRENT_BINARY_DIR}/images/${FNAME} COPY_ONLY)
  #ENDFOREACH(FNAME ${MEDOP_DOC_RST})

  INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${SPHINX_EXECUTABLE_TO_FIND} -c ${CMAKE_CURRENT_BINARY_DIR} -b html -d doctrees -D latex_paper_size=a4 ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/html)")

  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ DESTINATION share/doc/salome/dev/XMED)

ENDIF(SPHINX_STATUS)
