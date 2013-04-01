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

SET(top_builddir ${CMAKE_BINARY_DIR})
SET(top_srcdir ${CMAKE_SOURCE_DIR})
SET(srcdir ${CMAKE_CURRENT_SOURCE_DIR})
SET(builddir ${CMAKE_CURRENT_BINARY_DIR})
SET(datadir ${CMAKE_INSTALL_PREFIX}/share)
SET(docdir ${datadir}/doc/salome)

SET(indoxfiles Doxyfile_med_user.in static/header.html.in)
FOREACH(indoxfile ${indoxfiles})
  STRING(REGEX REPLACE ".in" "" baseindoxfile ${indoxfile})
  SET(input ${CMAKE_CURRENT_SOURCE_DIR}/${indoxfile})
  SET(output ${CMAKE_CURRENT_BINARY_DIR}/${baseindoxfile})
  CONFIGURE_FILE(${input} ${output})
  MESSAGE(STATUS "Creation of ${output}")
ENDFOREACH(indoxfile ${indoxfiles})
FILE(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/medcouplingexamples.doxy" input)
FILE(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/BuildPyExamplesFromCPP.py" pythondocexamplesgenerator)
FILE(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}" output)

ADD_CUSTOM_TARGET(usr_docs
  COMMAND ${PYTHON_EXECUTABLE} ${pythondocexamplesgenerator} ${input} ${output}
  COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile_med_user
  COMMAND ${PYTHON_EXECUTABLE} -c "import shutil, sys; shutil.rmtree(r'''${CMAKE_INSTALL_PREFIX}/share/doc/salome/gui/MED''', True); shutil.copytree(r'''${CMAKE_CURRENT_BINARY_DIR}/doc_ref_user/html''', r'''${CMAKE_INSTALL_PREFIX}/share/doc/salome/gui/MED'''); shutil.copy(r'''${CMAKE_CURRENT_SOURCE_DIR}/images/head.png''', r'''${CMAKE_INSTALL_PREFIX}/share/doc/salome/gui/MED''')"
  VERBATIM 
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}             
  )
