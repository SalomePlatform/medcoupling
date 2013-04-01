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
# Author : Anthony Geay (CEA/DEN)

INCLUDE(${SALOME_MACROS_DIR}/InstallAndCompilePythonFile.cmake)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(MEDLoader.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDLoader.i PROPERTIES SWIG_DEFINITIONS "-shadow")

SET (MEDLoader_SWIG_DPYS_FILES
    MEDLoaderCommon.i
    MEDLoaderTypemaps.i)

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/GaussPoints
)

SET (SWIG_MODULE_MEDLoader_EXTRA_DEPS ${MEDLoader_SWIG_DPYS_FILES}
    ${medloader_HEADERS_HXX}
    ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX}
    ${interpkernel_HEADERS_HXX} ${interpkernel_HEADERS_TXX})
    
SWIG_ADD_MODULE(MEDLoader python MEDLoader.i)
SWIG_LINK_LIBRARIES(MEDLoader ${PYTHON_LIBS} ${PLATFORM_LIBS} medloader medcoupling)
SET_TARGET_PROPERTIES(_MEDLoader PROPERTIES COMPILE_FLAGS "${PLATFORM_DEFINITIONS}")
IF(${MACHINE} STREQUAL WINDOWS)
  SET_TARGET_PROPERTIES(_MEDLoader PROPERTIES DEBUG_OUTPUT_NAME _MEDLoader_d)
ENDIF(${MACHINE} STREQUAL WINDOWS)
INSTALL(TARGETS _MEDLoader DESTINATION ${MED_salomelib_LIBS})

INSTALL(FILES MEDLoader.i MEDLoaderTypemaps.i MEDLoaderCommon.i DESTINATION ${MED_salomeinclude_HEADERS})

SET(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDLoader.py ${CMAKE_CURRENT_SOURCE_DIR}/libMEDLoader_Swig.py)
INSTALL_AND_COMPILE_PYTHON_FILE("${PYFILES_TO_INSTALL}" ${MED_salomescript_PYTHON})

INSTALL(FILES MEDLoaderDataForTest.py MEDLoaderTest.py MEDLoaderTest2.py MEDLoaderTest3.py DESTINATION ${MED_salomescript_PYTHON})

ADD_TEST(MEDLoaderTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDLoaderTest.py)
ADD_TEST(MEDLoaderTest2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDLoaderTest2.py)
ADD_TEST(MEDLoaderTest3 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDLoaderTest3.py)
ADD_TEST(MEDLoaderExamplesTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDLoaderExamplesTest.py)
ADD_TEST(SauvLoaderTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SauvLoaderTest.py)
