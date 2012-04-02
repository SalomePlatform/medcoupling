# Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/InstallAndCompilePythonFile.cmake)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(MEDLoader.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDLoader.i PROPERTIES SWIG_FLAGS "-shadow")

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDES_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  )

SWIG_ADD_MODULE(MEDLoader python MEDLoader.i)

SWIG_LINK_LIBRARIES(MEDLoader ${PYTHON_LIBRARIES} medloader)

INCLUDE_DIRECTORIES(
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/GaussPoints
  )

INSTALL(TARGETS _MEDLoader DESTINATION ${MED_salomepythondir})
SET(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDLoader.py ${CMAKE_CURRENT_SOURCE_DIR}/libMEDLoader_Swig.py)
INSTALL_AND_COMPILE_PYTHON_FILE("${PYFILES_TO_INSTALL}" ${MED_salomescript_PYTHON})

ADD_TEST(MEDLoaderTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDLoaderTest.py)
ADD_TEST(MEDLoaderTest2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDLoaderTest2.py)
ADD_TEST(MEDLoaderTest3 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDLoaderTest3.py)
ADD_TEST(SauvLoaderTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SauvLoaderTest.py)
