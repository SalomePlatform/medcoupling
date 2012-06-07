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

INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/InstallAndCompilePythonFile.cmake)
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES SWIG_FLAGS "-shadow")

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MPI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader/Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  )

SWIG_ADD_MODULE(ParaMEDMEM python ParaMEDMEM.i)
SWIG_LINK_LIBRARIES(ParaMEDMEM ${PYTHON_LIBRARIES} paramedmem medloader)

SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES SWIG_FLAGS "-shadow")

INSTALL(TARGETS _ParaMEDMEM DESTINATION ${MED_salomepythondir})
INSTALL(FILES ParaMEDMEM.i DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL_AND_COMPILE_PYTHON_FILE(${CMAKE_CURRENT_BINARY_DIR}/ParaMEDMEM.py ${MED_salomescript_PYTHON})
