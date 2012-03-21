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

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES SWIG_FLAGS "-shadow")

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_PATH}
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
  ${MPI_INCLUDES_DIR}
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
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/ParaMEDMEM.py DESTINATION ${MED_salomescript_PYTHON})

ADD_TEST(MEDCouplingBasicsTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingBasicsTest.py)
ADD_TEST(MEDCouplingRemapperTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingRemapperTest.py)
