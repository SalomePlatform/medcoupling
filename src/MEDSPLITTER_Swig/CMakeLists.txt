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

SET_SOURCE_FILES_PROPERTIES(libMEDSPLITTER_Swig.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(libMEDSPLITTER_Swig.i PROPERTIES SWIG_FLAGS "-shadow")

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDES_DIR}
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
  ${BOOST_INCLUDES_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDSPLITTER
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM_SWIG
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SWIG_ADD_MODULE(libMEDSPLITTER_Swig python libMEDSPLITTER_Swig.i)
SWIG_LINK_LIBRARIES(libMEDSPLITTER_Swig ${PYTHON_LIBRARIES} medsplittercpp)
SET_TARGET_PROPERTIES(_libMEDSPLITTER_Swig PROPERTIES COMPILE_FLAGS "${BOOST_FLAGS}")

INSTALL(TARGETS _libMEDSPLITTER_Swig DESTINATION ${MED_salomepythondir})
INSTALL(FILES libMEDSPLITTER_Swig.i DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/libMEDSPLITTER_Swig.py MEDSPLITTER.py test_MEDSPLITTER.py DESTINATION ${MED_salomescript_PYTHON})
INSTALL_AND_COMPILE_PYTHON_FILE(${CMAKE_CURRENT_SOURCE_DIR}/MEDSPLITTER.py ${MED_salomescript_PYTHON})

