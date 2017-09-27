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

ADD_DEFINITIONS(${PYTHON_DEFINITIONS} ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS})

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES CPLUSPLUS ON)
IF ("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES SWIG_FLAGS "-py3")
ELSE()
  SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES SWIG_DEFINITIONS "-shadow")
ENDIF()
SET(SWIG_MODULE_ParaMEDMEM_EXTRA_FLAGS "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")

SET (ParaMEDMEM_SWIG_DPYS_FILES
    ParaMEDMEM.typemap)

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${NUMPY_INCLUDE_DIR}
  ${MPI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  )

SET (SWIG_MODULE_ParaMEDMEM_EXTRA_DEPS ${ParaMEDMEM_SWIG_DPYS_FILES}
    ${paramedmem_HEADERS_HXX}
    ${medloader_HEADERS_HXX}
    ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX}
    ${interpkernel_HEADERS_HXX} ${interpkernel_HEADERS_TXX})

SWIG_ADD_MODULE(ParaMEDMEM python ParaMEDMEM.i)
SWIG_LINK_LIBRARIES(ParaMEDMEM ${PYTHON_LIBRARIES} paramedmem)
SWIG_CHECK_GENERATION(ParaMEDMEM)

INSTALL(TARGETS _ParaMEDMEM DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})
INSTALL(FILES ParaMEDMEM.i DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/ParaMEDMEM.py test_InterpKernelDEC.py test_NonCoincidentDEC.py test_StructuredCoincidentDEC.py DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})

SALOME_INSTALL_SCRIPTS(${CMAKE_CURRENT_BINARY_DIR}/ParaMEDMEM.py ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})
INSTALL(FILES test_InterpKernelDEC.py test_NonCoincidentDEC.py test_StructuredCoincidentDEC.py DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})
