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

FIND_PACKAGE(PythonLibs)

SET_SOURCE_FILES_PROPERTIES(MEDCouplingClient.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDCouplingClient.i PROPERTIES SWIG_FLAGS "-shadow")

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_PATH}
  ${OMNIORB_INCLUDES_DIR}
  ${KERNEL_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}/../../../idl
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCouplingCorba/Client
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  )

SWIG_ADD_MODULE(MEDCouplingClient python MEDCouplingClient.i)
SWIG_LINK_LIBRARIES(MEDCouplingClient ${PYTHON_LIBRARIES} medcouplingclient)
SET_TARGET_PROPERTIES(_MEDCouplingClient PROPERTIES COMPILE_FLAGS "${OMNIORB_FLAGS}")

INSTALL(TARGETS _MEDCouplingClient DESTINATION ${MED_salomepythondir})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/MEDCouplingClient.py DESTINATION ${MED_salomescript_PYTHON})
