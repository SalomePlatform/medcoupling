# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (CEA/DEN)

ADD_DEFINITIONS(${PYTHON_DEFINITIONS} ${OMNIORB_DEFINITIONS})

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(MEDCouplingClient.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDCouplingClient.i PROPERTIES SWIG_DEFINITIONS "-shadow")
SET(SWIG_MODULE_MEDCouplingClient_EXTRA_FLAGS ${SWIG_EXTRA_FLAGS_FOR_NUMPYANDSCIPY})

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${OMNIORB_INCLUDE_DIR} ${OMNIORB_INCLUDE_DIRS}
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
IF(WIN32)
  SET_TARGET_PROPERTIES(_MEDCouplingClient PROPERTIES DEBUG_OUTPUT_NAME _MEDCouplingClient_d)
ENDIF(WIN32)

INSTALL(TARGETS _MEDCouplingClient DESTINATION ${SALOME_INSTALL_PYTHON})
INSTALL(FILES MEDCouplingClient.i DESTINATION ${SALOME_INSTALL_HEADERS})
SET(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDCouplingClient.py)
INSTALL_AND_COMPILE_PYTHON_FILE("${PYFILES_TO_INSTALL}" ${SALOME_INSTALL_PYTHON})
