# Copyright (C) 2012-2014  CEA/DEN, EDF R&D
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

ADD_SUBDIRECTORY(Client)

ADD_DEFINITIONS(${PYTHON_DEFINITIONS} ${OMNIORB_DEFINITIONS})

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(MEDCouplingCorba.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDCouplingCorba.i PROPERTIES SWIG_DEFINITIONS "-shadow")
SET(SWIG_MODULE_MEDCouplingCorba_EXTRA_FLAGS ${SWIG_EXTRA_FLAGS_FOR_NUMPYANDSCIPY})

SET(MEDCouplingCorba_SWIG_DPYS_FILES
    MEDCouplingCorbaTypemaps.i)

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${OMNIORB_INCLUDE_DIR} ${OMNIORB_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}/../../idl
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCouplingCorba
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  )

# _ABR_ Ensure dependency mechanism on all SWIG files and headers
SET (SWIG_MODULE_MEDCouplingCorba_EXTRA_DEPS ${MEDCouplingCorba_SWIG_DPYS_FILES} 
    ${medcouplingcorba_HEADERS_HXX} 
    ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX}
    ${interpkernel_HEADERS_HXX} ${interpkernel_HEADERS_TXX})
    
SWIG_ADD_MODULE(MEDCouplingCorba python MEDCouplingCorba.i)
SWIG_LINK_LIBRARIES(MEDCouplingCorba ${PYTHON_LIBRARIES} medcouplingcorba ${PLATFORM_LIBS})
IF(WIN32)
  SET_TARGET_PROPERTIES(_MEDCouplingCorba PROPERTIES DEBUG_OUTPUT_NAME _MEDCouplingCorba_d)
ENDIF(WIN32)

INSTALL(TARGETS _MEDCouplingCorba DESTINATION ${SALOME_INSTALL_PYTHON})
INSTALL(FILES MEDCouplingCorba.i MEDCouplingCorbaTypemaps.i DESTINATION ${SALOME_INSTALL_HEADERS})
SET(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDCouplingCorba.py)
INSTALL_AND_COMPILE_PYTHON_FILE("${PYFILES_TO_INSTALL}" ${SALOME_INSTALL_PYTHON})
