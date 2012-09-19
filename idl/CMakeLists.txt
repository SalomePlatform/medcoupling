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
# Author : Anthony Geay (CEA/DEN)

INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/UseOMNIORB.cmake)

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIRS}
  ${KERNEL_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}
  )

SET(SalomeIDLMED_IDLSOURCES
  MED.idl
  MED_Gen.idl
  Compo1Py.idl
  MEDCouplingCorbaServant.idl
  ParaMEDCouplingCorbaServant.idl
  MEDDataManager.idl
  MEDCalculator.idl
  MEDEventListener.idl
  MEDOPFactory.idl
  )

IF(MPI_IS_OK)
  SET(SalomeIDLMED_IDLSOURCES ${SalomeIDLMED_IDLSOURCES} ParaMEDMEMComponent.idl)
ENDIF(MPI_IS_OK)

SET(SalomeIDLMEDTests_IDLSOURCES
  MEDCouplingCorbaServantTest.idl
  )

SET(IDL_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${KERNEL_ROOT_DIR}/idl/salome
  )

OMNIORB_ADD_MODULE(SalomeIDLMED "${SalomeIDLMED_IDLSOURCES}" "${IDL_INCLUDE_DIRS}")
TARGET_LINK_LIBRARIES(SalomeIDLMED ${CORBA_LIBS} ${KERNEL_LDFLAGS} ${SalomeIDLKernel})
SET_TARGET_PROPERTIES(SalomeIDLMED PROPERTIES COMPILE_FLAGS "${OMNIORB_DEFINITIONS}")

OMNIORB_ADD_MODULE(SalomeIDLMEDTests "${SalomeIDLMEDTests_IDLSOURCES}" "${IDL_INCLUDE_DIRS}")
TARGET_LINK_LIBRARIES(SalomeIDLMEDTests ${CORBA_LIBS} ${KERNEL_LDFLAGS} ${SalomeIDLKernel})
SET_TARGET_PROPERTIES(SalomeIDLMEDTests PROPERTIES COMPILE_FLAGS "${OMNIORB_DEFINITIONS}")

INSTALL(TARGETS SalomeIDLMED SalomeIDLMEDTests DESTINATION ${MED_salomelib_LIBS})
