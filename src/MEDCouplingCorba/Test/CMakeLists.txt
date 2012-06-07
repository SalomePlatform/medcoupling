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

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIRS}
  ${KERNEL_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}/../../../idl
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../Client
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  )

SET(TestMEDCouplingCorbaServ_SOURCES
  TestMEDCouplingCorbaServ.cxx
  MEDCouplingMeshFieldFactoryComponent.cxx
  MEDCouplingMeshFieldFactoryComponentServ.cxx
  )

ADD_EXECUTABLE(TestMEDCouplingCorbaServ ${TestMEDCouplingCorbaServ_SOURCES})
SET_TARGET_PROPERTIES(TestMEDCouplingCorbaServ PROPERTIES COMPILE_FLAGS "${OMNIORB_FLAGS}")
TARGET_LINK_LIBRARIES(TestMEDCouplingCorbaServ medcouplingcorba SalomeIDLMEDTests ${OMNIORB_LIBS})

SET(TestMEDCouplingCorbaClt_SOURCES
  TestMEDCouplingCorbaClt.cxx
  MEDCouplingMeshFieldFactoryComponent.cxx
  MEDCouplingMeshFieldFactoryComponentClt.cxx
  )

ADD_EXECUTABLE(TestMEDCouplingCorbaClt ${TestMEDCouplingCorbaClt_SOURCES})
SET_TARGET_PROPERTIES(TestMEDCouplingCorbaClt PROPERTIES COMPILE_FLAGS "${OMNIORB_FLAGS} ${CPPUNIT_FLAGS}")
TARGET_LINK_LIBRARIES(TestMEDCouplingCorbaClt medcouplingclient SalomeIDLMEDTests ${OMNIORB_LIBS} ${CPPUNIT_LIBS})

INSTALL(TARGETS TestMEDCouplingCorbaServ TestMEDCouplingCorbaClt DESTINATION ${MED_salomebin_BINS})
