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

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDES_DIR}
  ${MPI_INCLUDES_DIR}
  ${KERNEL_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}/../../idl
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDCouplingCorba
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCouplingCorba
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  )

SET(paramedmemcompo_SOURCES
  MPIMEDCouplingFieldDoubleServant.cxx
  ParaMEDMEMComponent_i.cxx
  )

ADD_LIBRARY(paramedmemcompo SHARED ${paramedmemcompo_SOURCES})
SET_TARGET_PROPERTIES(paramedmemcompo PROPERTIES COMPILE_FLAGS "${MPI_FLAGS} ${OMNIORB_FLAGS}")
TARGET_LINK_LIBRARIES(paramedmemcompo medcouplingcorba paramedmem)
INSTALL(TARGETS paramedmemcompo DESTINATION ${MED_salomelib_LIBS})

FILE(GLOB paramedmemcompo_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${paramedmemcompo_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
