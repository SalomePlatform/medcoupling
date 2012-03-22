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
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
  ${KERNEL_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}/../../idl
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(MEDMEMImpl_SOURCES
  MEDMEM_Med_i.cxx
  MEDMEM_Family_i.cxx
  MEDMEM_Field_i.cxx
  MEDMEM_Group_i.cxx
  MEDMEM_Mesh_i.cxx
  MEDMEM_Support_i.cxx
  MEDMEM_convert.cxx
  MEDMEM_GMesh_i.cxx
  )

ADD_LIBRARY(MEDMEMImpl SHARED ${MEDMEMImpl_SOURCES})
SET_TARGET_PROPERTIES(MEDMEMImpl PROPERTIES COMPILE_FLAGS "${OMNIORB_FLAGS} ${HDF5_FLAGS} ${MED3_FLAGS} ${XDR_FLAGS}")
TARGET_LINK_LIBRARIES(MEDMEMImpl medmem SalomeIDLMED)
INSTALL(TARGETS MEDMEMImpl DESTINATION ${MED_salomelib_LIBS})

FILE(GLOB MEDMEMImpl_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${MEDMEMImpl_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
