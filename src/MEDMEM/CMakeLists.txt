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
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
  ${XDR_INCLUDES_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(medmem_SOURCES
  MEDMEM_Init.cxx
  MEDMEM_RCBase.cxx
  MEDMEM_CellModel.cxx
  MEDMEM_Connectivity.cxx
  MEDMEM_Coordinate.cxx
  MEDMEM_DriversDef.cxx
  MEDMEM_DriverFactory.cxx
  MEDMEM_Family.cxx
  MEDMEM_Field.cxx
  MEDMEM_GenDriver.cxx
  MEDMEM_Group.cxx
  MEDMEM_Exception.cxx
  MEDMEM_MedMeshDriver.cxx
  MEDMEM_GibiMeshDriver.cxx
  MEDMEM_PorflowMeshDriver.cxx
  MEDMEM_SkyLineArray.cxx
  MEDMEM_Mesh.cxx
  MEDMEM_Meshing.cxx
  MEDMEM_Support.cxx
  MEDMEM_Unit.cxx
  MEDMEM_VtkMedDriver.cxx
  MEDMEM_VtkMeshDriver.cxx
  MEDMEM_Grid.cxx
  MEDMEM_MedVersion.cxx
  MEDMEM_TopLevel.cxx
  MEDMEM_ConnectZone.cxx
  MEDMEM_EnsightMedDriver.cxx
  MEDMEM_EnsightFieldDriver.cxx
  MEDMEM_EnsightUtils.cxx
  MEDMEM_EnsightMeshDriver.cxx
  MEDMEM_DriverTools.cxx
  MEDMEM_GaussLocalization.cxx
  MEDMEM_Extractor.cxx
  PointLocator.cxx
  MEDMEM_Remapper.cxx
  MEDMEM_GMesh.cxx
  MEDMEM_MeshFuse.cxx
  MEDMEM_MedFileBrowser.cxx
  )

ADD_LIBRARY(medmem SHARED ${medmem_SOURCES})
SET_TARGET_PROPERTIES(medmem PROPERTIES COMPILE_FLAGS "${HDF5_FLAGS} ${MED3_FLAGS} ${XDR_FLAGS}")
TARGET_LINK_LIBRARIES(medmem interpkernel ${MED3_LIBS_C_ONLY} ${HDF5_LIBS} ${XDR_LIBS})
INSTALL(TARGETS medmem DESTINATION ${MED_salomelib_LIBS})

FILE(GLOB medmem_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${medmem_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
FILE(GLOB medmem_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${medmem_HEADERS_TXX} DESTINATION ${MED_salomeinclude_HEADERS})
