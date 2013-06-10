# Copyright (C) 2012-2013  CEA/DEN, EDF R&D
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
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${XDR_INCLUDE_DIRS}
  ${CPPUNIT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  )

SET(MEDMEMCppTest_SOURCES
  MEDMEMTest.cxx
  MEDMEMTest_Array.cxx
  MEDMEMTest_AsciiFieldDriver.cxx
  MEDMEMTest_Connectivity.cxx
  MEDMEMTest_Coordinate.cxx
  MEDMEMTest_DesactivateFacesComputation.cxx
  MEDMEMTest_GetVolumeAbs.cxx
  MEDMEMTest_DriverFactory.cxx
  MEDMEMTest_Exception.cxx
  MEDMEMTest_Extractor.cxx
  MEDMEMTest_Family.cxx
  MEDMEMTest_Field.cxx
  MEDMEMTest_Formulae.cxx
  MEDMEMTest_GaussLocalization.cxx
  MEDMEMTest_GibiMeshDriver.cxx
  MEDMEMTest_Grid.cxx
  MEDMEMTest_Group.cxx
  MEDMEMTest_MedFileBrowser.cxx
  MEDMEMTest_MedFieldDriver.cxx
  MEDMEMTest_MedMeshDriver.cxx
  MEDMEMTest_MedVersion.cxx
  MEDMEMTest_MeshAndMeshing.cxx
  MEDMEMTest_ModulusArray.cxx
  MEDMEMTest_PartialDescendingConnectivity.cxx
  MEDMEMTest_PointerOf.cxx
  MEDMEMTest_PorflowMeshDriver.cxx
  MEDMEMTest_Remapper.cxx
  MEDMEMTest_SkyLineArray.cxx
  MEDMEMTest_Support.cxx
  MEDMEMTest_TopLevel.cxx
  MEDMEMTest_Unit.cxx
  MEDMEMTest_VtkFieldDriver.cxx
  MEDMEMTest_VtkMedDriver.cxx
  MEDMEMTest_VtkMeshDriver.cxx
  MEDMEMTest_nArray.cxx
  MEDMEMTest_Utils.hxx
  MEDMEMTest_Utils.cxx
  MEDMEMTest_MeshFuse.cxx
  )

SET(TestMEDMEM_SOURCES
  TestMEDMEM.cxx
  )

SET(MEDMEMTest_Array_fault_SOURCES
  MEDMEMTest_Array_fault.cxx
  )

SET(MEDMEMTest_Connectivity_fault_SOURCES
  MEDMEMTest_Connectivity_fault.cxx
  )

SET(MEDMEMTest_Coordinate_fault_SOURCES
  MEDMEMTest_Coordinate_fault.cxx
  )

SET(MEDMEMTest_Family_fault_SOURCES
  MEDMEMTest_Family_fault.cxx
  )

SET(MEDMEMTest_Field_fault_SOURCES
  MEDMEMTest_Field_fault.cxx
  )

SET(MEDMEMTest_Grid_fault_SOURCES
  MEDMEMTest_Grid_fault.cxx
  )

SET(MEDMEMTest_MeshAndMeshing_fault_SOURCES
  MEDMEMTest_MeshAndMeshing_fault.cxx
  )

SET(MEDMEMTest_Support_fault_SOURCES
  MEDMEMTest_Support_fault.cxx
  )

SET(MEDMEMTest_VtkFieldDriver_fault_SOURCES
  MEDMEMTest_VtkFieldDriver_fault.cxx
  )

SET(MEDMEMTest_VtkMeshDriver_fault_SOURCES
  MEDMEMTest_VtkMeshDriver_fault.cxx
  )

ADD_LIBRARY(MEDMEMCppTest SHARED ${MEDMEMCppTest_SOURCES})
SET_TARGET_PROPERTIES(MEDMEMCppTest PROPERTIES COMPILE_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${XDR_DEFINITIONS} ${CPPUNIT_DEFINITIONS} ${PLATFORM_DEFINITIONS}")
TARGET_LINK_LIBRARIES(MEDMEMCppTest medmem ${CPPUNIT_LIBS})

SET(ALLMEDMEMCPP_TESTS TestMEDMEM MEDMEMTest_Array_fault MEDMEMTest_Connectivity_fault MEDMEMTest_Coordinate_fault MEDMEMTest_Family_fault MEDMEMTest_Field_fault MEDMEMTest_Grid_fault MEDMEMTest_MeshAndMeshing_fault MEDMEMTest_Support_fault MEDMEMTest_VtkFieldDriver_fault MEDMEMTest_VtkMeshDriver_fault)
FOREACH(MEDMEMCPP_TEST ${ALLMEDMEMCPP_TESTS})
  ADD_EXECUTABLE(${MEDMEMCPP_TEST} ${${MEDMEMCPP_TEST}_SOURCES})
  SET_TARGET_PROPERTIES(${MEDMEMCPP_TEST} PROPERTIES COMPILE_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${XDR_DEFINITIONS} ${PLATFORM_DEFINITIONS}")
  TARGET_LINK_LIBRARIES(${MEDMEMCPP_TEST} MEDMEMCppTest ${PLATFORM_LIBS})
ENDFOREACH(MEDMEMCPP_TEST ${ALLMEDMEMCPP_TESTS})

INSTALL(TARGETS MEDMEMCppTest DESTINATION ${MED_salomelib_LIBS})
INSTALL(TARGETS ${ALLMEDMEMCPP_TESTS} DESTINATION ${MED_salomebin_BINS})
INSTALL(FILES testMEDMEM.py DESTINATION ${MED_salomebin_BINS})
INSTALL(FILES MEDMEMTest.hxx MEDMEMCppTestExport.hxx DESTINATION ${MED_salomeinclude_HEADERS})
