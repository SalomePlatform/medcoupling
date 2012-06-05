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
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
  ${XDR_INCLUDES_DIR}
  ${CPPUNIT_INCLUDES_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(dumpInterpolation_SOURCES dumpInterpolation.cxx)

SET(convertToPoly_SOURCES convertToPoly.cxx)

SET(ensight2med_SOURCES ensight2med.cxx)

SET(med2ensight_SOURCES med2ensight.cxx)

SET(med2vtk_SOURCES med2vtk.cxx)

SET(duplicateMED_SOURCES duplicateMED.cxx)

SET(medmembin_test_SOURCES med_test.cxx)

SET(duplicateMEDMESH_SOURCES duplicateMEDMESH.cxx)

SET(test_write_groups_poly_SOURCES test_write_groups_poly.cxx)

SET(test_MEDMEM_ModulusArray_SOURCES test_MEDMEM_ModulusArray.cxx)

SET(test_MEDMEM_Array_SOURCES test_MEDMEM_Array.cxx)

SET(test_MEDMEM_SkyLineArray_SOURCES test_MEDMEM_SkyLineArray.cxx)

SET(test_MEDMEM_CellModel_SOURCES test_MEDMEM_CellModel.cxx)

SET(testUPointerOf_SOURCES tests/testUPointerOf.cxx)

SET(testUCoordinate_SOURCES tests/testUCoordinate.cxx)

SET(testUUnit_SOURCES tests/testUUnit.cxx)

SET(testUGeoNameMeshEntities_SOURCES tests/testUGeoNameMeshEntities.cxx)

SET(testUMedException_SOURCES tests/testUMedException.cxx)

SET(testUModulusArray_SOURCES tests/testUModulusArray.cxx)

SET(testUSkyLineArray_SOURCES tests/testUSkyLineArray.cxx)

SET(testUArray_SOURCES tests/testUArray.cxx)

SET(testUCellModel_SOURCES tests/testUCellModel.cxx)

SET(readEntete_SOURCES tests/readEntete.cxx)

SET(readCoordinate_SOURCES tests/readCoordinate.cxx)

SET(test_copie_field__SOURCES test_copie_field_.cxx)

SET(test_copie_fieldT_SOURCES test_copie_fieldT.cxx)

SET(test_copie_coordinate_SOURCES test_copie_coordinate.cxx)

SET(test_copie_medarray_SOURCES test_copie_medarray.cxx)

SET(test_copie_connectivity_SOURCES test_copie_connectivity.cxx)

SET(test_copie_support_SOURCES test_copie_support.cxx)

SET(test_copie_family_SOURCES test_copie_family.cxx)

SET(test_copie_group_SOURCES test_copie_group.cxx)

SET(test_copie_mesh_SOURCES test_copie_mesh.cxx)

SET(test_affect_medarray_SOURCES test_affect_medarray.cxx)

SET(test_MEDMEM_Meshing_SOURCES test_MEDMEM_Meshing.cxx)

SET(test_MEDMEM_MeshingPoly_SOURCES test_MEDMEM_MeshingPoly.cxx)

SET(test_operation_fielddouble_SOURCES test_operation_fielddouble.cxx)

SET(test_operation_fieldint_SOURCES test_operation_fieldint.cxx)

SET(test_gibi_driver_SOURCES test_gibi_driver.cxx)

SET(test_porflow_driver_SOURCES test_porflow_driver.cxx)

SET(test_grid_SOURCES test_grid.cxx)

SET(test_MEDMEM_PolyDriverMedMeshRead_SOURCES test_MEDMEM_PolyDriverMedMeshRead.cxx)

SET(test_MEDMEM_PolyDriverMedMeshWrite_SOURCES test_MEDMEM_PolyDriverMedMeshWrite.cxx)

SET(test_MEDMEM_poly3D_SOURCES test_MEDMEM_poly3D.cxx)

SET(test_MEDMEM_nArray_SOURCES test_MEDMEM_nArray.cxx)

SET(test_MEDMEM_Meshing_poly_SOURCES test_MEDMEM_Meshing_poly.cxx)

SET(test_profil_MedFieldDriver_SOURCES test_profil_MedFieldDriver.cxx)

SET(test_profil_gauss_MedFieldDriver_SOURCES test_profil_gauss_MedFieldDriver.cxx)

SET(test_GaussLocalization_SOURCES test_GaussLocalization.cxx)

SET(testAnalFile_SOURCES testAnalFile.cxx)

SET(test_MEDMEM_MeshingFlica_SOURCES test_MEDMEM_MeshingFlica.cxx)

SET(testConvertPolygon_SOURCES testConvertPolygon.cxx)

SET(ALLMEDMEMBIN_TESTS
  dumpInterpolation
  convertToPoly
  ensight2med
  med2ensight
  med2vtk
  duplicateMED
  medmembin_test
  duplicateMEDMESH
  test_write_groups_poly
  test_MEDMEM_ModulusArray
  test_MEDMEM_Array
  test_MEDMEM_SkyLineArray
  test_MEDMEM_CellModel
  testUPointerOf
  testUCoordinate
  testUUnit
  testUGeoNameMeshEntities
  testUMedException
  testUModulusArray
  testUSkyLineArray
  testUArray
  testUCellModel
  readEntete
  readCoordinate
  test_copie_field_
  test_copie_fieldT
  test_copie_coordinate
  test_copie_medarray
  test_copie_connectivity
  test_copie_support
  test_copie_family
  test_copie_group
  test_copie_mesh
  test_affect_medarray
  test_MEDMEM_Meshing
  test_MEDMEM_MeshingPoly
  test_operation_fielddouble
  test_operation_fieldint
  test_gibi_driver
  test_porflow_driver
  test_grid
  test_MEDMEM_PolyDriverMedMeshRead
  test_MEDMEM_PolyDriverMedMeshWrite
  test_MEDMEM_poly3D
  test_MEDMEM_nArray
  test_MEDMEM_Meshing_poly
  test_profil_MedFieldDriver
  test_profil_gauss_MedFieldDriver
  test_GaussLocalization
  testAnalFile
  test_MEDMEM_MeshingFlica
  testConvertPolygon
  )

FOREACH(ALLMEDMEMBIN_TEST ${ALLMEDMEMBIN_TESTS})
  ADD_EXECUTABLE(${ALLMEDMEMBIN_TEST} ${${ALLMEDMEMBIN_TEST}_SOURCES})
  SET_TARGET_PROPERTIES(${ALLMEDMEMBIN_TEST} PROPERTIES COMPILE_FLAGS "${HDF5_FLAGS} ${MED3_FLAGS} ${XDR_FLAGS}")
  TARGET_LINK_LIBRARIES(${ALLMEDMEMBIN_TEST} medmem)
ENDFOREACH(ALLMEDMEMBIN_TEST ${ALLMEDMEMBIN_TESTS})

INSTALL(TARGETS ${ALLMEDMEMBIN_TESTS} DESTINATION ${MED_salomebin_BINS})

SET(MEDMEMBINTEST_PYOPT_TEST -c"import sys;import os;import subprocess;os.putenv(\"srcdir\",sys.argv[1]);sys.exit(subprocess.call([\"sh\",sys.argv[2]]))")
ADD_TEST(check_convertToPoly ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_convertToPoly)
ADD_TEST(check_duplicateMED ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_duplicateMED)
ADD_TEST(check_ensight2med ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_ensight2med)
ADD_TEST(check_med2ensight ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_med2ensight)
ADD_TEST(check_med2vtk ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_med2vtk)
ADD_TEST(check_readCoordinate ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_readCoordinate)
ADD_TEST(check_readEntete ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_readEntete)
ADD_TEST(check_testAnalFile ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_testAnalFile)
ADD_TEST(check_test_MEDMEM_PolyDriverMedMeshRead ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_MEDMEM_PolyDriverMedMeshRead)
ADD_TEST(check_test_MEDMEM_PolyDriverMedMeshWrite ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_MEDMEM_PolyDriverMedMeshWrite)
ADD_TEST(check_test_MEDMEM_poly3D ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_MEDMEM_poly3D)
ADD_TEST(check_test_copie_connectivity ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_copie_connectivity)
ADD_TEST(check_test_copie_family ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_copie_family)
ADD_TEST(check_test_copie_fieldT ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_copie_fieldT)
ADD_TEST(check_test_copie_field_ ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_copie_field_)
ADD_TEST(check_test_copie_group ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_copie_group)
ADD_TEST(check_test_copie_mesh ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_copie_mesh)
ADD_TEST(check_test_copie_support ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_copie_support)
ADD_TEST(check_test_gibi_driver ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_gibi_driver)
ADD_TEST(check_test_grid ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_grid)
ADD_TEST(check_test_operation_fielddouble ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_operation_fielddouble)
ADD_TEST(check_test_operation_fieldint ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_operation_fieldint)
ADD_TEST(check_test_porflow_driver ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_porflow_driver)
ADD_TEST(check_test_profil_MedFieldDriver ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_profil_MedFieldDriver)
ADD_TEST(check_test_profil_gauss_MedFieldDriver ${PYTHON_EXECUTABLE} ${MEDMEMBINTEST_PYOPT_TEST} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_profil_gauss_MedFieldDriver)
