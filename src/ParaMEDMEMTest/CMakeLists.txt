# Copyright (C) 2012-2026  CEA, EDF
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this library; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email :
# webmaster.salome@opencascade.com
#

add_subdirectory(MPIAccess) # maybe to delete one day ...

add_definitions(${MPI_DEFINITIONS} ${CPPUNIT_DEFINITIONS} ${HDF5_DEFINITIONS}
                ${MEDFILE_DEFINITIONS})

include_directories(
  ${MPI_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CPPUNIT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM/MPIAccess
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNELTest # for TestIKUtils.hxx
)

set(ParaMEDMEMTest_SOURCES
    ParaMEDMEMTest.cxx
    ParaMEDMEMTest_MPIProcessorGroup.cxx
    ParaMEDMEMTest_ByStringMPIProcessorGroup.cxx
    ParaMEDMEMTest_BlockTopology.cxx
    ParaMEDMEMTest_InterpKernelDEC.cxx
    ParaMEDMEMTest_StructuredCoincidentDEC.cxx
    ParaMEDMEMTest_ICoco.cxx
    ParaMEDMEMTest_Gauthier1.cxx
    ParaMEDMEMTest_FabienAPI.cxx
    ParaMEDMEMTest_NonCoincidentDEC.cxx
    ParaMEDMEMTest_OverlapDEC.cxx
    ParaMEDMEMTest_MEDLoader.cxx)

add_library(ParaMEDMEMTest ${ParaMEDMEMTest_SOURCES})
set_target_properties(ParaMEDMEMTest PROPERTIES COMPILE_FLAGS "")
target_link_libraries(
  ParaMEDMEMTest paramedmem paramedloader medcouplingremapper
  InterpKernelTestUtils ${CPPUNIT_LIBRARIES})
install(TARGETS ParaMEDMEMTest DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

set(TESTSParaMEDMEM)
set(TestParaMEDMEM_SOURCES TestParaMEDMEM.cxx)
list(APPEND TESTSParaMEDMEM TestParaMEDMEM)

set(test_perf_SOURCES test_perf.cxx)
list(APPEND TESTSParaMEDMEM test_perf)

if(MPI2_IS_OK)
  set(ParaMEDMEMTestMPI2_1_SOURCES MPI2Connector.cxx ParaMEDMEMTestMPI2_1.cxx)
  list(APPEND TESTSParaMEDMEM ParaMEDMEMTestMPI2_1)

  set(ParaMEDMEMTestMPI2_2_SOURCES MPI2Connector.cxx ParaMEDMEMTestMPI2_2.cxx)
  list(APPEND TESTSParaMEDMEM ParaMEDMEMTestMPI2_2)
endif(MPI2_IS_OK)

foreach(bintestparamem ${TESTSParaMEDMEM})
  add_executable(${bintestparamem} ${${bintestparamem}_SOURCES})
  target_link_libraries(${bintestparamem} ParaMEDMEMTest)
endforeach(bintestparamem ${TESTSParaMEDMEM})

# Now add CMake tests - test_perf, ParaMEDMEMTestMPI2_1 and ParaMEDMEMTestMPI2_2
# are left aside, as they are too specific
#
salome_generate_tests_environment(tests_env)

# -- some tests require 2, 3, 4 or 5 procs -- MPICH does not support
# --oversubscribe:
if(NOT ${MPIEXEC_EXECUTABLE} MATCHES "mpich")
  set(_oversub_opt "--oversubscribe")
endif()

add_test(NAME TestParaMEDMEM_Proc2 COMMAND ${MPIEXEC} -np 2 ${_oversub_opt}
                                           $<TARGET_FILE:TestParaMEDMEM>)
set_tests_properties(TestParaMEDMEM_Proc2 PROPERTIES ENVIRONMENT "${tests_env}")
add_test(NAME TestParaMEDMEM_Proc3 COMMAND ${MPIEXEC} -np 3 ${_oversub_opt}
                                           $<TARGET_FILE:TestParaMEDMEM>)
set_tests_properties(TestParaMEDMEM_Proc3 PROPERTIES ENVIRONMENT "${tests_env}")
add_test(NAME TestParaMEDMEM_Proc4 COMMAND ${MPIEXEC} -np 4 ${_oversub_opt}
                                           $<TARGET_FILE:TestParaMEDMEM>)
set_tests_properties(TestParaMEDMEM_Proc4 PROPERTIES ENVIRONMENT "${tests_env}")
add_test(NAME TestParaMEDMEM_Proc5 COMMAND ${MPIEXEC} -np 5 ${_oversub_opt}
                                           $<TARGET_FILE:TestParaMEDMEM>)
set_tests_properties(TestParaMEDMEM_Proc5 PROPERTIES ENVIRONMENT "${tests_env}")

# Installation rules
install(TARGETS ${TESTSParaMEDMEM} DESTINATION ${MEDCOUPLING_INSTALL_BINS})
set(COMMON_HEADERS_HXX MPIMainTest.hxx ParaMEDMEMTest.hxx MPI2Connector.hxx)
install(FILES ${COMMON_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# Application tests

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/ParaMEDMEMTest)
install(TARGETS ${TESTSParaMEDMEM} DESTINATION ${TEST_INSTALL_DIRECTORY})

# export MPIEXEC and _oversub_opt to CTestTestfile.cmake of salome test
# mechanism
configure_file(CTestTestfileInstall.cmake.in "CTestTestfileST.cmake" @ONLY)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfileST.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
