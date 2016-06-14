# Copyright (C) 2012-2016  CEA/DEN, EDF R&D
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

ADD_DEFINITIONS(${MPI_DEFINITIONS} ${CPPUNIT_DEFINITIONS})

INCLUDE_DIRECTORIES(
  ${MPI_INCLUDE_DIRS}
  ${CPPUNIT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNELTest   # for TestIKUtils.hxx
  )

SET(ParaMEDMEMTest_SOURCES
  ParaMEDMEMTest.cxx
  ParaMEDMEMTest_MPIProcessorGroup.cxx
  ParaMEDMEMTest_BlockTopology.cxx
  ParaMEDMEMTest_InterpKernelDEC.cxx
  ParaMEDMEMTest_StructuredCoincidentDEC.cxx
  ParaMEDMEMTest_ICoco.cxx
  ParaMEDMEMTest_Gauthier1.cxx
  ParaMEDMEMTest_FabienAPI.cxx
  ParaMEDMEMTest_NonCoincidentDEC.cxx
  ParaMEDMEMTest_OverlapDEC.cxx
  MPIAccessDECTest.cxx
  test_AllToAllDEC.cxx
  test_AllToAllvDEC.cxx
  test_AllToAllTimeDEC.cxx
  test_AllToAllvTimeDEC.cxx
  test_AllToAllvTimeDoubleDEC.cxx
  MPIAccessTest.cxx
  test_MPI_Access_Send_Recv.cxx
  test_MPI_Access_Cyclic_Send_Recv.cxx
  test_MPI_Access_SendRecv.cxx
  test_MPI_Access_ISend_IRecv.cxx
  test_MPI_Access_Cyclic_ISend_IRecv.cxx
  test_MPI_Access_ISendRecv.cxx
  test_MPI_Access_Probe.cxx
  test_MPI_Access_IProbe.cxx
  test_MPI_Access_Cancel.cxx
  test_MPI_Access_Send_Recv_Length.cxx
  test_MPI_Access_ISend_IRecv_Length.cxx
  test_MPI_Access_ISend_IRecv_Length_1.cxx
  test_MPI_Access_Time.cxx
  test_MPI_Access_Time_0.cxx
  test_MPI_Access_ISend_IRecv_BottleNeck.cxx
  )

ADD_LIBRARY(ParaMEDMEMTest ${ParaMEDMEMTest_SOURCES})
SET_TARGET_PROPERTIES(ParaMEDMEMTest PROPERTIES COMPILE_FLAGS "")
TARGET_LINK_LIBRARIES(ParaMEDMEMTest paramedmem paramedloader medcouplingremapper InterpKernelTestUtils ${CPPUNIT_LIBRARIES})
INSTALL(TARGETS ParaMEDMEMTest DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

SET(TESTSParaMEDMEM)
SET(TestParaMEDMEM_SOURCES
  TestParaMEDMEM.cxx
  )
SET(TESTSParaMEDMEM ${TESTSParaMEDMEM} TestParaMEDMEM)

SET(TestMPIAccessDEC_SOURCES
  TestMPIAccessDEC.cxx
  )
SET(TESTSParaMEDMEM ${TESTSParaMEDMEM} TestMPIAccessDEC)

SET(TestMPIAccess_SOURCES
  TestMPIAccess.cxx
  )
SET(TESTSParaMEDMEM ${TESTSParaMEDMEM} TestMPIAccess)

SET(test_perf_SOURCES
  test_perf.cxx
  )
SET(TESTSParaMEDMEM ${TESTSParaMEDMEM} test_perf)

IF(MPI2_IS_OK)
  SET(ParaMEDMEMTestMPI2_1_SOURCES
    MPI2Connector.cxx
    ParaMEDMEMTestMPI2_1.cxx
    )
  SET(TESTSParaMEDMEM ${TESTSParaMEDMEM} ParaMEDMEMTestMPI2_1)

  SET(ParaMEDMEMTestMPI2_2_SOURCES
    MPI2Connector.cxx
    ParaMEDMEMTestMPI2_2.cxx
    )
  SET(TESTSParaMEDMEM ${TESTSParaMEDMEM} ParaMEDMEMTestMPI2_2)
ENDIF(MPI2_IS_OK)

FOREACH(bintestparamem ${TESTSParaMEDMEM})
  ADD_EXECUTABLE(${bintestparamem} ${${bintestparamem}_SOURCES})
  TARGET_LINK_LIBRARIES(${bintestparamem} ParaMEDMEMTest)
ENDFOREACH(bintestparamem ${TESTSParaMEDMEM})

# Now add CMake tests - test_perf, ParaMEDMEMTestMPI2_1 and ParaMEDMEMTestMPI2_2
# are left aside, as they are too specific
#
SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

#  -- some tests require 2, 3, 4 or 5 procs --
ADD_TEST(NAME TestParaMEDMEM_Proc2 COMMAND ${MPIEXEC} -np 2 $<TARGET_FILE:TestParaMEDMEM>)
SET_TESTS_PROPERTIES(TestParaMEDMEM_Proc2 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(NAME TestParaMEDMEM_Proc3 COMMAND ${MPIEXEC} -np 3 $<TARGET_FILE:TestParaMEDMEM>)
SET_TESTS_PROPERTIES(TestParaMEDMEM_Proc3 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(NAME TestParaMEDMEM_Proc4 COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:TestParaMEDMEM>)
SET_TESTS_PROPERTIES(TestParaMEDMEM_Proc4 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(NAME TestParaMEDMEM_Proc5 COMMAND ${MPIEXEC} -np 5 $<TARGET_FILE:TestParaMEDMEM>)
SET_TESTS_PROPERTIES(TestParaMEDMEM_Proc5 PROPERTIES ENVIRONMENT "${tests_env}")

ADD_TEST(NAME TestMPIAccess_Proc2 COMMAND ${MPIEXEC} -np 2 $<TARGET_FILE:TestMPIAccess>)
SET_TESTS_PROPERTIES(TestMPIAccess_Proc2 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(NAME TestMPIAccess_Proc3 COMMAND ${MPIEXEC} -np 3 $<TARGET_FILE:TestMPIAccess>)
SET_TESTS_PROPERTIES(TestMPIAccess_Proc3 PROPERTIES ENVIRONMENT "${tests_env}")

ADD_TEST(NAME TestMPIAccessDEC_Proc4 COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:TestMPIAccessDEC>)
SET_TESTS_PROPERTIES(TestMPIAccessDEC_Proc4 PROPERTIES ENVIRONMENT "${tests_env}")

# Installation rules
INSTALL(TARGETS ${TESTSParaMEDMEM} DESTINATION ${MEDCOUPLING_INSTALL_BINS})
SET(COMMON_HEADERS_HXX
  MPIMainTest.hxx
  MPIAccessDECTest.hxx
  MPIAccessTest.hxx
  ParaMEDMEMTest.hxx
  MPI2Connector.hxx
)
INSTALL(FILES ${COMMON_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
