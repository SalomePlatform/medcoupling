# Copyright (C) 2012-2023  CEA, EDF
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
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ParaMEDMEM/MPIAccess
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
)

SET(MPIAccessTest_SOURCES
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

ADD_LIBRARY(MPIAccessTest ${MPIAccessTest_SOURCES})
SET_TARGET_PROPERTIES(MPIAccessTest PROPERTIES COMPILE_FLAGS "")
TARGET_LINK_LIBRARIES(MPIAccessTest paramedmem  ${CPPUNIT_LIBRARIES})
INSTALL(TARGETS MPIAccessTest DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

SET(TESTSMPIAccess)
SET(TestMPIAccessDEC_SOURCES
  TestMPIAccessDEC.cxx
  )
LIST(APPEND TESTSMPIAccess TestMPIAccessDEC)

SET(TestMPIAccess_SOURCES
  TestMPIAccess.cxx
  )
LIST(APPEND TESTSMPIAccess TestMPIAccess)

FOREACH(bintest ${TESTSMPIAccess})
  ADD_EXECUTABLE(${bintest} ${${bintest}_SOURCES})
  TARGET_LINK_LIBRARIES(${bintest} MPIAccessTest)
ENDFOREACH()

# Now add CMake tests 
SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

#  -- some tests require 2, 3, 4 or 5 procs --
# MPICH does not support --oversubscribe:
IF(NOT ${MPIEXEC_EXECUTABLE} MATCHES "mpich")
    SET(_oversub_opt "--oversubscribe")
ENDIF()

ADD_TEST(NAME TestMPIAccess_Proc2 COMMAND ${MPIEXEC} -np 2 ${_oversub_opt} $<TARGET_FILE:TestMPIAccess>)
SET_TESTS_PROPERTIES(TestMPIAccess_Proc2 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(NAME TestMPIAccess_Proc3 COMMAND ${MPIEXEC} -np 3 ${_oversub_opt} $<TARGET_FILE:TestMPIAccess>)
SET_TESTS_PROPERTIES(TestMPIAccess_Proc3 PROPERTIES ENVIRONMENT "${tests_env}")

ADD_TEST(NAME TestMPIAccessDEC_Proc4 COMMAND ${MPIEXEC} -np 4 ${_oversub_opt} $<TARGET_FILE:TestMPIAccessDEC>)
SET_TESTS_PROPERTIES(TestMPIAccessDEC_Proc4 PROPERTIES ENVIRONMENT "${tests_env}")

# Installation rules
INSTALL(TARGETS ${TESTSMPIAccess} DESTINATION ${MEDCOUPLING_INSTALL_BINS})
SET(COMMON_HEADERS_HXX
  MPIAccessDECTest.hxx
  MPIAccessTest.hxx
)
INSTALL(FILES ${COMMON_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
