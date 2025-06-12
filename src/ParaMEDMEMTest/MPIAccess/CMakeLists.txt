# Copyright (C) 2012-2025  CEA, EDF
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

add_definitions(${MPI_DEFINITIONS} ${CPPUNIT_DEFINITIONS})

include_directories(
  ${MPI_INCLUDE_DIRS}
  ${CPPUNIT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ParaMEDMEM/MPIAccess
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases)

set(MPIAccessTest_SOURCES
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
    test_MPI_Access_ISend_IRecv_BottleNeck.cxx)

add_library(MPIAccessTest ${MPIAccessTest_SOURCES})
set_target_properties(MPIAccessTest PROPERTIES COMPILE_FLAGS "")
target_link_libraries(MPIAccessTest paramedmem ${CPPUNIT_LIBRARIES})
install(TARGETS MPIAccessTest DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

set(TESTSMPIAccess)
set(TestMPIAccessDEC_SOURCES TestMPIAccessDEC.cxx)
list(APPEND TESTSMPIAccess TestMPIAccessDEC)

set(TestMPIAccess_SOURCES TestMPIAccess.cxx)
list(APPEND TESTSMPIAccess TestMPIAccess)

foreach(bintest ${TESTSMPIAccess})
  add_executable(${bintest} ${${bintest}_SOURCES})
  target_link_libraries(${bintest} MPIAccessTest)
endforeach()

# Now add CMake tests
salome_generate_tests_environment(tests_env)

# -- some tests require 2, 3, 4 or 5 procs -- MPICH does not support
# --oversubscribe:
if(NOT ${MPIEXEC_EXECUTABLE} MATCHES "mpich")
  set(_oversub_opt "--oversubscribe")
endif()

add_test(NAME TestMPIAccess_Proc2 COMMAND ${MPIEXEC} -np 2 ${_oversub_opt}
                                          $<TARGET_FILE:TestMPIAccess>)
set_tests_properties(TestMPIAccess_Proc2 PROPERTIES ENVIRONMENT "${tests_env}")
add_test(NAME TestMPIAccess_Proc3 COMMAND ${MPIEXEC} -np 3 ${_oversub_opt}
                                          $<TARGET_FILE:TestMPIAccess>)
set_tests_properties(TestMPIAccess_Proc3 PROPERTIES ENVIRONMENT "${tests_env}")

add_test(NAME TestMPIAccessDEC_Proc4 COMMAND ${MPIEXEC} -np 4 ${_oversub_opt}
                                             $<TARGET_FILE:TestMPIAccessDEC>)
set_tests_properties(TestMPIAccessDEC_Proc4 PROPERTIES ENVIRONMENT
                                                       "${tests_env}")

# Installation rules
install(TARGETS ${TESTSMPIAccess} DESTINATION ${MEDCOUPLING_INSTALL_BINS})
set(COMMON_HEADERS_HXX MPIAccessDECTest.hxx MPIAccessTest.hxx)
install(FILES ${COMMON_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
