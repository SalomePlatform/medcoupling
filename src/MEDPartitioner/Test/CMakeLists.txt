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

add_definitions(${BOOST_DEFINITIONS} ${CPPUNIT_DEFINITIONS})

include_directories(
  ${CPPUNIT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNELTest # for BasicMainTest.hxx
)

set(MEDPARTITIONERTest_HEADERS_HXX MEDPARTITIONERTest.hxx)

set(MEDPARTITIONERTest_SOURCES MEDPARTITIONERTest.cxx)

set(MEDPARTITIONERTest_LDFLAGS medpartitionercpp ${CPPUNIT_LIBRARIES})

if(MEDCOUPLING_USE_MPI)
  if(MEDCOUPLING_PARTITIONER_PARMETIS)
    set(MEDPARTITIONERTest_SOURCES ${MEDPARTITIONERTest_SOURCES}
                                   MEDPARTITIONERTestPara.cxx)
  endif(MEDCOUPLING_PARTITIONER_PARMETIS)
endif(MEDCOUPLING_USE_MPI)

add_library(MEDPARTITIONERTest ${MEDPARTITIONERTest_SOURCES})
target_link_libraries(MEDPARTITIONERTest ${MEDPARTITIONERTest_LDFLAGS})
install(TARGETS MEDPARTITIONERTest DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

add_executable(TestMEDPARTITIONER TestMEDPARTITIONER.cxx)
target_link_libraries(TestMEDPARTITIONER MEDPARTITIONERTest)
install(TARGETS TestMEDPARTITIONER DESTINATION ${MEDCOUPLING_INSTALL_BINS})

install(FILES ${MEDPARTITIONERTest_HEADERS_HXX}
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

salome_generate_tests_environment(tests_env)
add_test(TestMEDPARTITIONER ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDLoader/Test/MCTestLauncher.py
         ${CMAKE_CURRENT_BINARY_DIR}/TestMEDPARTITIONER)
set_tests_properties(TestMEDPARTITIONER PROPERTIES ENVIRONMENT "${tests_env}")

# Application tests

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDPartitioner)
install(TARGETS MEDPARTITIONERTest TestMEDPARTITIONER
        DESTINATION ${TEST_INSTALL_DIRECTORY})

install(
  FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
