# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
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

ADD_DEFINITIONS(${BOOST_DEFINITIONS} ${CPPUNIT_DEFINITIONS})

INCLUDE_DIRECTORIES(
  ${CPPUNIT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNELTest # for BasicMainTest.hxx
  )

SET(MEDPARTITIONERTest_HEADERS_HXX
  MEDPARTITIONERTest.hxx
  )

SET(MEDPARTITIONERTest_SOURCES
  MEDPARTITIONERTest.cxx
  )

SET(MEDPARTITIONERTest_LDFLAGS medpartitionercpp ${CPPUNIT_LIBRARIES})

IF(MEDTOOL_USE_MPI)
  IF(MEDTOOL_MED_PARTITIONER_PARMETIS)
    SET(MEDPARTITIONERTest_SOURCES ${MEDPARTITIONERTest_SOURCES} MEDPARTITIONERTestPara.cxx)
  ENDIF(MEDTOOL_MED_PARTITIONER_PARMETIS)
ENDIF(MEDTOOL_USE_MPI)

ADD_LIBRARY(MEDPARTITIONERTest SHARED ${MEDPARTITIONERTest_SOURCES})
TARGET_LINK_LIBRARIES(MEDPARTITIONERTest ${MEDPARTITIONERTest_LDFLAGS})
INSTALL(TARGETS MEDPARTITIONERTest DESTINATION ${MEDTOOL_INSTALL_LIBS})

ADD_EXECUTABLE(TestMEDPARTITIONER TestMEDPARTITIONER.cxx)
TARGET_LINK_LIBRARIES(TestMEDPARTITIONER MEDPARTITIONERTest)
INSTALL(TARGETS TestMEDPARTITIONER DESTINATION ${MEDTOOL_INSTALL_BINS})

INSTALL(FILES ${MEDPARTITIONERTest_HEADERS_HXX} DESTINATION ${MEDTOOL_INSTALL_HEADERS})

ADD_TEST(TestMEDPARTITIONER TestMEDPARTITIONER)

# Application tests

SET(TEST_INSTALL_DIRECTORY ${MEDTOOL_INSTALL_SCRIPT_SCRIPTS}/test/MEDCoupling/MEDPartitioner)
INSTALL(TARGETS MEDPARTITIONERTest TestMEDPARTITIONER DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
