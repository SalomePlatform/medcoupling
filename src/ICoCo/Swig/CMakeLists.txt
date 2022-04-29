# Copyright (C) 2012-2021  CEA/DEN, EDF R&D
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
# Author: Adrien Bruneton (CEA/DES)
#

#
# Only managing Python tests here. Python wrapping is done solely at the "medcoupling" level
# (see PyWrapping source directory)
#

IF(MEDCOUPLING_BUILD_PY_TESTS)
    INCLUDE(tests.set)
ENDIF()

SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH NOCHECK ${CMAKE_CURRENT_BINARY_DIR}/../../PyWrapping)
SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

FOREACH(test ${ALL_TESTS})
  GET_FILENAME_COMPONENT(testname ${test} NAME_WE)
  ADD_TEST(NAME ${testname}
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  SET_TESTS_PROPERTIES(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
ENDFOREACH()

# Application tests

SET(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/ICoCo_Swig)
INSTALL(FILES ${ALL_TESTS} DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
INSTALL(FILES tests.set DESTINATION ${TEST_INSTALL_DIRECTORY})
