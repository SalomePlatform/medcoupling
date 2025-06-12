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
# Author: Adrien Bruneton (CEA/DES)
#

#
# Only managing Python tests here. Python wrapping is done solely at the
# "medcoupling" level (see PyWrapping source directory)
#

if(MEDCOUPLING_BUILD_PY_TESTS)
  include(tests.set)
endif()

salome_accumulate_environment(PYTHONPATH NOCHECK
                              ${CMAKE_CURRENT_BINARY_DIR}/../../PyWrapping)
salome_generate_tests_environment(tests_env)

foreach(test ${ALL_TESTS})
  get_filename_component(testname ${test} NAME_WE)
  add_test(NAME ${testname} COMMAND ${PYTHON_EXECUTABLE}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  set_tests_properties(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
endforeach()

# Application tests

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/ICoCo_Swig)
install(FILES ${ALL_TESTS} DESTINATION ${TEST_INSTALL_DIRECTORY})

install(
  FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
install(FILES tests.set DESTINATION ${TEST_INSTALL_DIRECTORY})
