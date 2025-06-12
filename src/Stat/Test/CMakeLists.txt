# Copyright (C) 2024  CEA, EDF
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

salome_accumulate_environment(PYTHONPATH NOCHECK
                              ${CMAKE_CURRENT_BINARY_DIR}/../../Stat)
salome_generate_tests_environment(tests_env)

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/Stat)

add_test(
  NAME MEDCouplingIterativeStatistics
  COMMAND ${PYTHON_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/TestMEDCouplingIterativeStatistics.py)
set_tests_properties(MEDCouplingIterativeStatistics PROPERTIES ENVIRONMENT
                                                               "${tests_env}")

install(FILES TestMEDCouplingIterativeStatistics.py
        DESTINATION ${TEST_INSTALL_DIRECTORY})
install(
  FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
