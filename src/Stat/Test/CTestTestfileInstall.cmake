# Copyright (C) 2024-2026  CEA, EDF
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

set(TEST_NAMES TestMEDCouplingIterativeStatistics)

set(PYTHONPATH $ENV{PYTHONPATH})
cmake_path(APPEND PYTHONPATH "../../bin")
foreach(tfile ${TEST_NAMES})
  set(TEST_NAME ${COMPONENT_NAME}_${tfile})
  add_test(${TEST_NAME} python3 ${tfile}.py)
  set(TEST_ENVIRONMENT)
  set_tests_properties(
    ${TEST_NAME} PROPERTIES LABELS "${COMPONENT_NAME}" TIMEOUT ${TIMEOUT}
                            ENVIRONMENT "${PYTHONPATH}")
endforeach()
