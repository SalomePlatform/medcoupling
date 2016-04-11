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

INCLUDE(${SWIG_USE_FILE})

ADD_DEFINITIONS(${PYTHON_DEFINITIONS} ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS})

SET_SOURCE_FILES_PROPERTIES(MEDCoupling.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDCoupling.i PROPERTIES SWIG_DEFINITIONS "-shadow")
SET(SWIG_MODULE_MEDCoupling_EXTRA_FLAGS "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")

SET (MEDCoupling_SWIG_DPYS_FILES
    MEDCouplingCommon.i
    MEDCouplingMemArray.i
    MEDCouplingFieldDiscretization.i
    MEDCouplingFinalize.i
    MEDCouplingTypemaps.i)

SET (MC_pyTestFiles
    MEDCouplingBasicsTest1.py
    MEDCouplingBasicsTest2.py
    MEDCouplingBasicsTest3.py
    MEDCouplingBasicsTest4.py
    MEDCouplingBasicsTest5.py
    MEDCouplingRemapperTest.py
    MEDCouplingDataForTest.py
    MEDCouplingNumPyTest.py
    MEDCouplingPickleTest.py
)

SET (MC_Swig_interf
    MEDCoupling.i
    MEDCouplingCommon.i
    MEDCouplingRefCountObject.i
    MEDCouplingMemArray.i
    MEDCouplingFieldDiscretization.i
    MEDCouplingTimeDiscretization.i
    MEDCouplingFinalize.i
    MEDCouplingRemapper.i
    MEDCouplingTypemaps.i
    MEDCouplingDataArrayTypemaps.i
)

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${NUMPY_INCLUDE_DIR}
  ${PTHREAD_INCLUDE_DIR} # pthread dependancy due to python2.7 library
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  ${PROJECT_BINARY_DIR}/doc
  )

# _ABR_ Ensure dependency mechanism on all SWIG files and headers
SET (SWIG_MODULE_MEDCoupling_EXTRA_DEPS ${MEDCoupling_SWIG_DPYS_FILES}
    ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX}
    ${interpkernel_HEADERS_HXX} ${interpkernel_HEADERS_TXX})

# SWIG must run after the doc if we want to have the docstrings extracted from Doxygen
# into the Python module:
IF(MEDCOUPLING_BUILD_DOC)
    LIST(APPEND SWIG_MODULE_MEDCoupling_EXTRA_FLAGS -DWITH_DOCSTRINGS)
    LIST(APPEND SWIG_MODULE_MEDCoupling_EXTRA_DEPS
        ${PROJECT_BINARY_DIR}/doc/MEDCoupling_doc.i
        swig_ready)
ENDIF()

SWIG_ADD_MODULE(MEDCoupling python MEDCoupling.i)
SWIG_LINK_LIBRARIES(MEDCoupling ${PYTHON_LIBRARIES} ${PLATFORM_LIBS} medcoupling)

SET_SOURCE_FILES_PROPERTIES(MEDCouplingRemapper.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDCouplingRemapper.i PROPERTIES SWIG_DEFINITIONS "-shadow")
SET(SWIG_MODULE_MEDCouplingRemapper_EXTRA_FLAGS "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")

# _ABR_ Ensure dependency mechanism on all SWIG files and headers
SET (SWIG_MODULE_MEDCouplingRemapper_EXTRA_DEPS ${MEDCoupling_SWIG_DPYS_FILES}
    ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX}
    ${interpkernel_HEADERS_HXX} ${interpkernel_HEADERS_TXX})

SWIG_ADD_MODULE(MEDCouplingRemapper python MEDCouplingRemapper.i)
SWIG_LINK_LIBRARIES(MEDCouplingRemapper ${PYTHON_LIBRARIES} medcouplingremapper)

IF(WIN32)
  SET_TARGET_PROPERTIES(_MEDCouplingRemapper PROPERTIES DEBUG_OUTPUT_NAME _MEDCouplingRemapper_d)
  SET_TARGET_PROPERTIES(_MEDCoupling PROPERTIES DEBUG_OUTPUT_NAME _MEDCoupling_d)
ENDIF(WIN32)
INSTALL(TARGETS ${SWIG_MODULE_MEDCoupling_REAL_NAME} ${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME} DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

SET(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDCoupling.py ${CMAKE_CURRENT_BINARY_DIR}/MEDCouplingRemapper.py)
SALOME_INSTALL_SCRIPTS("${PYFILES_TO_INSTALL}" ${MEDCOUPLING_INSTALL_PYTHON})

INSTALL(FILES ${MC_Swig_interf} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
INSTALL(FILES ${MC_pyTestFiles} DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})
INSTALL(FILES MEDCouplingExamplesTest.py DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})

SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

ADD_TEST(MEDCouplingBasicsTest1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingBasicsTest1.py)
SET_TESTS_PROPERTIES(MEDCouplingBasicsTest1 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(MEDCouplingBasicsTest2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingBasicsTest2.py)
SET_TESTS_PROPERTIES(MEDCouplingBasicsTest2 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(MEDCouplingBasicsTest3 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingBasicsTest3.py)
SET_TESTS_PROPERTIES(MEDCouplingBasicsTest3 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(MEDCouplingBasicsTest4 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingBasicsTest4.py)
SET_TESTS_PROPERTIES(MEDCouplingBasicsTest4 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(MEDCouplingBasicsTest5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingBasicsTest5.py)
SET_TESTS_PROPERTIES(MEDCouplingBasicsTest5 PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(MEDCouplingExamplesTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingExamplesTest.py)
SET_TESTS_PROPERTIES(MEDCouplingExamplesTest PROPERTIES ENVIRONMENT "${tests_env}")
ADD_TEST(MEDCouplingRemapperTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingRemapperTest.py)
SET_TESTS_PROPERTIES(MEDCouplingRemapperTest PROPERTIES ENVIRONMENT "${tests_env}")

IF(NUMPY_FOUND)
  ADD_TEST(MEDCouplingNumPyTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingNumPyTest.py)
  SET_TESTS_PROPERTIES(MEDCouplingNumPyTest PROPERTIES ENVIRONMENT "${tests_env}")
  ADD_TEST(MEDCouplingPickleTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDCouplingPickleTest.py)
  SET_TESTS_PROPERTIES(MEDCouplingPickleTest PROPERTIES ENVIRONMENT "${tests_env}")
ENDIF(NUMPY_FOUND)

# Application tests

SET(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDCoupling_Swig)
LIST(APPEND MC_pyTestFiles MEDCouplingExamplesTest.py)
INSTALL(FILES ${MC_pyTestFiles} DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
