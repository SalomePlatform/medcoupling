# Copyright (C) 2012-2024  CEA, EDF
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

# Set MEDCOUPLING_SERVER_BASED_SEARCH=ON when building docs for putting on www.salome-platform.org
# MEDCOUPLING_SERVER_BASED_SEARCH is used in Doxyfile
#
IF(MEDCOUPLING_SERVER_BASED_SEARCH)
  SET(MEDCOUPLING_SERVER_BASED_SEARCH "YES")
ELSE()
  SET(MEDCOUPLING_SERVER_BASED_SEARCH "NO")
ENDIF()

configure_file(Doxyfile.in Doxyfile @ONLY)
configure_file(static/header.html.in static/header.html @ONLY)


IF(NOT MEDCOUPLING_ENABLE_PYTHON)
  ADD_CUSTOM_TARGET(doxygen_documentation
    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
    VERBATIM
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
ELSE()
  INCLUDE(doxy2swig/doxy2swig.cmake)

  # >>>>>>>>>>>>>>>>
  # For ease of maintenance, documentation for code examples is splitted in
  # several files. We here splice to a single file before running Doxygen.

  SET(EXAMPLE_FILES # files to concatenate: order is important!
    doxfiles/examples/medcouplingexamplesheader.doxy
    doxfiles/examples/medcouplingexamplesfields.doxy
    doxfiles/examples/medcouplingexamplesmeshes.doxy
    doxfiles/examples/medcouplingexamplesarrays.doxy
    doxfiles/examples/medcouplingexamplesother.doxy
    doxfiles/examples/medcouplingexamplesfooter.doxy
    )

  # This function adds IN_FILE contents to the end of OUT_FILE
  FUNCTION(concat IN_FILE OUT_FILE)
    FILE(READ ${IN_FILE} CONTENTS)
    FILE(APPEND ${OUT_FILE} ${CONTENTS})
  ENDFUNCTION()

  FILE(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}/tmp/medcouplingexamples.in" input)
  FILE(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/BuildPyExamplesFromCPP.py" pythondocexamplesgenerator)
  FILE(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}" output)

  # Prepare a temporary file to "concat" to:
  FILE(WRITE ${input} "")
  # Call the "concat" function for each example file
  FOREACH(EXAMPLE_FILE ${EXAMPLE_FILES})
    concat(${EXAMPLE_FILE} ${input})
  ENDFOREACH()
  # Note: The reason for writing to a temporary is so the real target file only
  # gets updated if its content has changed.
  # <<<<<<<<<<<<<<<<

  ADD_CUSTOM_TARGET(doxygen_documentation
    COMMAND ${PYTHON_EXECUTABLE} ${pythondocexamplesgenerator} ${input} ${output}
    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
    VERBATIM
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  # Create dummy target gathering the generation of all .i files:
  ADD_CUSTOM_TARGET(swig_ready)
  ADD_DEPENDENCIES(swig_ready doxygen_documentation)

  MEDCOUPLING_SWIG_DOCSTRING_CONFIGURE(doxygen_documentation swig_ready MEDCoupling)
  MEDCOUPLING_SWIG_DOCSTRING_CONFIGURE(doxygen_documentation swig_ready MEDLoader)

  SET(doxyfile_path ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  FILE(STRINGS ${doxyfile_path} enabled_sections REGEX "ENABLED_SECTIONS")
  IF(enabled_sections)
    FILE(READ ${doxyfile_path} doxyfile_content)
    STRING(REPLACE ${enabled_sections} "${enabled_sections} ENABLE_EXAMPLES" new_doxyfile ${doxyfile_content})
    FILE(WRITE ${doxyfile_path} ${new_doxyfile})
  ELSE()
    FILE(APPEND ${doxyfile_path}
      "#Temporary variable to enable python documentation sections\nENABLED_SECTIONS = ENABLE_EXAMPLES")
  ENDIF()

  # Swig generation to have docstrings correctly populated:
  INCLUDE(doxy2swig/doxy2swig.cmake)
ENDIF()

INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc_ref_dev/html/ DESTINATION ${MEDCOUPLING_INSTALL_DOC}/developer)

SET(MAKE_CLEAN_FILES doc_ref_dev tmp)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${MAKE_CLEAN_FILES}")
