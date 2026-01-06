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

# Set MEDCOUPLING_SERVER_BASED_SEARCH=ON when building docs for putting on
# www.salome-platform.org MEDCOUPLING_SERVER_BASED_SEARCH is used in Doxyfile
#
if(MEDCOUPLING_SERVER_BASED_SEARCH)
  set(MEDCOUPLING_SERVER_BASED_SEARCH "YES")
else()
  set(MEDCOUPLING_SERVER_BASED_SEARCH "NO")
endif()

configure_file(Doxyfile.in Doxyfile @ONLY)
configure_file(static/header.html.in static/header.html @ONLY)

if(NOT MEDCOUPLING_ENABLE_PYTHON)
  add_custom_target(
    doxygen_documentation
    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
    VERBATIM
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
else()
  include(doxy2swig/doxy2swig.cmake)

  # >>>>>>>>>>>>>>>>
  # For ease of maintenance, documentation for code examples is splitted in
  # several files. We here splice to a single file before running Doxygen.

  set(EXAMPLE_FILES # files to concatenate: order is important!
      doxfiles/examples/medcouplingexamplesheader.doxy
      doxfiles/examples/medcouplingexamplesfields.doxy
      doxfiles/examples/medcouplingexamplesmeshes.doxy
      doxfiles/examples/medcouplingexamplesarrays.doxy
      doxfiles/examples/medcouplingexamplesother.doxy
      doxfiles/examples/medcouplingexamplesfooter.doxy)

  # This function adds IN_FILE contents to the end of OUT_FILE
  function(concat IN_FILE OUT_FILE)
    file(READ ${IN_FILE} CONTENTS)
    file(APPEND ${OUT_FILE} ${CONTENTS})
  endfunction()

  file(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}/tmp/medcouplingexamples.in"
       input)
  file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/BuildPyExamplesFromCPP.py"
       pythondocexamplesgenerator)
  file(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}" output)

  # Prepare a temporary file to "concat" to:
  file(WRITE ${input} "")
  # Call the "concat" function for each example file
  foreach(EXAMPLE_FILE ${EXAMPLE_FILES})
    concat(${EXAMPLE_FILE} ${input})
  endforeach()
  # Note: The reason for writing to a temporary is so the real target file only
  # gets updated if its content has changed.
  # <<<<<<<<<<<<<<<<

  add_custom_target(
    doxygen_documentation
    COMMAND ${PYTHON_EXECUTABLE} ${pythondocexamplesgenerator} ${input}
            ${output}
    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
    VERBATIM
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Create dummy target gathering the generation of all .i files:
  add_custom_target(swig_ready)
  add_dependencies(swig_ready doxygen_documentation)

  medcoupling_swig_docstring_configure(doxygen_documentation swig_ready
                                       MEDCoupling)
  medcoupling_swig_docstring_configure(doxygen_documentation swig_ready
                                       MEDLoader)

  set(doxyfile_path ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  file(STRINGS ${doxyfile_path} enabled_sections REGEX "ENABLED_SECTIONS")
  if(enabled_sections)
    file(READ ${doxyfile_path} doxyfile_content)
    string(REPLACE ${enabled_sections} "${enabled_sections} ENABLE_EXAMPLES"
                   new_doxyfile ${doxyfile_content})
    file(WRITE ${doxyfile_path} ${new_doxyfile})
  else()
    file(
      APPEND ${doxyfile_path}
      "#Temporary variable to enable python documentation sections\nENABLED_SECTIONS = ENABLE_EXAMPLES"
    )
  endif()

  # Swig generation to have docstrings correctly populated:
  include(doxy2swig/doxy2swig.cmake)
endif()

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc_ref_dev/html/
        DESTINATION ${MEDCOUPLING_INSTALL_DOC}/developer)

set(MAKE_CLEAN_FILES doc_ref_dev tmp)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                    "${MAKE_CLEAN_FILES}")
