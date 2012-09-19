# Copyright (C) 2007-2012  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
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

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(libMEDMEM_Swig.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(libMEDMEM_Swig.i PROPERTIES SWIG_DEFINITIONS "-shadow")
SET_SOURCE_FILES_PROPERTIES(libMEDMEM_SwigPYTHON_wrap.cxx PROPERTIES COMPILE_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${XDR_DEFINITIONS} ${PLATFORM_DEFINITIONS}")

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${XDR_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SWIG_ADD_MODULE(libMEDMEM_Swig python libMEDMEM_Swig.i)
SWIG_LINK_LIBRARIES(libMEDMEM_Swig ${PYTHON_LIBS} medmem)
SET_TARGET_PROPERTIES(_libMEDMEM_Swig PROPERTIES COMPILE_FLAGS "${PLATFORM_DEFINITIONS}")

IF(${MACHINE} STREQUAL WINDOWS)
  SET_TARGET_PROPERTIES(_libMEDMEM_Swig PROPERTIES DEBUG_OUTPUT_NAME _libMEDMEM_Swig_d)
ENDIF(${MACHINE} STREQUAL WINDOWS)
INSTALL(TARGETS _libMEDMEM_Swig DESTINATION ${MED_salomepythondir})

INSTALL(FILES libMEDMEM_Swig.i my_typemap.i MEDMEM_SWIG_Templates.hxx DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL(FILES medutilities.py ${CMAKE_CURRENT_BINARY_DIR}/libMEDMEM_Swig.py DESTINATION ${MED_salomepythondir})
SET(SCRIPTS
  medmem.py
  med_test1.py
  med_test2.py
  med_test3.py
  med_test_grid.py
  med_test_skin.py
  medMeshing_test.py
  med_field_anal.py
  med_opfield_test.py
  med_opsupp_test.py
  test_gibi.py
  test_porflow.py
  testDriverAscii.py
  testGaussLocalization.py
  test_profil_MedFieldDriver.py
  testWriteAndFam.py
  testMedMemGeneral.py
  testMedMemCxxTests.py
  dumpMEDMEM.py
  ensightMedEnsight_test.py
  medEnsightMed_test.py
  medNumPy_test.py
  med2sauv 
  sauv2med 
  avs2med
)
INSTALL(FILES ${SCRIPTS} DESTINATION ${MED_salomescript_PYTHON}
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
