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

SET_SOURCE_FILES_PROPERTIES(libMedCorba_Swig.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(libMedCorba_Swig.i PROPERTIES SWIG_FLAGS "-shadow")

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDES_DIR}
  ${OMNIORB_INCLUDES_DIR}
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
  ${KERNEL_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}/../../idl
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MED
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM_I
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM_SWIG
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SWIG_ADD_MODULE(libMedCorba_Swig python libMedCorba_Swig.i)
SWIG_LINK_LIBRARIES(libMedCorba_Swig ${PYTHON_LIBRARIES} MEDMEMImpl MEDEngine)
SET_TARGET_PROPERTIES(_libMedCorba_Swig PROPERTIES COMPILE_FLAGS "${OMNIORB_FLAGS} ${HDF5_FLAGS} ${MED3_FLAGS} ${XDR_FLAGS}")

SET(MED_MedCorba_Swig_salomescript_DATA
  medcorba_test.py batchmode_medcorba_test.py batchmode_medcorba_test1.py
  )

INSTALL(TARGETS _libMedCorba_Swig DESTINATION ${MED_salomepythondir})
INSTALL(FILES libMedCorba_Swig.i DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/libMedCorba_Swig.py DESTINATION ${MED_salomescript_PYTHON})
INSTALL(FILES ${MED_MedCorba_Swig_salomescript_DATA} DESTINATION ${MED_salomescript_DATA})
