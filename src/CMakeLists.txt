# Copyright (C) 2015  CEA/DEN, EDF R&D
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

SET(MED_ENABLE_MICROMED ${SALOME_MED_MICROMED})

ADD_SUBDIRECTORY(INTERP_KERNEL)
IF(SALOME_BUILD_TESTS)
  ADD_SUBDIRECTORY(INTERP_KERNELTest)
ENDIF(SALOME_BUILD_TESTS)
ADD_SUBDIRECTORY(MEDCoupling)

IF(SALOME_MED_ENABLE_PYTHON)
  ADD_SUBDIRECTORY(MEDCoupling_Swig)
ENDIF(SALOME_MED_ENABLE_PYTHON)

IF(NOT SALOME_MED_MICROMED)
  ADD_SUBDIRECTORY(MEDLoader)
  IF(SALOME_MED_ENABLE_RENUMBER)
    ADD_SUBDIRECTORY(RENUMBER)
    IF(SALOME_MED_ENABLE_PYTHON)
      ADD_SUBDIRECTORY(RENUMBER_Swig)
    ENDIF(SALOME_MED_ENABLE_PYTHON)
  ENDIF(SALOME_MED_ENABLE_RENUMBER)
  IF(SALOME_MED_ENABLE_PARTITIONER)
    ADD_SUBDIRECTORY(MEDPartitioner)
    IF(SALOME_MED_ENABLE_PYTHON)
      ADD_SUBDIRECTORY(MEDPartitioner_Swig)
    ENDIF(SALOME_MED_ENABLE_PYTHON)
  ENDIF(SALOME_MED_ENABLE_PARTITIONER)
ENDIF(NOT SALOME_MED_MICROMED)

IF(SALOME_USE_MPI)
  ADD_SUBDIRECTORY(ParaMEDMEM)
  IF(NOT SALOME_MED_MICROMED)
    ADD_SUBDIRECTORY(ParaMEDLoader)
  ENDIF(NOT SALOME_MED_MICROMED)
  IF(SALOME_MED_ENABLE_PYTHON)
    ADD_SUBDIRECTORY(ParaMEDMEM_Swig)
  ENDIF(SALOME_MED_ENABLE_PYTHON)
  IF(SALOME_BUILD_TESTS)
    IF(NOT SALOME_MED_MICROMED)
      ADD_SUBDIRECTORY(ParaMEDMEMTest)
    ENDIF(NOT SALOME_MED_MICROMED)
  ENDIF(SALOME_BUILD_TESTS)
ENDIF(SALOME_USE_MPI)

# Application tests
INSTALL(FILES CTestTestfileInstall.cmake
  DESTINATION ${MEDTOOL_INSTALL_TESTS}
  RENAME CTestTestfile.cmake)

INSTALL(FILES CTestTestfileInstallMEDCoupling.cmake
  DESTINATION ${MEDTOOL_INSTALL_TESTS}/MEDCoupling
  RENAME CTestTestfile.cmake)
