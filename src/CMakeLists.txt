# Copyright (C) 2015-2016  CEA/DEN, EDF R&D
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

# InterpKERNEL
ADD_SUBDIRECTORY(INTERP_KERNEL)
IF(MEDCOUPLING_BUILD_TESTS)
  ADD_SUBDIRECTORY(INTERP_KERNELTest)
ENDIF(MEDCOUPLING_BUILD_TESTS)

# MEDCoupling
ADD_SUBDIRECTORY(MEDCoupling)
IF(MEDCOUPLING_ENABLE_PYTHON)
  ADD_SUBDIRECTORY(MEDCoupling_Swig)
ENDIF(MEDCOUPLING_ENABLE_PYTHON)

IF(NOT MEDCOUPLING_MICROMED)
  # MEDLoader
  ADD_SUBDIRECTORY(MEDLoader)

  # RENUMBER
  IF(MEDCOUPLING_ENABLE_RENUMBER)
    ADD_SUBDIRECTORY(RENUMBER)
    IF(MEDCOUPLING_ENABLE_PYTHON)
      ADD_SUBDIRECTORY(RENUMBER_Swig)
    ENDIF(MEDCOUPLING_ENABLE_PYTHON)
  ENDIF(MEDCOUPLING_ENABLE_RENUMBER)

  # MEDPartitioner
  IF(MEDCOUPLING_ENABLE_PARTITIONER)
    ADD_SUBDIRECTORY(MEDPartitioner)
    IF(MEDCOUPLING_ENABLE_PYTHON)
      ADD_SUBDIRECTORY(MEDPartitioner_Swig)
    ENDIF(MEDCOUPLING_ENABLE_PYTHON)
  ENDIF(MEDCOUPLING_ENABLE_PARTITIONER)
ENDIF(NOT MEDCOUPLING_MICROMED)

IF(MEDCOUPLING_USE_MPI)
  # ParaMEDMEM
  ADD_SUBDIRECTORY(ParaMEDMEM)

  # ParaMEDLoader
  IF(NOT MEDCOUPLING_MICROMED)
    ADD_SUBDIRECTORY(ParaMEDLoader)
  ENDIF(NOT MEDCOUPLING_MICROMED)
  IF(MEDCOUPLING_ENABLE_PYTHON)
    ADD_SUBDIRECTORY(ParaMEDMEM_Swig)
  ENDIF(MEDCOUPLING_ENABLE_PYTHON)
  IF(MEDCOUPLING_BUILD_TESTS)
    IF(NOT MEDCOUPLING_MICROMED)
      ADD_SUBDIRECTORY(ParaMEDMEMTest)
    ENDIF(NOT MEDCOUPLING_MICROMED)
  ENDIF(MEDCOUPLING_BUILD_TESTS)
ENDIF(MEDCOUPLING_USE_MPI)

# Application tests
INSTALL(FILES CTestTestfileInstall.cmake
  DESTINATION ${MEDCOUPLING_INSTALL_TESTS}
  RENAME CTestTestfile.cmake)
