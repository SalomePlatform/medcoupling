# Copyright (C) 2015-2025  CEA, EDF
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

if(MED_INT_IS_LONG)
  add_definitions(-DMED_INT_IS_LONG)
endif(MED_INT_IS_LONG)

# InterpKERNEL
add_subdirectory(INTERP_KERNEL)
if(MEDCOUPLING_BUILD_TESTS)
  add_subdirectory(INTERP_KERNELTest)
endif(MEDCOUPLING_BUILD_TESTS)

# MEDCoupling
add_subdirectory(MEDCoupling)
add_subdirectory(ICoCo)

# Shape recognition
if(MEDCOUPLING_ENABLE_SHAPERECOGN)
  add_subdirectory(ShapeRecogn)
endif(MEDCOUPLING_ENABLE_SHAPERECOGN)

if(MEDCOUPLING_ENABLE_PYTHON)
  add_subdirectory(MEDCoupling_Swig)

  find_package(SalomeBasicIterativeStatistics)
  if(BASICITERATIVESTATISTICS_FOUND)
    add_subdirectory(Stat)
  endif(BASICITERATIVESTATISTICS_FOUND)

  add_subdirectory(PyWrapping)

endif(MEDCOUPLING_ENABLE_PYTHON)

if(NOT MEDCOUPLING_MICROMED)
  # MEDLoader
  add_subdirectory(MEDLoader)

  # RENUMBER
  if(MEDCOUPLING_ENABLE_RENUMBER)
    add_subdirectory(RENUMBER)
    if(MEDCOUPLING_ENABLE_PYTHON)
      add_subdirectory(RENUMBER_Swig)
    endif(MEDCOUPLING_ENABLE_PYTHON)
  endif(MEDCOUPLING_ENABLE_RENUMBER)

  # MEDPartitioner
  if(MEDCOUPLING_ENABLE_PARTITIONER)
    add_subdirectory(MEDPartitioner)
    if(MEDCOUPLING_ENABLE_PYTHON)
      add_subdirectory(MEDPartitioner_Swig)
    endif(MEDCOUPLING_ENABLE_PYTHON)
  endif(MEDCOUPLING_ENABLE_PARTITIONER)
endif(NOT MEDCOUPLING_MICROMED)

if(MEDCOUPLING_USE_MPI)
  # ParaMEDMEM
  add_subdirectory(ParaMEDMEM)

  # ParaMEDLoader
  if(NOT MEDCOUPLING_MICROMED)
    add_subdirectory(ParaMEDLoader)
  endif(NOT MEDCOUPLING_MICROMED)
  if(MEDCOUPLING_ENABLE_PYTHON)
    add_subdirectory(ParaMEDMEM_Swig)
  endif(MEDCOUPLING_ENABLE_PYTHON)
  if(MEDCOUPLING_BUILD_TESTS)
    if(NOT MEDCOUPLING_MICROMED)
      add_subdirectory(ParaMEDMEMTest)
    endif(NOT MEDCOUPLING_MICROMED)
  endif(MEDCOUPLING_BUILD_TESTS)
endif(MEDCOUPLING_USE_MPI)

# Application tests
configure_file(CTestTestfileInstall.cmake.in "CTestTestfileST.cmake" @ONLY)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfileST.cmake
  DESTINATION ${MEDCOUPLING_INSTALL_TESTS}
  RENAME CTestTestfile.cmake)
