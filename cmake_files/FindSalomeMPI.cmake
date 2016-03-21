# Copyright (C) 2013-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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
# Author: Adrien Bruneton
#

# MPI detection for Salome
#
#  !! Please read the generic detection procedure in SalomeMacros.cmake !!
# 

SALOME_FIND_PACKAGE_AND_DETECT_CONFLICTS(MPI MPIEXEC 2)
MARK_AS_ADVANCED(MPI_EXTRA_LIBRARY MPI_LIBRARY)

SET(MPI_INCLUDE_DIRS ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})
SET(MPI_LIBRARIES ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

IF(MPI_FOUND) 
  # Detect if function MPI_Publish_name is provided by the external MPI library 
  # otherwise take ours.
  include(CheckSymbolExists)
  SET(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
  SET(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
  CHECK_SYMBOL_EXISTS(MPI_Publish_name mpi.h MPI2_IS_OK)
  SET(MPI_DEFINITIONS "${MPI_CXX_COMPILE_FLAGS}")
  IF(MPI2_IS_OK)
    MESSAGE(STATUS "Your mpi implementation is compatible with mpi2 ... adding -DHAVE_MPI2")
    SET(MPI_DEFINITIONS "${MPI_CXX_COMPILE_FLAGS} -DHAVE_MPI2")
  ENDIF(MPI2_IS_OK)

  SALOME_ACCUMULATE_HEADERS(MPI_INCLUDE_DIRS)
  SALOME_ACCUMULATE_ENVIRONMENT(LD_LIBRARY_PATH ${MPI_LIBRARIES})
ENDIF()
