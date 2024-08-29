# Copyright (C) 2012-2024  CEA, EDF
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

# Its REALLY misleading to enable MPI and compile para_medfile and not having
# parallel loading and writing abilities if this option is OFF. It might be
# simpler for everyone if it just fails to compile when hdf5 is not parallel.
# if(HDF5_ENABLE_PARALLEL OR HDF5_IS_PARALLEL)
#   add_definitions("-DHDF5_IS_PARALLEL")
# endif(HDF5_ENABLE_PARALLEL OR HDF5_IS_PARALLEL)

set(paramedloader_SOURCES ParaMEDLoader.cxx ParaMEDFileMesh.cxx)

add_library(paramedloader ${paramedloader_SOURCES})

target_include_directories(
  paramedloader PUBLIC
  "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/para_medfile>"
)

target_compile_definitions(paramedloader PUBLIC HDF5_IS_PARALLEL)
if(NOT DEFINED MSVC)
  target_compile_options(paramedloader PRIVATE -Wsign-compare -Wconversion)
endif()

target_link_libraries(paramedloader paramedmem medloader)
install(
  TARGETS paramedloader
  EXPORT ${PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

file(GLOB paramedloader_HEADERS_HXX "${PROJECT_SOURCE_DIR}/include/para_medfile/*")
install(FILES ${paramedloader_HEADERS_HXX}
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
