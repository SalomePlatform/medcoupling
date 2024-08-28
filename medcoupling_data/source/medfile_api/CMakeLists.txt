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
#
# Author : Anthony Geay (CEA/DEN)

set(MEDLOADER_SOURCES
    CrackAlgo.cxx
    MEDLoader.cxx
    MEDLoaderBase.cxx
    MEDLoaderTraits.cxx
    MEDFileUtilities.cxx
    MEDFileMesh.cxx
    MEDFileMeshElt.cxx
    MEDFileBasis.cxx
    MEDFileMeshLL.cxx
    MEDFileField.cxx
    MEDFileFieldInternal.cxx
    MEDFileFieldGlobs.cxx
    MEDFileField1TS.cxx
    MEDFileFieldMultiTS.cxx
    MEDFileJoint.cxx
    MEDFileEquivalence.cxx
    MEDFileParameter.cxx
    MEDFileData.cxx
    MEDFileFieldOverView.cxx
    MEDFileMeshReadSelector.cxx
    MEDFileMeshSupport.cxx
    MEDFileStructureElement.cxx
    MEDFileEntities.cxx
    MEDFileBlowStrEltUp.cxx
    SauvMedConvertor.cxx
    SauvReader.cxx
    SauvWriter.cxx
    MEDMESHConverterUtilities.cxx
    MeshFormatReader.cxx
    MeshFormatWriter.cxx
    libmesh5.cxx)

# add_definitions(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${XDR_DEFINITIONS})

add_library(medloader ${MEDLOADER_SOURCES})
add_library(MEDCoupling::medloader ALIAS medloader)

if(NOT DEFINED MSVC)
    target_compile_options(medloader PRIVATE -Wsign-compare -Wconversion)
endif()

# if(WIN32)
#   if(MEDCOUPLING_MEDLOADER_USE_XDR)
#     add_definitions(-DNOMINMAX)
#     set(MEDLOADER_XDR_INCLUDE_DIRS ${XDR_INCLUDE_DIRS})
#     set(MEDLOADER_XDR_LIBRARIES ${XDR_LIBRARIES})
#   endif(MEDCOUPLING_MEDLOADER_USE_XDR)
# else(WIN32)
#   set(MEDLOADER_XDR_INCLUDE_DIRS ${XDR_INCLUDE_DIRS})
#   set(MEDLOADER_XDR_LIBRARIES ${XDR_LIBRARIES})
# endif(WIN32)

target_compile_definitions(medloader PUBLIC MED_INT_IS_LONG)

target_include_directories(
  medloader ${warning_guard}
  PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/medfile_api>"
)

set_target_properties(
  medloader PROPERTIES
  COMPILE_FLAGS ""
  INSTALL_RPATH ${MEDCOUPLING_INSTALL_LIBS_DEPS_RPATH}
)

# TODO: remove this dependency by splitting medloader into
# compound_mesh (not dependant and stays there) and
# medreader/writer (moves to medcoupling_io)

# find_package(HDF5 REQUIRED)
# find_package(medfile REQUIRED)
# 
# if(NOT TARGET HDF5::HDF5)
#   if(BUILD_SHARED_LIBS)
#     add_library(HDF5::HDF5 ALIAS hdf5::hdf5-shared)
#   else()
#     add_library(HDF5::HDF5 ALIAS hdf5::hdf5-static)
#   endif()
# endif()

target_link_libraries(medloader PUBLIC medcouplingcpp medfile::medfile HDF5::HDF5)

install(
  TARGETS medloader
  # RUNTIME_DEPENDENCIES DIRECTORIES ${MEDCOUPLING_INSTALL_LIBS}
  EXPORT ${CMAKE_PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS}
  INCLUDES DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

file(GLOB medloader_HEADERS "${PROJECT_SOURCE_DIR}/include/medfile_api/*")
install(FILES ${medloader_HEADERS} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
