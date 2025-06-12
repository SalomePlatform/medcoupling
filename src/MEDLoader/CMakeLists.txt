# Copyright (C) 2012-2025  CEA, EDF
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

add_definitions(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${XDR_DEFINITIONS})

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

if(MEDCOUPLING_ENABLE_PYTHON)
  add_subdirectory(Swig)
endif(MEDCOUPLING_ENABLE_PYTHON)

if(MEDCOUPLING_BUILD_TESTS)
  add_subdirectory(Test)
endif(MEDCOUPLING_BUILD_TESTS)

set(MEDLOADER_XDR_INCLUDE_DIRS)
set(MEDLOADER_XDR_LIBRARIES)

if(WIN32)
  if(MEDCOUPLING_MEDLOADER_USE_XDR)
    add_definitions(-DNOMINMAX)
    set(MEDLOADER_XDR_INCLUDE_DIRS ${XDR_INCLUDE_DIRS})
    set(MEDLOADER_XDR_LIBRARIES ${XDR_LIBRARIES})
  endif(MEDCOUPLING_MEDLOADER_USE_XDR)
else(WIN32)
  set(MEDLOADER_XDR_INCLUDE_DIRS ${XDR_INCLUDE_DIRS})
  set(MEDLOADER_XDR_LIBRARIES ${XDR_LIBRARIES})
endif(WIN32)

include_directories(
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MEDLOADER_XDR_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints)

set(medloader_SOURCES
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

add_library(medloader ${medloader_SOURCES})
set_target_properties(medloader PROPERTIES COMPILE_FLAGS "")
target_link_libraries(medloader medcouplingcpp ${MEDFILE_C_LIBRARIES}
                      ${HDF5_LIBRARIES} ${MEDLOADER_XDR_LIBRARIES})
install(
  TARGETS medloader
  EXPORT ${PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

file(GLOB medloader_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
file(GLOB medloader_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
install(FILES ${medloader_HEADERS_HXX} ${medloader_HEADERS_TXX}
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
set(medloader_HEADERS_HXX PARENT_SCOPE)
