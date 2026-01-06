# Copyright (C) 2024-2026  CEA, EDF
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

add_definitions(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS})

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

if(MEDCOUPLING_ENABLE_PYTHON)
  add_subdirectory(Swig)
endif(MEDCOUPLING_ENABLE_PYTHON)

if(MEDCOUPLING_BUILD_TESTS)
  add_subdirectory(Test)
endif(MEDCOUPLING_BUILD_TESTS)

include_directories(
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${LAPACKE_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases)

set(shaperecogn_SOURCES
    MathOps.cxx
    Nodes.cxx
    NodesBuilder.cxx
    Areas.cxx
    AreasBuilder.cxx
    ShapeRecognMeshBuilder.cxx
    ShapeRecognMesh.cxx
    PrimitiveType.cxx)

add_library(shaperecogn ${shaperecogn_SOURCES})
set_target_properties(shaperecogn PROPERTIES COMPILE_FLAGS "")
target_link_libraries(shaperecogn medcouplingcpp ${LAPACK_LIBRARIES})
install(
  TARGETS shaperecogn
  EXPORT ${PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

file(GLOB shaperecogn_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")

# To allow usage as SWIG dependencies:
set(shaperecogn_HEADERS_HXX PARENT_SCOPE)
