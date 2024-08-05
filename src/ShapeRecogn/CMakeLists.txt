# Copyright (C) 2012-2024  CEA, EDF
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

ADD_DEFINITIONS(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS})

IF (NOT DEFINED MSVC)
  ADD_DEFINITIONS(-Wsign-compare -Wconversion)
ENDIF()

IF(MEDCOUPLING_ENABLE_PYTHON)
  ADD_SUBDIRECTORY(Swig)
ENDIF(MEDCOUPLING_ENABLE_PYTHON)

IF(MEDCOUPLING_BUILD_TESTS)
  ADD_SUBDIRECTORY(Test)
ENDIF(MEDCOUPLING_BUILD_TESTS)

INCLUDE_DIRECTORIES(
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${LAPACKE_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
)

SET(shaperecogn_SOURCES
  MathOps.cxx
  Nodes.cxx
  NodesBuilder.cxx
  Areas.cxx
  AreasBuilder.cxx
  ShapeRecognMeshBuilder.cxx
  ShapeRecognMesh.cxx
)

ADD_LIBRARY(shaperecogn ${shaperecogn_SOURCES})
SET_TARGET_PROPERTIES(shaperecogn PROPERTIES COMPILE_FLAGS "")
TARGET_LINK_LIBRARIES(shaperecogn medcouplingcpp medloader ${MEDFILE_C_LIBRARIES} ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
INSTALL(TARGETS shaperecogn EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

FILE(GLOB shaperecogn_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${shaperecogn_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
SET(shaperecogn_HEADERS_HXX PARENT_SCOPE)
