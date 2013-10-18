# Copyright (C) 2012-2013  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (CEA/DEN)

ADD_DEFINITIONS(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS})

IF(SALOME_MED_ENABLE_PYTHON)
  ADD_SUBDIRECTORY(Swig)
ENDIF(SALOME_MED_ENABLE_PYTHON)

IF(SALOME_BUILD_TESTS)
  ADD_SUBDIRECTORY(Test)
ENDIF(SALOME_BUILD_TESTS)

INCLUDE_DIRECTORIES(
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(medloader_SOURCES
  MEDLoader.cxx
  MEDLoaderBase.cxx
  MEDFileUtilities.cxx
  MEDFileMesh.cxx
  MEDFileMeshElt.cxx
  MEDFileBasis.cxx
  MEDFileMeshLL.cxx
  MEDFileField.cxx
  MEDFileParameter.cxx
  MEDFileData.cxx
  MEDFileFieldOverView.cxx
  MEDFileMeshReadSelector.cxx
  SauvMedConvertor.cxx
  SauvReader.cxx
  SauvWriter.cxx
  )

ADD_LIBRARY(medloader SHARED ${medloader_SOURCES})
SET_TARGET_PROPERTIES(medloader PROPERTIES COMPILE_FLAGS "${XDR_DEFINITIONS}")
TARGET_LINK_LIBRARIES(medloader medcoupling ${MEDFILE_C_LIBRARIES} ${HDF5_LIBRARIES})
INSTALL(TARGETS medloader EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${SALOME_INSTALL_LIBS})

FILE(GLOB medloader_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${medloader_HEADERS_HXX} DESTINATION ${SALOME_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
SET(medloader_HEADERS_HXX PARENT_SCOPE)
