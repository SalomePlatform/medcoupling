# Copyright (C) 2012-2016  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (CEA/DEN)

ADD_DEFINITIONS(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${XDR_DEFINITIONS})

IF(MEDCOUPLING_ENABLE_PYTHON)
  ADD_SUBDIRECTORY(Swig)
ENDIF(MEDCOUPLING_ENABLE_PYTHON)

IF(MEDCOUPLING_BUILD_TESTS)
  ADD_SUBDIRECTORY(Test)
ENDIF(MEDCOUPLING_BUILD_TESTS)

SET(MEDLOADER_XDR_INCLUDE_DIRS)
SET(MEDLOADER_XDR_LIBRARIES)

IF(WIN32)
  IF(MEDCOUPLING_MEDLOADER_USE_XDR)
    ADD_DEFINITIONS(-DNOMINMAX)
    SET(MEDLOADER_XDR_INCLUDE_DIRS ${XDR_INCLUDE_DIRS})
    SET(MEDLOADER_XDR_LIBRARIES ${XDR_LIBRARIES})
  ENDIF(MEDCOUPLING_MEDLOADER_USE_XDR)
ELSE(WIN32)
  SET(MEDLOADER_XDR_INCLUDE_DIRS ${XDR_INCLUDE_DIRS})
ENDIF(WIN32)

INCLUDE_DIRECTORIES(
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MEDLOADER_XDR_INCLUDE_DIRS}
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
  MEDLoaderTraits.cxx
  MEDFileUtilities.cxx
  MEDFileMesh.cxx
  MEDFileMeshElt.cxx
  MEDFileBasis.cxx
  MEDFileMeshLL.cxx
  MEDFileField.cxx
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
  )

ADD_LIBRARY(medloader ${medloader_SOURCES})
SET_TARGET_PROPERTIES(medloader PROPERTIES COMPILE_FLAGS "")
TARGET_LINK_LIBRARIES(medloader medcoupling ${MEDFILE_C_LIBRARIES} ${HDF5_LIBRARIES} ${MEDLOADER_XDR_LIBRARIES})
INSTALL(TARGETS medloader EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

FILE(GLOB medloader_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
FILE(GLOB medloader_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${medloader_HEADERS_HXX} ${medloader_HEADERS_TXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
SET(medloader_HEADERS_HXX PARENT_SCOPE)
