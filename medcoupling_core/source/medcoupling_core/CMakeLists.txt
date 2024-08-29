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

# KILLER WARNINGS TO HELP INT64 MANAGEMENT ADD_DEFINITIONS(-Wextra -Wsynth
# -Wno-unused-parameter -pedantic -Wall -Wshadow  -Wpointer-arith -Wcast-qual
# -Wcast-align -Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor
# -Woverloaded-virtual -Wreturn-type  -Wno-endif-labels  -Wsign-compare
# -Wmissing-format-attribute -Wno-multichar -Wno-deprecated-declarations
# -Wpacked -Wredundant-decls -Wlong-long -Wdisabled-optimization
# -Wunknown-pragmas -Wundef  -Wreorder)

set(MEDCOUPLING_SOURCES
    amr/MEDCouplingCartesianAMRMesh.cxx
    amr/MEDCouplingAMRAttribute.cxx
    array/MEDCouplingMemArray.cxx
    array/MEDCouplingMemArrayFloat.cxx
    array/MEDCouplingMemArrayChar.cxx
    array/MEDCouplingMap.cxx
    array/MEDCouplingTraits.cxx
    array/MEDCouplingTimeLabel.cxx
    array/MEDCouplingRefCountObject.cxx
    array/MEDCouplingMatrix.cxx
    array/MEDCouplingSkyLineArray.cxx
    array/MEDCouplingPartDefinition.cxx
    mesh/MEDCouplingUMesh.cxx
    mesh/MEDCouplingUMesh_internal.cxx
    mesh/MEDCouplingUMesh_intersection.cxx
    mesh/MEDCoupling1GTUMesh.cxx
    mesh/MEDCouplingCMesh.cxx
    mesh/MEDCouplingIMesh.cxx
    mesh/MEDCouplingCurveLinearMesh.cxx
    mesh/MEDCouplingStructuredMesh.cxx
    mesh/MEDCouplingPointSet.cxx
    mesh/MEDCouplingMappedExtrudedMesh.cxx
    mesh/MEDCouplingMesh.cxx
    field/MEDCouplingField.cxx
    field/MEDCouplingFieldFloat.cxx
    field/MEDCouplingFieldDouble.cxx
    field/MEDCouplingFieldInt32.cxx
    field/MEDCouplingFieldInt64.cxx
    field/MEDCouplingTimeDiscretization.cxx
    field/MEDCouplingFieldDiscretization.cxx
    field/MEDCouplingFieldDiscretizationOnNodesFE.cxx
    field/MEDCouplingFieldTemplate.cxx
    field/MEDCouplingGaussLocalization.cxx
    field/MEDCouplingNatureOfField.cxx
    field/MEDCouplingMultiFields.cxx
    field/MEDCouplingDefinitionTime.cxx
    field/MEDCouplingFieldOverTime.cxx
    field/MEDCouplingVoronoi.cxx)

add_library(medcouplingcpp ${MEDCOUPLING_SOURCES})
add_library(MEDCoupling::medcoupling ALIAS medcouplingcpp)

if(NOT DEFINED MSVC)
  target_compile_options(medcouplingcpp PRIVATE -Wsign-compare -Wconversion)
endif()

target_include_directories(
  medcouplingcpp
  PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/medcoupling_core>"
)

set_target_properties(
  medcouplingcpp PROPERTIES
  INSTALL_RPATH ${MEDCOUPLING_INSTALL_LIBS_DEPS_RPATH}
  OUTPUT_NAME "medcoupling")

if(MEDCOUPLING_USE_64BIT_IDS)
  target_compile_definitions(medcouplingcpp PUBLIC "MEDCOULPING_USE_64BIT_IDS")
endif()

target_link_libraries(medcouplingcpp PUBLIC interpkernel)

install(
  TARGETS medcouplingcpp
  EXPORT ${CMAKE_PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS}
  INCLUDES DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# TODO: remove the remapper from medcoupling_data to put it in medcoupling_tools
# TODO: Needs to seperate para_medmem in two parts, para_remapper and para_polymesh

############
# Remapper #
############

set(MEDCOUPLINGREMAPPER_SOURCES remapper/MEDCouplingRemapper.cxx)
add_library(medcouplingremapper ${MEDCOUPLINGREMAPPER_SOURCES})
target_link_libraries(medcouplingremapper PUBLIC medcouplingcpp)

set_target_properties(
  medcouplingremapper PROPERTIES
  INSTALL_RPATH ${MEDCOUPLING_INSTALL_LIBS_DEPS_RPATH})

if(MEDCOUPLING_USE_64BIT_IDS)
  target_compile_definitions(medcouplingremapper PUBLIC "MEDCOULPING_USE_64BIT_IDS")
endif()

install(
  TARGETS medcouplingremapper
  EXPORT ${CMAKE_PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS}
  INCLUDES DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

file(GLOB medcoupling_HEADERS "${PROJECT_SOURCE_DIR}/include/polymesh/*")
install(FILES ${medcoupling_HEADERS} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
