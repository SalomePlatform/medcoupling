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
# Author : Adrien Bruneton (CEA/DES)

set(ICOCO_SOURCES ICoCoField.cpp
                  ICoCoMEDDoubleField.cxx
                  ICoCoMEDIntField.cxx)

add_library(medicoco ${ICOCO_SOURCES})
add_library(MEDCoupling::medicoco ALIAS medicoco)

target_include_directories(
  medicoco
  PUBLIC "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/icoco>")

target_link_libraries(medicoco PUBLIC medcouplingcpp)

set_target_properties(
  medicoco PROPERTIES
  INSTALL_RPATH ${MEDCOUPLING_INSTALL_LIBS_DEPS_RPATH}
)

if(NOT DEFINED MSVC)
  target_compile_options(medicoco PRIVATE -Wsign-compare -Wconversion)
endif()

install(
  TARGETS medicoco
  EXPORT ${CMAKE_PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS}
  INCLUDES DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

file(GLOB ICOCO_HEADERS "${PROJECT_SOURCE_DIR}/include/icoco/*")
install(FILES ${ICOCO_HEADERS} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
