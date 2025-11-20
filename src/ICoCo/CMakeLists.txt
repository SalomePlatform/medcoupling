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
# Author : Adrien Bruneton (CEA/DES)

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

if(MEDCOUPLING_ENABLE_PYTHON)
  add_subdirectory(Swig)
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases)

set(icoco_SOURCES ICoCoField.cpp # [ABN] Yes, .cpp, this is imposed by ICoCo.
                  ICoCoMEDDoubleField.cxx ICoCoMEDIntField.cxx)

add_library(medicoco ${icoco_SOURCES})
set_target_properties(medicoco PROPERTIES OUTPUT_NAME "medicoco")
target_link_libraries(medicoco medcouplingcpp)
install(
  TARGETS medicoco
  EXPORT ${PROJECT_NAME}TargetGroup
  ARCHIVE DESTINATION ${MEDCOUPLING_INSTALL_LIBS}
  RUNTIME DESTINATION ${MEDCOUPLING_INSTALL_BINS}
  LIBRARY DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

file(GLOB icoco_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx"
     "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${icoco_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
