# Copyright (C) 2012-2023  CEA/DEN, EDF R&D
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
# Author : Adrien Bruneton (CEA/DES)

IF (NOT DEFINED MSVC)
  ADD_DEFINITIONS(-Wsign-compare -Wconversion)
ENDIF()

IF(MEDCOUPLING_ENABLE_PYTHON)
  ADD_SUBDIRECTORY(Swig)
ENDIF()


INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  )

SET(icoco_SOURCES
  ICoCoField.cpp      # [ABN] Yes, .cpp, this is imposed by ICoCo.
  ICoCoMEDDoubleField.cxx
  ICoCoMEDIntField.cxx
)

ADD_LIBRARY(medicoco ${icoco_SOURCES})
SET_TARGET_PROPERTIES(medicoco PROPERTIES OUTPUT_NAME "medicoco")
TARGET_LINK_LIBRARIES(medicoco medcouplingcpp)
INSTALL(TARGETS medicoco EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

FILE(GLOB icoco_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx"
                                  "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${icoco_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
