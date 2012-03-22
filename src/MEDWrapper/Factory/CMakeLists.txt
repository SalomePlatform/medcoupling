# Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

INCLUDE_DIRECTORIES(
  ${HDF5_INCLUDES_DIR}
  ${BOOST_INCLUDES_DIR}
  ${MED3_INCLUDES_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../Base
  ${CMAKE_CURRENT_SOURCE_DIR}/../V2_1/Wrapper
  ${CMAKE_CURRENT_SOURCE_DIR}/../V2_2
  )

IF(MED_ENABLE_KERNEL)
  INCLUDE_DIRECTORIES(${KERNEL_ROOT_DIR}/include/salome)
ELSE(MED_ENABLE_KERNEL)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/adm_local_without_kernel/unix)
ENDIF(MED_ENABLE_KERNEL)  

SET(MEDWrapper_SOURCES
  MED_Factory.cxx
  )

SET(mprint_version_SOURCES
  mprint_version.cxx
  )

ADD_LIBRARY(MEDWrapper SHARED ${MEDWrapper_SOURCES})
SET_TARGET_PROPERTIES(MEDWrapper PROPERTIES COMPILE_FLAGS "${BOOST_FLAGS} ${MED3_FLAGS}")
TARGET_LINK_LIBRARIES(MEDWrapper MEDWrapperBase MEDWrapper_V2_1 MEDWrapper_V2_2)
INSTALL(TARGETS MEDWrapper DESTINATION ${MED_salomelib_LIBS})

ADD_EXECUTABLE(mprint_version ${mprint_version_SOURCES})
SET_TARGET_PROPERTIES(mprint_version PROPERTIES COMPILE_FLAGS "${BOOST_FLAGS} ${MED3_FLAGS}")
TARGET_LINK_LIBRARIES(mprint_version MEDWrapper)

FILE(GLOB MEDWrapper_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${MEDWrapper_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL(TARGETS mprint_version DESTINATION ${MED_salomebin_BINS})