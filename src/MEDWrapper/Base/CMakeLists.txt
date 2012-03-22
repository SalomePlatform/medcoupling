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
  ${KERNEL_ROOT_DIR}/include/salome
  )

SET(MEDWrapperBase_SOURCES
  MED_Structures.cxx
  MED_Wrapper.cxx
  MED_Algorithm.cxx
  MED_GaussUtils.cxx
  MED_CoordUtils.cxx
  MED_Utilities.cxx
  MED_GaussDef.cxx
  )

ADD_LIBRARY(MEDWrapperBase SHARED ${MEDWrapperBase_SOURCES})
SET_TARGET_PROPERTIES(MEDWrapperBase PROPERTIES COMPILE_FLAGS "${HDF5_FLAGS} ${BOOST_FLAGS}")
TARGET_LINK_LIBRARIES(MEDWrapperBase ${BOOST_LIB_THREAD} ${BOOST_LIB_DATE_TIME})
INSTALL(TARGETS MEDWrapperBase DESTINATION ${MED_salomelib_LIBS})

FILE(GLOB MEDWrapperBase_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${MEDWrapperBase_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
