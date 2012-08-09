# Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

ADD_SUBDIRECTORY(Test)

INCLUDE_DIRECTORIES(
  ${PARMETIS_INCLUDE_DIRS}
  ${METIS_INCLUDE_DIRS}
  ${SCOTCH_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIRS}
  ${METIS_INCLUDE_DIRS}
  ${LIBXML_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MPI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(medsplittercpp_HEADERS_HXX
  MEDSPLITTER_Topology.hxx
  MEDSPLITTER_MESHCollection.hxx
  MEDSPLITTER_MESHCollection.H
  MEDSPLITTER_MESHCollectionDriver.H
  MEDSPLITTER_MESHCollectionMedXMLDriver.H
  MEDSPLITTER_MESHCollectionMedAsciiDriver.H
  MEDSPLITTER_MESHCollectionDriver.hxx
  MEDSPLITTER_MESHCollectionMedXMLDriver.hxx
  MEDSPLITTER_MESHCollectionMedAsciiDriver.hxx
  MEDSPLITTER_ParallelTopology.hxx
  MEDSPLITTER_FaceModel.hxx
  MEDSPLITTER_Graph.hxx
  MEDSPLITTER_UserGraph.hxx
  MEDSPLITTER_API.hxx
  MEDSPLITTER_SequentialTopology.hxx
  MEDSPLITTER_utils.hxx
  MEDSPLITTER.hxx
  MEDSPLITTER_ParaDomainSelector.hxx
  MEDSPLITTER_MeshSendReceive.hxx
  MEDSPLITTER_JointExchangeData.hxx
  )

SET(medsplittercpp_SOURCES
  MEDSPLITTER_MESHCollection.cxx
  MEDSPLITTER_MESHCollectionDriver.cxx
  MEDSPLITTER_MESHCollectionMedXMLDriver.cxx
  MEDSPLITTER_MESHCollectionMedAsciiDriver.cxx
  MEDSPLITTER_ParallelTopology.cxx
  MEDSPLITTER_Graph.cxx
  MEDSPLITTER_UserGraph.cxx
  MEDSPLITTER_API.cxx
  MEDSPLITTER_ParaDomainSelector.cxx
  MEDSPLITTER_MeshSendReceive.cxx
  MEDSPLITTER_JointExchangeData.cxx
  )

SET(medsplittercpp_DEFINITIONS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${BOOST_DEFINITIONS} ${LIBXML_DEFINITIONS}")
SET(medsplittercpp_LDFLAGS medmem ${LIBXML_LIBS})

IF(MED_ENABLE_SCOTCH)
  SET(medsplittercpp_SOURCES ${medsplittercpp_SOURCES} MEDSPLITTER_SCOTCHGraph.cxx)
  SET(medsplittercpp_HEADERS_HXX ${medsplittercpp_HEADERS_HXX} MEDSPLITTER_SCOTCHGraph.hxx)
  SET(medsplittercpp_DEFINITIONS "${medsplittercpp_DEFINITIONS} ${SCOTCH_DEFINITIONS}")
  SET(medsplittercpp_LDFLAGS ${medsplittercpp_LDFLAGS} ${SCOTCH_LIBS})
ENDIF(MED_ENABLE_SCOTCH)

IF(MPI_IS_OK)
  IF(MED_ENABLE_PARMETIS)
    SET(medsplittercpp_SOURCES ${medsplittercpp_SOURCES} MEDSPLITTER_METISGraph.cxx)
    SET(medsplittercpp_HEADERS_HXX ${medsplittercpp_HEADERS_HXX} MEDSPLITTER_METISGraph.hxx)
    SET(medsplittercpp_DEFINITIONS "${medsplittercpp_DEFINITIONS} ${MPI_DEFINITIONS} ${PARMETIS_DEFINITIONS}")
    SET(medsplittercpp_LDFLAGS ${medsplittercpp_LDFLAGS} ${PARMETIS_LIBS} ${MPI_LIBS})
  ENDIF(MED_ENABLE_PARMETIS)
ELSE(MPI_IS_OK)
  IF(MED_ENABLE_METIS)
    SET(medsplittercpp_SOURCES ${medsplittercpp_SOURCES} MEDSPLITTER_METISGraph.cxx)
    SET(medsplittercpp_HEADERS_HXX ${medsplittercpp_HEADERS_HXX} MEDSPLITTER_METISGraph.hxx)
    SET(medsplittercpp_DEFINITIONS "${medsplittercpp_DEFINITIONS} ${METIS_DEFINITIONS}")
    SET(medsplittercpp_LDFLAGS ${medsplittercpp_LDFLAGS} ${METIS_LIBS})
  ENDIF(MED_ENABLE_METIS)
ENDIF(MPI_IS_OK)


IF(MPI_IS_OK)
  ADD_EXECUTABLE(medsplitter_para medsplitter_para.cxx)
  SET_TARGET_PROPERTIES(medsplitter_para PROPERTIES COMPILE_FLAGS "${medsplittercpp_DEFINITIONS}")
  TARGET_LINK_LIBRARIES(medsplitter_para medsplittercpp)
  INSTALL(TARGETS medsplitter_para DESTINATION ${MED_salomebin_BINS})
ELSE(MPI_IS_OK)
  ADD_EXECUTABLE(medsplitter medsplitter.cxx)
  SET_TARGET_PROPERTIES(medsplitter PROPERTIES COMPILE_FLAGS "${medsplittercpp_DEFINITIONS}")
  TARGET_LINK_LIBRARIES(medsplitter medsplittercpp)
  INSTALL(TARGETS medsplitter DESTINATION ${MED_salomebin_BINS})
ENDIF(MPI_IS_OK)

ADD_LIBRARY(medsplittercpp SHARED ${medsplittercpp_SOURCES})
SET_TARGET_PROPERTIES(medsplittercpp PROPERTIES COMPILE_FLAGS "${medsplittercpp_DEFINITIONS}")
TARGET_LINK_LIBRARIES(medsplittercpp ${medsplittercpp_LDFLAGS})
INSTALL(TARGETS medsplittercpp DESTINATION ${MED_salomelib_LIBS})

INSTALL(FILES ${medsplittercpp_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
