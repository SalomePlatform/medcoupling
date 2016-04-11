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

ADD_DEFINITIONS(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${LIBXML2_DEFINITIONS})

INCLUDE_DIRECTORIES(
  ${LIBXML2_INCLUDE_DIR}
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

IF(MEDCOUPLING_PARTITIONER_METIS)
  ADD_DEFINITIONS(${METIS_DEFINITIONS})
  IF(MEDCOUPLING_METIS_V5)
    ADD_DEFINITIONS("-DMED_ENABLE_METIS_V5")
  ENDIF(MEDCOUPLING_METIS_V5)
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIRS})
ENDIF(MEDCOUPLING_PARTITIONER_METIS)


########
# VERY IMPORTANT PUT METIS BEFORE SCOTCH because
# metis.h is also in SCOTCH install dir !!!
########
IF(MEDCOUPLING_PARTITIONER_SCOTCH)
  ADD_DEFINITIONS(${SCOTCH_DEFINITIONS})
  INCLUDE_DIRECTORIES(${SCOTCH_INCLUDE_DIRS})
ENDIF(MEDCOUPLING_PARTITIONER_SCOTCH)

IF(MEDCOUPLING_PARTITIONER_PARMETIS)
  ADD_DEFINITIONS(${PARMETIS_DEFINITIONS})
  INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIRS})
ENDIF(MEDCOUPLING_PARTITIONER_PARMETIS)

IF(MEDCOUPLING_USE_MPI)
  ADD_DEFINITIONS(${MPI_DEFINITIONS})
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_DIRS})
ENDIF(MEDCOUPLING_USE_MPI)

IF(MEDCOUPLING_BUILD_TESTS)
  ADD_SUBDIRECTORY(Test)
ENDIF(MEDCOUPLING_BUILD_TESTS)

SET(medpartitionercpp_HEADERS_HXX
  MEDPARTITIONER_MeshCollection.hxx
  MEDPARTITIONER_MeshCollectionDriver.hxx
  MEDPARTITIONER_MeshCollectionMedXmlDriver.hxx
  MEDPARTITIONER_MeshCollectionMedAsciiDriver.hxx
  MEDPARTITIONER_ParallelTopology.hxx
  MEDPARTITIONER_JointFinder.hxx
  MEDPARTITIONER_Graph.hxx
  MEDPARTITIONER_UserGraph.hxx
  MEDPARTITIONER_Utils.hxx
  MEDPARTITIONER.hxx
  MEDPARTITIONER_ParaDomainSelector.hxx
  MEDPARTITIONER_ConnectZone.hxx
  MEDPARTITIONER_Topology.hxx
  MEDPARTITIONER_MEDPartitioner.hxx
  )

SET(medpartitionercpp_SOURCES
  MEDPARTITIONER_MeshCollection.cxx
  MEDPARTITIONER_MeshCollectionDriver.cxx
  MEDPARTITIONER_MeshCollectionMedXmlDriver.cxx
  MEDPARTITIONER_MeshCollectionMedAsciiDriver.cxx
  MEDPARTITIONER_ParallelTopology.cxx
  MEDPARTITIONER_Graph.cxx
  MEDPARTITIONER_UserGraph.cxx
  MEDPARTITIONER_Utils.cxx
  MEDPARTITIONER_ParaDomainSelector.cxx
  MEDPARTITIONER_ConnectZone.cxx
  MEDPARTITIONER_metis.c
  MEDPARTITIONER_MEDPartitioner.cxx
  )

SET(medpartitionercpp_LDFLAGS
  ${MEDFILE_C_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${STDLIB}
  ${LIBXML2_LIBRARIES}
  interpkernel
  medcoupling
  medloader
)
IF(MEDCOUPLING_PARTITIONER_PARMETIS)
  SET(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX} MEDPARTITIONER_ParMetisGraph.hxx)
  SET(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES} MEDPARTITIONER_ParMetisGraph.cxx MEDPARTITIONER_MetisGraph.cxx)
  SET(medpartitionercpp_DEFINITIONS "${medpartitionercpp_DEFINITIONS} ${PARMETIS_DEFINITIONS}")
  SET(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${PARMETIS_LIBRARIES})
ENDIF(MEDCOUPLING_PARTITIONER_PARMETIS)
IF(MEDCOUPLING_PARTITIONER_METIS)
  SET(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX} MEDPARTITIONER_MetisGraph.hxx)
  SET(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES} MEDPARTITIONER_MetisGraph.cxx)
  SET(medpartitionercpp_DEFINITIONS "${medpartitionercpp_DEFINITIONS} ${METIS_DEFINITIONS}")
  SET(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${METIS_LIBRARIES})
ENDIF(MEDCOUPLING_PARTITIONER_METIS)
IF(MEDCOUPLING_PARTITIONER_SCOTCH)
  SET(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX} MEDPARTITIONER_ScotchGraph.hxx)
  SET(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES} MEDPARTITIONER_ScotchGraph.cxx)
  SET(medpartitionercpp_DEFINITIONS "${medpartitionercpp_DEFINITIONS} ${SCOTCH_DEFINITIONS}")
  SET(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${SCOTCH_LIBRARIES})
ENDIF(MEDCOUPLING_PARTITIONER_SCOTCH)

IF(${MEDCOUPLING_USE_MPI})
  SET(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES} MEDPARTITIONER_UtilsPara.cxx MEDPARTITIONER_JointFinder.cxx)
  ADD_EXECUTABLE(medpartitioner_para medpartitioner_para.cxx)
  SET(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${MPI_LIBRARIES})
  SET_TARGET_PROPERTIES(medpartitioner_para PROPERTIES COMPILE_FLAGS "${medpartitionercpp_DEFINITIONS}")
  TARGET_LINK_LIBRARIES(medpartitioner_para medpartitionercpp ${medpartitionercpp_LDFLAGS})
  INSTALL(TARGETS medpartitioner_para DESTINATION ${MEDCOUPLING_INSTALL_BINS})
ENDIF(${MEDCOUPLING_USE_MPI})

ADD_DEFINITIONS(${medpartitionercpp_DEFINITIONS})

ADD_LIBRARY(medpartitionercpp ${medpartitionercpp_SOURCES})
TARGET_LINK_LIBRARIES(medpartitionercpp ${medpartitionercpp_LDFLAGS} ${PLATFORM_LIBS} ${PTHREAD_LIBS})
INSTALL(TARGETS medpartitionercpp DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

ADD_EXECUTABLE(medpartitioner medpartitioner.cxx)
TARGET_LINK_LIBRARIES(medpartitioner medpartitionercpp ${medpartitionercpp_LDFLAGS})
INSTALL(TARGETS medpartitioner DESTINATION ${MEDCOUPLING_INSTALL_BINS})

INSTALL(FILES ${medpartitionercpp_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
