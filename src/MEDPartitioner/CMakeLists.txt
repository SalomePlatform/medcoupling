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

add_definitions(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS}
                ${LIBXML2_DEFINITIONS})

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

include_directories(
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
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints)

if(MEDCOUPLING_PARTITIONER_METIS)
  add_definitions(${METIS_DEFINITIONS})
  if(MEDCOUPLING_METIS_V5)
    add_definitions("-DMED_ENABLE_METIS_V5")
  endif(MEDCOUPLING_METIS_V5)
  include_directories(${METIS_INCLUDE_DIRS})
endif(MEDCOUPLING_PARTITIONER_METIS)

#
# VERY IMPORTANT PUT METIS BEFORE SCOTCH because metis.h is also in SCOTCH
# install dir !!!
#
if(MEDCOUPLING_PARTITIONER_PTSCOTCH)
  add_definitions(${PTSCOTCH_DEFINITIONS})
  include_directories(${PTSCOTCH_INCLUDE_DIRS})
endif(MEDCOUPLING_PARTITIONER_PTSCOTCH)

if(MEDCOUPLING_PARTITIONER_SCOTCH)
  add_definitions(${SCOTCH_DEFINITIONS})
  include_directories(${SCOTCH_INCLUDE_DIRS})
endif(MEDCOUPLING_PARTITIONER_SCOTCH)

if(MEDCOUPLING_PARTITIONER_PARMETIS)
  add_definitions(${PARMETIS_DEFINITIONS})
  include_directories(${PARMETIS_INCLUDE_DIRS})
endif(MEDCOUPLING_PARTITIONER_PARMETIS)

if(MEDCOUPLING_USE_MPI)
  add_definitions(${MPI_DEFINITIONS})
  include_directories(${MPI_INCLUDE_DIRS})
endif(MEDCOUPLING_USE_MPI)

if(MEDCOUPLING_BUILD_TESTS)
  add_subdirectory(Test)
endif(MEDCOUPLING_BUILD_TESTS)

set(medpartitionercpp_HEADERS_HXX
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
    MEDPARTITIONER_MEDPartitioner.hxx)

set(medpartitionercpp_SOURCES
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
    MEDPARTITIONER_MEDPartitioner.cxx)

set(medpartitionercpp_LDFLAGS
    ${MEDFILE_C_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${STDLIB}
    ${LIBXML2_LIBRARIES}
    interpkernel
    medcouplingcpp
    medloader
    -lpthread)
if(MEDCOUPLING_PARTITIONER_PARMETIS)
  set(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX}
                                    MEDPARTITIONER_ParMetisGraph.hxx)
  set(medpartitionercpp_SOURCES
      ${medpartitionercpp_SOURCES} MEDPARTITIONER_ParMetisGraph.cxx
      MEDPARTITIONER_MetisGraph.cxx)
  set(medpartitionercpp_DEFINITIONS
      "${medpartitionercpp_DEFINITIONS} ${PARMETIS_DEFINITIONS}")
  set(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS}
                                ${PARMETIS_LIBRARIES})
endif(MEDCOUPLING_PARTITIONER_PARMETIS)
if(MEDCOUPLING_PARTITIONER_METIS)
  set(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX}
                                    MEDPARTITIONER_MetisGraph.hxx)
  set(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES}
                                MEDPARTITIONER_MetisGraph.cxx)
  set(medpartitionercpp_DEFINITIONS
      "${medpartitionercpp_DEFINITIONS} ${METIS_DEFINITIONS}")
  set(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${METIS_LIBRARIES})
endif(MEDCOUPLING_PARTITIONER_METIS)
if(MEDCOUPLING_PARTITIONER_SCOTCH)
  set(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX}
                                    MEDPARTITIONER_ScotchGraph.hxx)
  set(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES}
                                MEDPARTITIONER_ScotchGraph.cxx)
  set(medpartitionercpp_DEFINITIONS
      "${medpartitionercpp_DEFINITIONS} ${SCOTCH_DEFINITIONS}")
  set(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS}
                                ${SCOTCH_LIBRARIES})
endif(MEDCOUPLING_PARTITIONER_SCOTCH)
if(MEDCOUPLING_PARTITIONER_PTSCOTCH)
  set(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX}
                                    MEDPARTITIONER_PTScotchGraph.hxx)
  set(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES}
                                MEDPARTITIONER_PTScotchGraph.cxx)
  set(medpartitionercpp_DEFINITIONS
      "${medpartitionercpp_DEFINITIONS} ${PTSCOTCH_DEFINITIONS}")
  set(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS}
                                ${PTSCOTCH_LIBRARIES})
endif(MEDCOUPLING_PARTITIONER_PTSCOTCH)

if(${MEDCOUPLING_USE_MPI})
  set(medpartitionercpp_SOURCES
      ${medpartitionercpp_SOURCES} MEDPARTITIONER_UtilsPara.cxx
      MEDPARTITIONER_JointFinder.cxx)
  add_executable(medpartitioner_para medpartitioner_para.cxx)
  set(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${MPI_LIBRARIES})
  set_target_properties(
    medpartitioner_para PROPERTIES COMPILE_FLAGS
                                   "${medpartitionercpp_DEFINITIONS}")
  target_link_libraries(medpartitioner_para medpartitionercpp
                        ${medpartitionercpp_LDFLAGS})
  install(TARGETS medpartitioner_para DESTINATION ${MEDCOUPLING_INSTALL_BINS})
endif(${MEDCOUPLING_USE_MPI})

add_definitions(${medpartitionercpp_DEFINITIONS})

add_library(medpartitionercpp ${medpartitionercpp_SOURCES})
target_link_libraries(medpartitionercpp ${medpartitionercpp_LDFLAGS}
                      ${PLATFORM_LIBS} ${PTHREAD_LIBS})
install(TARGETS medpartitionercpp DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

add_executable(medpartitioner medpartitioner.cxx)
target_link_libraries(medpartitioner medpartitionercpp
                      ${medpartitionercpp_LDFLAGS})
install(TARGETS medpartitioner DESTINATION ${MEDCOUPLING_INSTALL_BINS})

install(FILES ${medpartitionercpp_HEADERS_HXX}
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
