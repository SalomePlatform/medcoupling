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

IF(CPPUNIT_IS_OK)
  ADD_SUBDIRECTORY(Test)
ENDIF(CPPUNIT_IS_OK)

INCLUDE_DIRECTORIES(
  ${PARMETIS_INCLUDE_DIRS}
  ${METIS_INCLUDE_DIRS}
  ${SCOTCH_INCLUDE_DIRS}
  ${METIS_INCLUDE_DIRS}
  ${LIBXML_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MPI_INCLUDE_DIRS}
  ${PTHREADS_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

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
  MEDPARTITIONER_SkyLineArray.hxx
  MEDPARTITIONER_Topology.hxx
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
  MEDPARTITIONER_SkyLineArray.cxx
  )

SET(medpartitionercpp_DEFINITIONS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${LIBXML_DEFINITIONS} ${MPI_DEFINITIONS}")
SET(medpartitionercpp_LDFLAGS
  ${MED3_LIBS_C_ONLY}
  ${HDF5_LIBS}
  ${STDLIB}
  ${LIBXML_LIBS}
  ${MPI_LIBS}
  interpkernel
  medcoupling
  medloader
)
IF(MED_ENABLE_PARMETIS)
  SET(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX} MEDPARTITIONER_MetisGraph.hxx)
  SET(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES} MEDPARTITIONER_ParMetisGraph.cxx)
  SET(medpartitionercpp_DEFINITIONS "${medpartitionercpp_DEFINITIONS} ${PARMETIS_DEFINITIONS}")
  SET(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${PARMETIS_LIBS})
ELSE(MED_ENABLE_PARMETIS)  
  IF(MED_ENABLE_METIS)
    SET(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX} MEDPARTITIONER_MetisGraph.hxx)  
    SET(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES} MEDPARTITIONER_MetisGraph.cxx)
    SET(medpartitionercpp_DEFINITIONS "${medpartitionercpp_DEFINITIONS} ${METIS_DEFINITIONS}")
    SET(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${METIS_LIBS})
  ENDIF(MED_ENABLE_METIS)
  IF(MED_ENABLE_SCOTCH)
    SET(medpartitionercpp_HEADERS_HXX ${medpartitionercpp_HEADERS_HXX} MEDPARTITIONER_ScotchGraph.hxx)
    SET(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES} MEDPARTITIONER_ScotchGraph.cxx)
    SET(medpartitionercpp_DEFINITIONS "${medpartitionercpp_DEFINITIONS} ${SCOTCH_DEFINITIONS}")
    SET(medpartitionercpp_LDFLAGS ${medpartitionercpp_LDFLAGS} ${SCOTCH_LIBS})
  ENDIF(MED_ENABLE_SCOTCH)    
ENDIF(MED_ENABLE_PARMETIS)  
  
IF(MPI_IS_OK)
  SET(medpartitionercpp_SOURCES ${medpartitionercpp_SOURCES}
    MEDPARTITIONER_UtilsPara.cxx
    MEDPARTITIONER_JointFinder.cxx
  )
ENDIF(MPI_IS_OK)			  
  
ADD_LIBRARY(medpartitionercpp SHARED ${medpartitionercpp_SOURCES})
SET_TARGET_PROPERTIES(medpartitionercpp PROPERTIES COMPILE_FLAGS "${PLATFORM_DEFINITIONS} ${medpartitionercpp_DEFINITIONS}")
TARGET_LINK_LIBRARIES(medpartitionercpp ${medpartitionercpp_LDFLAGS} ${PLATFORM_LIBS} ${PTHREAD_LIBS})
INSTALL(TARGETS medpartitionercpp DESTINATION ${MED_salomelib_LIBS})
  
IF(MED_ENABLE_PARMETIS)
  ADD_EXECUTABLE(medpartitioner_para medpartitioner_para.cxx)
  SET_TARGET_PROPERTIES(medpartitioner_para PROPERTIES COMPILE_FLAGS "${medpartitionercpp_DEFINITIONS}")
  TARGET_LINK_LIBRARIES(medpartitioner_para medpartitionercpp ${medpartitionercpp_LDFLAGS})
  INSTALL(TARGETS medpartitioner_para DESTINATION ${MED_salomebin_BINS})
ELSE(MED_ENABLE_PARMETIS)
  ADD_EXECUTABLE(medpartitioner medpartitioner.cxx)
  SET_TARGET_PROPERTIES(medpartitioner PROPERTIES COMPILE_FLAGS "${PLATFORM_DEFINITIONS} ${medpartitionercpp_DEFINITIONS}")
  TARGET_LINK_LIBRARIES(medpartitioner medpartitionercpp ${medpartitionercpp_LDFLAGS})
  INSTALL(TARGETS medpartitioner DESTINATION ${MED_salomebin_BINS})
ENDIF(MED_ENABLE_PARMETIS)

INSTALL(FILES ${medpartitionercpp_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
