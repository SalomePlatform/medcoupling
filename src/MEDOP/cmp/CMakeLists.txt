# Copyright (C) 2012-2013  CEA/DEN, EDF R&D
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
  ${OMNIORB_INCLUDE_DIR} ${OMNIORB_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIRS}
  ${KERNEL_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/idl
  ${CMAKE_SOURCE_DIR}/src/MEDCoupling
  ${CMAKE_SOURCE_DIR}/src/INTERP_KERNEL/Bases
  ${CMAKE_SOURCE_DIR}/src/INTERP_KERNEL
  ${CMAKE_SOURCE_DIR}/src/MEDLoader
)

SET(MEDOPFactoryEngine_SOURCES
  MEDDataManager_i.cxx
  MEDCalculator_i.cxx
  MEDOPFactory_i.cxx
)
SET(COMMON_LIBS
  medloader
  SalomeIDLMED
  ${TOOLSDS}
  ${SalomeHDFPersist}
  ${SalomeContainer}
  ${SalomeCommunication}
  ${SalomeKernelHelpers}
  ${SalomeLifeCycleCORBA}
  ${SALOMELocalTrace}
  ${SALOMEBasics} 
  ${SalomeGenericObj}
  medcoupling
  interpkernel
  ${MED3_LIBS_C_ONLY}
  ${HDF5_LIBS}
  ${OMNIORB_LIBS}
  ${PLATFORM_LIBS}
)
SET(COMMON_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${PLATFORM_DEFINITIONS} ${OMNIORB_DEFINITIONS}")
# This undefines the macros MIN and MAX which are specified in the windows headers 
IF(WINDOWS)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -DNOMINMAX")
ENDIF(WINDOWS)

ADD_LIBRARY(MEDOPFactoryEngine SHARED ${MEDOPFactoryEngine_SOURCES})
SET_TARGET_PROPERTIES(MEDOPFactoryEngine PROPERTIES COMPILE_FLAGS "${COMMON_FLAGS}")
TARGET_LINK_LIBRARIES(MEDOPFactoryEngine SalomeIDLMED ${COMMON_LIBS})

INSTALL(TARGETS MEDOPFactoryEngine DESTINATION ${MED_salomelib_LIBS})

FILE(GLOB MEDOPFactoryEngine_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${MEDOPFactoryEngine_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})

INSTALL(FILES test_medop_components.py PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ DESTINATION ${MED_salomebin_BINS}/xmed)
