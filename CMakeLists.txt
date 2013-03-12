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

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.7 FATAL_ERROR)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

ENABLE_TESTING()
SET(VERSION "7.0.0")
SET(VERSION_DEV "1")
SET(WITH_MEDMEMGUI "0")

#Defining all options first
option(SALOME_USE_MPI "Use MPI containers" OFF)

SET(KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR})
IF(KERNEL_ROOT_DIR)
  FILE(TO_CMAKE_PATH $ENV{KERNEL_ROOT_DIR} KERNEL_ROOT_DIR)
  SET(MED_ENABLE_KERNEL ON)
  
  # directories for packages detection
  SET(SALOME_MACROS_DIR ${KERNEL_ROOT_DIR}/salome_adm/cmake_files)
  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SALOME_MACROS_DIR})

  find_package(PLATFORM)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindPYTHON.cmake)
  if(SALOME_USE_MPI)
    INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindMPI.cmake)
  endif(SALOME_USE_MPI)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindOMNIORB.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindPTHREADS.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindHDF5.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindBOOST.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindLIBXML2.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindSWIG.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindCPPUNIT.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindDOXYGEN.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindSPHINX.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindSalomeLIBBATCH.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/deprecated/FindKERNEL.cmake)

  SET(GUI_ROOT_DIR $ENV{GUI_ROOT_DIR})
  FILE(TO_CMAKE_PATH $ENV{GUI_ROOT_DIR} GUI_ROOT_DIR)
  IF(GUI_ROOT_DIR)
  SET(MED_ENABLE_GUI ON)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/deprecated/FindCAS.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/deprecated/FindQT4.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/deprecated/FindOPENGL.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/FindSalomeVTK.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/deprecated/FindQWT.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/deprecated/FindSIPPYQT.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/deprecated/FindGUI.cmake)
  ENDIF(GUI_ROOT_DIR)
ELSE(KERNEL_ROOT_DIR)
  SET(SALOME_MACROS_DIR ${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files)
  
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindPLATFORM.cmake)
  if(SALOME_USE_MPI)
    INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindMPI.cmake)
  endif(SALOME_USE_MPI)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindBOOST.cmake)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindCPPUNIT.cmake)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindHDF5.cmake)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindPYTHON.cmake)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindLIBXML2.cmake)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindPTHREADS.cmake)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindSWIG.cmake)

ENDIF(KERNEL_ROOT_DIR)

INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindPARMETIS.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindMETIS.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindSCOTCH.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindSPLITTER.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindRENUMBER.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindMEDFILE.cmake)
IF(WIN32)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindXDR.cmake)
ENDIF(WIN32)

EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; t=sys.argv[-1].split(\".\") ; t[:]=(int(elt) for elt in t) ; sys.stdout.write(\"0x%02x%02x%02x\"%tuple(t))" ${VERSION}
  OUTPUT_VARIABLE XVERSION)

SET(MED_salomebin_BINS bin/salome)

SET(MED_salomelib_LIBS lib/salome)

SET(MED_salomeidl_IDLS idl/salome)
SET(KERNEL_salomeidl_IDLS ${MED_salomeidl_IDLS})

SET(MED_salomeinclude_HEADERS include/salome)

SET(MED_salomeres_DATA share/salome/resources/med)

SET(MED_salomeres_SCRIPTS share/salome/resources/med)

SET(MED_salomescript_SCRIPTS bin/salome)

SET(MED_salomescript_DATA bin/salome)

SET(MED_salomescript_PYTHON bin/salome)
SET(KERNEL_salomescript_PYTHON ${MED_salomescript_PYTHON})

SET(MED_appliskel_SCRIPTS bin/salome/appliskel)

SET(MED_appliskel_PYTHON bin/salome/appliskel)

SET(MED_pythondir lib/python${PYTHON_VERSION}/site-packages)

SET(MED_salomepythondir ${MED_pythondir}/salome)
SET(KERNEL_salomepythondir ${MED_salomepythondir})

SET(MED_sharedpkgpython_PYTHON ${MED_salomepythondir}/shared_modules)

SET(MED_admlocalcmake adm_local/cmake_files)

IF(MED_ENABLE_KERNEL)
  ADD_SUBDIRECTORY(idl)
ELSE(MED_ENABLE_KERNEL)
  ADD_SUBDIRECTORY(adm_local_without_kernel)
ENDIF(MED_ENABLE_KERNEL)
ADD_SUBDIRECTORY(adm_local)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(resources)
ADD_SUBDIRECTORY(bin)

SET(input ${CMAKE_CURRENT_SOURCE_DIR}/MED_version.h.in)
SET(output ${CMAKE_CURRENT_BINARY_DIR}/MED_version.h)
MESSAGE(STATUS "Creation of ${output}")
CONFIGURE_FILE(${input} ${output} @ONLY)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/MED_version.h DESTINATION ${MED_salomeinclude_HEADERS})
