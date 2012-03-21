# Copyright (C) 2007-2011  CEA/DEN, EDF R&D, OPEN CASCADE
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

SET(KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR})

STRING(REPLACE "\\" "/" KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR})

IF(KERNEL_ROOT_DIR)
  SET(MED_ENABLE_KERNEL ON)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindPLATFORM.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindPYTHON.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindOMNIORB.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindPTHREADS.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindHDF5.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindBOOST.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindLIBXML2.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindSWIG.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindCPPUNIT.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindDOXYGEN.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindMPI.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindLIBBATCH.cmake)
  INCLUDE(${KERNEL_ROOT_DIR}/salome_adm/cmake_files/FindKERNEL.cmake)
ELSE(KERNEL_ROOT_DIR)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindPLATFORM.cmake)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local_without_kernel/cmake_files/FindMPI.cmake)
ENDIF(KERNEL_ROOT_DIR)

SET(GUI_ROOT_DIR $ENV{GUI_ROOT_DIR})

IF(GUI_ROOT_DIR)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/FindCAS.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/FindQT4.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/FindOPENGL.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/FindVTK.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/FindQWT.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/FindSIPPYQT.cmake)
  INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/FindGUI.cmake)
ENDIF(GUI_ROOT_DIR)

INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindPARMETIS.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindMETIS.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindSCOTCH.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindSPLITTER.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindRENUMBER.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindMEDFILE.cmake)
IF(WINDOWS)
  INCLUDE(${CMAKE_SOURCE_DIR}/adm_local/cmake_files/FindXDR.cmake)
ENDIF(WINDOWS)

SET(MED_salomebin_BINS bin/salome)

SET(MED_salomelib_LIBS lib/salome)

SET(MED_salomeidl_IDLS idl/salome)

SET(MED_salomeinclude_HEADERS include/salome)

SET(MED_salomeres_DATA share/salome/resources/MED)

SET(MED_salomeres_SCRIPTS share/salome/resources/MED)

SET(MED_salomescript_SCRIPTS bin/salome)

SET(MED_salomescript_DATA bin/salome)

SET(MED_salomescript_PYTHON bin/salome)

SET(MED_appliskel_SCRIPTS bin/salome/appliskel)

SET(MED_appliskel_PYTHON bin/salome/appliskel)

SET(MED_pythondir lib/python${PYTHON_VERSION}/site-packages)

SET(MED_salomepythondir ${MED_pythondir}/salome)

SET(MED_sharedpkgpython_PYTHON ${MED_salomepythondir}/shared_modules)

###############################################

IF(MED_ENABLE_KERNEL)
IF(MED_ENABLE_GUI)
  SET(ACLOCAL_AMFLAGS
      -I
      adm_local/unix/config_files
      -I
      ${KERNEL_ROOT_DIR}/salome_adm/unix/config_files
      -I
      ${GUI_ROOT_DIR}/adm_local/unix/config_files
  )
ELSE(MED_ENABLE_GUI)
  SET(ACLOCAL_AMFLAGS
      -I
      adm_local/unix/config_files
      -I
      ${KERNEL_ROOT_DIR}/salome_adm/unix/config_files
  )
ENDIF(MED_ENABLE_GUI)
ELSE(MED_ENABLE_KERNEL)
IF(MED_ENABLE_GUI)
  SET(ACLOCAL_AMFLAGS
      -I
      adm_local/unix/config_files
      -I
      adm_local_without_kernel/unix/config_files
      -I
      ${GUI_ROOT_DIR}/adm_local/unix/config_files
  )
ELSE(MED_ENABLE_GUI)
  SET(ACLOCAL_AMFLAGS
      -I
      adm_local/unix/config_files
      -I
      adm_local_without_kernel/unix/config_files
  )
ENDIF(MED_ENABLE_GUI)
ENDIF(MED_ENABLE_KERNEL)

IF(MED_ENABLE_KERNEL)
  SET(IDLDIR
      idl
  )
ENDIF(MED_ENABLE_KERNEL)

SET(SUBDIRS
    adm_local
)

IF(NOT MED_ENABLE_KERNEL)
SET(SUBDIRS
    ${SUBDIRS}
    adm_local_without_kernel
)
ENDIF(NOT MED_ENABLE_KERNEL)

SET(SUBDIRS
    ${SUBDIRS}
    ${IDLDIR}
    src
    doc
    resources
    bin
)

SET(DIST_SUBDIRS
    adm_local
    adm_local_without_kernel
    idl
    src
    doc
    resources
    bin
)

SET(DISTCLEANFILES
    a.out
    aclocal.m4
    configure
    local-install.sh
    hack_libtool
    adm_local/unix/config_files/config.guess
    adm_local/unix/config_files/config.sub
    adm_local/unix/config_files/depcomp
    adm_local/unix/config_files/install-sh
    adm_local/unix/config_files/libtool.m4
    adm_local/unix/config_files/ltmain.sh
    adm_local/unix/config_files/ltoptions.m4
    adm_local/unix/config_files/ltsugar.m4
    adm_local/unix/config_files/ltversion.m4
    adm_local/unix/config_files/lt~obsolete.m4
    adm_local/unix/config_files/missing
    adm_local/unix/config_files/py-compile
)

SET(salomeinclude_DATA
    MED_version.h
)

SET(EXTRA_DIST
    ${EXTRA_DIST}
    build_configure
    clean_configure
    LICENCE
    INSTALL
    README
    build_cmake
    build_cmake.bat
)

SET(input ${CMAKE_CURRENT_SOURCE_DIR}/Makefile.in)
SET(output ${CMAKE_CURRENT_BINARY_DIR}/Makefile)
MESSAGE(STATUS "Creation of ${output}")
CONFIGURE_FILE(${input} ${output})
SET(input ${CMAKE_CURRENT_SOURCE_DIR}/MED_version.h.in)
SET(output ${CMAKE_CURRENT_BINARY_DIR}/MED_version.h)
MESSAGE(STATUS "Creation of ${output}")
CONFIGURE_FILE(${input} ${output})

FOREACH(dir ${SUBDIRS})
  IF(NOT dir STREQUAL .)
    ADD_SUBDIRECTORY(${dir})
  ENDIF(NOT dir STREQUAL .)
ENDFOREACH(dir ${SUBDIRS})
            
FOREACH(f ${salomeinclude_DATA})
  SET(DEST include/salome)
  STRING(COMPARE EQUAL ${f} SALOMEconfig.h.in test_SALOMEconfig.h.in)
  IF(test_SALOMEconfig.h.in)
    INSTALL(FILES SALOMEconfig.ref.in DESTINATION ${DEST} RENAME SALOMEconfig.h.in)
  ELSE(test_SALOMEconfig.h.in)
    SET(dummy dummy-NOTFOUND)
    MARK_AS_ADVANCED(dummy)
    FIND_FILE(dummy ${f} PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)
    IF(dummy)
      GET_FILENAME_COMPONENT(ext ${f} EXT)
      IF(ext STREQUAL .py)
        IF(DEST STREQUAL bin/salome)
          SET(PERMS)
          SET(PERMS ${PERMS} OWNER_READ OWNER_WRITE OWNER_EXECUTE)
          SET(PERMS ${PERMS} GROUP_READ GROUP_EXECUTE)
          SET(PERMS ${PERMS} WORLD_READ WORLD_EXECUTE)
          INSTALL(FILES ${f} DESTINATION ${DEST} PERMISSIONS ${PERMS})
        ELSE(DEST STREQUAL bin/salome)
          INSTALL(FILES ${f} DESTINATION ${DEST})
        ENDIF(DEST STREQUAL bin/salome)
      ELSE(ext STREQUAL .py)
        INSTALL(FILES ${f} DESTINATION ${DEST})
      ENDIF(ext STREQUAL .py)  
    ELSE(dummy)
      GET_FILENAME_COMPONENT(ext ${f} EXT)
      IF(ext STREQUAL .qm)
        STRING(REGEX REPLACE .qm .ts input ${f})
        SET(input ${CMAKE_CURRENT_SOURCE_DIR}/${input})
        SET(output ${CMAKE_CURRENT_BINARY_DIR}/${f})
        EXECUTE_PROCESS(COMMAND ${QT_LRELEASE_EXECUTABLE} ${input} -qm ${output})
      ENDIF(ext STREQUAL .qm)
      INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${f} DESTINATION ${DEST})
    ENDIF(dummy)
  ENDIF(test_SALOMEconfig.h.in)
  GET_FILENAME_COMPONENT(ext ${f} EXT)
  IF(ext STREQUAL .py)
    INSTALL(CODE "SET(PYTHON_FILE ${f})")
    INSTALL(CODE "SET(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})")
    INSTALL(CODE "SET(DEST ${DEST})")
    INSTALL(CODE "SET(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})")
    #
    STRING(REPLACE "\\" "/" KERNEL_ROOT_DIR ${KERNEL_ROOT_DIR})
    INSTALL(SCRIPT ${KERNEL_ROOT_DIR}/salome_adm/cmake_files/install_and_compile_python_file.cmake)
  ENDIF(ext STREQUAL .py)
ENDFOREACH(f ${salomeinclude_DATA})
