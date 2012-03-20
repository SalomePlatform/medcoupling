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
  ${OMNIORB_INCLUDES_DIR}
  $ENV{KERNEL_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}
  )

SET(IDL_FILES
  MED.idl
  MED_Gen.idl
  Compo1Py.idl
  MEDCouplingCorbaServant.idl
  MEDCouplingCorbaServantTest.idl
  ParaMEDCouplingCorbaServant.idl
  )

SET(IDLCXXFLAGS -I${CMAKE_CURRENT_SOURCE_DIR} -I$ENV{KERNEL_ROOT_DIR}/idl/salome)

MACRO( CompileIDLPythonAndInstall OMNIORB_IDL_PYTHON IDLPYFLAGS IDL_FILE FULLDIR )
  FILE(MAKE_DIRECTORY ${FULLDIR})
  STRING(REPLACE " " ";" IDLPYFLAGS2 ${IDLPYFLAGS})
  MESSAGE(STATUS "Compiling ${IDL_FILE} into ${FULLDIR}")
  EXECUTE_PROCESS(COMMAND ${OMNIORB_IDL_PYTHON} ${IDLPYFLAGS2} -C${FULLDIR} ${IDL_FILE})
ENDMACRO( CompileIDLPythonAndInstall )

INSTALL(CODE "MACRO( CompileIDLPythonAndInstall OMNIORB_IDL_PYTHON IDLPYFLAGS IDL_FILE FULLDIR )
  FILE(MAKE_DIRECTORY \${FULLDIR})
  STRING(REPLACE \" \" \";\" IDLPYFLAGS2 \${IDLPYFLAGS})
  MESSAGE(STATUS \"Compiling \${IDL_FILE} into \${FULLDIR}\")
  EXECUTE_PROCESS(COMMAND \${OMNIORB_IDL_PYTHON} \${IDLPYFLAGS2} -C\${FULLDIR} \${IDL_FILE})
ENDMACRO( CompileIDLPythonAndInstall )")

FOREACH(input ${IDL_FILES})
  STRING(REGEX REPLACE ".idl" "" base ${input})
  SET(src ${CMAKE_CURRENT_BINARY_DIR}/${base}SK.cc)
  SET(outputs ${src})
  SET(dynsrc ${CMAKE_CURRENT_BINARY_DIR}/${base}DynSK.cc)
  SET(outputs ${outputs} ${dynsrc})
  SET(inc ${CMAKE_CURRENT_BINARY_DIR}/${base}.hh)
  SET(outputs ${outputs} ${inc})
  SET(input ${CMAKE_CURRENT_SOURCE_DIR}/${input})
  SET(flags ${IDLCXXFLAGS} ${OMNIORB_IDLCXXFLAGS})
  STRING(REGEX MATCH "-bcxx" ISBCXX ${flags})
  IF(NOT ISBCXX)
    SET(flags -bcxx ${flags})
  ENDIF(NOT ISBCXX)
  ADD_CUSTOM_COMMAND(OUTPUT ${outputs}
    COMMAND ${OMNIORB_IDL} ${flags} ${input}
    MAIN_DEPENDENCY ${input})
  
  INSTALL(FILES ${input} DESTINATION idl/salome)
  SET(IDL_HEADER ${CMAKE_CURRENT_BINARY_DIR}/${base}.hh)
  INSTALL(FILES ${IDL_HEADER} DESTINATION include/salome)
  INSTALL(CODE "SET(IDL_FILE ${input})")
  INSTALL(CODE "SET(DIR ${salomepythondir})")
  IF(WINDOWS)
    INSTALL(CODE "SET(DIR bin/salome)")
  ENDIF(WINDOWS)
  INSTALL(CODE "SET(CMAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})")
  INSTALL(CODE "SET(OMNIORB_IDL_PYTHON ${OMNIORB_IDL_PYTHON})")
  SET(flags)
  FOREACH(f ${IDLCXXFLAGS})
    SET(flags "${flags} ${f}")
  ENDFOREACH(f ${IDLCXXFLAGS})
  STRING(REGEX MATCH "-bpython" ISBPYTHON ${flags})
  IF(NOT ISBPYTHON)
    SET(flags "-bpython ${flags}")
  ENDIF(NOT ISBPYTHON)
  SET(IDLPYFLAGS ${flags})
  STRING(REPLACE "\\" "/" IDLPYFLAGS ${IDLPYFLAGS})
  INSTALL(CODE "SET(IDLPYFLAGS ${IDLPYFLAGS})")
  STRING(REPLACE "\\" "/" KERNEL_ROOT_DIR ${KERNEL_ROOT_DIR})
  INSTALL(CODE "CompileIDLPythonAndInstall( \"${OMNIORB_IDL_PYTHON}\" \"${IDLPYFLAGS}\" \"${input}\" \"${CMAKE_INSTALL_PREFIX}/${MED_salomepythondir}\" )")
ENDFOREACH(input ${IDL_FILES})

SET(SalomeIDLMED_SOURCES
  MEDSK.cc
  MED_GenSK.cc
  MEDCouplingCorbaServantSK.cc
  Compo1PySK.cc
  ParaMEDCouplingCorbaServantSK.cc
  )

ADD_LIBRARY(SalomeIDLMED SHARED ${SalomeIDLMED_SOURCES})
TARGET_LINK_LIBRARIES(SalomeIDLMED ${CORBA_LIBS} ${KERNEL_LDFLAGS} ${SalomeIDLKernel})
SET_TARGET_PROPERTIES(SalomeIDLMED PROPERTIES COMPILE_FLAGS "${OMNIORB_FLAGS}")

SET(SalomeIDLMEDTests_SOURCES
  MEDCouplingCorbaServantTestSK.cc
)

ADD_LIBRARY(SalomeIDLMEDTests SHARED ${SalomeIDLMEDTests_SOURCES})
TARGET_LINK_LIBRARIES(SalomeIDLMEDTests ${CORBA_LIBS} ${KERNEL_LDFLAGS} ${SalomeIDLKernel})
SET_TARGET_PROPERTIES(SalomeIDLMEDTests PROPERTIES COMPILE_FLAGS "${OMNIORB_FLAGS}")

#INSTALL(TARGETS medcoupling DESTINATION ${MED_salomelib_LIBS})

#FILE(GLOB medcoupling_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
#FILE(GLOB medcoupling_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
#INSTALL(FILES ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX} DESTINATION ${MED_salomeinclude_HEADERS})
