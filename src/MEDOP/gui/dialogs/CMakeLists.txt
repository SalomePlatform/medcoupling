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

INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/UseQT4EXT.cmake)

MACRO (QT4_WRAP_UI_XX outfiles )
  QT4_EXTRACT_OPTIONS(ui_files ui_options ${ARGN})
  FOREACH (it ${ui_files})
    GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
    GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
    SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/ui_${outfile}.hxx)
    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
      COMMAND ${QT_UIC_EXECUTABLE}
      ARGS ${ui_options} -o ${outfile} ${infile}
      MAIN_DEPENDENCY ${infile})
    SET(${outfiles} ${${outfiles}} ${outfile})
  ENDFOREACH (it)
ENDMACRO (QT4_WRAP_UI_XX)

SET(MEDOPGUI_dialogs_FORMS
  GenericDialog.ui
  DlgAlias.ui
  DlgUseInWorkspace.ui
  DlgChangeUnderlyingMesh.ui
  DlgImageToMed.ui
)
QT4_WRAP_UI_XX(MEDOPGUI_dialogs_FORMS_HEADERS ${MEDOPGUI_dialogs_FORMS})

SET(MEDOPGUI_dialogs_FILES
  GenericDialog.hxx
  DlgAlias.hxx
  DlgUseInWorkspace.hxx
  DlgChangeUnderlyingMesh.hxx
  DlgImageToMed.hxx
)
QT4_WRAP_CPP(MEDOPGUI_dialogs_HEADERS_MOC ${MEDOPGUI_dialogs_FILES})

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIRS}
  ${CAS_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${KERNEL_ROOT_DIR}/include/salome
  ${GUI_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

SET(MEDOPGUI_dialogs_SOURCES
  GenericDialog.cxx
  DlgAlias.cxx
  DlgUseInWorkspace.cxx
  DlgChangeUnderlyingMesh.cxx
  DlgImageToMed.cxx
)

ADD_LIBRARY(MEDOPGUI_dialogs SHARED ${MEDOPGUI_dialogs_SOURCES} ${MEDOPGUI_dialogs_HEADERS_MOC} ${MEDOPGUI_dialogs_FORMS_HEADERS})
SET_TARGET_PROPERTIES(MEDOPGUI_dialogs PROPERTIES COMPILE_FLAGS "${CAS_DEFINITIONS} ${BOOST_DEFINITIONS} ${PLATFORM_DEFINITIONS} ${OMNIORB_DEFINITIONS}")
TARGET_LINK_LIBRARIES(MEDOPGUI_dialogs ${CAS_KERNEL} ${QT_LIBS} ${SalomeIDLKernel} ${SalomeGuiHelpers} ${OMNIORB_LIBS} ${qtx} ${suit} ${PLATFORM_LIBS} ${SalomeLifeCycleCORBA} ${SalomeKernelHelpers} ${SalomeGuiHelpers})
INSTALL(TARGETS MEDOPGUI_dialogs DESTINATION ${MED_salomelib_LIBS})

#
# ---------------------------------------------------------
# Unit test programs
# ---------------------------------------------------------
#
ADD_EXECUTABLE(dlgTester dlgTester.cxx)
SET_TARGET_PROPERTIES(dlgTester PROPERTIES COMPILE_FLAGS "${PLATFORM_DEFINITIONS}")
TARGET_LINK_LIBRARIES(dlgTester MEDOPGUI_dialogs ${QT_LIBS} ${PLATFORM_LIBS})
INSTALL(TARGETS dlgTester DESTINATION ${MED_salomebin_BINS})
#ADD_TEST(dlgTester dlgTester)

FILE(GLOB COMMON_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${COMMON_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})

