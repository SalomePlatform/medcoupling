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

SET(MEDOPGUI_dialogs_SOURCES
  GenericDialog.cxx
  DlgAlias.cxx
  DlgUseInWorkspace.cxx
  DlgChangeUnderlyingMesh.cxx
  DlgImageToMed.cxx
)

SET(MEDOPGUI_dialogs_HEADERS GenericDialog.hxx DlgAlias.hxx DlgUseInWorkspace.hxx DlgChangeUnderlyingMesh.hxx DlgImageToMed.hxx)
SET(MEDOPGUI_dialogs_FORMS DlgAlias.ui DlgChangeUnderlyingMesh.ui DlgImageToMed.ui DlgUseInWorkspace.ui GenericDialog.ui)
#SET(MEDGUI_RESOURCES images.qrc)

QT4_WRAP_CPP(MEDOPGUI_dialogs_HEADERS_MOC ${MEDOPGUI_dialogs_HEADERS})
QT4_WRAP_UI_XX(MEDOPGUI_dialogs_FORMS_HEADERS ${MEDOPGUI_dialogs_FORMS})
#QT4_ADD_RESOURCES(MEDGUI_RESOURCES_RCC ${MEDGUI_RESOURCES})

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIRS}
  ${CAS_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIRS}
  ${KERNEL_ROOT_DIR}/include/salome
  ${GUI_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  )

ADD_LIBRARY(MEDOPGUI_dialogs SHARED ${MEDOPGUI_dialogs_SOURCES} ${MEDOPGUI_dialogs_HEADERS_MOC} ${MEDOPGUI_dialogs_FORMS_HEADERS})
SET_TARGET_PROPERTIES(MEDOPGUI_dialogs PROPERTIES COMPILE_FLAGS "${CAS_DEFINITIONS} ${BOOST_DEFINITIONS}")
TARGET_LINK_LIBRARIES(MEDOPGUI_dialogs medcalculator ${QT_MT_LIBS} ${CAS_KERNEL} ${BOOST_LIBS} ${GUI_LDFLAGS} ${SalomeApp} ${CAM} ${LightApp} ${qtx} ${suit} ${SalomeObject})
INSTALL(TARGETS MEDOPGUI_dialogs DESTINATION ${MED_salomelib_LIBS})

ADD_EXECUTABLE(dlgTester dlgTester.cxx)
SET_TARGET_PROPERTIES(dlgTester PROPERTIES COMPILE_FLAGS "${CAS_DEFINITIONS} ${BOOST_DEFINITIONS}")
TARGET_LINK_LIBRARIES(dlgTester MEDOPGUI_dialogs)
#ADD_TEST(dlgTester dlgTester)