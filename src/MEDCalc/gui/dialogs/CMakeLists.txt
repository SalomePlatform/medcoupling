# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
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

ADD_DEFINITIONS(${CAS_DEFINITIONS} ${BOOST_DEFINITIONS} ${OMNIORB_DEFINITIONS})

INCLUDE(UseQt4Ext)

SET(MEDCALCGUI_dialogs_FORMS
  GenericDialog.ui
  DlgAlias.ui
  DlgUseInWorkspace.ui
  DlgChangeUnderlyingMesh.ui
  DlgInterpolateField.ui
  DlgImageToMed.ui
)
QT4_WRAP_UI(MEDCALCGUI_dialogs_FORMS_HEADERS ${MEDCALCGUI_dialogs_FORMS})

SET(MEDCALCGUI_dialogs_FILES
  GenericDialog.hxx
  DlgAlias.hxx
  DlgUseInWorkspace.hxx
  DlgChangeUnderlyingMesh.hxx
  DlgInterpolateField.hxx
  DlgImageToMed.hxx
)
QT4_WRAP_CPP(MEDCALCGUI_dialogs_HEADERS_MOC ${MEDCALCGUI_dialogs_FILES})

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIR} ${OMNIORB_INCLUDE_DIRS}
  ${CAS_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${GUI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${QT4_INCLUDE_DIRECTORY}
)

SET(MEDCALCGUI_dialogs_SOURCES
  GenericDialog.cxx
  DlgAlias.cxx
  DlgUseInWorkspace.cxx
  DlgChangeUnderlyingMesh.cxx
  DlgInterpolateField.cxx
  DlgImageToMed.cxx
)

ADD_LIBRARY(MEDCALCGUI_dialogs SHARED ${MEDCALCGUI_dialogs_SOURCES} ${MEDCALCGUI_dialogs_HEADERS_MOC} ${MEDCALCGUI_dialogs_FORMS_HEADERS})
TARGET_LINK_LIBRARIES(MEDCALCGUI_dialogs ${CAS_KERNEL} ${QT_LIBRARIES} ${KERNEL_SalomeIDLKernel}
    ${GUI_SalomeGuiHelpers} ${OMNIORB_LIBRARIES} ${GUI_qtx} ${GUI_suit} ${PLATFORM_LIBRARIES}
    ${KERNEL_SalomeLifeCycleCORBA} ${KERNEL_SalomeKernelHelpers} ${GUI_SalomeGuiHelpers})
INSTALL(TARGETS MEDCALCGUI_dialogs DESTINATION ${SALOME_INSTALL_LIBS})

#
# ---------------------------------------------------------
# Unit test programs
# ---------------------------------------------------------
#
ADD_EXECUTABLE(dlgTester dlgTester.cxx)
TARGET_LINK_LIBRARIES(dlgTester MEDCALCGUI_dialogs ${QT_LIBS} ${PLATFORM_LIBS})
INSTALL(TARGETS dlgTester DESTINATION ${SALOME_INSTALL_BINS})
#ADD_TEST(dlgTester dlgTester)

FILE(GLOB COMMON_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${COMMON_HEADERS_HXX} DESTINATION ${SALOME_INSTALL_HEADERS})
