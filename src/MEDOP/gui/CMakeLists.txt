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

ADD_DEFINITIONS(${MEDFILE_DEFINITIONS} ${XDR_DEFINITIONS} ${CAS_DEFINITIONS} ${BOOST_DEFINITIONS} ${OMNIORB_DEFINITIONS})
INCLUDE(${QT_USE_FILE})
SET(QT_LIBS "${QT_LIBRARIES}")
INCLUDE(UseQt4Ext)

ADD_SUBDIRECTORY(dialogs)


SET(MEDOPGUI_SOURCES
  MEDOPFactoryClient.cxx
  XmedConsoleDriver.cxx
  WorkspaceController.cxx
  XmedDataModel.cxx
  MEDEventListener_i.cxx
  MEDOPModule.cxx
  DatasourceController.cxx
  factory.cxx
)

SET(MEDOPGUI_HEADERS MEDOPModule.hxx MEDEventListener_i.hxx WorkspaceController.hxx DatasourceController.hxx)

QT4_WRAP_CPP(MEDOPGUI_HEADERS_MOC ${MEDOPGUI_HEADERS})

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIR} ${OMNIORB_INCLUDE_DIRS}
  ${CAS_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR}
  ${GUI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/dialogs
  ${CMAKE_CURRENT_BINARY_DIR}/dialogs
  ${CMAKE_CURRENT_SOURCE_DIR}/../cmp
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDGUI
  ${PROJECT_BINARY_DIR}/idl
  ${PROJECT_BINARY_DIR}
)

SET(MEDOPGUITS_SOURCES
  MED_images.ts
  MED_msg_en.ts
  MED_msg_fr.ts
  MED_msg_ja.ts
  MEDOP_msg_en.ts
  MEDOP_msg_fr.ts
  MEDOP_msg_ja.ts
)

SET(COMMON_FLAGS
  ${CAS_KERNEL}
  ${QT_MT_LIBRARIES}
  ${OMNIORB_LIBRARIES}
  ${PLATFORM_LIBRARIES}
  ${BOOST_LIBRARIES}
  SalomeIDLMED
  MEDOPGUI_dialogs
  MEDOPFactoryEngine
  ${GUI_qtx}
  ${GUI_suit}
  ${GUI_SalomeObject}
  ${KERNEL_SalomeLifeCycleCORBA}
  ${KERNEL_SalomeKernelHelpers}
  ${GUI_SalomeApp}
  ${GUI_SalomeGuiHelpers}
  ${GUI_SalomeTreeData}
  ${KERNEL_OpUtil}
  ${GUI_CAM}
  ${GUI_LightApp}
  ${GUI_PyConsole}
)

ADD_LIBRARY(MEDOPGUI SHARED ${MEDOPGUI_SOURCES} ${MEDOPGUI_HEADERS_MOC})
IF(WITH_MEDMEMGUI)
  SET_TARGET_PROPERTIES(MEDOPGUI PROPERTIES COMPILE_FLAGS "-D__WITH_MEDMEMGUI__")
ENDIF(WITH_MEDMEMGUI)
TARGET_LINK_LIBRARIES(MEDOPGUI ${COMMON_FLAGS})

INSTALL(TARGETS MEDOPGUI DESTINATION ${SALOME_INSTALL_LIBS})
QT4_INSTALL_TS_RESOURCES("${MEDOPGUITS_SOURCES}" "${SALOME_MED_INSTALL_RES_DATA}")

FILE(GLOB COMMON_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${COMMON_HEADERS_HXX} DESTINATION ${SALOME_INSTALL_HEADERS})
