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

ADD_SUBDIRECTORY(dialogs)

INCLUDE(${GUI_ROOT_DIR}/adm_local/cmake_files/UseQT4EXT.cmake)

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
#SET(MEDGUI_FORMS MEDGUIFileContentDial.ui MEDGUISelectComponents.ui)
#SET(MEDGUI_RESOURCES images.qrc)

QT4_WRAP_CPP(MEDGUI_HEADERS_MOC ${MEDOPGUI_HEADERS})
#QT4_WRAP_UI(MEDGUI_FORMS_HEADERS ${MEDGUI_FORMS})
#QT4_ADD_RESOURCES(MEDGUI_RESOURCES_RCC ${MEDGUI_RESOURCES})

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIRS}
  ${CAS_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIRS}
  ${KERNEL_ROOT_DIR}/include/salome
  ${GUI_ROOT_DIR}/include/salome
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/../../../idl
  ${CMAKE_CURRENT_SOURCE_DIR}/../cmp
  ${CMAKE_CURRENT_SOURCE_DIR}/dialogs
  ${CMAKE_CURRENT_BINARY_DIR}/dialogs
  ${CMAKE_CURRENT_BINARY_DIR}/../../..
  )

SET(MEDOPGUITS_SOURCES
  MEDOP_msg_en.ts
  MEDOP_msg_fr.ts
  )
#-D__WITH_MEDMEMGUI__ libMEDGUI.la
ADD_LIBRARY(MEDOPGUI SHARED ${MEDOPGUI_SOURCES} ${MEDGUI_HEADERS_MOC})
SET_TARGET_PROPERTIES(MEDOPGUI PROPERTIES COMPILE_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${XDR_DEFINITIONS} ${CAS_DEFINITIONS} ${BOOST_DEFINITIONS}")
TARGET_LINK_LIBRARIES(MEDOPGUI MEDOPGUI_dialogs MEDOPFactoryEngine SalomeIDLMED ${QT_MT_LIBS} ${CAS_KERNEL} ${BOOST_LIBS} ${GUI_LDFLAGS} ${SalomeApp} ${CAM} ${LightApp} ${qtx} ${suit} ${SalomeObject})
INSTALL(TARGETS MEDOPGUI DESTINATION ${MED_salomelib_LIBS})
QT4_INSTALL_TS_RESOURCES("${MEDOPGUITS_SOURCES}" "${MED_salomeres_DATA}")
