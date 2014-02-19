# Copyright (C) 2012-2014  CEA/DEN, EDF R&D, OPEN CASCADE
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

ADD_DEFINITIONS(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${XDR_DEFINITIONS} ${CAS_DEFINITIONS} ${BOOST_DEFINITIONS} ${OMNIORB_DEFINITIONS})

INCLUDE(${QT_USE_FILE})
INCLUDE(UseQt4Ext)

SET(MEDGUI_SOURCES
  MEDGUIFileContentDial.cxx
  MEDGUISelection.cxx
  MEDGUISelectComponents.cxx
  MEDGUIDataBaseDockWidget.cxx
)

SET(MEDGUI_HEADERS MEDGUIFileContentDial.h MEDGUISelectComponents.h)
SET(MEDGUI_FORMS MEDGUIFileContentDial.ui MEDGUISelectComponents.ui)
#SET(MEDGUI_RESOURCES images.qrc)

QT4_WRAP_CPP(MEDGUI_HEADERS_MOC ${MEDGUI_HEADERS})
QT4_WRAP_UI(MEDGUI_FORMS_HEADERS ${MEDGUI_FORMS})
#QT4_ADD_RESOURCES(MEDGUI_RESOURCES_RCC ${MEDGUI_RESOURCES})

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIR} ${OMNIORB_INCLUDE_DIRS}
  ${CAS_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIRS}
  ${GUI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/../../idl
  ${CMAKE_CURRENT_BINARY_DIR}/../..
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCalculator
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(MEDGUITS_SOURCES
  MED_images.ts
  MED_msg_fr.ts
  MED_msg_en.ts
  MED_msg_ja.ts
)

ADD_LIBRARY(MEDGUI SHARED ${MEDGUI_SOURCES} ${MEDGUI_HEADERS_MOC} ${MEDGUI_FORMS_HEADERS})
TARGET_LINK_LIBRARIES(MEDGUI medcalculator ${QT_MT_LIBS} ${CAS_KERNEL} ${BOOST_LIBRARIES} ${GUI_LDFLAGS} ${SalomeApp} ${CAM} ${LightApp} ${qtx} ${suit} ${SalomeObject} ${PLATFORM_LIBS} ${SalomeLifeCycleCORBA})
INSTALL(TARGETS MEDGUI DESTINATION ${SALOME_INSTALL_LIBS})
QT4_INSTALL_TS_RESOURCES("${MEDGUITS_SOURCES}" "${SALOME_MED_INSTALL_RES_DATA}")
