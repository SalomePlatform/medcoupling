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
# Author : Anthony Geay (CEA/DEN)

INCLUDE(UseOmniORB)

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIR} ${OMNIORB_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
)

SET(SalomeIDLMED_IDLSOURCES
  MEDCouplingCorbaServant.idl
  ParaMEDCouplingCorbaServant.idl
  MEDDataManager.idl
  MEDCalculator.idl
  MEDEventListener.idl
  MEDPresentationManager.idl
  MEDOPFactory.idl
)

IF(SALOME_USE_MPI)
  SET(SalomeIDLMED_IDLSOURCES ${SalomeIDLMED_IDLSOURCES} ParaMEDMEMComponent.idl)
ENDIF(SALOME_USE_MPI)

SET(SalomeIDLMEDTests_IDLSOURCES
  MEDCouplingCorbaServantTest.idl
)

SET(_idl_include_dirs
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${KERNEL_ROOT_DIR}/idl/salome
)

SET(_idl_link_flags
  ${KERNEL_SalomeIDLKernel}
)

OMNIORB_ADD_MODULE(SalomeIDLMED      "${SalomeIDLMED_IDLSOURCES}"      "${_idl_include_dirs}" "${_idl_link_flags}")
OMNIORB_ADD_MODULE(SalomeIDLMEDTests "${SalomeIDLMEDTests_IDLSOURCES}" "${_idl_include_dirs}" "SalomeIDLMED")

INSTALL(TARGETS SalomeIDLMED EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${SALOME_INSTALL_LIBS})
INSTALL(TARGETS SalomeIDLMEDTests                              DESTINATION ${SALOME_INSTALL_LIBS})
