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

INCLUDE_DIRECTORIES(
  ${HDF5_INCLUDES_DIR}
  )

SET(TOOLS_HDFI_SOURCES
  MEDattrFermer.cxx
  MEDattrNumEcrire.cxx
  MEDattrNumLire.cxx
  MEDattrOuvrir.cxx
  MEDattrStringEcrire.cxx
  MEDattrStringLire.cxx
  MEDdatagroupCreer.cxx
  MEDdatagroupFermer.cxx
  MEDdatagroupOuvrir.cxx
  MEDdatasetFermer.cxx
  MEDdatasetNumEcrire.cxx
  MEDdatasetNumLire.cxx
  MEDdatasetOuvrir.cxx
  MEDdatasetStringEcrire.cxx
  MEDdatasetStringLire.cxx
  MEDfichierCreer.cxx
  MEDfichierFermer.cxx
  MEDfichierOuvrir.cxx
  MEDindiceInfo.cxx
  MEDindiceNum.cxx
  MEDmodeErreurVerrouiller.cxx
  MEDnObjets.cxx
  MEDobjetIdentifer.cxx
  )

SET(TOOLS_MISC_SOURCES
  MED1cstring.cxx
  MEDnomDataset.cxx
  MED2cstring.cxx
  MEDnomEntite.cxx
  MEDGeometrieElement.cxx
  MEDnomGeometrie.cxx
  MEDcstringFree.cxx
  MEDparametresGeometrie.cxx
  MEDfstring.cxx
  )

SET(API_CI_SOURCES
  MEDchampCr.cxx
  MEDchampEcr.cxx
  MEDchampInfo.cxx
  MEDchampLire.cxx
  MEDconnEcr.cxx
  MEDconnLire.cxx
  MEDcoordEcr.cxx
  MEDcoordLire.cxx
  MEDdimLire.cxx
  MEDelementsEcr.cxx
  MEDelementsLire.cxx
  MEDequivCr.cxx
  MEDequivEcr.cxx
  MEDequivInfo.cxx
  MEDequivLire.cxx
  MEDfam2groA.cxx
  MEDfam2groB.cxx
  MEDfamCr.cxx
  MEDfamEcr.cxx
  MEDfamInfo.cxx
  MEDfamLire.cxx
  MEDfamMaaCr.cxx
  MEDfamMaaInfo.cxx
  MEDfamMaaLire.cxx
  MEDfermer.cxx
  MEDfichDesEcr.cxx
  MEDfichEntete.cxx
  MEDgro2famA.cxx
  MEDgro2famB.cxx
  MEDlFichDes.cxx
  MEDmaaCr.cxx
  MEDmaaInfo.cxx
  MEDnChamp.cxx
  MEDnCorres.cxx
  MEDnEntMaa.cxx
  MEDnEntites.cxx
  MEDnEquiv.cxx
  MEDnFam.cxx
  MEDnMaa.cxx
  MEDnPasdetemps.cxx
  MEDnProfil.cxx
  MEDnVal.cxx
  MEDnValProfil.cxx
  MEDnbnoisEcr.cxx
  MEDnbnoisLire.cxx
  MEDnbnomaEcr.cxx
  MEDnbnomaLire.cxx
  MEDnbnosoEcr.cxx
  MEDnbnosoLire.cxx
  MEDnoeudsEcr.cxx
  MEDnoeudsLire.cxx
  MEDnomEcr.cxx
  MEDnomLire.cxx
  MEDnumEcr.cxx
  MEDnumLire.cxx
  MEDouvrir.cxx
  MEDpasdetempsInfo.cxx
  MEDprofilEcr.cxx
  MEDprofilLire.cxx
  MEDprofilInfo.cxx
  MEDunvCr.cxx
  MEDunvLire.cxx
  MEDformatConforme.cxx
  MEDversionConforme.cxx
  MEDversionDonner.cxx
  MEDversionLire.cxx
  MEDbodyFittedEcr.cxx
  MEDbodyFittedLire.cxx
  MEDfamGridEcr.cxx
  MEDfamGridLire.cxx
  MEDgridCr.cxx
  MEDgridEcr.cxx
  MEDgridInfo.cxx
  MEDgridLire.cxx
  MEDnGrid.cxx
  )

SET(med_V2_1_SOURCES
  ${TOOLS_HDFI_SOURCES}
  ${TOOLS_MISC_SOURCES}
  ${API_CI_SOURCES}
  )


SET(mdump_V2_1_SOURCES
  mdump_V2_1.cxx
  )

SET(test1_V2_1_SOURCES
  test1_V2_1.cxx
  )

ADD_LIBRARY(med_V2_1 SHARED ${med_V2_1_SOURCES})
SET_TARGET_PROPERTIES(med_V2_1 PROPERTIES COMPILE_FLAGS "-D${MACHINE} ${HDF5_FLAGS}")
TARGET_LINK_LIBRARIES(med_V2_1 ${HDF5_LIBS})
INSTALL(TARGETS med_V2_1 DESTINATION ${MED_salomelib_LIBS})

ADD_EXECUTABLE(mdump_V2_1 ${mdump_V2_1_SOURCES})
SET_TARGET_PROPERTIES(mdump_V2_1 PROPERTIES COMPILE_FLAGS "-D${MACHINE} ${HDF5_FLAGS}")
TARGET_LINK_LIBRARIES(mdump_V2_1 med_V2_1)
ADD_EXECUTABLE(test1_V2_1 ${mdump_V2_1_SOURCES})
SET_TARGET_PROPERTIES(test1_V2_1 PROPERTIES COMPILE_FLAGS "-D${MACHINE} ${HDF5_FLAGS}")
TARGET_LINK_LIBRARIES(test1_V2_1 med_V2_1)
INSTALL(TARGETS mdump_V2_1 test1_V2_1 DESTINATION ${MED_salomebin_BINS})

FILE(GLOB med_V2_1_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${med_V2_1_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
