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

SET(MED_RESOURCES_FILES
  #boitenew.cnc
  #boitenew.inp
  #boitenew.xyz
  #Case1.cnc
  #Case1.inp
  #Case1.xyz
  #cube.cnc
  #cube.inp
  #cube.xyz
  #test3.cnc
  #test3.inp
  #test3.xyz
  #titi.cnc
  #titi.inp
  #titi.xyz
  #carre_en_quad4.med
  #carre_en_quad4_seg2.med
  #cas_defaut_domaine_fluide.med
  #ChampsDarcy.med
  #cube_hexa8.med
  #cube_hexa8_quad4.med
  #darcy_1.1_res.med
  #darcy_1.3_resCASTEM.med
  #darcy_1.3_resPORFLOW.med
  #darcy_1.3_resTRACES.med
  #darcy2_Castem_EFMH.med
  #darcy2_Castem_qua_EFMH.med
  #darcy2_Castem_qua_VF.med
  #Darcy3_3D_H_10x10x10_2.med
  #Darcy3_3D_H_10x10x10.sauve
  #dx200_dy1_avec_2couches.sauve
  #elle_2D_QT_10x10.sauve
  #elle_2D_QT_2x2.sauve
  #elle_2D_QT_40x40.sauve
  #elle_2D_QT_4x4.sauve
  #elle_3D_HPr_10x10x10_2.med
  #elle_3D_HPr_10x10x10.sauve
  #elle_3D_HPr_2x2x2_2.med
  #elle_3D_HPr_2x2x2.sauve
  #elle_3D_HPr_4x4x4_2.med
  #elle_3D_HPr_4x4x4.sauve
  #extendedtransport53_triangles.med
  #geomMesh_nomorereadable21.med
  #geomMesh22.med
  #H_CastCast_EFMH_I129_COUPLEX1.med
  #H_CastCast_VF_I129_COUPLEX1.med
  #H_CastCast_VF_Se79_COUPLEX1.med
  #H_CastPorf_I129_COUPLEX1.med
  #H_CastPorf_Se79_COUPLEX1.med
  #H_PorfCast_EFMH_I129_COUPLEX1.med
  #H_PorfCast_EFMH_Se79_COUPLEX1.med
  #H_PorfPorf_I129_COUPLEX1.med
  #H_Traces_I129_COUPLEX1.med
  #H_Traces_Se79_COUPLEX1.med
  # inclusion_2d_raf.sauve
  # inclusion_2d.sauve
  # mail_ktest1-3-hexa.sauve
  # mail_ktest1-3-tetra.sauve
  # mail_ktest3-1.sauve
  # mail_ktest3-2.sauve
  # maill.00.med
  # maill.00_nofield.med
  # maill.00_nomesh.med
  # maill.00_without_seg2.med
  # maill.0.med
  # maillage_5_5_5.med
  # maillage_andra2_100elts.sauve
  # maillage_cas2_2d.sauve
  # maillage_cas4_234elts.sauve
  # maillage_CHEMVAL_100elts.sauve
  # maillage_CHEMVAL_40elts.sauve
  # maillage_chemvalIV_cas1_100elts.sauve
  # maillage_chemvalIV_cas1_40elts.med
  # maillage_chemvalIV_cas1_40elts.sauve
  # maillage_UniSegFam.med
  # mail_test1-1-qua.sauve
  # mail_test1-1-tri.sauve
  # mail_test1-2-qua.sauve
  # mail_test1-2-tri.sauve
  # mail-test1-4-1.sauve
  # mail-test1-4-2.sauve
  # mesh_nomorereadable21.med
  # mesh.med
  # Mistrat.med
  #Old_ChampsDarcy.med
  #Old_darcy_1.1_res.med
  #Old_darcy_1.3_resCASTEM.med
  #Old_darcy_1.3_resPORFLOW.med
  #Old_darcy_1.3_resTRACES.med
  #Old_darcy2_Castem_EFMH.med
  #Old_darcy2_Castem_qua_EFMH.med
  #Old_darcy2_Castem_qua_VF.med
  #Old_H_CastCast_EFMH_I129_COUPLEX1.med
  #Old_H_CastCast_VF_I129_COUPLEX1.med
  #Old_H_CastCast_VF_Se79_COUPLEX1.med
  #Old_H_CastPorf_I129_COUPLEX1.med
  #Old_H_CastPorf_Se79_COUPLEX1.med
  #Old_H_PorfCast_EFMH_I129_COUPLEX1.med
  #Old_H_PorfCast_EFMH_Se79_COUPLEX1.med
  #Old_H_PorfPorf_I129_COUPLEX1.med
  #Old_H_PorfPorf_Se79_COUPLEX1.med
  #Old_H_Traces_I129_COUPLEX1.med
  #Old_H_Traces_Se79_COUPLEX1.med
  #Old_maillage_chemvalIV_cas1_40elts.med
  #pointe_nomorereadable21.med
  pointe.med
  #poly3D.med
  #polyedres.med
  #polygones.med
  #recoll_bord.med
  #test19.med
  # test_2D.med
  # trio_2D.med
  # TimeStamps.med
  # zzzz121b.med
  # zzzz121b_without_tr6.med
  UnitTetra.med
  GeneralTetra.med
  NudgedSimpler.med
  NudgedTetra.med
  CornerTetra.med
  SimpleIncludedTetra.med
  SimpleIncludingTetra.med
  Test2D.med
  Test2Dpoly.med
  Test3D.med
  Test3Dpoly.med
  UnitTetraDegenT.med
  DegenEdgeXY.med
  DegenFaceXYZ.med
  DegenTranslatedInPlane.med
  ComplexIncludedTetra.med
  ComplexIncludingTetra.med
  HalfstripOnly.med
  HalfstripOnly2.med
  SimpleHalfstripOnly.med
  GenTetra1.med
  GenTetra2.med
  TrickyTetra1.med
  LargeUnitTetra.med
  # LargeInconsistentTetra.med
  DividedUnitTetraSimpler.med
  DividedUnitTetra.med
  NudgedDividedUnitTetra.med
  NudgedDividedUnitTetraSimpler.med
  DividedGenTetra1.med
  DividedGenTetra2.med
  Box1.med
  Box2.med
  Box3.med
  Box1Moderate.med
  Box2Moderate.med
  BoxModSmall1.med
  BoxModSmall2.med
  BoxEvenSmaller1.med
  TinyBox.med
  BoxHexa1.med
  BoxHexa2.med
  MovedHexaBox1.med
  MovedHexaBox2.med
  BoxTetra2.med
  square1.med
  # square1_split
  # square1_split1.med
  # square1_split2.med
  # square1_split3.med
  square2.med
  # square2_split
  # square2_split1.med
  # square2_split2.med
  # testStructCart3D.med
  # Mesh3D_10_2d1.med
  # Mesh3D_10_2d2.med
  # Mesh3D_11.med
  Pol1.fig
  Pol2.fig
  Pol3.fig
  Pol4.fig
  # blow5_ascii.case
  # blow5_ascii.geo
  # blow5_ascii_cd_displacement
  # blow5_ascii_cd_thickness
  # blow5_ascii_pd_displacement
  # blow5_ascii_pd_thickness
  #test_2D.sauve
  allPillesTest.sauv
  BDC-714.sauv
  portico_3subs.sauv
  agitateur.med
  )

# IF(MED_ENABLE_GUI)
#   FILE(GLOB MED_RESOURCES_PNG "${CMAKE_CURRENT_SOURCE_DIR}/*.png")
#   SET(MED_RESOURCES_FILES
#     ${MED_RESOURCES_FILES}
#     ModuleMed.png
#     Data.png
#     Infos.png
#     Structure.png
#     explore_med_file.png
#     field_selection.png
#     med_mesh.png
#     med_field.png
#     mesh_selection.png
#     )
# ENDIF(MED_ENABLE_GUI)

FOREACH(resfile ${MED_RESOURCES_FILES})
  CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/${resfile}" "${CMAKE_CURRENT_BINARY_DIR}/${resfile}" COPYONLY)
ENDFOREACH(resfile)

INSTALL(FILES ${MED_RESOURCES_FILES} DESTINATION ${MEDTOOL_INSTALL_RES_DATA})

#MESSAGE(STATUS "Creation of ${CMAKE_CURRENT_BINARY_DIR}/MEDCatalog.xml")
#CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/MEDCatalog.xml.in ${CMAKE_CURRENT_BINARY_DIR}/MEDCatalog.xml @ONLY)
#MESSAGE(STATUS "Creation of ${CMAKE_CURRENT_BINARY_DIR}/SalomeApp.xml")
#CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/SalomeApp.xml.in ${CMAKE_CURRENT_BINARY_DIR}/SalomeApp.xml @ONLY)

#INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/MEDCatalog.xml ${CMAKE_CURRENT_BINARY_DIR}/SalomeApp.xml DESTINATION ${SALOME_MED_INSTALL_RES_DATA})
