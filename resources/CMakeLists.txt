# Copyright (C) 2012-2016  CEA/DEN, EDF R&D
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

SET(MED_test_fig_files
  Pol1.fig
  Pol2.fig
  Pol3.fig
  Pol4.fig
)

SET(MED_other_FILES
  pointe.med
  poly3D.med
  polygones.med
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
  square1_split
  square1_split1.med
  square1_split2.med
  square1_split3.med
  square2.med
  square2_split
  square2_split1.med
  square2_split2.med
  Mesh3D_10_2d1.med
  Mesh3D_10_2d2.med
  Mesh3D_11.med
  allPillesTest.sauv
  BDC-714.sauv
  portico_3subs.sauv
  agitateur.med
  test_MED_MAIL.sauv
  )

SET(MED_RESOURCES_FILES ${MED_test_fig_files})

IF(NOT MEDCOUPLING_MICROMED)
    FOREACH(resfile ${MED_other_FILES})
        LIST(APPEND MED_RESOURCES_FILES ${resfile})    
    ENDFOREACH()
ENDIF()

# MEDFiles are always copied to build dir (for tests) but not always for install
IF(MEDCOUPLING_BUILD_TESTS)
  FOREACH(resfile ${MED_RESOURCES_FILES})
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/${resfile}" "${CMAKE_CURRENT_BINARY_DIR}/${resfile}" COPYONLY)
  ENDFOREACH(resfile)
ENDIF()

IF(MEDCOUPLING_WITH_FILE_EXAMPLES)
    INSTALL(FILES ${MED_RESOURCES_FILES} DESTINATION ${MEDCOUPLING_INSTALL_RES_DATA})
ENDIF()
