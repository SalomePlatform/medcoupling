// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

// File      : UnitTetraIntersectionBaryTest.cxx
// Created   : Thu Dec 11 15:54:41 2008
// Author    : Edward AGAPOV (eap)
//
#include "UnitTetraIntersectionBaryTest.hxx"

#include "UnitTetraIntersectionBary.hxx"
#include "TetraAffineTransform.hxx"
#include "InterpolationUtils.hxx"
#include "SplitterTetra.txx"

#include <iostream>

using namespace INTERP_KERNEL;

namespace INTERP_TEST
{
  void fill_UnitTetraIntersectionBary(UnitTetraIntersectionBary& bary, double nodes[][3])
  {
    int faceConn[4][3] = { { 0, 1, 2 },// inverse order
                           { 0, 3, 1 },
                           { 1, 3, 2 },
                           { 3, 0, 2 } };
//     int faceConn[4][3] = { { 0, 2, 1 },
//                            { 0, 1, 3 },
//                            { 1, 2, 3 },
//                            { 3, 2, 0 } };
    bary.init(true);
    for ( int i = 0; i < 4; ++i ) {
      int* faceNodes = faceConn[ i ];
      TransformedTriangle tri(nodes[faceNodes[0]], nodes[faceNodes[1]], nodes[faceNodes[2]]);
      tri.calculateIntersectionVolume();
      bary.addSide( tri );
    }
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_1()
  {
    // cutting tetra coincides with the unit one
    double nodes[4][3] = { { 0.0, 0.0, 0.0 },
                           { 1.0, 0.0, 0.0 },
                           { 0.0, 1.0, 0.0 },
                           { 0.0, 0.0, 1.0 } };
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.166667, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[1], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[2], 1e-5);
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_2()
  {
    // cutting tetra fully include the unit one
    double nodes[4][3] = { {-0.1,-0.1,-0.1 },
                           { 1.5,-0.1,-0.1 },
                           {-0.1, 1.5,-0.1 },
                           {-0.1,-0.1, 1.5 } };
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.166667, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[1], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[2], 1e-5);
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_3()
  {
    // cutting tetra is same as the unit one but moved up by 0.5
    double nodes[4][3] = { { 0.0, 0.0, 0.5 },
                           { 1.0, 0.0, 0.5 },
                           { 0.0, 1.0, 0.5 },
                           { 0.0, 0.0, 1.5 } };
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.020833333333333332, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.125, baryCenter[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.125, baryCenter[1], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.625, baryCenter[2], 1e-5);
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_4()
  {
    // same as previous but no cutting sides lay on the sides of unit tetra
    double nodes[4][3] = { {-0.2,-0.2, 0.5 },
                           { 1.0, 0.0, 0.5 },
                           { 0.0, 1.0, 0.5 },
                           { 0.0, 0.0, 2.0 } };
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.020833333333333332, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.125, baryCenter[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.125, baryCenter[1], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.625, baryCenter[2], 1e-5);
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_5()
  {
    // cutting tetra is similar and parallel to the UT but moved (-0.1,-0.1,-0.1)
    double nodes[4][3] = { {-0.1,-0.1,-0.1 },
                           { 1.1,-0.1,-0.1 },
                           {-0.1, 1.1,-0.1 },
                           {-0.1,-0.1, 1.1 } };
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1215, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.225, baryCenter[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.225, baryCenter[1], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.225, baryCenter[2], 1e-5);
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_6()
  {
    // cutting tetra is deeped into the UT with one corner
    double nodes[4][3] = { { 0.2, 0.2, 0.2 },
                           { 1.0, 0.2, 0.2 },
                           { 0.9, 1.0, 0.2 },
                           { 0.9, 9.0, 1.0 } };
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.000441855, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.353463 , baryCenter[0], 1e-5 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.33877  , baryCenter[1], 1e-5 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.207767 , baryCenter[2], 1e-5 );
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_7()
  {
    // cutting tetra passes through the UT with one corner
    double nodes[4][3] = { {-0.2, 0.2, 0.2 },
                           { 1.0, 0.2, 0.2 },
                           { 0.9, 1.0, 0.2 },
                           { 0.9, 0.9, 1.0 } };
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0103501, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.215578 , baryCenter[0], 1e-5 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.341363 , baryCenter[1], 1e-5 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.263903 , baryCenter[2], 1e-5 );
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_8()
  {
    // cutting tetra passes through the UT with one edge
    double nodes[4][3] = { { 0.5, 0.2, -0.2 }, // O
                           {-0.5,-0.2, -0.2 }, // OX
                           { 1.0,-0.5, -0.2 }, // OY
                           { 0.5, 0.2,  1.5 } };//OZ
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0349217, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.332275  , baryCenter[0], 1e-2 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0565892 , baryCenter[1], 1e-3 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.308713  , baryCenter[2], 1e-2 );
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_9()
  {
    // cutting tetra touches the UT at an edge, intersection volume == 0
    double nodes[4][3] = { { 1.0, 0.0, 0.0 }, // 0
                           {-1.0, 2.0, 2.0 }, // OX
                           {-1.0,-2.0, 2.0 }, // OY
                           { 1.0, 0.0, 2.0 } };//OZ
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( !ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, vol, 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( -1. , baryCenter[0], 1e-5 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( -1. , baryCenter[1], 1e-5 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( -1. , baryCenter[2], 1e-5 );
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_10()
  {
    // cutting tetra fully includes theUT and touches it at an edge
    double nodes[4][3] = { { 1.0, 0.0, 0.0 }, // 0
                           {-1.0,-4.0, 2.0 }, // OX
                           {-1.0, 4.0, 2.0 }, // OY
                           { 1.0, 0.0,-2.0 } };//OZ
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.166667, vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[1], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, baryCenter[2], 1e-5);
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_11()
  {
    // cutting tetra intersects the UT and touches it at an edge
    double nodes[4][3] = { { 1.0, 0.0, 0.0 }, // 0
                           {-1.0,-4.0, 2.0 }, // OX
                           {-1.0, 4.0, 2.0 }, // OY
                           {-1.0, 0.0,-1.0 } };//OZ
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.15873 , vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.250000, baryCenter[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.230952, baryCenter[1], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.260714, baryCenter[2], 1e-5);
  }
  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_12()
  {
    // cutting tetra has one corner inside the UT and one its side passes through an UT edge
    double nodes[4][3] = { { 0.25, 0.25, 0.25 }, // 0
                           { 1.75,-0.25,-0.25 }, // OX
                           { 0.5 , 0.25, 0.25 }, // OY
                           { 0.5 , 0   , 0.5  } };//OZ
    UnitTetraIntersectionBary bary;
    fill_UnitTetraIntersectionBary(bary,nodes);
    double baryCenter[3];
    bool ok    = bary.getBary( baryCenter );
    double vol = bary.getVolume();
    CPPUNIT_ASSERT( ok );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.005208 , vol, 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.562500, baryCenter[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.125000, baryCenter[1], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.250000, baryCenter[2], 1e-5);
  }

  struct __MESH_DUMMY
  {
    typedef int MyConnType;
  };

  void UnitTetraIntersectionBaryTest::test_UnitTetraIntersectionBary_13()
  {
    double T[] = {
      66.6666666666666714,133.333333333333343,66.6666666666666714,
      100,200,100,
      100,100,100,
      200,200,0 };
    
    double S[] = {
      100,166.666666666666657,66.6666666666666714,
      100,150,50,
      75,150,75,
      100,100,100};

    int conn[4] = { 0,1,2,3 };
    
    const double* tnodes[4]={ T, T+3, T+6, T+9 };
    const double* snodes[4]={ S, S+3, S+6, S+9 };
    
    __MESH_DUMMY dummyMesh;
    SplitterTetra<__MESH_DUMMY> src( dummyMesh, snodes, conn );
    double volume = src.intersectTetra( tnodes );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6944.4444444444443,volume,1e-9);
  }

  void UnitTetraIntersectionBaryTest::test_TetraAffineTransform_reverseApply()
  {
    double nodes[12] = { -4.0, 9.0, 3.0, 
                         11.0, 0.0, 2.0, 
                         0.0, 0.0, 0.0, 
                         2.0, 1.0,10.0 };
    //    double pSrc[3] = { -4.0, 9.0, 3.0 };
    double pSrc[3] = { 40., -20., 100. };
    double pDest[] = {1,1,1};
    TetraAffineTransform a(nodes);
    a.apply( pDest, pSrc );
    a.reverseApply( pDest, pDest );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( pSrc[0], pDest[0], 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( pSrc[1], pDest[1], 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( pSrc[2], pDest[2], 1e-12);
  }

  void UnitTetraIntersectionBaryTest::test_barycentric_coords()
  {
    // compute barycentric coordinates
    double nodes[4][3] = { {11.0, 0.0, 2.0 },
                           {-4.0, 9.0, 3.0 },
                           { 0.0, 0.0, 0.0 }, 
                           { 6.0, 1.0,10.0 }};
    std::vector<const double*> n (4);
    n[0] = &nodes[0][0];
    n[1] = &nodes[1][0];
    n[2] = &nodes[2][0];
    n[3] = &nodes[3][0];
    double p  [3] = { 2, 2, 5 }, bc[4];
    barycentric_coords(n, p, bc);
    double bcSum = 0;
    double p2 [3] = { 0,0,0 };
    for ( int i = 0; i < 4; ++i ) {
      bcSum += bc[i];
      p2[0] += bc[i]*n[i][0];
      p2[1] += bc[i]*n[i][1];
      p2[2] += bc[i]*n[i][2];
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., bcSum, 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( p[0], p2[0], 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( p[1], p2[1], 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( p[2], p2[2], 1e-12);
  }

  /* Conventions:
  *   - for HEXA8, point 5 is taken to be the origin (see med file ref connec):
  *          0 ------ 3
            /|       /|
           / |      / |
          1 ------ 2  |
          |  |     |  |
          |  |     |  |
          |  4-----|- 7
          | /      | /
          5 ------ 6
   */
  void UnitTetraIntersectionBaryTest::test_cuboid_mapped_coords_3D()
  {
    double nodes[8][3] = { { 0.0, 2.0, 4.0 }, //0
                           { 0.0, 0.0, 4.0 },
                           { 1.0, 0.0, 4.0 },
                           { 1.0, 2.0, 4.0 },
                           { 0.0, 2.0, 0.0 }, // 4
                           { 0.0, 0.0, 0.0 },
                           { 1.0, 0.0, 0.0 },
                           { 1.0, 2.0, 0.0 }
    };
    // Translate cube:
    for (int i=0; i < 8; ++i)
      for (int j=0; j < 3; ++j)
        nodes[i][j] += 15.0;

    std::vector<const double*> n (8);
    for (int i=0; i<8; i++)
      n[i] = &nodes[i][0];

    {
        // middle point
        double p[3] = { 15.5, 16.0, 17.0 }, bc[3];
        cuboid_mapped_coords(n, p, bc);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, bc[0], 1e-12);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, bc[1], 1e-12);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, bc[2], 1e-12);
    }
    {
      // point 1
      double p[3] = { 15.0, 15.0, 19.0 }, bc[3];
      cuboid_mapped_coords(n, p, bc);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, bc[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, bc[1], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[2], 1e-12);
    }
    {
      // point 7
      double p[3] = { 16.0, 17.0, 15.0 }, bc[3];
      cuboid_mapped_coords(n, p, bc);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[1], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, bc[2], 1e-12);
    }
    {
      // point 3
      double p[3] = { 16.0, 17.0, 19.0 }, bc[3];
      cuboid_mapped_coords(n, p, bc);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[1], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[2], 1e-12);
    }
    {
      // point outside
      double p[3] = { 2.0, 16.0, 18.0 }, bc[3];
      CPPUNIT_ASSERT_THROW(cuboid_mapped_coords(n, p, bc), INTERP_KERNEL::Exception);
    }

  }

  /* Convention
      - for QUAD4, point 0 is taken to be the origin (again see med file ref connec):

         1------2
         |      |
         |      |
         0------3
  */
  void UnitTetraIntersectionBaryTest::test_quad_mapped_coords_2D()
  {

    double nodes[4][2] = { { 0.0, 0.0 },
                           { 0.0, 1.0 },
                           { 2.0, 3.0 },
                           { 1.0, 0.0 } };

    // Translate quad4:
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 2; ++j)
        nodes[i][j] += 15.0;

    std::vector<const double*> n (4);
    for (int i=0; i<4; i++)
      n[i] = &nodes[i][0];

    {
      // middle point
      double p[2] = { 15.75, 16.0 }, bc[2];
      quad_mapped_coords(n, p, bc);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, bc[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, bc[1], 1e-12);
    }

    {
      // middle point of seg
      double p[2] = { 15.5, 15.0 }, bc[2];
      quad_mapped_coords(n, p, bc);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, bc[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, bc[1], 1e-12);
    }

    {
      // point 1
      double p[2] = { 15.0, 16.0 }, bc[2];
      quad_mapped_coords(n, p, bc);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, bc[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[1], 1e-12);
    }
    {
      // point 2
      double p[2] = { 17.0, 18.0 }, bc[2];
      quad_mapped_coords(n, p, bc);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[1], 1e-12);
    }
    {
      // point 3
      double p[2] = { 16.0, 15.0 }, bc[2];
      quad_mapped_coords(n, p, bc);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, bc[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, bc[1], 1e-12);
    }
    {
      // point outside
      double p[2] = { 18.0, 18.0 }, bc[2];
      CPPUNIT_ASSERT_THROW(quad_mapped_coords(n, p, bc), INTERP_KERNEL::Exception);
    }
  }


}
