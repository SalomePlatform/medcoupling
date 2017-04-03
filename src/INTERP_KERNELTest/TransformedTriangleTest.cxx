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

#include "TransformedTriangleTest.hxx"

#include <iostream>

using namespace INTERP_KERNEL;

namespace INTERP_TEST
{

  /**
   * Creates the TransformedTriangle objects used by the tests.
   *
   */
  void TransformedTriangleTest::setUp() 
  {
    // tri1 -> no unstable double products - no changes brought about by preCalculateDoubleProducts
    //         this allows the testing of calcUnstableT
    // tri2 -> unstable double products - for testing calcStableC / preCalculateDoubleProducts

    // triangle to test unstable C and T calculations
    p1[0] = -1.5 ; p1[1] = 0.5; p1[2] = 0.5;
    q1[0] = 2.0 ; q1[1] = 0.4; q1[2] = 0.6;
    r1[0] = 1.0 ; r1[1] = 2.4; r1[2] = 1.2;
    hp1 = 1 - p1[0] - p1[1] - p1[2];
    hq1 = 1 - q1[0] - q1[1] - q1[2];
    hr1 = 1 - r1[0] - r1[1] - r1[2]; 
    Hp1 = 1 - p1[0] - p1[1];
    Hq1 = 1 - q1[0] - q1[1];
    Hr1 = 1 - r1[0] - r1[1];

    //  std::cout <<std::endl<< "constructing tri1..." << std::endl;
    tri1 = new TransformedTriangle(p1, q1, r1);
 

    // triangle to test stable C calculation
    const double err = 1.5e-3;
    
    p2[0] = 0.000000000084654984189118; p2[1] = -0.000000000000000027536546231654231688873; p2[2] = 0.0000000000000001649875466831349431;
    q2[0] = -p2[0] +err; q2[1] = -p2[1] + err; q2[2] = -p2[2] +err;
    r2[0] = 2.01 ; r2[1] = 1.8; r2[2] = 0.92;
    
    hp2 = 1 - p2[0] - p2[1] - p2[2];
    hq2 = 1 - q2[0] - q2[1] - q2[2];
    hr2 = 1 - r2[0] - r2[1] - r2[2]; 
    Hp2 = 1 - p2[0] - p2[1];
    Hq2 = 1 - q2[0] - q2[1];
    Hr2 = 1 - r2[0] - r2[1];
    
    tri2 = new TransformedTriangle(p2, q2, r2);
  
  

  }

  /**
   * Liberates the transformed triangle objects used by the test suite
   * 
   */
  void TransformedTriangleTest::tearDown() 
  {
    delete tri1;
    delete tri2;
  }

  /// Tests that _coords has correct values after construction of object is finished
  /// \brief Status : pass
  void TransformedTriangleTest::test_constructor() {
    // test that _coords has correct values after constructor is called

    double good_values1[15] = 
      {
        p1[0], p1[1], p1[2], hp1, Hp1,
        q1[0], q1[1], q1[2], hq1, Hq1,
        r1[0], r1[1], r1[2], hr1, Hr1
      };

    double good_values2[15] = 
      {
        p2[0], p2[1], p2[2], hp2, Hp2,
        q2[0], q2[1], q2[2], hq2, Hq2,
        r2[0], r2[1], r2[2], hr2, Hr2
      };

  
    for(int i = 0 ; i < 15 ; ++i)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(good_values1[i], tri1->_coords[i], ERR_TOL);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(good_values2[i], tri2->_coords[i], ERR_TOL);
      }

    CPPUNIT_ASSERT_EQUAL(true, tri1->_is_double_products_calculated);
    CPPUNIT_ASSERT_EQUAL(true, tri2->_is_double_products_calculated);
  }

  /// Tests the calculation of double products (without the corrections)
  /// \brief Status : pass
  void TransformedTriangleTest::test_calcUnstableC() {
    typedef TransformedTriangle::TriSegment TriSegment;

    // test that the correct c-values are calculated
  
    double correct_c_vals[24] = 
      { 
        p1[0] * q1[1] - p1[1] * q1[0], 
        p1[1] * q1[2] - p1[2] * q1[1], 
        p1[2] * q1[0] - p1[0] * q1[2],
        p1[0] * hq1 - hp1 * q1[0],
        p1[1] * hq1 - hp1 * q1[1],
        p1[2] * hq1 - hp1 * q1[2],
        Hp1 * q1[0] - p1[0] * Hq1,
        p1[1] * Hq1 - Hp1 * q1[1],
        q1[0] * r1[1] - q1[1] * r1[0], 
        q1[1] * r1[2] - q1[2] * r1[1], 
        q1[2] * r1[0] - q1[0] * r1[2],
        q1[0] * hr1 - hq1 * r1[0],
        q1[1] * hr1 - hq1 * r1[1],
        q1[2] * hr1 - hq1 * r1[2],
        Hq1 * r1[0] - q1[0] * Hr1,
        q1[1] * Hr1 - Hq1 * r1[1],
        r1[0]*p1[1]-r1[1]*p1[0], 
        r1[1]*p1[2]-r1[2]*p1[1], 
        r1[2]*p1[0]-r1[0]*p1[2],
        r1[0] * hp1 - hr1 * p1[0],
        r1[1] * hp1 - hr1 * p1[1],
        r1[2] * hp1 - hr1 * p1[2],
        Hr1 * p1[0] - r1[0] * Hp1,
        r1[1] * Hp1 - Hr1 * p1[1]
      };

    double c_vals[3 * 8];
    for(TriSegment seg = TransformedTriangle::PQ ; seg <= TransformedTriangle::RP ; seg = TriSegment(seg + 1)) {
      
      c_vals[seg*8 + 0] = tri1->calcUnstableC(seg, TransformedTriangle::C_XY);
      c_vals[seg*8 + 1] = tri1->calcUnstableC(seg, TransformedTriangle::C_YZ);
      c_vals[seg*8 + 2] = tri1->calcUnstableC(seg, TransformedTriangle::C_ZX);
      c_vals[seg*8 + 3] = tri1->calcUnstableC(seg, TransformedTriangle::C_XH);
      c_vals[seg*8 + 4] = tri1->calcUnstableC(seg, TransformedTriangle::C_YH);
      c_vals[seg*8 + 5] = tri1->calcUnstableC(seg, TransformedTriangle::C_ZH);
      c_vals[seg*8 + 6] = tri1->calcUnstableC(seg, TransformedTriangle::C_01);
      c_vals[seg*8 + 7] = tri1->calcUnstableC(seg, TransformedTriangle::C_10);

    }
    
    for(int i = 0 ; i < 3*8 ; ++i) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL( correct_c_vals[i], c_vals[i], ERR_TOL );
    }


  }

  /// Tests the calculation of triple products (without corrections)
  /// \brief Status : pass
  void TransformedTriangleTest::test_calcUnstableT()
  {
    typedef TransformedTriangle::TetraCorner TetraCorner;

    // correct values calculated by determinants (Grandy, [15])
    const double correct_t_vals[4] = 
      {
        p1[0]*(q1[1]*r1[2] - q1[2]*r1[1]) -
        q1[0]*(p1[1]*r1[2] - p1[2]*r1[1]) +
        r1[0]*(p1[1]*q1[2] - p1[2]*q1[1]),

        -(hp1*(q1[1]*r1[2] - q1[2]*r1[1]) -
          hq1*(p1[1]*r1[2] - p1[2]*r1[1]) +
          hr1*(p1[1]*q1[2] - p1[2]*q1[1])),

        -(p1[0]*(hq1*r1[2] - q1[2]*hr1) -
          q1[0]*(hp1*r1[2] - p1[2]*hr1) +
          r1[0]*(hp1*q1[2] - p1[2]*hq1)),
    
        -(p1[0]*(q1[1]*hr1 - r1[1]*hq1) -
          q1[0]*(p1[1]*hr1 - r1[1]*hp1) +
          r1[0]*(p1[1]*hq1 - q1[1]*hp1))
      };
    

    // test that triple products are correctly calculated
    for(TetraCorner corner = TransformedTriangle::O ; corner <= TransformedTriangle::Z ; corner = TetraCorner(corner + 1)) 
      {
      
        for(int row = 1 ; row < 4 ; ++row)
          {
            const double t = tri1->calcTByDevelopingRow(corner, row, false);
            //    std::cout << std::endl  << " Corner = " << corner  << " Row = " << row << " got: " << t << 
            //  " expected: " << correct_t_vals[corner]<< std::endl;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(correct_t_vals[corner], t, ERR_TOL);    
          }
      }
  }

  /// Tests the consistency correction
  /// \brief Status : fails because it is not significant - the consistency correction is not brought into play
  void TransformedTriangleTest::test_calcStableC_Consistency()
  {

    typedef TransformedTriangle::TriSegment TriSegment;
    typedef TransformedTriangle::TetraCorner TetraCorner;

    // grandy, eq 14
    double correct_c_vals[24] = 
      { 
        p2[0] * q2[1] - p2[1] * q2[0], 
        p2[1] * q2[2] - p2[2] * q2[1], 
        p2[2] * q2[0] - p2[0] * q2[2],
        p2[0] * hq2 - hp2 * q2[0],
        p2[1] * hq2 - hp2 * q2[1],
        p2[2] * hq2 - hp2 * q2[2],
        Hp2 * q2[0] - p2[0] * Hq2,
        p2[1] * Hq2 - Hp2 * q2[1],
        q2[0] * r2[1] - q2[1] * r2[0], 
        q2[1] * r2[2] - q2[2] * r2[1], 
        q2[2] * r2[0] - q2[0] * r2[2],
        q2[0] * hr2 - hq2 * r2[0],
        q2[1] * hr2 - hq2 * r2[1],
        q2[2] * hr2 - hq2 * r2[2],
        Hq2 * r2[0] - q2[0] * Hr2,
        q2[1] * Hr2 - Hq2 * r2[1],
        r2[0]*p2[1]-r2[1]*p2[0], 
        r2[1]*p2[2]-r2[2]*p2[1], 
        r2[2]*p2[0]-r2[0]*p2[2],
        r2[0] * hp2 - hr2 * p2[0],
        r2[1] * hp2 - hr2 * p2[1],
        r2[2] * hp2 - hr2 * p2[2],
        Hr2 * p2[0] - r2[0] * Hp2,
        r2[1] * Hp2 - Hr2 * p2[1]
      };


    // number of inconsistent cases found : 
    // should be (at least) 1 for the test to be meaningful
    int num_cases = 0; 

    // find unstable products to check for consistency (Grandy [46])  
    for(TriSegment seg = TransformedTriangle::PQ ; seg <= TransformedTriangle::RP ; seg = TriSegment(seg + 1)) 
      {
        const double c_xy = tri2->calcUnstableC(seg, TransformedTriangle::C_XY);
        const double c_yz = tri2->calcUnstableC(seg, TransformedTriangle::C_YZ);
        const double c_zx = tri2->calcUnstableC(seg, TransformedTriangle::C_ZX);
        const double c_xh = tri2->calcUnstableC(seg, TransformedTriangle::C_XH);
        const double c_yh = tri2->calcUnstableC(seg, TransformedTriangle::C_YH);
        const double c_zh = tri2->calcUnstableC(seg, TransformedTriangle::C_ZH);
      
        const int num_zero = (c_yz*c_xh == 0.0 ? 1 : 0) + (c_zx*c_yh == 0.0 ? 1 : 0) + (c_xy*c_zh == 0.0 ? 1 : 0);
        const int num_neg = (c_yz*c_xh < 0.0 ? 1 : 0) + (c_zx*c_yh < 0.0 ? 1 : 0) + (c_xy*c_zh < 0.0 ? 1 : 0);
      
        if((num_zero == 1 && num_neg != 1) || num_zero == 2 || (num_neg == 0 && num_zero !=3) || num_neg == 3 )
          {
            ++num_cases;
  
            double min_dist = -1.0; // initialised first time through loop
            TetraCorner min_corner = TransformedTriangle::O;
  
            for(TetraCorner corner = TransformedTriangle::O ; corner <= TransformedTriangle::Z ; corner = TetraCorner(corner + 1))
              {
                // calculate distance from each corner of tetraeder to the segment
                // formula : ( (Q-P) x (P - corner) )^2 / norm(Q-P)^2
      
                const double ptP[3] = { tri2->_coords[5*seg], tri2->_coords[5*seg + 1], tri2->_coords[5*seg + 2] };
                const double ptQ[3] = { tri2->_coords[5*( (seg+1) % 3)], tri2->_coords[5*( (seg+1) % 3) + 1], tri2->_coords[5*( (seg+1) % 3) + 2] };
                const double ptCorner[3] = { 
                  corner == TransformedTriangle::X ? 1.0 : 0.0,
                  corner == TransformedTriangle::Y ? 1.0 : 0.0,
                  corner == TransformedTriangle::Z ? 1.0 : 0.0,
                };

                const double diff_21[3] = { ptQ[0] - ptP[0], ptQ[1] - ptP[1], ptQ[2] - ptP[2] };
                const double diff_1_corner[3] = { ptP[0] - ptCorner[0], ptP[1] - ptCorner[1], ptP[2] - ptCorner[2] };
      
                const double cross[3] = { 
                  diff_21[1]*diff_1_corner[2] - diff_21[2]*diff_1_corner[1],  
                  diff_21[2]*diff_1_corner[0] - diff_21[0]*diff_1_corner[2],
                  diff_21[0]*diff_1_corner[1] - diff_21[1]*diff_1_corner[0]
                };

                const double cross_sq = cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2];

                const double norm_pq = diff_21[0]*diff_21[0] + diff_21[1]*diff_21[1] + diff_21[2]*diff_21[2];

                if(corner == TransformedTriangle::O || (cross_sq / norm_pq) < min_dist)
                  {
                    min_dist = cross_sq / norm_pq;
                    min_corner = corner;
                  }
              }
  
            // now check if the corresponding double products are zero
            static const DoubleProduct DOUBLE_PRODUCTS[12] = 
              {
                TransformedTriangle::C_YZ, TransformedTriangle::C_XY, TransformedTriangle::C_ZX, // O
                TransformedTriangle::C_ZH, TransformedTriangle::C_YZ, TransformedTriangle::C_YH, // X
                TransformedTriangle::C_ZH, TransformedTriangle::C_ZX, TransformedTriangle::C_XH, // Y
                TransformedTriangle::C_XY, TransformedTriangle::C_YH, TransformedTriangle::C_XH  // Z
              };

            for(int i = 0; i < 3 ; ++i) 
              {
                DoubleProduct dp = DOUBLE_PRODUCTS[3*min_corner + i];
                //        std::cout << std::endl << "in test inconsistent (seg,dp) :(" << seg <<", " << dp << ")" << std::endl;
                CPPUNIT_ASSERT_EQUAL(0.0, tri2->calcStableC(seg, dp));
                correct_c_vals[8*seg + dp] = 0.0;
              }
          }

      }
  
    if(num_cases < 1)
      {
        CPPUNIT_FAIL("Consistency test not pertinent");
      }

    //  std::cout << std::endl << "Number of geometric inconsistencies : " << num_cases << std::endl; 
    
    // check that all other double products have right value too
    double c_vals[8*3];

    for(TriSegment seg = TransformedTriangle::PQ ; seg <= TransformedTriangle::RP ; seg = TriSegment(seg + 1)) {
      
      c_vals[seg*8 + 0] = tri2->calcStableC(seg, TransformedTriangle::C_XY);
      c_vals[seg*8 + 1] = tri2->calcStableC(seg, TransformedTriangle::C_YZ);
      c_vals[seg*8 + 2] = tri2->calcStableC(seg, TransformedTriangle::C_ZX);
      c_vals[seg*8 + 3] = tri2->calcStableC(seg, TransformedTriangle::C_XH);
      c_vals[seg*8 + 4] = tri2->calcStableC(seg, TransformedTriangle::C_YH);
      c_vals[seg*8 + 5] = tri2->calcStableC(seg, TransformedTriangle::C_ZH);
      c_vals[seg*8 + 6] = tri2->calcStableC(seg, TransformedTriangle::C_01);
      c_vals[seg*8 + 7] = tri2->calcStableC(seg, TransformedTriangle::C_10);

    }

    for(int i = 0 ; i < 24 ; ++i)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(correct_c_vals[i], c_vals[i], ERR_TOL);
      }
  }

}
