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

#include "BBTreeTest.hxx"
#include <iostream>
#include <vector>
#include "DirectedBoundingBox.hxx"

namespace INTERP_TEST
{


  void BBTreeTest::setUp() 
  {
  }

 
  void BBTreeTest::tearDown() 
  {
  }

  /**
   * Test that creates a tree in 2D and check that 
   * the results are correct in three
   * cases :
   * a non matching search
   * a standard case
   * a bbox overlapping the bboxes of the tree
   */
  void BBTreeTest::test_BBTree() {
    //bbox tree creation
    const int N=10;
    double* bbox=new double[4*N*N];
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
        {
          bbox[4*(i*N+j)]=i;
          bbox[4*(i*N+j)+1]=i+1;
          bbox[4*(i*N+j)+2]=j;
          bbox[4*(i*N+j)+3]=j+1;
        }
    BBTree<2> tree(bbox,0,0,N*N);
    std::vector <int> elems; 
  
    //box outside the tree
    double bbox1[4]={-2.0, -1.0, 0.0, 1.0};
    tree.getIntersectingElems(bbox1,elems);
    CPPUNIT_ASSERT_EQUAL(0,(int)elems.size());
    elems.clear();
  
    //box intersecting 4 tree elems
    double bbox2[4]={2.5, 3.5, 0.5, 1.5};
    tree.getIntersectingElems(bbox2,elems);
    CPPUNIT_ASSERT_EQUAL(4,(int)elems.size());
    elems.clear();
  
    //box exactly superimposed to two tree elems
    double bbox3[4]={5.0,6.0,7.0,9.0};
    tree.getIntersectingElems(bbox3,elems);
    CPPUNIT_ASSERT_EQUAL(2,(int)elems.size());
    elems.clear();

    double xx[2]={1.0,1.0};
    tree.getElementsAroundPoint(xx,elems);
    CPPUNIT_ASSERT_EQUAL(4,(int)elems.size());

    delete[] bbox;
  }

  void BBTreeTest::test_DirectedBB_3D()
  {
    // a rectangle 1x2 extruded along vector (10,0,10)
    const int nbP = 8, dim = 3;
    double coords[nbP*dim] =
      {
        0,0,0,    2,0,0,   2,1,0,   0,1,0, 
        10,0,10, 12,0,10, 12,1,10, 10,1,10
      };
    INTERP_KERNEL::DirectedBoundingBox bb( coords, nbP, dim);
    bb.enlarge( 1e-12 );

    // corners of extrusion are IN
    for ( int i = 0; i < nbP*dim; i+=dim )
      CPPUNIT_ASSERT( !bb.isOut( coords + i ));

    // points near corners of extrusion are OUT
    double p[nbP*dim] = 
      {
        0,0,3,  6,0,3,   5,1,2,   0,1,2, 
        8,0,9, 11,0,8, 11,0.5,8, 8,0.5,9
      };
    for ( int i = 0; i < nbP*dim; i+=dim )
      CPPUNIT_ASSERT( bb.isOut( p + i ));

    // the extrusions  shifted by 3 in XOY plane are OUT
    double shifted_X[nbP*dim] =
      {
        3,0,0,    5,0,0,   5,1,0,   3,1,0, 
        13,0,10, 15,0,10, 15,1,10, 13,1,10
      };
    double shifted_x[nbP*dim] =
      {
        -3,0,0, -1,0,0, -1,1,0, -3,1,0, 
        7,0,10, 9,0,10, 9,1,10, 7,1,10
      };
    double shifted_Y[nbP*dim] =
      {
        0,3,0,    2,3,0,   2,4,0,   0,4,0, 
        10,3,10, 12,3,10, 12,4,10, 10,4,10
      };
    double shifted_y[nbP*dim] =
      {
        0,-3,0,    2,-3,0,   2,-2,0,   0,-2,0, 
        10,-3,10, 12,-3,10, 12,-2,10, 10,-2,10
      };
    INTERP_KERNEL::DirectedBoundingBox shiftedBB_x( shifted_x, nbP, dim);
    INTERP_KERNEL::DirectedBoundingBox shiftedBB_X( shifted_X, nbP, dim);
    INTERP_KERNEL::DirectedBoundingBox shiftedBB_y( shifted_y, nbP, dim);
    INTERP_KERNEL::DirectedBoundingBox shiftedBB_Y( shifted_Y, nbP, dim);

    CPPUNIT_ASSERT( bb.isDisjointWith( shiftedBB_x ));
    CPPUNIT_ASSERT( bb.isDisjointWith( shiftedBB_X ));
    CPPUNIT_ASSERT( bb.isDisjointWith( shiftedBB_y ));
    CPPUNIT_ASSERT( bb.isDisjointWith( shiftedBB_Y ));

    // intersecting box is IN
    double inters_coords[nbP*dim] =
      {
        0,0,0,    2,0,0,   2,1,0,   0,1,0, 
        0,0,2,    2,0,2,   2,1,2,   0,1,2
      };
    INTERP_KERNEL::DirectedBoundingBox ibb( inters_coords, nbP, dim);
    CPPUNIT_ASSERT( !bb.isDisjointWith( ibb ));

    // overlapping non-directed BB
    double overlappingBB[2*dim] =
      {
        5,6, 0, 1, -5,4
      };
    CPPUNIT_ASSERT( !bb.isDisjointWith( overlappingBB ));

    // non-overlapping non-directed BB
    double nonoverlappingBB_1[2*dim] =
      {
        5,6, 0,1, -5,2
      };
    CPPUNIT_ASSERT( bb.isDisjointWith( nonoverlappingBB_1 ));
    double nonoverlappingBB_2[2*dim] =
      {
        5,6, 0,1, 7,20
      };
    CPPUNIT_ASSERT( bb.isDisjointWith( nonoverlappingBB_2 ));
  }

  void BBTreeTest::test_DirectedBB_2D()
  {
    // a segment of length 2 extruded along vector (10,10)
    const int nbP = 4, dim = 2;
    double coords[nbP*dim] =
      {
        0,0,    2,0,
        10,10, 12,10,
      };
    INTERP_KERNEL::DirectedBoundingBox bb( coords, nbP, dim);
    bb.enlarge( 1e-12 );

    // corners of extrusion are IN
    for ( int i = 0; i < nbP*dim; i+=dim )
      CPPUNIT_ASSERT( !bb.isOut( coords + i ));

    // points near corners of extrusion are OUT
    double p[nbP*dim] = 
      {
        1,2,  4,1,
        11,8, 8,9,
      };
    for ( int i = 0; i < nbP*dim; i+=dim )
      CPPUNIT_ASSERT( bb.isOut( p + i ));

    // the extrusions shifted by 3 along OX are OUT
    double shifted_X[nbP*dim] =
      {
        3,0,    5,0,
        13,10, 15,10,
      };
    double shifted_x[nbP*dim] =
      {
        -3,0, -1,0,
        7,10, 9,10,
      };
    INTERP_KERNEL::DirectedBoundingBox shiftedBB_x( shifted_x, nbP, dim);
    INTERP_KERNEL::DirectedBoundingBox shiftedBB_X( shifted_X, nbP, dim);

    CPPUNIT_ASSERT( bb.isDisjointWith( shiftedBB_x ));
    CPPUNIT_ASSERT( bb.isDisjointWith( shiftedBB_X ));

    // intersecting box is IN
    double inters_coords[nbP*dim] =
      {
        0,0,    2,0, 
        0,2,    2,2
      };
    INTERP_KERNEL::DirectedBoundingBox ibb( inters_coords, nbP, dim);
    CPPUNIT_ASSERT( !bb.isDisjointWith( ibb ));

    // overlapping non-directed BB
    double overlappingBB[2*dim] =
      {
        5,6, -5,4
      };
    CPPUNIT_ASSERT( !bb.isDisjointWith( overlappingBB ));

    // non-overlapping non-directed BB
    double nonoverlappingBB_1[2*dim] =
      {
        5,6, -5,2
      };
    CPPUNIT_ASSERT( bb.isDisjointWith( nonoverlappingBB_1 ));
    double nonoverlappingBB_2[2*dim] =
      {
        5,6, 7,20
      };
    CPPUNIT_ASSERT( bb.isDisjointWith( nonoverlappingBB_2 ));
  }

  void BBTreeTest::test_DirectedBB_1D()
  {
    const int nbP = 2, dim = 1;
    double coords[nbP*dim] =
      {
        0, 10
      };
    INTERP_KERNEL::DirectedBoundingBox bb( coords, nbP, dim);
    bb.enlarge( 1e-12 );

    // coords are IN
    for ( int i = 0; i < nbP*dim; i+=dim )
      CPPUNIT_ASSERT( !bb.isOut( coords + i ));

    // points near ends are OUT
    double p[nbP*dim] = 
      {
        -0.0001, 10.1
      };
    for ( int i = 0; i < nbP*dim; i+=dim )
      CPPUNIT_ASSERT( bb.isOut( p + i ));

    // shifted boxes are OUT
    double shifted_X[nbP*dim] =
      {
        10.1, 11
      };
    double shifted_x[nbP*dim] =
      {
        -3.0, -0.001
      };
    INTERP_KERNEL::DirectedBoundingBox shiftedBB_x( shifted_x, nbP, dim);
    INTERP_KERNEL::DirectedBoundingBox shiftedBB_X( shifted_X, nbP, dim);

    CPPUNIT_ASSERT( bb.isDisjointWith( shiftedBB_x ));
    CPPUNIT_ASSERT( bb.isDisjointWith( shiftedBB_X ));

    // intersecting box is IN
    double inters_coords[nbP*dim] =
      {
        -2,2
      };
    INTERP_KERNEL::DirectedBoundingBox ibb( inters_coords, nbP, dim);
    CPPUNIT_ASSERT( !bb.isDisjointWith( ibb ));

    // overlapping non-directed BB
    double overlappingBB[2*dim] =
      {
        -5,4
      };
    CPPUNIT_ASSERT( !bb.isDisjointWith( overlappingBB ));

    // non-overlapping non-directed BB
    double nonoverlappingBB_1[2*dim] =
      {
        -5,-2
      };
    CPPUNIT_ASSERT( bb.isDisjointWith( nonoverlappingBB_1 ));
    double nonoverlappingBB_2[2*dim] =
      {
        11,16
      };
    CPPUNIT_ASSERT( bb.isDisjointWith( nonoverlappingBB_2 ));
  }

}
