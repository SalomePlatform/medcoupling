//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "BBTreeTest.hxx"
#include <iostream>
#include <vector>
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


}
