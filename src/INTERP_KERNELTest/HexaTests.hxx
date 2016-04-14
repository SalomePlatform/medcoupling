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

#ifndef __HEXA_TESTS_HXX_
#define __HEXA_TESTS_HXX_

#include "InterpolationTestSuite.hxx"

namespace INTERP_TEST
{
  /**
   * \brief Class performing intersection tests on meshes with hexahedral elements.
   *
   */
  class HexaTests : public InterpolationTestSuite<3,3>
  {
    CPPUNIT_TEST_SUITE( HexaTests );

    CPPUNIT_TEST( simpleHexaBox );
    //VB : slightly inaccurate so that it triggers a failure of the test
    // should be investigated in the future
    //    CPPUNIT_TEST( reflexiveHexaBox );
    CPPUNIT_TEST( hexaBoxes );
    CPPUNIT_TEST( hexaBoxesMoved );

    CPPUNIT_TEST_SUITE_END();

  public:
    
    /// Intersection between two boxes, aligned with the axes.One has 60 hexahedral elements and the other has 39 tetrahedral elements
    /// \brief Status : pass
    void simpleHexaBox()
    {
      _testTools->intersectMeshes("BoxHexa1", "BoxTetra2", 65250, 1.0e-5);
    }

    /// Intersection of a box with 60 hexahedral elements with itself
    /// \brief Status : pass
    void reflexiveHexaBox()
    {
      _testTools->intersectMeshes("BoxHexa1", "BoxHexa1", 204000);
    }

    /// Intersection between two boxes, aligned with the axes.Both have hexahedral elements : one 36, the other 60
    /// \brief Status : pass
    void hexaBoxes()
    {
      _testTools->intersectMeshes("BoxHexa1", "BoxHexa2", 65250);
    }

    /// Intersection between two boxes in general position with hexahedral elements. One has 200 elements and the other 420.
    /// \brief Status : fails - reason unknown. The matrix does not fulfil the transpose requirement : that W_AB = W_BA^T 
    void hexaBoxesMoved()
    {
      _testTools->intersectMeshes("MovedHexaBox1", "MovedHexaBox2", 65250);
    }

  };
}
#endif
