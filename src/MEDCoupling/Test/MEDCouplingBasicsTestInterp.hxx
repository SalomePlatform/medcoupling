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
// Author : Anthony Geay (CEA/DEN)

#ifndef __MEDCOUPLINGBASICSTESTINTERP_HXX__
#define __MEDCOUPLINGBASICSTESTINTERP_HXX__

#include "MEDCouplingBasicsTest.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingMultiFields;

  class MEDCouplingBasicsTestInterp : public MEDCouplingBasicsTest
  {
    CPPUNIT_TEST_SUITE(MEDCouplingBasicsTestInterp);
    CPPUNIT_TEST( test2DInterpP0P0_1 );
    CPPUNIT_TEST( test2DInterpP0P0PL_1 );
    CPPUNIT_TEST( test2DInterpP0P0PL_2 );
    CPPUNIT_TEST( test2DInterpP0P0PL_3 );
    CPPUNIT_TEST( test2DInterpP0P0PL_4 );
    CPPUNIT_TEST( test2DInterpP0P1_1 );
    CPPUNIT_TEST( test2DInterpP0P1PL_1 );
    CPPUNIT_TEST( test2DInterpP0P1PL_2 );
    CPPUNIT_TEST( test2DInterpP1P0_1 );
    CPPUNIT_TEST( test2DInterpP1P0PL_1 );
    CPPUNIT_TEST( test2DInterpP1P0PL_2 );
    CPPUNIT_TEST( test2DInterpP1P1_1 );
    CPPUNIT_TEST( test2DInterpP1P1PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P0_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P0PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P1_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P1PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P0_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P0PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P1_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P1PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P0_2 );
    CPPUNIT_TEST( test3DSurfInterpP0P0_3 );

    CPPUNIT_TEST( testInterpolationCC );
    CPPUNIT_TEST( testInterpolationCU1D );
    CPPUNIT_TEST( testInterpolationCU2D );
    CPPUNIT_TEST( testInterpolationCU3D );

    CPPUNIT_TEST( test3DInterpP0P0_1 );
    CPPUNIT_TEST( test3DInterpP0P0PL_1 );
    CPPUNIT_TEST( test3DInterpP0P0PL_2 );
    CPPUNIT_TEST( test3DInterpP0P0PL_3 );
    CPPUNIT_TEST( test3DInterpP0P0PL_4 );
    CPPUNIT_TEST( test3DInterpP0P1_1 );
    CPPUNIT_TEST( test3DInterpP0P1PL_1 );
    CPPUNIT_TEST( test3DInterpP1P0_1 );
    CPPUNIT_TEST( test3DInterpP1P0PL_1 );
    CPPUNIT_TEST( test3DInterpP1P1_1 );
    CPPUNIT_TEST( test3DInterpP1P1PL_1 );
    CPPUNIT_TEST( test3DInterpP0P0Empty );
    CPPUNIT_TEST( test2DInterpP0IntegralUniform );
    CPPUNIT_TEST( test3DSurfInterpP0IntegralUniform );
    CPPUNIT_TEST( test3DInterpP0IntegralUniform );
    CPPUNIT_TEST( test2DInterpP1IntegralUniform );
    CPPUNIT_TEST( test3DInterpP1IntegralUniform );
    CPPUNIT_TEST( test2DInterpP1P0Bary_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P0Bary_1 );
    CPPUNIT_TEST( test3DInterpP1P0Bary_1 );
    CPPUNIT_TEST( test3DTo1DInterpP0P0PL_1 );

    CPPUNIT_TEST( test2D1DBasicInterpP0P0 );
    CPPUNIT_TEST( test2D1DSegQuadInterpP0P0_1 );
    CPPUNIT_TEST( test2D1DSegQuadInterpP0P0_2 );
    CPPUNIT_TEST( test2D1DSegQuadInterpP0P0_3 );
    CPPUNIT_TEST( test2D1DSegQuadInterpP0P0_4 );
    CPPUNIT_TEST( test2D1DSegQuadInterpP0P0_5 );
    CPPUNIT_TEST( test2D1DSegQuadInterpP0P0_6 );
    CPPUNIT_TEST( test2D1DSegTriInterpP0P0_1 );
    CPPUNIT_TEST( test2D1DSegTriInterpP0P0_2 );
    CPPUNIT_TEST( test2D1DSegTriInterpP0P0_3 );
    CPPUNIT_TEST( test2D1DSegTriInterpP0P0_4 );
    CPPUNIT_TEST( test2D1DSegTriInterpP0P0_5 );
    CPPUNIT_TEST( test2D1DSegTriInterpP0P0_6 );
    CPPUNIT_TEST( test3D2DBasicInterpP0P0 );
    CPPUNIT_TEST( test3D2DQuadHexaInterpP0P0_1 );
    CPPUNIT_TEST( test3D2DQuadHexaInterpP0P0_2 );
    CPPUNIT_TEST( test3D2DQuadHexaInterpP0P0_3 );
    CPPUNIT_TEST( test3D2DQuadHexaInterpP0P0_4 );
    CPPUNIT_TEST( test3D2DQuadHexaInterpP0P0_5 );
    CPPUNIT_TEST( test3D2DQuadHexaInterpP0P0_6 );
    CPPUNIT_TEST( test3D2DTriHexaInterpP0P0_1 );
    CPPUNIT_TEST( test3D2DTriHexaInterpP0P0_2 );
    CPPUNIT_TEST( test3D2DTriHexaInterpP0P0_3 );
    CPPUNIT_TEST( test3D2DTriHexaInterpP0P0_4 );
    CPPUNIT_TEST( test3D2DTriHexaInterpP0P0_5 );
    CPPUNIT_TEST( test3D2DTriHexaInterpP0P0_6 );
    CPPUNIT_TEST( test3D2DQuadTetraInterpP0P0_1 );
    CPPUNIT_TEST( test3D2DQuadTetraInterpP0P0_2 );
    CPPUNIT_TEST( test3D2DQuadTetraInterpP0P0_3 );
    CPPUNIT_TEST( test3D2DQuadTetraInterpP0P0_4 );
    CPPUNIT_TEST( test3D2DQuadTetraInterpP0P0_5 );
    CPPUNIT_TEST( test3D2DQuadTetraInterpP0P0_6 );
    CPPUNIT_TEST( test3D2DTriTetraInterpP0P0_1 );
    CPPUNIT_TEST( test3D2DTriTetraInterpP0P0_2 );
    CPPUNIT_TEST( test3D2DTriTetraInterpP0P0_3 );
    CPPUNIT_TEST( test3D2DTriTetraInterpP0P0_4 );
    CPPUNIT_TEST( test3D2DTriTetraInterpP0P0_5 );
    CPPUNIT_TEST( test3D2DTriTetraInterpP0P0_6 );

    CPPUNIT_TEST( test1DInterp_1 );
    CPPUNIT_TEST( test2DCurveInterpP0P0_1 );
    CPPUNIT_TEST( test2DCurveInterpP0P0_2 );
    CPPUNIT_TEST( test2DCurveInterpP0P1_1 );
    CPPUNIT_TEST( test2DCurveInterpP1P0_1 );
    CPPUNIT_TEST( test2DCurveInterpP1P1_1 );
    CPPUNIT_TEST_SUITE_END();
  public:
    void test2DInterpP0P0_1();
    void test2DInterpP0P0PL_1();
    void test2DInterpP0P0PL_2();
    void test2DInterpP0P0PL_3();
    void test2DInterpP0P0PL_4();
    void test2DInterpP0P1_1();
    void test2DInterpP0P1PL_1();
    void test2DInterpP0P1PL_2();
    void test2DInterpP1P0_1();
    void test2DInterpP1P0PL_1();
    void test2DInterpP1P0PL_2();
    void test2DInterpP1P1_1();
    void test2DInterpP1P1PL_1();
    void test3DSurfInterpP0P0_1();
    void test3DSurfInterpP0P0PL_1();
    void test3DSurfInterpP0P1_1();
    void test3DSurfInterpP0P1PL_1();
    void test3DSurfInterpP1P0_1();
    void test3DSurfInterpP1P0PL_1();
    void test3DSurfInterpP1P1_1();
    void test3DSurfInterpP1P1PL_1();
    void test3DSurfInterpP0P0_2();
    void test3DSurfInterpP0P0_3();
    void test3DInterpP0P0_1();
    void test3DInterpP0P0PL_1();
    void test3DInterpP0P0PL_2();
    void test3DInterpP0P0PL_3();
    void test3DInterpP0P0PL_4();
    void test3DInterpP0P1_1();
    void test3DInterpP0P1PL_1();
    void test3DInterpP1P0_1();
    void test3DInterpP1P0PL_1();
    void test3DInterpP1P1_1();
    void test3DInterpP1P1PL_1();

    void testInterpolationCC();
    void testInterpolationCU1D();
    void testInterpolationCU2D();
    void testInterpolationCU3D();

    void test3DInterpP0P0Empty();
    void test2DInterpP0IntegralUniform();
    void test3DSurfInterpP0IntegralUniform();
    void test3DInterpP0IntegralUniform();
    void test2DInterpP1IntegralUniform();
    void test3DInterpP1IntegralUniform();
    void test2DInterpP1P0Bary_1();
    void test3DSurfInterpP1P0Bary_1();
    void test3DInterpP1P0Bary_1();
    void test3DTo1DInterpP0P0PL_1();

    void test2D1DBasicInterpP0P0();
    void test2D1DSegQuadInterpP0P0_1();
    void test2D1DSegQuadInterpP0P0_2();
    void test2D1DSegQuadInterpP0P0_3();
    void test2D1DSegQuadInterpP0P0_4();
    void test2D1DSegQuadInterpP0P0_5();
    void test2D1DSegQuadInterpP0P0_6();
    void test2D1DSegTriInterpP0P0_1();
    void test2D1DSegTriInterpP0P0_2();
    void test2D1DSegTriInterpP0P0_3();
    void test2D1DSegTriInterpP0P0_4();
    void test2D1DSegTriInterpP0P0_5();
    void test2D1DSegTriInterpP0P0_6();
    void test3D2DBasicInterpP0P0();
    void test3D2DQuadHexaInterpP0P0_1();
    void test3D2DQuadHexaInterpP0P0_2();
    void test3D2DQuadHexaInterpP0P0_3();
    void test3D2DQuadHexaInterpP0P0_4();
    void test3D2DQuadHexaInterpP0P0_5();
    void test3D2DQuadHexaInterpP0P0_6();
    void test3D2DTriHexaInterpP0P0_1();
    void test3D2DTriHexaInterpP0P0_2();
    void test3D2DTriHexaInterpP0P0_3();
    void test3D2DTriHexaInterpP0P0_4();
    void test3D2DTriHexaInterpP0P0_5();
    void test3D2DTriHexaInterpP0P0_6();
    void test3D2DQuadTetraInterpP0P0_1();
    void test3D2DQuadTetraInterpP0P0_2();
    void test3D2DQuadTetraInterpP0P0_3();
    void test3D2DQuadTetraInterpP0P0_4();
    void test3D2DQuadTetraInterpP0P0_5();
    void test3D2DQuadTetraInterpP0P0_6();
    void test3D2DTriTetraInterpP0P0_1();
    void test3D2DTriTetraInterpP0P0_2();
    void test3D2DTriTetraInterpP0P0_3();
    void test3D2DTriTetraInterpP0P0_4();
    void test3D2DTriTetraInterpP0P0_5();
    void test3D2DTriTetraInterpP0P0_6();

    void test1DInterp_1();
    void test2DCurveInterpP0P0_1();
    void test2DCurveInterpP0P0_2();
    void test2DCurveInterpP0P1_1();
    void test2DCurveInterpP1P0_1();
    void test2DCurveInterpP1P1_1();
  };
}

#endif
