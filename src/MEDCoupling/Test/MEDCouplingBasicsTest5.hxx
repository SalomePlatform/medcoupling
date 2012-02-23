// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __MEDCOUPLINGBASICSTEST5_HXX__
#define __MEDCOUPLINGBASICSTEST5_HXX__

#include "MEDCouplingBasicsTest.hxx"

#include <map>
#include <vector>

namespace ParaMEDMEM
{
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingMultiFields;

  class MEDCouplingBasicsTest5 : public MEDCouplingBasicsTest
  {
    CPPUNIT_TEST_SUITE(MEDCouplingBasicsTest5);
    CPPUNIT_TEST( testUMeshTessellate2D1 );
    CPPUNIT_TEST( testIntersect2DMeshesTmp4 );
    CPPUNIT_TEST( testGetCellIdsCrossingPlane1 );
    CPPUNIT_TEST( testBuildSlice3D1 );
    CPPUNIT_TEST( testBuildSlice3DSurf1 );
    CPPUNIT_TEST( testDataArrayDoubleAdvSetting1 );
    CPPUNIT_TEST( testDataArrayIntAdvSetting1 );
    CPPUNIT_TEST_SUITE_END();
  public:
    void testUMeshTessellate2D1();
    void testIntersect2DMeshesTmp4();
    void testGetCellIdsCrossingPlane1();
    void testBuildSlice3D1();
    void testBuildSlice3DSurf1();
    void testDataArrayDoubleAdvSetting1();
    void testDataArrayIntAdvSetting1();
  };
}

#endif
