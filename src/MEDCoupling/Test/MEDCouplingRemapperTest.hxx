//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef __MEDCOUPLINGREMAPPERTEST_HXX__
#define __MEDCOUPLINGREMAPPERTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include <map>
#include <vector>

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;

  class MEDCouplingRemapperTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(MEDCouplingRemapperTest);
    CPPUNIT_TEST( test2DInterpP0P0_1 );
    CPPUNIT_TEST( test2DInterpP0P0R_1 );
    CPPUNIT_TEST( test2DInterpMultiMethods );
    CPPUNIT_TEST( testMultiDimCombi );
    CPPUNIT_TEST( testNatureOfField );
    CPPUNIT_TEST( testExtruded );
    CPPUNIT_TEST_SUITE_END();
  public:
    void test2DInterpP0P0_1();
    void test2DInterpP0P0R_1();
    void test2DInterpMultiMethods();
    void testMultiDimCombi();
    void testNatureOfField();
    void testExtruded();
  private:
    MEDCouplingUMesh *build3DSourceMesh_2();
    MEDCouplingUMesh *build3DTargetMesh_2();
    MEDCouplingUMesh *build1DTargetMesh_1();
    MEDCouplingUMesh *build1DTargetMesh_2();
    MEDCouplingUMesh *build2DSourceMesh_1();
    MEDCouplingUMesh *build2DTargetMesh_1();
    MEDCouplingUMesh *build2DTargetMeshPerm_1();
    MEDCouplingUMesh *build2DTargetMesh_2();
    MEDCouplingUMesh *build2DTargetMesh_3();
    MEDCouplingUMesh *build3DSurfSourceMesh_1();
    MEDCouplingUMesh *build3DSurfSourceMesh_2();
    MEDCouplingUMesh *build3DSurfTargetMesh_1();
    MEDCouplingUMesh *build3DSurfTargetMeshPerm_1();
    MEDCouplingUMesh *build3DSurfTargetMesh_2();
    MEDCouplingUMesh *build3DSourceMesh_1();
    MEDCouplingUMesh *build3DTargetMesh_1();
    MEDCouplingUMesh *build3DExtrudedUMesh_1(MEDCouplingUMesh *&mesh2D);
    double sumAll(const std::vector< std::map<int,double> >& matrix);
  };
}

#endif
