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

#include "UnitTetra3D2DIntersectionTest.hxx"

#include "TetraAffineTransform.hxx"
#include "InterpolationUtils.hxx"
#include "SplitterTetra.txx"

#include <iostream>

using namespace INTERP_KERNEL;

namespace INTERP_TEST
{
  struct __MESH_DUMMY
  {
    typedef int MyConnType;
    static const int MY_SPACEDIM=3;
  };

  static SplitterTetra<__MESH_DUMMY>* buildSplitterTetra()
  {
    const int conn[4] = { 0,1,2,3 };

    const double targetCoords[] = { -20., 0.,10.,
                                    -20.,10.,10.,
                                    -12., 0.,10.,
                                    -20., 0.,18. };

    const double* tetraCoords[]={ targetCoords, targetCoords+3, targetCoords+6, targetCoords+9 };

    __MESH_DUMMY dummyMesh;
    SplitterTetra<__MESH_DUMMY>* targetTetra = new SplitterTetra<__MESH_DUMMY>( dummyMesh, tetraCoords, conn );
    return targetTetra;
  }

  void UnitTetra3D2DIntersectionTest::test_UnitTetra3D2DIntersection_1()
  {
    const int conn[4] = { 0,1,2 };

    const double sourceCoords[] = { -20., 0., 10.,
                                    -12., 0., 10.,
                                    -20.,10., 10. };

    SplitterTetra<__MESH_DUMMY>* targetTetra = buildSplitterTetra();
    const double dimCaracteristic = 1.;
    const double precision = 1.e-12;
    std::multiset<TriangleFaceKey> listOfTetraFacesTreated;
    std::set<TriangleFaceKey> listOfTetraFacesColinear;

    const double* sourceTriCoords[] = { sourceCoords, sourceCoords+3, sourceCoords+6 };
    double surface = targetTetra->intersectSourceFace(NORM_TRI3,
                                                      3,
                                                      conn,
                                                      sourceTriCoords,
                                                      dimCaracteristic,
                                                      precision,
                                                      listOfTetraFacesTreated,
                                                      listOfTetraFacesColinear);
    delete targetTetra;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(40.,surface,precision);

    CPPUNIT_ASSERT_EQUAL(4,(int)listOfTetraFacesTreated.size());
    std::multiset<TriangleFaceKey> correctListOfTetraFacesTreated;
    TriangleFaceKey key1 = TriangleFaceKey(0, 1, 2);
    correctListOfTetraFacesTreated.insert(key1);
    TriangleFaceKey key2 = TriangleFaceKey(0, 1, 3);
    correctListOfTetraFacesTreated.insert(key2);
    TriangleFaceKey key3 = TriangleFaceKey(0, 2, 3);
    correctListOfTetraFacesTreated.insert(key3);
    TriangleFaceKey key4 = TriangleFaceKey(1, 2, 3);
    correctListOfTetraFacesTreated.insert(key4);
    CPPUNIT_ASSERT(correctListOfTetraFacesTreated == listOfTetraFacesTreated);

    CPPUNIT_ASSERT_EQUAL(1,(int)listOfTetraFacesColinear.size());
    std::set<TriangleFaceKey> correctListOfTetraFacesColinear;
    correctListOfTetraFacesColinear.insert(key1);
    CPPUNIT_ASSERT(correctListOfTetraFacesColinear == listOfTetraFacesColinear);

  }

  void UnitTetra3D2DIntersectionTest::test_UnitTetra3D2DIntersection_2()
  {
    const int conn[4] = { 0,1,2,3 };

    const double sourceCoords[] = { -20., 0., 10.,
                                    -12., 0., 10.,
                                    -12.,10., 10.,
                                    -20.,10., 10. };

    SplitterTetra<__MESH_DUMMY>* targetTetra = buildSplitterTetra();
    const double dimCaracteristic = 1.;
    const double precision = 1.e-12;
    std::multiset<TriangleFaceKey> listOfTetraFacesTreated;
    std::set<TriangleFaceKey> listOfTetraFacesColinear;

    const double* sourceQuadCoords[] = { sourceCoords, sourceCoords+3, sourceCoords+6, sourceCoords+9 };
    double surface = targetTetra->intersectSourceFace(NORM_QUAD4,
                                                      4,
                                                      conn,
                                                      sourceQuadCoords,
                                                      dimCaracteristic,
                                                      precision,
                                                      listOfTetraFacesTreated,
                                                      listOfTetraFacesColinear);
    delete targetTetra;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(40.,surface,precision);

    CPPUNIT_ASSERT_EQUAL(4,(int)listOfTetraFacesTreated.size());
    std::multiset<TriangleFaceKey> correctListOfTetraFacesTreated;
    TriangleFaceKey key1 = TriangleFaceKey(0, 1, 2);
    correctListOfTetraFacesTreated.insert(key1);
    TriangleFaceKey key2 = TriangleFaceKey(0, 1, 3);
    correctListOfTetraFacesTreated.insert(key2);
    TriangleFaceKey key3 = TriangleFaceKey(0, 2, 3);
    correctListOfTetraFacesTreated.insert(key3);
    TriangleFaceKey key4 = TriangleFaceKey(1, 2, 3);
    correctListOfTetraFacesTreated.insert(key4);
    CPPUNIT_ASSERT(correctListOfTetraFacesTreated == listOfTetraFacesTreated);

    CPPUNIT_ASSERT_EQUAL(1,(int)listOfTetraFacesColinear.size());
    std::set<TriangleFaceKey> correctListOfTetraFacesColinear;
    correctListOfTetraFacesColinear.insert(key1);
    CPPUNIT_ASSERT(correctListOfTetraFacesColinear == listOfTetraFacesColinear);

 }

  void UnitTetra3D2DIntersectionTest::test_UnitTetra3D2DIntersection_3()
  {
    const int conn[4] = { 0,1,2 };

    const double sourceCoords[] = { -20., 0., 16.,
                                    -18., 0., 16.,
                                    -20.,2.5, 16. };

    SplitterTetra<__MESH_DUMMY>* targetTetra = buildSplitterTetra();
    const double dimCaracteristic = 1.;
    const double precision = 1.e-12;
    std::multiset<TriangleFaceKey> listOfTetraFacesTreated;
    std::set<TriangleFaceKey> listOfTetraFacesColinear;

    const double* sourceTri2Coords[] = { sourceCoords, sourceCoords+3, sourceCoords+6 };
    double surface = targetTetra->intersectSourceFace(NORM_TRI3,
                                                      3,
                                                      conn,
                                                      sourceTri2Coords,
                                                      dimCaracteristic,
                                                      precision,
                                                      listOfTetraFacesTreated,
                                                      listOfTetraFacesColinear);
    delete targetTetra;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5,surface,precision);

    CPPUNIT_ASSERT_EQUAL(0,(int)listOfTetraFacesTreated.size());

    CPPUNIT_ASSERT_EQUAL(0,(int)listOfTetraFacesColinear.size());
 }

}
