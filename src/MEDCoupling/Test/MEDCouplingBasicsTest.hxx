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

#ifndef __MEDCOUPLINGBASICSTEST_HXX__
#define __MEDCOUPLINGBASICSTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include <map>
#include <vector>

namespace MEDCoupling
{
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingMultiFields;

  class MEDCouplingBasicsTest : public CppUnit::TestFixture
  {
  public:
    static MEDCouplingUMesh *build3DSourceMesh_2();
    static MEDCouplingUMesh *build3DTargetMesh_2();
    static MEDCouplingUMesh *build1DTargetMesh_1();
    static MEDCouplingUMesh *build2DSourceMesh_1();
    static MEDCouplingUMesh *build2DTargetMesh_1();
    static MEDCouplingUMesh *build2DTargetMeshPerm_1();
    static MEDCouplingUMesh *build2DTargetMesh_2();
    static MEDCouplingUMesh *buildCU1DMesh_U();
    static MEDCouplingUMesh *buildCU2DMesh_U();
    static MEDCouplingUMesh *buildCU3DMesh_U();
    static MEDCouplingUMesh *build3DSurfSourceMesh_1();
    static MEDCouplingUMesh *build3DSurfSourceMesh_2();
    static MEDCouplingUMesh *build3DSurfTargetMesh_1();
    static MEDCouplingUMesh *build3DSurfTargetMeshPerm_1();
    static MEDCouplingUMesh *build3DSurfTargetMesh_2();
    static MEDCouplingUMesh *build3DSourceMesh_1();
    static MEDCouplingUMesh *build3DTargetMesh_1();
    static MEDCouplingUMesh *build2DTargetMeshMergeNode_1();
    static MEDCouplingUMesh *build3DTargetMeshMergeNode_1();
    static MEDCouplingUMesh *build3DExtrudedUMesh_1(MEDCouplingUMesh *&mesh2D);
    static void build3DExtrudedUMesh_2(MEDCouplingUMesh *&meshN, MEDCouplingUMesh *&meshTT, MEDCouplingUMesh *&meshTF);
    static MEDCouplingUMesh *build2DTargetMeshMerged_1();
    static MEDCouplingUMesh *build2DCurveMesh(double dx, double dy);
    static MEDCouplingUMesh *build1DMesh(double dx);
    static MEDCouplingUMesh *build1DSourceMesh_2();
    static MEDCouplingUMesh *build1DTargetMesh_2();
    static MEDCouplingUMesh *build2DCurveSourceMesh_2();
    static MEDCouplingUMesh *build2DCurveTargetMesh_2();
    static MEDCouplingUMesh *build1DTargetMesh_3();
    static MEDCouplingUMesh *build2DCurveTargetMesh_3();
    static MEDCouplingUMesh *build2DTargetMesh_3();
    static MEDCouplingUMesh *build3DTargetMesh_3();
    static MEDCouplingUMesh *build2DTargetMesh_4();
    static MEDCouplingUMesh *build1DMultiTypes_1();
    static MEDCouplingUMesh *build2DMultiTypes_1();
    static MEDCouplingUMesh *build3DMultiTypes_1();
    static MEDCouplingUMesh *buildHexa8Mesh_1();
    static MEDCouplingUMesh *buildPointe_1(MEDCouplingUMesh *&m1);

    static MEDCouplingUMesh *build2D1DSourceMesh();
    static MEDCouplingUMesh *build2D1DTargetMesh();
    static MEDCouplingUMesh *build2D1DSegSourceMesh(const double shiftX = 0.,
                                                    const double inclinationX = 0.);
    static MEDCouplingUMesh *build2D1DQuadTargetMesh(const double inclinaisonX = 0.);
    static MEDCouplingUMesh *build2D1DTriTargetMesh(const double inclinaisonX = 0.);
    static MEDCouplingUMesh *build3D2DSourceMesh();
    static MEDCouplingUMesh *build3D2DTargetMesh();
    static MEDCouplingUMesh* build3D2DQuadSourceMesh(const double shiftX = 0.,
                                                     const double inclinationX = 0.);
    static MEDCouplingUMesh* build3D2DTriSourceMesh(const double shiftX = 0.,
                                                    const double inclinationX = 0.);
    static MEDCouplingUMesh* build3D2DTetraTargetMesh(const double inclinaisonX = 0.);
    static MEDCouplingUMesh* build3D2DHexaTargetMesh(const double inclinaisonX = 0.);

    static DataArrayDouble *buildCoordsForMultiTypes_1();
    static MEDCouplingMultiFields *buildMultiFields_1();
    static std::vector<MEDCouplingFieldDouble *> buildMultiFields_2();
    static double sumAll(const std::vector< std::map<int,double> >& matrix);
  protected:
    static int countNonZero(const std::vector< std::map<int,double> >& matrix);

    static void test2D1DMeshesIntersection(MEDCouplingUMesh *sourceMesh,
                                           MEDCouplingUMesh *targetMesh,
                                           const double correctSurf,
                                           const int correctDuplicateFacesNbr,
                                           const int correctTotalIntersectFacesNbr = -1);
    static void test3D2DMeshesIntersection(MEDCouplingUMesh *sourceMesh,
                                           MEDCouplingUMesh *targetMesh,
                                           const double correctSurf,
                                           const int correctDuplicateFacesNbr,
                                           const int correctTotalIntersectFacesNbr = -1);
  };
}

#endif
