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

#include "MEDCouplingBasicsTestInterp.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3DSurf.hxx"
#include "Interpolation3D.txx"
#include "Interpolation2D1D.txx"
#include "Interpolation2D3D.txx"
#include "InterpolationCC.txx"
#include "InterpolationCU.txx"
#include "Interpolation2DCurve.hxx"
#include "Interpolation1D.txx"

#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "MEDCouplingNormalizedCartesianMesh.txx"

#include <cmath>
#include <functional>

using namespace MEDCoupling;

typedef std::vector<std::map<int,double> > IntersectionMatrix;

void MEDCouplingBasicsTestInterp::test2DInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[3]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Convex, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<3;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP0P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  //
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP0P0PL_2()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int> cellsIds(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  //
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP0P0PL_3()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  //
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP0P0PL_4()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  cellsIds.resize(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  //
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP0P1_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
      CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333329,res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[3][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[5][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333329,res[6][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[7][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[8][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[8][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25,sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP0P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP0P1PL_2()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  cellsIds.resize(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP1P0_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[3][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333333,res[1][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333333,res[2][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666667,res[3][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[2][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[3][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[4][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP1P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[1][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[2][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[4][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP1P0PL_2()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int >cellsIds(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[1][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[2][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[4][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP1P1_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_2();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
      CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333334,res[0][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05416666666666665,res[1][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666666,res[1][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333334,res[2][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05416666666666665,res[3][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666668,res[3][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1416666666666666,res[4][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02499999999999999,res[4][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02499999999999999,res[4][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09999999999999999,res[4][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666666,res[5][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09583333333333333,res[5][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333333,res[6][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666667,res[7][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09583333333333331,res[7][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.04166666666666668,res[8][3],1.e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP1P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[0][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[2][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,res[4][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,res[4][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[8][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(25.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[3]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP0P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP0P1_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
      CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333329*sqrt(2.),res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[3][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[5][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333329*sqrt(2.),res[6][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[7][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[8][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[8][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25*sqrt(2.),sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP0P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP1P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[3][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333333*sqrt(2.),res[1][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333333*sqrt(2.),res[2][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666667*sqrt(2.),res[3][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[2][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[3][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[4][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP1P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[1][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[2][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[4][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP1P1_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_2();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
      CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333334*sqrt(2.),res[0][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05416666666666665*sqrt(2.),res[1][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666666*sqrt(2.),res[1][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333334*sqrt(2.),res[2][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05416666666666665*sqrt(2.),res[3][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666668*sqrt(2.),res[3][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1416666666666666*sqrt(2.),res[4][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02499999999999999*sqrt(2.),res[4][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02499999999999999*sqrt(2.),res[4][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09999999999999999*sqrt(2.),res[4][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666666*sqrt(2.),res[5][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09583333333333333*sqrt(2.),res[5][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333333*sqrt(2.),res[6][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666667*sqrt(2.),res[7][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09583333333333331*sqrt(2.),res[7][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.04166666666666668*sqrt(2.),res[8][3],1.e-12);
      res.clear();
    }
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP1P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[0][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[2][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,res[4][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,res[4][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[8][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(25.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP0P0_2()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMeshPerm_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::Triangulation);
  {
    myInterpolator.setOrientation(2);
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
    res.clear();
  }
  {
    myInterpolator.setOrientation(0);
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.125*sqrt(2.),res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75*sqrt(2.),sumAll(res),1e-12);
    res.clear();
  }
  {
    myInterpolator.setOrientation(1);
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.875*sqrt(2.),sumAll(res),1e-12);
    res.clear();
  }
  {
    myInterpolator.setOrientation(-1);
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),sumAll(res),1e-12);
    res.clear();
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

/*!
 * Test of precision option implemented by Fabien that represents distance of "barycenter" to the other cell.
 */
void MEDCouplingBasicsTestInterp::test3DSurfInterpP0P0_3()
{
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  double vecTrans[3]={0.,0.,1.e-10};
  double vec[3]={0.,-1.,0.};
  double pt[3]={-0.3,-0.3,5.e-11};
  const int N=32;
  const double deltaA=M_PI/N;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::Triangulation);
  myInterpolator.setMaxDistance3DSurfIntersect(1e-9);
  for(int i=0;i<N;i++)
    {
      res.clear();
      MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_2();
      sourceMesh->rotate(pt,vec,i*deltaA);
      MEDCouplingUMesh *targetMesh=build3DSurfSourceMesh_2();
      targetMesh->translate(vecTrans);
      targetMesh->rotate(pt,vec,i*deltaA);
      MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
      MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
      CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
      sourceMesh->decrRef();
      targetMesh->decrRef();
    }
  //
  myInterpolator.setMaxDistance3DSurfIntersect(1e-11);
  for(int i=0;i<N;i++)
    {
      res.clear();
      MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_2();
      sourceMesh->rotate(pt,vec,i*deltaA);
      MEDCouplingUMesh *targetMesh=build3DSurfSourceMesh_2();
      targetMesh->translate(vecTrans);
      targetMesh->rotate(pt,vec,i*deltaA);
      MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
      MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
      CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,sumAll(res),1e-12);
      sourceMesh->decrRef();
      targetMesh->decrRef();
    }
  //
  res.clear();
  myInterpolator.setMaxDistance3DSurfIntersect(-1.);//unactivate fabien lookup
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_2();
  MEDCouplingUMesh *targetMesh=build3DSurfSourceMesh_2();
  targetMesh->translate(vecTrans);
  myInterpolator.setBoundingBoxAdjustment(1e-11);
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper0(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper0(targetMesh);
  myInterpolator.interpolateMeshes(sourceWrapper0,targetWrapper0,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,sumAll(res),1e-12);
  sourceMesh->decrRef();
  targetMesh->decrRef();
  //
  res.clear();
  sourceMesh=build3DSurfSourceMesh_2();
  targetMesh=build3DSurfSourceMesh_2();
  targetMesh->translate(vecTrans);
  myInterpolator.setBoundingBoxAdjustment(1e-9);
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper1(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper1(targetMesh);
  myInterpolator.interpolateMeshes(sourceWrapper1,targetWrapper1,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  sourceMesh->decrRef();
  targetMesh->decrRef();
  //keeping the same bbox adj == 1.e-11 but trying rotation
  res.clear();
  sourceMesh=build3DSurfSourceMesh_2();
  sourceMesh->rotate(pt,vec,M_PI/4.);
  targetMesh=build3DSurfSourceMesh_2();
  targetMesh->translate(vecTrans);
  targetMesh->rotate(pt,vec,M_PI/4.);
  myInterpolator.setBoundingBoxAdjustment(1e-11);
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper2(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper2(targetMesh);
  myInterpolator.interpolateMeshes(sourceWrapper2,targetWrapper2,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.e6,sumAll(res),1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[0][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(41666.66666666667,res[0][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[0][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[0][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[0][10],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(41666.66666666667,res[1][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[1][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[1][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[2][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[2][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[2][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[2][9],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[2][11],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(395833.3333333333,res[3][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[3][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333331,res[3][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[3][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(395833.3333333333,res[3][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[4][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[4][4],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[4][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[4][9],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[4][10],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[5][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333331,res[5][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[5][4],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(395833.3333333333,res[5][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(395833.3333333333,res[5][10],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[6][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(250000,res[6][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(541666.6666666667,res[6][9],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[6][11],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333331,res[7][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(479166.6666666667,res[7][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(333333.3333333333,res[7][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(624999.9999999997,res[7][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(479166.6666666667,res[7][4],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(479166.6666666667,res[7][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333333,res[7][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333331,res[7][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333333,res[7][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333333,res[7][9],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333331,res[7][10],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(479166.6666666667,res[7][11],1e-7);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][9],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0P0PL_2()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  std::vector<int> cellsIds(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][9],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0P0PL_3()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][9],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0P0PL_4()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  cellsIds.resize(targetMesh->getNumberOfCells());
  for(int j=0;j<targetMesh->getNumberOfCells();j++)
    cellsIds[j]=j;
  targetMesh->convertToPolyTypes(&cellsIds[0],&cellsIds[0]+cellsIds.size());
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][10],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][9],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0P1_1()
{
  MEDCouplingUMesh *sourceMesh=build3DTargetMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSourceMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
    CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(244444.4444444445,res[0][4],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[0][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(291666.6666666666,res[0][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[0][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[1][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(140277.7777777778,res[1][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(119444.4444444444,res[1][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[1][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(119444.4444444444,res[1][4],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[1][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(26388.88888888889,res[1][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(348611.1111111111,res[2][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888888,res[2][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(244444.4444444444,res[3][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333334,res[3][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(291666.6666666666,res[3][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[3][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(536111.111111111,res[4][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(297222.2222222221,res[4][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(223611.1111111111,res[5][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[5][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[5][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(26388.88888888892,res[5][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[6][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(536111.1111111109,res[7][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(297222.2222222221,res[7][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.1111111111,res[8][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.11111111111,res[8][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666666,res[8][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.11111111111,res[8][4],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[8][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[8][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1466666.666666668,res[8][7],1e-7);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DTargetMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSourceMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
    CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][2],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.,sumAll(res),1e-12);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP1P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
    CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[0][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(140277.7777777778,res[1][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(223611.1111111111,res[1][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.1111111111,res[1][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(119444.4444444444,res[2][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(244444.4444444445,res[2][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.11111111111,res[2][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[3][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[3][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[3][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(536111.1111111109,res[3][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[3][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(244444.4444444445,res[4][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(119444.4444444445,res[4][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.11111111111,res[4][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[5][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[5][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(536111.1111111109,res[5][4],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[5][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666666,res[5][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(291666.6666666666,res[6][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(26388.88888888889,res[6][1],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(348611.1111111112,res[6][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(291666.6666666667,res[6][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666666,res[6][8],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[7][0],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[7][2],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[7][3],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(297222.2222222221,res[7][4],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(26388.88888888892,res[7][5],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[7][6],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(297222.2222222222,res[7][7],1e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1466666.666666668,res[7][8],1e-7);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP1P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( int i = 0; i < 4; ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    res.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
    CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.75,res[0][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25,res[0][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][7],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[5][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[5][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[6][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[6][2],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[6][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[6][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25,res[7][6],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.75,res[7][8],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP1P1_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_2();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_2();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  double res3D[8][28]= {{124999.999883775978, 245370.370390364464, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 203703.703634892299, 187500.000094145857, 0.0, 0.0, 4629.6296266718, 0.0, 215277.777751402784, 209722.222322299582, 0.0, 0.0, 0.0, 0.0, 104166.666590829205, 121296.296368812196, 0.0, 250000.000003472145},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120370.370368827047, 0.0, 0.0, 38888.888897777797, 0.0, 0.0, 45370.3703701697596, 0.0, 0.0, 45370.3703701697596, 83333.3333263888926, 0.0},
                        {0.0, 0.0, 0.0, 97222.2222222221753, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 97222.2222222221608, 0.0, 97222.2222222222044, 41666.6666666666642, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 277777.777787084982, 199074.074074073927, 0.0, 0.0, 0.0, 4629.62962962962774, 0.0, 321759.259254934732, 83333.3333333333139, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4629.62962667180363, 0.0, 0.0, 251388.88888319055, 194444.444454861077, 0.0, 79629.6296194135939, 250000.000003472145, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 85185.1851851851534, 4629.62962962962774, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 118518.518518518511, 0.0, 41666.6666666666642, 83333.3333333333285, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 324074.07407629228, 0.0, 0.0, 0.0, 247685.185185184964, 6481.48148148147993, 0.0, 173611.11111196311, 0.0, 164814.814814814832, 0.0, 4629.62962962962865, 208333.33333418527, 0.0, 83333.3333333333285, 203703.703697273799, 249999.999999999767, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {125000.000000000015, 423611.111111110775, 134259.259259259241, 194444.444444444351, 164814.814814814745, 164351.851851851825, 203703.703703703592, 249999.999999999825, 0.0, 0.0, 0.0, 0.0, 6481.48148148147902, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 118518.518518518453, 0.0, 4629.62962962962956, 83333.3333333333139, 85185.1851851851825, 41666.6666666666642, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  int i=0;
  double sum = 0;
  //cout.precision(18);
  for(std::vector<std::map<int,double> >::const_iterator iter1=res.begin();iter1!=res.end();iter1++,i++)
    {
      //cout<< "res3D[" <<i<< "][]={";
      for(int j=0;j<28;j++)
        {
          std::map<int,double>::const_iterator iter2=(*iter1).find(j);
          if(iter2!=(*iter1).end())
            {
              //cout<< iter2->second<< ", ";
              sum += iter2->second;
              CPPUNIT_ASSERT_DOUBLES_EQUAL(res3D[i][j],(*iter2).second,1.e-5);
            }
          else
            {
              //cout << "0.0, ";
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res3D[i][j],1e-14);
            }
        }
      //cout << "}" << endl;
    }
  //cout << "Sum = " << sum << endl;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000,sum,1.e-5);
  //clean-up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP1P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_2();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_2();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20.,res[0][24],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[1][26],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][21],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.,res[3][23],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][14],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.,res[5][17],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.,res[6][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(97.,sumAll(res),1e-12);
  //clean-up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0P0Empty()
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(2);
  sourceMesh->allocateCells(0);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(0,0);
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(0);
  targetMesh->finishInsertingCells();
  myCoords=DataArrayDouble::New();
  myCoords->alloc(0,2);
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::testInterpolationCC()
{
  double arr1[3] = { 0/2., 1/2., 2/2. };
  double arr2[4] = { 0/3, 1/3., 2/3., 3/3. };
  MEDCouplingCMesh* mesh[2];
  for ( int i = 0; i < 2; ++i )
    {
      const double* arr = i ? arr1 : arr2;
      const int nb_coord = i ? 3 : 4;
      DataArrayDouble* coords = DataArrayDouble::New();
      coords->useArray( arr, /*ownership=*/false, CPP_DEALLOC, nb_coord, 1 );

      mesh[i] = MEDCouplingCMesh::New();
      mesh[i]->setCoords( coords, coords, coords );
      coords->decrRef();
    }
  MEDCouplingNormalizedCartesianMesh<3> targetWrapper(mesh[1]);
  MEDCouplingNormalizedCartesianMesh<3> sourceWrapper(mesh[0]);
  CPPUNIT_ASSERT_EQUAL( 27,int( sourceWrapper.getNumberOfElements()));
  CPPUNIT_ASSERT_EQUAL( 3, int( sourceWrapper.nbCellsAlongAxis(0)));
  CPPUNIT_ASSERT_EQUAL( 3, int( sourceWrapper.nbCellsAlongAxis(1)));
  CPPUNIT_ASSERT_EQUAL( 3, int( sourceWrapper.nbCellsAlongAxis(2)));
  CPPUNIT_ASSERT_THROW( sourceWrapper.nbCellsAlongAxis(3), INTERP_KERNEL::Exception);

  INTERP_KERNEL::InterpolationCC myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");

  CPPUNIT_ASSERT_EQUAL(8,int( res.size()));
  CPPUNIT_ASSERT_EQUAL(8,int( res[0].size()));
  const double precis = 1e-7;
  std::set<double> vals;
  double sum = 0;
  for ( int i = 0; i < (int)res.size(); ++i )
    for ( std::map<int,double>::iterator s_v = res[i].begin(); s_v != res[i].end(); ++s_v)
      {
        sum += s_v->second;
        double vvv;
#ifdef WIN32
        double vv = s_v->second / precis;
        if(vv>=0.0)
          {
            vvv = floor(vv+0.5);
          }
        else
          {
            vvv = ceil(vv-0.5);
          }
#else
        vvv = round( s_v->second / precis );
#endif
        vals.insert( precis * vvv );
      }
  //cout << "tgt: " << i << " src: " << s_v->first << " - w: " << s_v->second << endl;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, sum, precis );

  std::set<double>::iterator v = vals.begin();
  CPPUNIT_ASSERT_EQUAL( 4, int( vals.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00462963, *v++, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00925926, *v++, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.01851850, *v++, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703700, *v++, precis );

  mesh[0]->decrRef();
  mesh[1]->decrRef();
}

void MEDCouplingBasicsTestInterp::testInterpolationCU1D()
{
  MEDCouplingCMesh* meshC = MEDCouplingCMesh::New();
  DataArrayDouble* coords = DataArrayDouble::New();
  double arr[4] = { -1/3., 1/3., 2/3., 4/3. };
  coords->useArray( arr, /*ownership=*/false, CPP_DEALLOC, 4, 1 );
  meshC->setCoords( coords );
  coords->decrRef();

  MEDCouplingUMesh * meshU = buildCU1DMesh_U();

  MEDCouplingNormalizedCartesianMesh<1>      sourceWrapper(meshC);
  MEDCouplingNormalizedUnstructuredMesh<1,1> targetWrapper(meshU);
  INTERP_KERNEL::InterpolationCU myInterpolator;
  std::vector<std::map<int,double> > res;
  const double precis = 1e-13;
  myInterpolator.setPrecision(precis);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");

//   std::cout.precision(18);
//   for ( int i = 0; i < (int)res.size(); ++i )
//     for ( std::map<int,double>::iterator s_v = res[i].begin(); s_v != res[i].end(); ++s_v)
//     {
//       std::cout << "CPPUNIT_ASSERT_DOUBLES_EQUAL( "<<s_v->second<<" ,res["<<i<<"]["<<s_v->first<<"],precis);"<<std::endl;
//     }

  double sum = sumAll(res);
  CPPUNIT_ASSERT_EQUAL(3,int( res.size()));
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, sum, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[1][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.33333333333333 ,res[1][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.08333333333333 ,res[1][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25000000000000 ,res[2][2],precis);

  meshC->decrRef();
  meshU->decrRef();
}

void MEDCouplingBasicsTestInterp::testInterpolationCU2D()
{
  MEDCouplingCMesh* meshC = MEDCouplingCMesh::New();
  DataArrayDouble* coords = DataArrayDouble::New();
  double arr[4] = { -1/3., 1/3., 2/3., 4/3. };
  coords->useArray( arr, /*ownership=*/false, CPP_DEALLOC, 4, 1 );
  meshC->setCoords( coords, coords );
  coords->decrRef();

  MEDCouplingUMesh * meshU = buildCU2DMesh_U();

  MEDCouplingNormalizedCartesianMesh<2>      sourceWrapper(meshC);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(meshU);
  INTERP_KERNEL::InterpolationCU myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");

  const double precis = 1e-7;
  double sum = sumAll(res);
  CPPUNIT_ASSERT_EQUAL(5,int( res.size()));
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, sum, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1111111 ,res[0][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[0][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[0][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0277778 ,res[0][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[1][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0277778 ,res[1][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1111111 ,res[1][6],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[1][7],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0277778 ,res[2][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[2][5],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[2][7],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1111111 ,res[2][8],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0416667 ,res[3][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0138889 ,res[3][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0277778 ,res[3][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0416667 ,res[3][5],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0138889 ,res[4][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0972222 ,res[4][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0138889 ,res[4][5],precis);

  std::vector<std::map<int,double> > resRev;
  myInterpolator.interpolateMeshesRev(targetWrapper,sourceWrapper,resRev,"P0P0");

  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[0][0] ,resRev[0][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[0][1] ,resRev[1][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[3][1] ,resRev[1][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[4][1] ,resRev[1][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[3][2] ,resRev[2][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[4][2] ,resRev[2][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[0][3] ,resRev[3][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[1][3] ,resRev[3][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[0][4] ,resRev[4][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[1][4] ,resRev[4][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[2][4] ,resRev[4][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[3][4] ,resRev[4][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[2][5] ,resRev[5][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[3][5] ,resRev[5][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[4][5] ,resRev[5][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[1][6] ,resRev[6][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[1][7] ,resRev[7][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[2][7] ,resRev[7][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( res[2][8] ,resRev[8][2],precis);

  meshC->decrRef();
  meshU->decrRef();
}

void MEDCouplingBasicsTestInterp::testInterpolationCU3D()
{
  MEDCouplingCMesh* meshC = MEDCouplingCMesh::New();
  DataArrayDouble* coords = DataArrayDouble::New();
  double arr[4] = { -1/3., 1/3., 2/3., 4/3. };
  coords->useArray( arr, /*ownership=*/false, CPP_DEALLOC, 4, 1 );
  meshC->setCoords( coords, coords, coords );
  coords->decrRef();

  MEDCouplingUMesh * meshU = buildCU3DMesh_U();

  MEDCouplingNormalizedCartesianMesh<3>      sourceWrapper(meshC);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(meshU);
  INTERP_KERNEL::InterpolationCU myInterpolator;
  std::vector<std::map<int,double> > res;
  const double precis = 1e-13;
  myInterpolator.setPrecision(precis);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");

  double sum = sumAll(res);
  CPPUNIT_ASSERT_EQUAL(8,int( res.size()));
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, sum, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.02700000000000 ,res[0][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00299999999999 ,res[1][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.02999999999999 ,res[1][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03000000000000 ,res[1][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00300000000000 ,res[2][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.02999999999999 ,res[2][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.02999999999999 ,res[2][6],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00033333333333 ,res[3][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[3][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[3][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[3][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[3][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[3][5],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[3][6],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[3][7],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[3][8],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00299999999999 ,res[4][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.02999999999999 ,res[4][9],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03000000000000 ,res[4][18],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00033333333333 ,res[5][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[5][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[5][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[5][9],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[5][10],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[5][11],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[5][18],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[5][19],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[5][20],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00033333333333 ,res[6][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[6][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[6][6],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[6][9],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[6][12],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[6][15],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00333333333333 ,res[6][18],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[6][21],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03333333333333 ,res[6][24],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.7037037037e-05 ,res[7][0],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00037037037037 ,res[7][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00037037037037 ,res[7][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00037037037037 ,res[7][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][5],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00037037037037 ,res[7][6],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][7],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][8],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00037037037037 ,res[7][9],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][10],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][11],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][12],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703703703703 ,res[7][13],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703703703703 ,res[7][14],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][15],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703703703703 ,res[7][16],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703703703703 ,res[7][17],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00037037037037 ,res[7][18],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][19],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][20],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][21],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703703703703 ,res[7][22],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703703703703 ,res[7][23],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00370370370370 ,res[7][24],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703703703703 ,res[7][25],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703703703703 ,res[7][26],precis);


  meshC->decrRef();
  meshU->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP0IntegralUniform()
{
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(5,myInterpolator.toIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.fromIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  res.clear();
  targetMesh->decrRef();
  //
  targetMesh=build2DTargetMeshPerm_1();
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper2(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator2;
  CPPUNIT_ASSERT(myInterpolator2.getMeasureAbsStatus());
  CPPUNIT_ASSERT_EQUAL(5,myInterpolator2.toIntegralUniform(targetWrapper2,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  res.clear();
  myInterpolator2.setMeasureAbsStatus(false);
  CPPUNIT_ASSERT(!myInterpolator2.getMeasureAbsStatus());
  CPPUNIT_ASSERT_EQUAL(5,myInterpolator2.toIntegralUniform(targetWrapper2,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.125,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75,sumAll(res),1e-12);
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP0IntegralUniform()
{
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  std::vector<std::map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(5,myInterpolator.toIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.fromIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP0IntegralUniform()
{
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  INTERP_KERNEL::Interpolation3D myInterpolator;
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  std::vector<std::map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(8,myInterpolator.toIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000.,res[0][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[0][1],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[0][2],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[0][3],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[0][4],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[0][5],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[0][6],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3375000.,res[0][7],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000.,sumAll(res),1e-6);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.fromIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000.,res[0][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[1][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[2][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[3][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[4][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[5][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[6][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3375000.,res[7][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000.,sumAll(res),1e-6);
  res.clear();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP1IntegralUniform()
{
  MEDCouplingUMesh *targetMesh=build2DSourceMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(4,myInterpolator.toIntegralUniform(targetWrapper,res,"P1"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333331,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333331,res[0][3],1e-12);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.fromIntegralUniform(targetWrapper,res,"P1"));
  CPPUNIT_ASSERT_EQUAL(4,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333331,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333331,res[3][0],1e-12);
  res.clear();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DInterpP1IntegralUniform()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(sourceMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(9,myInterpolator.toIntegralUniform(targetWrapper,res,"P1"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][1],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(500000.,res[0][2],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][3],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][4],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(500000.,res[0][5],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][6],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][7],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2000000.,res[0][8],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000.,sumAll(res),1e-6);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.fromIntegralUniform(targetWrapper,res,"P1"));
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[1][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(500000.,res[2][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[3][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[4][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(500000.,res[5][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[6][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[7][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2000000.,res[8][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000.,sumAll(res),1e-6);
  sourceMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DInterpP1P0Bary_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Barycentric,INTERP_KERNEL::BarycentricGeo2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666669,res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[0][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[0][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625,res[1][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[1][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625,res[2][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[2][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625,res[3][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[3][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625,res[3][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[4][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[4][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DSurfInterpP1P0Bary_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  std::vector<std::map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Barycentric,INTERP_KERNEL::BarycentricGeo2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666669*sqrt(2.),res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[0][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[0][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625*sqrt(2.),res[1][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[1][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625*sqrt(2.),res[2][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[2][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625*sqrt(2.),res[3][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[3][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625*sqrt(2.),res[3][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[4][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[4][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

#include <iomanip>
void MEDCouplingBasicsTestInterp::test3DInterpP1P0Bary_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_2();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_2();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::Barycentric);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());

  double res3D[5][28]={{104166.66658918398, 885416.666685817763, 135416.666666666541, 36458.3333333335031, 31249.9999999999018, 145833.333333333256, 41666.6666666667516, 124999.999999999971, 177083.333326388849, 0.0, 31249.9999999999636, 0.0, 41666.666620792399, 159722.22229009436, 0.0, 0.0, 41666.6666631944681, 125000, 43499.2283723790752, 164351.851924000395, 36458.3333372396883, 0.0, 0.0, 125000.000001736029, 34722.2221800900952, 13599.5370788455439, 0.0, 167438.27159690368},
                       {0.0, 41666.6664479170649, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 125000.000161457952, 0.0, 0.0, 0.0, 0.0, 111111.11112005508, 0.0, 0.0, 291666.666656249959, 41666.6666666666933, 6944.4444415638809, 270833.333520485845, 0.0, 0.0, 124999.999989583303, 41666.6665798612958, 20833.3333186342825, 145833.333354303701, 83333.3333263888198, 27777.7777501651799},
                       {0.0, 93750.0000000000728, 125000.000000000058, 0.0, 0.0, 72916.666666666526, 291666.666666666628, 41666.6666666667152, 197916.66666666657, 166666.666666666802, 218750.000000000116, 41666.6666666665697, 0.0, 0.0, 0.0, 0.0, 0.0, 41666.6666666666861, 0.0, 0.0, 0.0, 0.0, 0.0, 41666.6666666666642, 0.0, 0.0, 0.0, 0.0},
                       {72916.6666484848247, 82465.2777799315081, 0.0, 0.0, 217447.916666666686, 197916.666666666802, 0.0, 41666.6666666666715, 0.0, 0.0, 0.0, 0.0, 290364.583310396119, 125000.000018181803, 41666.6666666666351, 166666.666666666599, 0.0, 41666.6666666665551, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 27777.7777734705051, 0.0, 0.0, 27777.7778028684952},
                       {72916.6666461071727, 172309.027782170655, 70312.5000000000437, 253906.250000000029, 0.0, 0.0, 0.0, 41666.666666666657, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 258246.527775988478, 71180.5555571812583, 253906.250006944378, 41666.6666666666861, 0.0, 41666.6666649305407, 20833.3333186342534, 6944.44445267237552, 0.0, 27777.7777953707919}};

  double sum = 0;
  int i=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=res.begin();iter1!=res.end();iter1++,i++)
    {
      for(int j=0;j<28;j++)
        {
          std::map<int,double>::const_iterator iter2=(*iter1).find(j);
          if(iter2!=(*iter1).end())
            {
              sum += iter2->second;
              CPPUNIT_ASSERT_DOUBLES_EQUAL(res3D[i][j],(*iter2).second,1.e-5);
            }
          else
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res3D[i][j],1e-14);
            }
        }
    }
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000,sum,1.e-5);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3DTo1DInterpP0P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DTargetMesh_1();
  MEDCouplingUMesh *targetMesh=build1DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  std::vector<std::map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test1DInterp_1()
{
  //      c1   c0    c2    - pay attention to cell order!
  // S: o---o------o---o
  // T:   o---o------o---o
  //      n0  n1     n2  n3
  //
  // ---+---+------+---+---> X
  //    0.  1.     3.  4.   
  MEDCouplingUMesh *sourceMesh=build1DMesh(0);
  MEDCouplingUMesh *targetMesh=build1DMesh(0.5);
  //
  MEDCouplingNormalizedUnstructuredMesh<1,1> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<1,1> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation1D myInterpolator;
  const double precis = 1e-13;
  myInterpolator.setPrecision(precis);

  // P0P0
  std::vector<std::map<int,double> > res;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL( 3, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5, res[0][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[0][2], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[1][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[1][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[2][2], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.5, sumAll(res), precis);

  // P1P0
  res.clear();
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL( 3, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[0][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5, res[0][2], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, res[1][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[2][3], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.5, sumAll(res), precis);

  // P0P1
  res.clear();
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");

  CPPUNIT_ASSERT_EQUAL( 4, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[0][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.5, res[1][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[2][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, res[2][2], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.5, sumAll(res), precis);

  // P1P1
  res.clear();
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL( 4, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[0][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, res[1][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[1][2], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, res[2][2], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[2][3], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.5, sumAll(res), precis);

  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DCurveInterpP0P0_1()
{
  // coincident meshes
  MEDCouplingUMesh *sourceMesh=build2DCurveMesh(0,0);
  MEDCouplingUMesh *targetMesh=build2DCurveMesh(0,0);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,1> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,1> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2DCurve myInterpolator;
  const double precis = 1e-13;
  myInterpolator.setPrecision(precis);
  std::vector<std::map<int,double> > res;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");

  CPPUNIT_ASSERT_EQUAL( 2, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( sqrt(2.),res[0][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., res[1][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.+sqrt(2.), sumAll(res), precis);

  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DCurveInterpP0P0_2()
{
  // equal meshes shifted one from another along X by 0.5
  MEDCouplingUMesh *sourceMesh=build2DCurveMesh(0.5,0);
  MEDCouplingUMesh *targetMesh=build2DCurveMesh(0,0);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,1> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,1> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2DCurve myInterpolator;
  const double precis = 1e-13;
  myInterpolator.setPrecision(precis);
  myInterpolator.setMedianPlane(1.);// median line on target
  std::vector<std::map<int,double> > res;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");

  double tolInters = myInterpolator.getBoundingBoxAdjustmentAbs() * sqrt(2.);
  CPPUNIT_ASSERT_EQUAL( 2, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0,res[0][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( tolInters,res[0][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, res[1][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5+tolInters, sumAll(res), precis);

  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DCurveInterpP0P1_1()
{
  // coincident meshes
  MEDCouplingUMesh *sourceMesh=build2DCurveMesh(0,0);
  MEDCouplingUMesh *targetMesh=build2DCurveMesh(0,0);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,1> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,1> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2DCurve myInterpolator;
  const double precis = 1e-13;
  myInterpolator.setPrecision(precis);
  std::vector<std::map<int,double> > res;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");

  const double len1 = 1., len0 = sqrt(2.);
  CPPUNIT_ASSERT_EQUAL( 3, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len1, res[0][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len0, res[1][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len1, res[1][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len0, res[2][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( len0+len1, sumAll(res), precis);

  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DCurveInterpP1P0_1()
{
  // coincident meshes
  MEDCouplingUMesh *sourceMesh=build2DCurveMesh(0,0);
  MEDCouplingUMesh *targetMesh=build2DCurveMesh(0,0);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,1> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,1> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2DCurve myInterpolator;
  const double precis = 1e-13;
  myInterpolator.setPrecision(precis);
  std::vector<std::map<int,double> > res;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");

  const double len1 = 1., len0 = sqrt(2.);
  CPPUNIT_ASSERT_EQUAL( 2, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len1, res[1][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len0, res[0][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len1, res[1][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len0, res[0][2], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( len0+len1, sumAll(res), precis);

  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2DCurveInterpP1P1_1()
{
  // coincident meshes
  MEDCouplingUMesh *sourceMesh=build2DCurveMesh(0,0);
  MEDCouplingUMesh *targetMesh=build2DCurveMesh(0,0);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,1> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,1> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2DCurve myInterpolator;
  const double precis = 1e-13;
  myInterpolator.setPrecision(precis);
  std::vector<std::map<int,double> > res;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");

  const double len1 = 1., len0 = sqrt(2.);
  CPPUNIT_ASSERT_EQUAL( 3, int( res.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len1, res[0][0], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*(len0+len1), res[1][1], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*len0, res[2][2], precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( len0+len1, sumAll(res), precis);

  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2D1DBasicInterpP0P0()
{
  MEDCouplingUMesh *sourceMesh=build2D1DSourceMesh();
  MEDCouplingUMesh *targetMesh=build2D1DTargetMesh();

  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D1D myInterpolator;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::Geometric2D);
  std::vector<std::map<int,double> > matrix;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,matrix,"P0P0");

  CPPUNIT_ASSERT_EQUAL(2,(int)matrix.size());

  CPPUNIT_ASSERT_DOUBLES_EQUAL(3., matrix[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., matrix[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., matrix[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8., matrix[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., matrix[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5., matrix[0][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6., matrix[0][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., matrix[0][7],1e-12);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., matrix[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., matrix[1][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., matrix[1][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4., matrix[1][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5., matrix[1][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0., matrix[1][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6., matrix[1][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3., matrix[1][7],1e-12);

  INTERP_KERNEL::Interpolation2D3D::DuplicateFacesType duplicateFaces = myInterpolator.retrieveDuplicateFaces();
  CPPUNIT_ASSERT_EQUAL(1,(int)duplicateFaces.size());

  INTERP_KERNEL::Interpolation2D3D::DuplicateFacesType correctDuplicateFaces;
  std::set<int> face6;
  face6.insert(0);
  face6.insert(1);
  correctDuplicateFaces[6] = face6;

  CPPUNIT_ASSERT(correctDuplicateFaces == duplicateFaces);

  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test2D1DSegQuadInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh();
  MEDCouplingUMesh *targetMesh=build2D1DQuadTargetMesh();
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 16., 0, 4);
}

void MEDCouplingBasicsTestInterp::test2D1DSegQuadInterpP0P0_2()
{
  const double shiftX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build2D1DQuadTargetMesh();
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 2. * 16., 4, 2 * 4);
}

void MEDCouplingBasicsTestInterp::test2D1DSegQuadInterpP0P0_3()
{
  const double shiftX = 1.5;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build2D1DQuadTargetMesh(inclinationX);
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 20., 0, 4);
}

void MEDCouplingBasicsTestInterp::test2D1DSegQuadInterpP0P0_4()
{
  const double shiftX = 3.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build2D1DQuadTargetMesh(inclinationX);
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 2. * 20., 4, 2 * 4);
}

void MEDCouplingBasicsTestInterp::test2D1DSegQuadInterpP0P0_5()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build2D1DQuadTargetMesh(inclinationX);
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 12., 0, 3);
}

void MEDCouplingBasicsTestInterp::test2D1DSegQuadInterpP0P0_6()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build2D1DQuadTargetMesh();
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 10., 0, 2);
}

void MEDCouplingBasicsTestInterp::test2D1DSegTriInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh();
  MEDCouplingUMesh *targetMesh=build2D1DTriTargetMesh();
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 16., 0, 4);
}

void MEDCouplingBasicsTestInterp::test2D1DSegTriInterpP0P0_2()
{
  const double shiftX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build2D1DTriTargetMesh();
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 2. * 16., 4, 2 * 4);
}

void MEDCouplingBasicsTestInterp::test2D1DSegTriInterpP0P0_3()
{
  const double shiftX = 1.5;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build2D1DTriTargetMesh(inclinationX);
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 20., 0, 8);
}

void MEDCouplingBasicsTestInterp::test2D1DSegTriInterpP0P0_4()
{
  const double shiftX = 3.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build2D1DTriTargetMesh(inclinationX);
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 2. * 20., 4, 8);
}

void MEDCouplingBasicsTestInterp::test2D1DSegTriInterpP0P0_5()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build2D1DTriTargetMesh(inclinationX);
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 12., 0, 6);
}

void MEDCouplingBasicsTestInterp::test2D1DSegTriInterpP0P0_6()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build2D1DSegSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build2D1DTriTargetMesh();
  test2D1DMeshesIntersection(sourceMesh, targetMesh, 20., 2, 4);
}

void MEDCouplingBasicsTestInterp::test3D2DBasicInterpP0P0()
{
  MEDCouplingUMesh *sourceMesh=build3D2DSourceMesh();
  MEDCouplingUMesh *targetMesh=build3D2DTargetMesh();

  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D3D myInterpolator;
  myInterpolator.setPrecision(1e-12);
  std::vector<std::map<int,double> > matrix;
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( size_t i = 0; i < sizeof(sp)/sizeof(sp[0]); ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    matrix.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,matrix,"P0P0");

    CPPUNIT_ASSERT_EQUAL(3,(int)matrix.size());

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[0][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(40.        ,matrix[0][2],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.         ,matrix[0][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5        ,matrix[0][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[0][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(32.        ,matrix[0][6],1e-12);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.*sqrt(3.),matrix[1][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[1][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(40.        ,matrix[1][2],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(80.        ,matrix[1][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[1][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(80.        ,matrix[1][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(80.        ,matrix[1][6],1e-12);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(32.        ,matrix[2][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[2][2],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[2][3],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.         ,matrix[2][4],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(80.        ,matrix[2][5],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(112.       ,matrix[2][6],1e-12);

    INTERP_KERNEL::Interpolation2D3D::DuplicateFacesType duplicateFaces = myInterpolator.retrieveDuplicateFaces();
    CPPUNIT_ASSERT_EQUAL(3,(int)duplicateFaces.size());

    INTERP_KERNEL::Interpolation2D3D::DuplicateFacesType correctDuplicateFaces;
    std::set<int> face2;
    face2.insert(0);
    face2.insert(1);
    correctDuplicateFaces[2] = face2;
    std::set<int> face5;
    face5.insert(1);
    face5.insert(2);
    correctDuplicateFaces[5] = face5;
    std::set<int> face6;
    face6.insert(0);
    face6.insert(1);
    face6.insert(2);
    correctDuplicateFaces[6] = face6;

    CPPUNIT_ASSERT(correctDuplicateFaces == duplicateFaces);
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTestInterp::test3D2DQuadHexaInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh();
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 240., 0, 20);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadHexaInterpP0P0_2()
{
  const double shiftX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 2. * 240., 20, 2 * 20);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadHexaInterpP0P0_3()
{
  const double shiftX = 1.5;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 300., 0, 20);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadHexaInterpP0P0_4()
{
  const double shiftX = 3.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 2. * 300., 20, 2 * 20);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadHexaInterpP0P0_5()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 180., 0, 15);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadHexaInterpP0P0_6()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 150., 0, 10);
}

void MEDCouplingBasicsTestInterp::test3D2DTriHexaInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh();
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 240., 0, 40);
}

void MEDCouplingBasicsTestInterp::test3D2DTriHexaInterpP0P0_2()
{
  const double shiftX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 2. * 240., 40, 2 * 40);
}

void MEDCouplingBasicsTestInterp::test3D2DTriHexaInterpP0P0_3()
{
  const double shiftX = 1.5;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 300., 0, 40);
}

void MEDCouplingBasicsTestInterp::test3D2DTriHexaInterpP0P0_4()
{
  const double shiftX = 3.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 2. * 300., 40, 2 * 40);
}

void MEDCouplingBasicsTestInterp::test3D2DTriHexaInterpP0P0_5()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 180., 0, 30);
}

void MEDCouplingBasicsTestInterp::test3D2DTriHexaInterpP0P0_6()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DHexaTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 150., 0, 20);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadTetraInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh();
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 240., 20, 40);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadTetraInterpP0P0_2()
{
  const double shiftX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 2. * 240., 20, 2 * 40);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadTetraInterpP0P0_3()
{
  const double shiftX = 1.5;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 300., 0, 100);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadTetraInterpP0P0_4()
{
  const double shiftX = 3.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 2. * 300., 20, 2 * 40);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadTetraInterpP0P0_5()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 180., 0, 45);
}

void MEDCouplingBasicsTestInterp::test3D2DQuadTetraInterpP0P0_6()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DQuadSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 150., 0, 30);
}

void MEDCouplingBasicsTestInterp::test3D2DTriTetraInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh();
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 240., 0, 40);
}

void MEDCouplingBasicsTestInterp::test3D2DTriTetraInterpP0P0_2()
{
  const double shiftX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 2. * 240., 40, 40 + 80);
}

void MEDCouplingBasicsTestInterp::test3D2DTriTetraInterpP0P0_3()
{
  const double shiftX = 1.5;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 300., 0);
}

void MEDCouplingBasicsTestInterp::test3D2DTriTetraInterpP0P0_4()
{
  const double shiftX = 3.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 2. * 300., 40, 40 + 80);
}

void MEDCouplingBasicsTestInterp::test3D2DTriTetraInterpP0P0_5()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh(inclinationX);
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 180., 0);
}

void MEDCouplingBasicsTestInterp::test3D2DTriTetraInterpP0P0_6()
{
  const double shiftX = 9.;
  const double inclinationX = 3.;
  MEDCouplingUMesh *sourceMesh=build3D2DTriSourceMesh(shiftX, inclinationX);
  MEDCouplingUMesh *targetMesh=build3D2DTetraTargetMesh();
  test3D2DMeshesIntersection(sourceMesh, targetMesh, 150., 0);
}

