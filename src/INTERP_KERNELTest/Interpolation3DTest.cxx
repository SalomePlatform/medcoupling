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

#include "Interpolation3DTest.hxx"

#include "MEDFileMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "VectorUtils.hxx"

// levels :
// 1 - titles and volume results
// 2 - symmetry / diagonal results and intersection matrix output
// 3 - empty
// 4 - empty
// 5 - misc
#include "Log.hxx"


//#define VOL_PREC 1.0e-6

using namespace MEDCoupling;
using namespace INTERP_KERNEL;

double Interpolation3DTest::sumRow(const IntersectionMatrix& m, int i) const
{
  double vol = 0.0;
  for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
    {
      if(iter->count(i) != 0.0)
        {
          std::map<int, double>::const_iterator iter2 = iter->find(i);
          vol += iter2->second;
        }
    }
  return vol;
}

double Interpolation3DTest::sumCol(const IntersectionMatrix& m, int i) const
{
  double vol = 0.0;
  const std::map<int, double>& col = m[i];
  for(std::map<int, double>::const_iterator iter = col.begin() ; iter != col.end() ; ++iter)
    {
      vol += std::fabs(iter->second);
    }
  return vol;
}


void Interpolation3DTest::getVolumes(MEDCoupling::MEDCouplingUMesh& mesh, double *tab) const
{
  MCAuto<MEDCouplingFieldDouble> vol=mesh->getMeasureField(true);
  std::copy(vol->getArray()->begin(),vol->getArray()->end(),tab);
}

double Interpolation3DTest::sumVolume(const IntersectionMatrix& m) const
{

  std::vector<double> volumes;
  for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
    {
      for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
        {
          volumes.push_back(iter2->second);
          //    vol += std::abs(iter2->second);
        }
    }

  // sum in ascending order to avoid rounding errors

  sort(volumes.begin(), volumes.end());
  const double vol = accumulate(volumes.begin(), volumes.end(), 0.0);

  return vol;
}

bool Interpolation3DTest::testVolumes(const IntersectionMatrix& m,  MEDCouplingUMesh& sMesh,  MEDCouplingUMesh& tMesh) const
{
  bool ok = true;

  // source elements
  double* sVol = new double[sMesh.getNumberOfCells()];
  getVolumes(sMesh, sVol);

  for(int i = 0; i < sMesh.getNumberOfCells(); ++i)
    {
      const double sum_row = sumRow(m, i+1);
      if(!epsilonEqualRelative(sum_row, sVol[i], VOL_PREC))
        {
          LOG(1, "Source volume inconsistent : vol of cell " << i << " = " << sVol[i] << " but the row sum is " << sum_row );
          ok = false;
        }
      LOG(1, "diff = " <<sum_row - sVol[i] );
    }

  // target elements
  double* tVol = new double[tMesh.getNumberOfCells()];
  getVolumes(tMesh, tVol);
  for(int i = 0; i < tMesh.getNumberOfCells(); ++i)
    {
      const double sum_col = sumCol(m, i);
      if(!epsilonEqualRelative(sum_col, tVol[i], VOL_PREC))
        {
          LOG(1, "Target volume inconsistent : vol of cell " << i << " = " << tVol[i] << " but the col sum is " << sum_col);
          ok = false;
        }
      LOG(1, "diff = " <<sum_col - tVol[i] );
    }
  delete[] sVol;
  delete[] tVol;

  return ok;
}

bool Interpolation3DTest::areCompatitable(const IntersectionMatrix& m1, const IntersectionMatrix& m2) const
{
  bool compatitable = true;
  int i = 0;
  for(IntersectionMatrix::const_iterator iter = m1.begin() ; iter != m1.end() ; ++iter)
    {
      for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
        {
          int j = iter2->first;
          if(m2.at(j-1).count(i+1) == 0)
            {
              if(!epsilonEqual(iter2->second, 0.0, VOL_PREC))
                {
                  LOG(2, "V1( " << i << ", " << j << ") exists, but V2( " << j - 1 << ", " << i + 1 << ") " << " does not " );
                  LOG(2, "(" << i << ", " << j << ") fails");
                  compatitable = false;
                }
            }
        }
      ++i;
    }
  if(!compatitable)
    {
      LOG(1, "*** matrices are not compatitable");
    }
  return compatitable;
}

bool Interpolation3DTest::testSymmetric(const IntersectionMatrix& m1, const IntersectionMatrix& m2) const
{

  int i = 0;
  bool isSymmetric = true;

  LOG(1, "Checking symmetry src - target" );
  isSymmetric = isSymmetric & areCompatitable(m1, m2) ;
  LOG(1, "Checking symmetry target - src" );
  isSymmetric = isSymmetric & areCompatitable(m2, m1);

  for(IntersectionMatrix::const_iterator iter = m1.begin() ; iter != m1.end() ; ++iter)
    {
      for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
        {
          int j = iter2->first;
          const double v1 = iter2->second;
          //if(m2[j - 1].count(i+1) > 0)
          //  {
          std::map<int, double> theMap =  m2.at(j-1);
          const double v2 = theMap[i + 1];
          if(v1 != v2)
            {
              LOG(2, "V1( " << i << ", " << j << ") = " << v1 << " which is different from V2( " << j - 1 << ", " << i + 1 << ") = " << v2 << " | diff = " << v1 - v2 );
              if(!epsilonEqualRelative(v1, v2, VOL_PREC))
                {
                  LOG(2, "(" << i << ", " << j << ") fails");
                  isSymmetric = false;
                }
            }
        }
      ++i;
    }
  if(!isSymmetric)
    {
      LOG(1, "*** matrices are not symmetric");
    }
  return isSymmetric;
}

bool Interpolation3DTest::testDiagonal(const IntersectionMatrix& m) const
{
  LOG(1, "Checking if matrix is diagonal" );
  int i = 1;
  bool isDiagonal = true;
  for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
    {
      for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
        {
          int j = iter2->first;
          const double vol = iter2->second;
          if(vol != 0.0 && (i != j))
            {
              LOG(2, "V( " << i - 1 << ", " << j << ") = " << vol << " which is not zero" );
              if(!epsilonEqual(vol, 0.0, VOL_PREC))
                {
                  LOG(2, "(" << i << ", " << j << ") fails");
                  isDiagonal = false;
                }
            }
        }
      ++i;
    }
  if(!isDiagonal)
    {
      LOG(1, "*** matrix is not diagonal");
    }
  return isDiagonal;
}

void Interpolation3DTest::dumpIntersectionMatrix(const IntersectionMatrix& m) const
{
  int i = 0;
  std::cout << "Intersection matrix is " << std::endl;
  for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
    {
      for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
        {

          std::cout << "V(" << i << ", " << iter2->first << ") = " << iter2->second << std::endl;

        }
      ++i;
    }
  std::cout << "Sum of volumes = " << sumVolume(m) << std::endl;
}

void Interpolation3DTest::setUp()
{
  interpolator = new Interpolation3D();
}

void Interpolation3DTest::tearDown()
{
  delete interpolator;
}

void Interpolation3DTest::calcIntersectionMatrix(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, IntersectionMatrix& m) const
{
  string dataDir = "";
  if ( getenv("MEDCOUPLING_ROOT_DIR") ) {
    dataDir = getenv("MEDCOUPLING_ROOT_DIR");
    dataDir += "/share/resources/med/";
  }
  else {
    dataDir = get_current_dir_name();
    dataDir += "/../../resources/";
  }

  LOG(1, std::endl << "=== -> intersecting src = " << mesh1 << ", target = " << mesh2 );

  LOG(5, "Loading " << mesh1 << " from " << mesh1path);
  MESH sMesh(MED_DRIVER, dataDir+mesh1path, mesh1);

  LOG(5, "Loading " << mesh2 << " from " << mesh2path);
  MESH tMesh(MED_DRIVER, dataDir+mesh2path, mesh2);

  m = interpolator->interpolateMeshes(sMesh, tMesh);

  // if reflexive, check volumes
  if(strcmp(mesh1path,mesh2path) == 0)
    {
      const bool row_and_col_sums_ok = testVolumes(m, sMesh, tMesh);
      CPPUNIT_ASSERT_EQUAL_MESSAGE("Row or column sums incorrect", true, row_and_col_sums_ok);
    }

  LOG(1, "Intersection calculation done. " << std::endl );

}

void Interpolation3DTest::intersectMeshes(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, const double correctVol, const double prec, bool doubleTest) const
{
  LOG(1, std::endl << std::endl << "=============================" );

  using std::string;
  const string path1 = string(mesh1path) + string(mesh1);
  const string path2 = string(mesh2path) + string(mesh2);

  const bool isTestReflexive = (path1.compare(path2) == 0);

  IntersectionMatrix matrix1;
  calcIntersectionMatrix(mesh1path, mesh1, mesh2path, mesh2, matrix1);

#if LOG_LEVEL >= 2
  dumpIntersectionMatrix(matrix1);
#endif

  std::cout.precision(16);

  const double vol1 = sumVolume(matrix1);

  if(!doubleTest)
    {
      LOG(1, "vol =  " << vol1 <<"  correctVol = " << correctVol );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(correctVol, vol1, prec * std::max(correctVol, vol1));

      if(isTestReflexive)
        {
          CPPUNIT_ASSERT_EQUAL_MESSAGE("Reflexive test failed", true, testDiagonal(matrix1));
        }
    }
  else
    {

      IntersectionMatrix matrix2;
      calcIntersectionMatrix(mesh2path, mesh2, mesh1path, mesh1, matrix2);

#if LOG_LEVEL >= 2
      dumpIntersectionMatrix(matrix2);
#endif

      const double vol2 = sumVolume(matrix2);

      LOG(1, "vol1 =  " << vol1 << ", vol2 = " << vol2 << ", correctVol = " << correctVol );

      CPPUNIT_ASSERT_DOUBLES_EQUAL(correctVol, vol1, prec * std::max(vol1, correctVol));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(correctVol, vol2, prec * std::max(vol2, correctVol));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(vol1, vol2, prec * std::max(vol1, vol2));
      CPPUNIT_ASSERT_EQUAL_MESSAGE("Symmetry test failed", true, testSymmetric(matrix1, matrix2));
    }

}



