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

#ifndef _TESTING_UTILS_HXX_
#define _TESTING_UTILS_HXX_

#include "Interpolation3D.hxx"
#include "MEDMEM_Mesh.hxx"

#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>

#include "VectorUtils.hxx"

// levels : 
// 1 - titles and volume results
// 2 - symmetry / diagonal results and intersection matrix output
// 3 - empty
// 4 - empty
// 5 - misc
#include "Log.hxx"

using namespace MEDMEM;
using namespace INTERP_KERNEL;
using namespace MED_EN;


double sumVolume(const IntersectionMatrix& m) 
{
  
  vector<double> volumes;
  for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
    {
      for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
        {
          volumes.push_back(iter2->second);
          //    vol += std::fabs(iter2->second);
        }
    }
  
  // sum in ascending order to avoid rounding errors

  sort(volumes.begin(), volumes.end());
  const double vol = accumulate(volumes.begin(), volumes.end(), 0.0);

  return vol;
}

#if 0

bool areCompatitable(const IntersectionMatrix& m1, const IntersectionMatrix& m2)
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
      
bool testSymmetric(const IntersectionMatrix& m1, const IntersectionMatrix& m2)
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

bool testDiagonal(const IntersectionMatrix& m)
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

#endif

void dumpIntersectionMatrix(const IntersectionMatrix& m) 
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

std::pair<int,int> countNumberOfMatrixEntries(const IntersectionMatrix& m)
{
  
  int numElems = 0;
  int numNonZero = 0;
  for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
    {
      numElems += iter->size();
      for(map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
        {
          if(!epsilonEqual(iter2->second, 0.0, VOL_PREC))
            {
              ++numNonZero;
            }
        }
    }
  return std::make_pair(numElems, numNonZero);
}


void calcIntersectionMatrix(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, IntersectionMatrix& m) 
{
  const std::string dataBaseDir = getenv("MED_ROOT_DIR");
  const std::string dataDir = dataBaseDir + "/share/salome/resources/med/";

  LOG(1, std::endl << "=== -> intersecting src = " << mesh1 << ", target = " << mesh2 );

  LOG(5, "Loading " << mesh1 << " from " << mesh1path);
  const MESH sMesh(MED_DRIVER, dataDir+mesh1path, mesh1);
  const int numSrcElems = sMesh.getNumberOfElements(MED_CELL, MED_ALL_ELEMENTS);
  LOG(1, "Source mesh has " << numSrcElems << " elements");


  LOG(5, "Loading " << mesh2 << " from " << mesh2path);
  const MESH tMesh(MED_DRIVER, dataDir+mesh2path, mesh2);
  const int numTargetElems = tMesh.getNumberOfElements(MED_CELL, MED_ALL_ELEMENTS);

  LOG(1, "Target mesh has " << numTargetElems << " elements");

  Interpolation3D* interpolator = new Interpolation3D();

  m = interpolator->interpolateMeshes(sMesh, tMesh);

  std::pair<int, int> eff = countNumberOfMatrixEntries(m);
  LOG(1, eff.first << " of " << numTargetElems * numSrcElems << " intersections calculated : ratio = " 
      << double(eff.first) / double(numTargetElems * numSrcElems));
  LOG(1, eff.second << " non-zero elements of " << eff.first << " total : filter efficiency = " 
      << double(eff.second) / double(eff.first));

  delete interpolator;

  LOG(1, "Intersection calculation done. " << std::endl );
  
}








#if 0
void intersectMeshes(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, const double correctVol, const double prec, bool doubleTest) 
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
      // CPPUNIT_ASSERT_DOUBLES_EQUAL(correctVol, vol1, prec * std::max(correctVol, vol1));
     
      if(isTestReflexive)
        {
          // CPPUNIT_ASSERT_EQUAL_MESSAGE("Reflexive test failed", true, testDiagonal(matrix1));
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

      // CPPUNIT_ASSERT_DOUBLES_EQUAL(correctVol, vol1, prec * std::max(vol1, correctVol));
      // CPPUNIT_ASSERT_DOUBLES_EQUAL(correctVol, vol2, prec * std::max(vol2, correctVol));
      // CPPUNIT_ASSERT_DOUBLES_EQUAL(vol1, vol2, prec * std::max(vol1, vol2));
      // CPPUNIT_ASSERT_EQUAL_MESSAGE("Symmetry test failed", true, testSymmetric(matrix1, matrix2));
    }

}



#endif


#endif
