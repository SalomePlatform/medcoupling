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
#include "TestInterpKernelUtils.hxx"

#include "MeshTestToolkit.hxx"

#include "MEDFileMesh.hxx"

#include "MEDCouplingNormalizedUnstructuredMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "MEDCouplingFieldDouble.hxx"

#include "Interpolation3DSurf.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3D.txx"

#include <map>
#include <cmath>
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm>


// levels :
// 1 - titles and volume results
// 2 - symmetry / diagonal results and intersection matrix output
// 3 - empty
// 4 - empty
// 5 - misc
#include "Log.hxx"

#include <cppunit/extensions/HelperMacros.h>

//#define VOL_PREC 1.0e-6
using namespace MEDCoupling;
using namespace INTERP_KERNEL;

namespace INTERP_TEST
{
  /**
   * Calculates the sum of a row of an intersection matrix
   *
   * @param m  an intersection matrix
   * @param i  the index of the row (1 <= i <= no. rows)
   * @return the sum of the values of row i
   *
   */
  template <int SPACEDIM, int MESHDIM>
  double MeshTestToolkit<SPACEDIM,MESHDIM>::sumRow(const IntersectionMatrix& m, int i) const
  {
    double vol = 0.0;
    for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
      {
        if(iter->count(i) != 0.0)
          {
            std::map<int, double>::const_iterator iter2 = iter->find(i);
            vol += fabs(iter2->second);
          }
      }
    return vol;
  }

  /**
   * Calculates the sum of a column of an intersection matrix
   *
   * @param m  an intersection matrix
   * @param i  the index of the column (0 <= i <= no. rows - 1)
   * @return the sum of the values of column i
   *
   */
  template <int SPACEDIM, int MESHDIM>
  double MeshTestToolkit<SPACEDIM,MESHDIM>::sumCol(const IntersectionMatrix& m, int i) const
  {
    double vol = 0.0;
    const std::map<int, double>& col = m[i];
    for(std::map<int, double>::const_iterator iter = col.begin() ; iter != col.end() ; ++iter)
      {
        vol += fabs(iter->second);
      }
    return vol;
  }

  /**
   * Gets the volumes of the elements in a mesh.
   *
   * @param mesh   the mesh
   * @param tab    pointer to double[no. elements of mesh] array in which to store the volumes
   */
  template <int SPACEDIM, int MESHDIM>
  void MeshTestToolkit<SPACEDIM,MESHDIM>::getVolumes(MEDCoupling::MEDCouplingUMesh& mesh, double *tab) const
  {
    MCAuto<MEDCouplingFieldDouble> vol=mesh.getMeasureField(true);
    std::copy(vol->getArray()->begin(),vol->getArray()->end(),tab);
  }

  /**
   * Sums all the elements (volumes) of an intersection matrix
   *
   * @param m  the intersection matrix
   * @return   the sum of the elements of m
   */

  template <int SPACEDIM, int MESHDIM>
  double MeshTestToolkit<SPACEDIM,MESHDIM>::sumVolume(const IntersectionMatrix& m) const
  {
    std::vector<double> volumes;
    for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
      {
        for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
          {
            volumes.push_back(fabs(iter2->second));
          }
      }

    // sum in ascending order to avoid rounding errors

    sort(volumes.begin(), volumes.end());
    const double vol = accumulate(volumes.begin(), volumes.end(), 0.0);

    return vol;
  }

  /**
   * Verifies if for a given intersection matrix the sum of each row is equal to the volumes
   * of the corresponding source elements and the sum of each column is equal to the volumes
   * of the corresponding target elements. This will be true as long as the meshes correspond
   * to the same geometry. The equalities are in the "epsilon-sense", making sure the relative
   * error is small enough.
   *
   * @param  m       the intersection matrix
   * @param  sMesh   the source mesh
   * @param  tMesh   the target mesh
   * @return true if the condition is verified, false if not.
   */
  template <int SPACEDIM, int MESHDIM>
  bool MeshTestToolkit<SPACEDIM,MESHDIM>::testVolumes(const IntersectionMatrix& m,  MEDCoupling::MEDCouplingUMesh& sMesh,  MEDCoupling::MEDCouplingUMesh& tMesh) const
  {
    bool ok = true;

    // source elements
    double* sVol = new double[sMesh.getNumberOfCells()];
    getVolumes(sMesh, sVol);

    for(int i = 0; i < sMesh.getNumberOfCells(); ++i)
      {
        const double sum_row = sumRow(m, i);
        if(!epsilonEqualRelative(sum_row, fabs(sVol[i]), _precision))
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
        if(!epsilonEqualRelative(sum_col,fabs(tVol[i]), _precision))
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

  /**
   * Verifies that two intersection matrices have the necessary elements to be able to be each others' transposes.
   *
   * @param m1  the first intersection matrix
   * @param m2  the second intersection matrix
   *
   * @return true if for each element (i,j) of m1, the element (j,i) exists in m2, false if not.
   */
  template <int SPACEDIM, int MESHDIM>
  bool MeshTestToolkit<SPACEDIM,MESHDIM>::areCompatitable(const IntersectionMatrix& m1, const IntersectionMatrix& m2) const
  {
    bool compatitable = true;
    int i = 0;
    for(IntersectionMatrix::const_iterator iter = m1.begin() ; iter != m1.end() ; ++iter)
      {
        for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
          {
            int j = iter2->first;
            if(m2.at(j).count(i) == 0)
              {
                if(!epsilonEqual(iter2->second, 0.0, _precision))
                  {
                    LOG(2, "V1( " << i << ", " << j << ") exists, but V2( " << j << ", " << i << ") " << " does not " );
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

  /**
   * Tests if two intersection matrices are each others' transposes.
   *
   * @param m1  the first intersection matrix
   * @param m2  the second intersection matrix
   * @return true if m1 = m2^T, false if not.
   */
  template <int SPACEDIM, int MESHDIM>
  bool MeshTestToolkit<SPACEDIM,MESHDIM>::testTranspose(const IntersectionMatrix& m1, const IntersectionMatrix& m2) const
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
            const double v1 = fabs(iter2->second);
            //if(m2[j - 1].count(i+1) > 0)
            //  {
            std::map<int, double> theMap =  m2.at(j);
            const double v2 = fabs(theMap[i]);
            if(v1 != v2)
              {
                LOG(2, "V1( " << i << ", " << j << ") = " << v1 << " which is different from V2( " << j << ", " << i << ") = " << v2 << " | diff = " << v1 - v2 );
                if(!epsilonEqualRelative(v1, v2, _precision))
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

  /**
   * Tests if an intersection matrix is diagonal.
   *
   * @param m   the intersection matrix
   * @return true if m is diagonal; false if not
   *
   */
  template <int SPACEDIM, int MESHDIM>
  bool MeshTestToolkit<SPACEDIM,MESHDIM>::testDiagonal(const IntersectionMatrix& m) const
  {
    LOG(1, "Checking if matrix is diagonal" );
    int i = 0;
    bool isDiagonal = true;
    for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
      {
        for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
          {
            int j = iter2->first;
            const double vol = iter2->second;
            if(vol != 0.0 && (i != j))
              {
                LOG(2, "V( " << i << ", " << j << ") = " << vol << " which is not zero" );
                if(!epsilonEqual(vol, 0.0, _precision))
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

  /**
   * Outputs the intersection matrix as a list of all its elements to std::cout.
   *
   * @param m the intersection matrix to output
   */
  template <int SPACEDIM, int MESHDIM>
  void MeshTestToolkit<SPACEDIM,MESHDIM>::dumpIntersectionMatrix(const IntersectionMatrix& m) const
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

  /**
   * Calculates the intersection matrix for two meshes.
   * If the source and target meshes are the same, a CppUnit assertion raised if testVolumes() returns false.
   *
   * @param  mesh1path   the path to the file containing the source mesh, relative to {$MEDCOUPLING_ROOT_DIR}/share/resources/med/
   * @param  mesh1       the name of the source mesh
   * @param  mesh2path   the path to the file containing the target mesh, relative to {$MEDCOUPLING_ROOT_DIR}/share/resources/med/
   * @param  mesh2       the name of the target mesh
   * @param  m           intersection matrix in which to store the result of the intersection
   */
  template <int SPACEDIM, int MESHDIM>
  void MeshTestToolkit<SPACEDIM,MESHDIM>::calcIntersectionMatrix(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, IntersectionMatrix& m) const
  {
    LOG(1, std::endl << "=== -> intersecting src = " << mesh1path << ", target = " << mesh2path );

    LOG(5, "Loading " << mesh1 << " from " << mesh1path);
    MCAuto<MEDFileUMesh> sMeshML=MEDFileUMesh::New(INTERP_TEST::getResourceFile(mesh1path).c_str(),mesh1);
    MCAuto<MEDCouplingUMesh> sMesh=sMeshML->getMeshAtLevel(0);

    LOG(5, "Loading " << mesh2 << " from " << mesh2path);
    MCAuto<MEDFileUMesh> tMeshML=MEDFileUMesh::New(INTERP_TEST::getResourceFile(mesh2path).c_str(),mesh2);
    MCAuto<MEDCouplingUMesh> tMesh=tMeshML->getMeshAtLevel(0);

    MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM> sMesh_wrapper(sMesh);
    MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM> tMesh_wrapper(tMesh);

    if (SPACEDIM==2 && MESHDIM==2)
      {
        Interpolation2D interpolator;
        interpolator.setOptions(_precision, LOG_LEVEL, _intersectionType,1);
        interpolator.interpolateMeshes(sMesh_wrapper, tMesh_wrapper,m,"P0P0");
      }
    else if (SPACEDIM==3 && MESHDIM==2)
      {
        Interpolation3DSurf interpolator;
        interpolator.setOptions(_precision,LOG_LEVEL, 0.5,_intersectionType,false,1);
        interpolator.interpolateMeshes(sMesh_wrapper, tMesh_wrapper,m,"P0P0");
      }
    else if (SPACEDIM==3 && MESHDIM==3)
      {
        Interpolation3D interpolator;
        interpolator.interpolateMeshes(sMesh_wrapper, tMesh_wrapper,m,"P0P0");
      }
    else
      {
        throw INTERP_KERNEL::Exception("Wrong dimensions");
      }
    // if reflexive, check volumes
    if(strcmp(mesh1path,mesh2path) == 0)
      {
        const bool row_and_col_sums_ok = testVolumes(m, *sMesh, *tMesh);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Row or column sums incorrect", true, row_and_col_sums_ok);
        const bool is_diagonal =testDiagonal(m);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Self intersection matrix is not diagonal", true, is_diagonal);
      }

    LOG(1, "Intersection calculation done. " << std::endl );
  }

  /**
   * Tests the intersection algorithm for two meshes.
   * Depending on the nature of the meshes, different tests will be performed. The sum of the elements will
   * be compared to the given total volume of the intersection in all cases. If the two meshes are the same, then
   * it will be confirmed that the intersection matrix is diagonal, otherwise the intersection matrices will be
   * calculated once which each mesh as source mesh, and it will be verified that the they are each others' transpose.
   *
   * @param  mesh1path   the path to the file containing the source mesh, relative to {$MEDCOUPLING_ROOT_DIR}/share/resources/med/
   * @param  mesh1       the name of the source mesh
   * @param  mesh2path   the path to the file containing the target mesh, relative to {$MEDCOUPLING_ROOT_DIR}/share/resources/med/
   * @param  mesh2       the name of the target mesh
   * @param  correctVol  the total volume of the intersection of the two meshes
   * @param  prec        maximum relative error to be tolerated in volume comparisions
   * @param  doubleTest  if false, only the test with mesh 1 as the source mesh and mesh 2 as the target mesh will be performed
   *
   */
  template <int SPACEDIM, int MESHDIM>
  void MeshTestToolkit<SPACEDIM,MESHDIM>::intersectMeshes(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, const double correctVol, const double prec, bool doubleTest) const
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

        CPPUNIT_ASSERT_EQUAL_MESSAGE("Symmetry test failed", true, testTranspose(matrix1, matrix2));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(correctVol, vol1, prec * std::max(vol1, correctVol));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(correctVol, vol2, prec * std::max(vol2, correctVol));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(vol1, vol2, prec * std::max(vol1, vol2));
      }
  }

  /**
   * Utility method used to facilitate the call to intersect meshes.
   * It calls intersectMeshes, using "mesh1.med" as file name for the mesh with name "mesh1" and
   * "mesh2.med" as file name for the mesh with name "mesh2". The rest of the arguments are passed
   * along as they are.
   *
   * @param  mesh1       the name of the source mesh
   * @param  mesh2       the name of the target mesh
   * @param  correctVol  the total volume of the intersection of the two meshes
   * @param  prec        maximum relative error to be tolerated in volume comparisions
   * @param  doubleTest  if false, only the test with mesh 1 as the source mesh and mesh 2 as the target mesh will be performed
   *
   */
  template <int SPACEDIM, int MESHDIM>
  void MeshTestToolkit<SPACEDIM,MESHDIM>::intersectMeshes(const char* mesh1, const char* mesh2, const double correctVol, const double prec, bool doubleTest) const
  {
    const std::string path1 = std::string(mesh1) + std::string(".med");
    std::cout << "here :" << path1 << std::endl;
    const std::string path2 = std::string(mesh2) + std::string(".med");

    intersectMeshes(path1.c_str(), mesh1, path2.c_str(), mesh2, correctVol, prec, doubleTest);
  }
}
