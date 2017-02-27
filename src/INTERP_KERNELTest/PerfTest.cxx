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

#include "Interpolation3D.hxx"
#include "Interpolation3D.txx"
#include "MeshTestToolkit.txx"
#include "Log.hxx"
#include "VectorUtils.hxx"
#include "TestInterpKernelUtils.hxx"

#include "MEDCouplingNormalizedUnstructuredMesh.hxx"

#include <cassert>
#include <string>

/**
 * \file PerfTest.cxx
 * Test program which takes two meshes and calculates their intersection matrix.
 *
 * USAGE : PerfTest mesh1 mesh2
 *         where mesh1 and mesh2 are the names of two meshes located in
 *         the files mesh1.med, mesh2.med in {$MEDCOUPLING_ROOT_DIR}/share/resources/med/
 *
 */

namespace INTERP_TEST
{
  /**
   * \brief Specialization of MeshTestToolkit for the purposes of performance testing.
   *
   */
  class PerfTestToolkit : public MeshTestToolkit<3,3>
  {

  public:

    /**
     * Calculates the intersection matrix for two meshes.
     * Outputs the names of the meshes intersected, the number of elements in each mesh,
     * the number of matrix elements and the number of non-zero matrix elements, etc.
     * These values help to determine how well the filtering algorithm is working.
     *
     * @param  mesh1path   the path to the file containing the source mesh, relative to {$MEDCOUPLING_ROOT_DIR}/share/resources/med/
     * @param  mesh1       the name of the source mesh
     * @param  mesh2path   the path to the file containing the target mesh, relative to {$MEDCOUPLING_ROOT_DIR}/share/resources/med/
     * @param  mesh2       the name of the target mesh
     * @param  m           intersection matrix in which to store the result of the intersection
     */
    void calcIntersectionMatrix(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, IntersectionMatrix& m)
    {
      LOG(1, std::endl << "=== -> intersecting src = " << mesh1 << ", target = " << mesh2 );

      LOG(5, "Loading " << mesh1 << " from " << mesh1path);
      MCAuto<MEDFileUMesh> sMeshML=MEDFileUMesh::New(INTERP_TEST::getResourceFile(mesh1path).c_str(),mesh1);
      MCAuto<MEDCouplingUMesh> sMesh=sMeshML->getMeshAtLevel(0);


      LOG(5, "Loading " << mesh2 << " from " << mesh2path);
      MCAuto<MEDFileUMesh> tMeshML=MEDFileUMesh::New(INTERP_TEST::getResourceFile(mesh2path).c_str(),mesh2);
      MCAuto<MEDCouplingUMesh> tMesh=tMeshML->getMeshAtLevel(0);

      MEDCouplingNormalizedUnstructuredMesh<3,3> sMesh_wrapper(sMesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> tMesh_wrapper(tMesh);

      Interpolation3D interpolator;
      interpolator.interpolateMeshes(sMesh_wrapper, tMesh_wrapper,m,"P0P0");

      std::pair<int, int> eff = countNumberOfMatrixEntries(m);
//      LOG(1, eff.first << " of " << numTargetElems * numSrcElems << " intersections calculated : ratio = "
//          << double(eff.first) / double(numTargetElems * numSrcElems));
      LOG(1, eff.second << " non-zero elements of " << eff.first << " total : filter efficiency = "
          << double(eff.second) / double(eff.first));

      LOG(1, "Intersection calculation done. " << std::endl );

    }

    /**
     * Counts the number of elements in an intersection matrix, and the number of these which are non-zero.
     *
     * @param m  the intersection matrix
     * @return  pair<int, int> containing as its first element the number of elements in m and as its second element the
     *                         number these which are non-zero
     */
    std::pair<int,int> countNumberOfMatrixEntries(const IntersectionMatrix& m)
    {

      int numElems = 0;
      int numNonZero = 0;
      for(IntersectionMatrix::const_iterator iter = m.begin() ; iter != m.end() ; ++iter)
        {
          numElems += iter->size();
          for(std::map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
            {
              if(!INTERP_KERNEL::epsilonEqual(iter2->second, 0.0, VOL_PREC))
                {
                  ++numNonZero;
                }
            }
        }
      return std::make_pair(numElems, numNonZero);
    }

  };
}

/**
 * Main method of the program.
 * Intersects the meshes and outputs some information about the calculation as well as the
 * intersection matrix on std::cout.
 *
 * @param argc  number of arguments given to the program (should be 3, the user giving 2 mesh names)
 * @param argv  vector to the arguments as strings.
 */
int main(int argc, char** argv)
{
  using INTERP_TEST::PerfTestToolkit;

  assert(argc == 3);

  // load meshes
  const std::string mesh1 = argv[1];
  const std::string mesh2 = argv[2];

  const std::string mesh1path = mesh1 + ".med";
  const std::string mesh2path = mesh2 + ".med";

  IntersectionMatrix m;

  PerfTestToolkit testTools;

  testTools.calcIntersectionMatrix(mesh1path.c_str(), mesh1.c_str(), mesh2path.c_str(), mesh2.c_str(), m);

  testTools.dumpIntersectionMatrix(m);

  return 0;

}

