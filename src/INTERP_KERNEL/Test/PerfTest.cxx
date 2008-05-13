#include "Interpolation3D.hxx"
#include "Interpolation3D.txx"
#include "MeshTestToolkit.txx"
#include "Log.hxx"
#include "VectorUtils.hxx"

#include "MEDMEM_Mesh.hxx"
#include "MEDNormalizedUnstructuredMesh.hxx"

#include <cassert>
#include <string>

using namespace MEDMEM;
using namespace MED_EN;

/**
 * \file PerfTest.cxx
 * Test program which takes two meshes and calculates their intersection matrix. 
 * 
 * USAGE : PerfTest mesh1 mesh2 
 *         where mesh1 and mesh2 are the names of two meshes located in
 *         the files mesh1.med, mesh2.med in {$MED_ROOT_DIR}/share/salome/resources/med/
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
     * @param  mesh1path   the path to the file containing the source mesh, relative to {$MED_ROOT_DIR}/share/salome/resources/med/
     * @param  mesh1       the name of the source mesh
     * @param  mesh2path   the path to the file containing the target mesh, relative to {$MED_ROOT_DIR}/share/salome/resources/med/
     * @param  mesh2       the name of the target mesh
     * @param  m           intersection matrix in which to store the result of the intersection
     */
    void calcIntersectionMatrix(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, IntersectionMatrix& m) 
    {
      const string dataBaseDir = getenv("MED_ROOT_DIR");
      const string dataDir = dataBaseDir + "/share/salome/resources/med/";

      LOG(1, std::endl << "=== -> intersecting src = " << mesh1 << ", target = " << mesh2 );
      
      LOG(5, "Loading " << mesh1 << " from " << mesh1path);
      const MESH sMesh(MED_DRIVER, dataDir+mesh1path, mesh1);
      const int numSrcElems = sMesh.getNumberOfElements(MED_CELL, MED_ALL_ELEMENTS);
      LOG(1, "Source mesh has " << numSrcElems << " elements");
    
    
      LOG(5, "Loading " << mesh2 << " from " << mesh2path);
      const MESH tMesh(MED_DRIVER, dataDir+mesh2path, mesh2);
      const int numTargetElems = tMesh.getNumberOfElements(MED_CELL, MED_ALL_ELEMENTS);
    
      LOG(1, "Target mesh has " << numTargetElems << " elements");
			
			MEDNormalizedUnstructuredMesh<3,3> sMesh_wrapper(&sMesh);
			MEDNormalizedUnstructuredMesh<3,3> tMesh_wrapper(&tMesh);
			
			Interpolation3D interpolator;
      interpolator.interpolateMeshes(sMesh_wrapper, tMesh_wrapper,m);
    
      std::pair<int, int> eff = countNumberOfMatrixEntries(m);
      LOG(1, eff.first << " of " << numTargetElems * numSrcElems << " intersections calculated : ratio = " 
	  << double(eff.first) / double(numTargetElems * numSrcElems));
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
	  for(map<int, double>::const_iterator iter2 = iter->begin() ; iter2 != iter->end() ; ++iter2)
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
  const string mesh1 = argv[1];
  const string mesh2 = argv[2];

  const string mesh1path = mesh1 + ".med";
  const string mesh2path = mesh2 + ".med";

  IntersectionMatrix m;

  PerfTestToolkit testTools;

  testTools.calcIntersectionMatrix(mesh1path.c_str(), mesh1.c_str(), mesh2path.c_str(), mesh2.c_str(), m);

  testTools.dumpIntersectionMatrix(m);
    
  return 0;

}

