#ifndef __INTERSECTORTETRA_TXX__
#define __INTERSECTORTETRA_TXX__

#include "IntersectorTetra.hxx"

#include "TetraAffineTransform.hxx"
#include "TransformedTriangle.hxx"
#include "MeshUtils.hxx"
#include "VectorUtils.hxx"
#include "Log.hxx"

#include <cmath>
#include <cassert>
#include <string>
#include <sstream>

using namespace MEDMEM;
using namespace MED_EN;

/// Smallest volume of the intersecting elements in the transformed space that will be returned as non-zero. 
/// Since the scale is always the same in the transformed space (the target tetrahedron is unitary), this number is independent of the scale of the meshes.
#define SPARSE_TRUNCATION_LIMIT 1.0e-14

namespace INTERP_KERNEL
{
  
  /**
   * Constructor creating object from target cell global number 
   * 
   * @param srcMesh     mesh containing the source elements
   * @param targetMesh  mesh containing the target elements
   * @param targetCell  global number of the target cell
   *
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::IntersectorTetra(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh,
                                                                                  const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& targetMesh,
                                                                                  ConnType targetCell)
    : _srcMesh(srcMesh), _t(0)
  {   
    // get array of corners of target tetraedron
    const double* tetraCorners[4];
    for(int i = 0 ; i < 4 ; ++i)
      tetraCorners[i] = getCoordsOfNode(i, targetCell, targetMesh);
    // create the affine transform
    createAffineTransform(tetraCorners);
  }

  /**
   * Constructor creating object from the four corners of the tetrahedron.
   *
   * @param srcMesh       mesh containing the source elements
   * @param tetraCorners  array of four pointers to double[3] arrays containing the coordinates of the
   *                      corners of the tetrahedron
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::IntersectorTetra(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh,
                                                                                  const double** tetraCorners)
    : _srcMesh(srcMesh), _t(0)
  {
    // create the affine transform
    createAffineTransform(tetraCorners);
  }

  /**
   * Destructor
   *
   * Deletes _t and the coordinates (double[3] vectors) in _nodes
   *
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::~IntersectorTetra()
  {
    delete _t;
    for(hash_map< int, double* >::iterator iter = _nodes.begin(); iter != _nodes.end() ; ++iter)
      delete[] iter->second;
	
  }
  
  /**
   * Calculates the volume of intersection of an element in the source mesh and the target element.
   * It first calculates the transformation that takes the target tetrahedron into the unit tetrahedron. After that, the 
   * faces of the source element are triangulated and the calculated transformation is applied 
   * to each triangle. The algorithm of Grandy, implemented in INTERP_KERNEL::TransformedTriangle is used
   * to calculate the contribution to the volume from each triangle. The volume returned is the sum of these contributions
   * divided by the determinant of the transformation.
   *
   * The class will cache the intermediary calculations of transformed nodes of source cells and volumes associated 
   * with triangulated faces to avoid having to recalculate these.
   *
   * @param element      global number of the source element (1 <= srcCell < # source cells)
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  double IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::intersectSourceCell(ConnType element)
  {
        
    //{ could be done on outside?
    // check if we have planar tetra element
    if(_t->determinant() == 0.0)
      {
        // tetra is planar
        LOG(2, "Planar tetra -- volume 0");
        return 0.0;
      }

    // get type of cell
    unsigned char nbOfNodes4Type = _srcMesh.getNumberOfNodesOfElement(element);
    // hack tony for the moment
    
    medGeometryElement type;
    if(nbOfNodes4Type==4)
      type=MED_EN::MED_TETRA4;
    else if(nbOfNodes4Type==8)
      type=MED_EN::MED_HEXA8;
    else
      std::cerr << "ca merde sec" << std::endl;
    ////const medGeometryElement type = _srcMesh.getNumberOfNodesOfElement(element);
    
    // get cell model for the element
    const CELLMODEL& cellModel= CELLMODEL_Map::retrieveCellModel(type);

    // halfspace filtering
    bool isOutside[8] = {true, true, true, true, true, true, true, true};
    bool isTargetOutside = false;

    // calculate the coordinates of the nodes
    for(int i = 0;i<nbOfNodes4Type;++i)
      {
        // we could store mapping local -> global numbers too, but not sure it is worth it
        const int globalNodeNum = getGlobalNumberOfNode(i, element, _srcMesh);
        if(_nodes.find(globalNodeNum) == _nodes.end()) 
          {
            calculateNode(globalNodeNum);
          }

        checkIsOutside(_nodes[globalNodeNum], isOutside);       
      }

    // halfspace filtering check
    // NB : might not be beneficial for caching of triangles
    for(int i = 0; i < 8; ++i)
      {
        if(isOutside[i])
          {
            isTargetOutside = true;
          }
      }

    double totalVolume = 0.0;
    
    if(!isTargetOutside)
      {
        for(int i = 1 ; i <= cellModel.getNumberOfConstituents(1) ; ++i)
          {
            const medGeometryElement faceType = cellModel.getConstituentType(1, i);
            const CELLMODEL&  faceModel= CELLMODEL_Map::retrieveCellModel(faceType);
       
            assert(faceModel.getDimension() == 2);

            int faceNodes[faceModel.getNumberOfNodes()];

            // get the nodes of the face
            for(int j = 1; j <= faceModel.getNumberOfNodes(); ++j)
              {
                const int locNodeNum = cellModel.getNodeConstituent(1, i, j);
                assert(locNodeNum >= 1);
                assert(locNodeNum <= cellModel.getNumberOfNodes());

                faceNodes[j-1] = getGlobalNumberOfNode(locNodeNum-1, element, _srcMesh);
              }

            switch(faceType)
              {
              case MED_EN::MED_TRIA3:
                {
                  // create the face key
                  TriangleFaceKey key = TriangleFaceKey(faceNodes[0], faceNodes[1], faceNodes[2]);
             
                  // calculate the triangle if needed
                  if(_volumes.find(key) == _volumes.end())
                    {
                      TransformedTriangle tri(_nodes[faceNodes[0]], _nodes[faceNodes[1]], _nodes[faceNodes[2]]);
                      calculateVolume(tri, key);
                      totalVolume += _volumes[key];
                    } else {    
                    // count negative as face has reversed orientation
                    totalVolume -= _volumes[key];
                  }
                }
                break;

              case MED_EN::MED_QUAD4:

                // simple triangulation of faces along a diagonal :
                // 
                // 2 ------ 3
                // |      / |
                // |     /  |
                // |    /   |
                // |   /    |
                // |  /     |
                // | /      |
                // 1 ------ 4
                //
                //? not sure if this always works 
                {
                  // calculate the triangles if needed

                  // local nodes 1, 2, 3
                  TriangleFaceKey key1 = TriangleFaceKey(faceNodes[0], faceNodes[1], faceNodes[2]);
                  if(_volumes.find(key1) == _volumes.end())
                    {
                      TransformedTriangle tri(_nodes[faceNodes[0]], _nodes[faceNodes[1]], _nodes[faceNodes[2]]);
                      calculateVolume(tri, key1);
                      totalVolume += _volumes[key1];
                    } else {
                    // count negative as face has reversed orientation
                    totalVolume -= _volumes[key1];
                  }

                  // local nodes 1, 3, 4
                  TriangleFaceKey key2 = TriangleFaceKey(faceNodes[0], faceNodes[2], faceNodes[3]);
                  if(_volumes.find(key2) == _volumes.end())
                    {
                      TransformedTriangle tri(_nodes[faceNodes[0]], _nodes[faceNodes[2]], _nodes[faceNodes[3]]);
                      calculateVolume(tri, key2);
                      totalVolume += _volumes[key2];
                    }
                  else
                    { 
                      // count negative as face has reversed orientation
                      totalVolume -= _volumes[key2];
                    }
                }
                break;

              default:
                std::cout << "+++ Error : Only elements with triangular and quadratilateral faces are supported at the moment." << std::endl;
                assert(false);
              }
          }
      }

    // reset if it is very small to keep the matrix sparse
    // is this a good idea?
    if(epsilonEqual(totalVolume, 0.0, SPARSE_TRUNCATION_LIMIT))
      {
        totalVolume = 0.0;
      }
    
    LOG(2, "Volume = " << totalVolume << ", det= " << _t->determinant());

    // NB : fault in article, Grandy, [8] : it is the determinant of the inverse transformation 
    // that should be used (which is equivalent to dividing by the determinant)
    return std::fabs(1.0 / _t->determinant() * totalVolume) ;
  }
}

#endif
