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
#ifndef __SPLITTERTETRA_TXX__
#define __SPLITTERTETRA_TXX__

#include "SplitterTetra.hxx"

#include "TetraAffineTransform.hxx"
#include "TransformedTriangle.hxx"
#include "MeshUtils.hxx"
#include "VectorUtils.hxx"
#include "CellModel.hxx"
#include "Log.hxx"
#include "UnitTetraIntersectionBary.hxx"

#include <cmath>
#include <cassert>
#include <string>
#include <sstream>

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
  /*template<class MyMeshType>
  SplitterTetra<MyMeshType>::SplitterTetra(const MyMeshType& srcMesh, const MyMeshType& targetMesh, typename MyMeshType::MyConnType targetCell)
    : _src_mesh(srcMesh), _t(0)
  {   
    // get array of corners of target tetraedron
    const double* tetraCorners[4];
    for(int i = 0 ; i < 4 ; ++i)
      tetraCorners[i] = getCoordsOfNode(i, targetCell, targetMesh);
    // create the affine transform
    createAffineTransform(tetraCorners);
    }*/

  /*!
   * output is expected to be allocated with 24*sizeof(void*) in order to store the 24 tetras.
   * These tetras have to be deallocated.
   */
  template<class MyMeshType>
  void SplitterTetra<MyMeshType>::splitIntoDualCells(SplitterTetra<MyMeshType> **output)
  {
    double tmp[12];
    const double *tmp2[4]={tmp,tmp+3,tmp+6,tmp+9};
    typename MyMeshType::MyConnType conn[4]={-1,-1,-1,-1};
    for(int i=0;i<24;i++)
      {
        splitMySelfForDual(tmp,i,conn[0]);
        output[i]=new SplitterTetra<MyMeshType>(_src_mesh,tmp2,conn);
      }
  }

  /**
   * Constructor creating object from the four corners of the tetrahedron.
   *
   * @param srcMesh       mesh containing the source elements
   * @param tetraCorners  array of four pointers to double[3] arrays containing the coordinates of the
   *                      corners of the tetrahedron
   */
  template<class MyMeshType>
  SplitterTetra<MyMeshType>::SplitterTetra(const MyMeshType& srcMesh, const double** tetraCorners, const typename MyMeshType::MyConnType *nodesId)
    : _t(0), _src_mesh(srcMesh)
  {
    std::copy(nodesId,nodesId+4,_conn);
    _coords[0]=tetraCorners[0][0]; _coords[1]=tetraCorners[0][1]; _coords[2]=tetraCorners[0][2];
    _coords[3]=tetraCorners[1][0]; _coords[4]=tetraCorners[1][1]; _coords[5]=tetraCorners[1][2];
    _coords[6]=tetraCorners[2][0]; _coords[7]=tetraCorners[2][1]; _coords[8]=tetraCorners[2][2];
    _coords[9]=tetraCorners[3][0]; _coords[10]=tetraCorners[3][1]; _coords[11]=tetraCorners[3][2];
    // create the affine transform
    createAffineTransform(tetraCorners);
  }

  /*!
   * This contructor is used to build part of 1/24th dual cell of tetraCorners.
   * @param i is in 0..23 included.
   * @param nodeId is the id of first node of this in target mesh in C mode.
   */
  /*template<class MyMeshType>
  SplitterTetra<MyMeshType>::SplitterTetra(const MyMeshType& srcMesh, const double** tetraCorners, int i, typename MyMeshType::MyConnType nodeId)
    : _t(0), _src_mesh(srcMesh), _conn(nodeId)
  {
    double *newCoords[4];
    splitMySelfForDual(tetraCorners,newCoords,i);
    createAffineTransform(newCoords);
    }*/

  /**
   * Destructor
   *
   * Deletes _t and the coordinates (double[3] vectors) in _nodes
   *
   */
  template<class MyMeshType>
  SplitterTetra<MyMeshType>::~SplitterTetra()
  {
    delete _t;
    for(hash_map< int, double* >::iterator iter = _nodes.begin(); iter != _nodes.end() ; ++iter)
      delete[] iter->second;
  }

  /*!
   * \Forget already calculated triangles, which is crucial for calculation of barycenter of intersection
   */
  template<class MyMeshType>
  void SplitterTetra<MyMeshType>::clearVolumesCache()
  {
    _volumes.clear();
  }

  /*!
   * This method destroys the 4 pointers pointed by tetraCorners[0],tetraCorners[1],tetraCorners[2] and tetraCorners[3]
   * @param i is in 0..23 included.
   * @param output is expected to be sized of 12 in order to.
   */
  template<class MyMeshType>
  void SplitterTetra<MyMeshType>::splitMySelfForDual(double* output, int i, typename MyMeshType::MyConnType& nodeId)
  {
    double *tmp[4];
    int offset=i/6;
    nodeId=_conn[offset];
    tmp[0]=_coords+3*offset; tmp[1]=_coords+((offset+1)%4)*3; tmp[2]=_coords+((offset+2)%4)*3; tmp[3]=_coords+((offset+3)%4)*3;
    int caseToTreat=i%6;
    int case1=caseToTreat/2;
    int case2=caseToTreat%2;
    const int tab[3][2]={{1,2},{3,2},{1,3}};
    const int *curTab=tab[case1];
    double pt0[3]; pt0[0]=(tmp[curTab[case2]][0]+tmp[0][0])/2.; pt0[1]=(tmp[curTab[case2]][1]+tmp[0][1])/2.; pt0[2]=(tmp[curTab[case2]][2]+tmp[0][2])/2.;
    double pt1[3]; pt1[0]=(tmp[0][0]+tmp[curTab[0]][0]+tmp[curTab[1]][0])/3.; pt1[1]=(tmp[0][1]+tmp[curTab[0]][1]+tmp[curTab[1]][1])/3.; pt1[2]=(tmp[0][2]+tmp[curTab[0]][2]+tmp[curTab[1]][2])/3.;
    double pt2[3]; pt2[0]=(tmp[0][0]+tmp[1][0]+tmp[2][0]+tmp[3][0])/4.; pt2[1]=(tmp[0][1]+tmp[1][1]+tmp[2][1]+tmp[3][1])/4.; pt2[2]=(tmp[0][2]+tmp[1][2]+tmp[2][2]+tmp[3][2])/4.;
    std::copy(pt1,pt1+3,output+case2*3);
    std::copy(pt0,pt0+3,output+(abs(case2-1))*3);
    std::copy(pt2,pt2+3,output+2*3);
    std::copy(tmp[0],tmp[0]+3,output+3*3);
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
   * @param element      global number of the source element in C mode.
   */
  template<class MyMeshType>
  double SplitterTetra<MyMeshType>::intersectSourceCell(typename MyMeshType::MyConnType element,
                                                        double*                         baryCentre)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    const NumberingPolicy numPol=MyMeshType::My_numPol;
    //{ could be done on outside?
    // check if we have planar tetra element
    if(_t->determinant() == 0.0)
      {
        // tetra is planar
        LOG(2, "Planar tetra -- volume 0");
        return 0.0;
      }

    // get type of cell
    NormalizedCellType normCellType=_src_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(element));
    const CellModel& cellModelCell=CellModel::getCellModel(normCellType);
    unsigned nbOfNodes4Type=cellModelCell.isDynamic() ? _src_mesh.getNumberOfNodesOfElement(OTT<ConnType,numPol>::indFC(element)) : cellModelCell.getNumberOfNodes();
    // halfspace filtering
    bool isOutside[8] = {true, true, true, true, true, true, true, true};
    bool isTargetOutside = false;

    // calculate the coordinates of the nodes
    int *cellNodes=new int[nbOfNodes4Type];
    for(int i = 0;i<(int)nbOfNodes4Type;++i)
      {
        // we could store mapping local -> global numbers too, but not sure it is worth it
        const int globalNodeNum = getGlobalNumberOfNode(i, OTT<ConnType,numPol>::indFC(element), _src_mesh);
        cellNodes[i]=globalNodeNum;
        if(_nodes.find(globalNodeNum) == _nodes.end()) 
          {
            //for(hash_map< int , double* >::iterator iter3=_nodes.begin();iter3!=_nodes.end();iter3++)
            //  std::cout << (*iter3).first << " ";
            //std::cout << std::endl << "*** " << globalNodeNum << std::endl;
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
        /// calculator of intersection barycentre
        UnitTetraIntersectionBary baryCalculator( _t->determinant() < 0.);

        // get nb of sons of a cell
        const ConnType* rawCellConn = _src_mesh.getConnectivityPtr() + OTT<ConnType,numPol>::conn2C( _src_mesh.getConnectivityIndexPtr()[ element ]);
        const int rawNbCellNodes = _src_mesh.getConnectivityIndexPtr()[ element+1 ] - _src_mesh.getConnectivityIndexPtr()[ element ];
        unsigned nbOfSons = cellModelCell.getNumberOfSons2(rawCellConn, rawNbCellNodes);

        for(unsigned ii = 0 ; ii < nbOfSons; ++ii)
          {
            // get sons connectivity
            NormalizedCellType faceType;
            int *faceNodes, nbFaceNodes;
            if ( cellModelCell.isDynamic() )
              {
                faceNodes=new int[nbOfNodes4Type];
                nbFaceNodes = cellModelCell.fillSonCellNodalConnectivity2(ii,rawCellConn,rawNbCellNodes,faceNodes,faceType);
                for ( int i = 0; i < nbFaceNodes; ++i )
                  faceNodes[i] = OTT<ConnType,numPol>::coo2C(faceNodes[i]);
              }
            else
              {
                faceType = cellModelCell.getSonType(ii);
                const CellModel& faceModel=CellModel::getCellModel(faceType);
                assert(faceModel.getDimension() == 2);
                faceNodes=new int[faceModel.getNumberOfNodes()];      
                cellModelCell.fillSonCellNodalConnectivity(ii,cellNodes,faceNodes);
              }
            // intersect a son with the unit tetra
            switch(faceType)
              {
              case NORM_TRI3:
                {
                  // create the face key
                  TriangleFaceKey key = TriangleFaceKey(faceNodes[0], faceNodes[1], faceNodes[2]);

                  // calculate the triangle if needed
                  if(_volumes.find(key) == _volumes.end())
                    {
                      TransformedTriangle tri(_nodes[faceNodes[0]], _nodes[faceNodes[1]], _nodes[faceNodes[2]]);
                      calculateVolume(tri, key);
                      totalVolume += _volumes[key];
                      if ( baryCentre )
                        baryCalculator.addSide( tri );
                    } else {    
                      // count negative as face has reversed orientation
                      totalVolume -= _volumes[key];
                    }
                }
                break;

              case NORM_QUAD4:

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

              case NORM_POLYGON:
                {
                  int nbTria = nbFaceNodes - 2; // split polygon into nbTria triangles
                  for ( int iTri = 0; iTri < nbTria; ++iTri )
                    {
                      TriangleFaceKey key = TriangleFaceKey(faceNodes[0], faceNodes[1+iTri], faceNodes[2+iTri]);
                      if(_volumes.find(key) == _volumes.end())
                        {
                          TransformedTriangle tri(_nodes[faceNodes[0]], _nodes[faceNodes[1+iTri]], _nodes[faceNodes[2+iTri]]);
                          calculateVolume(tri, key);
                          totalVolume += _volumes[key];
                        }
                      else
                        {
                          totalVolume -= _volumes[key];
                        }
                    }
                }
                break;

              default:
                std::cout << "+++ Error : Only elements with triangular and quadratilateral faces are supported at the moment." << std::endl;
                assert(false);
              }
            delete [] faceNodes;
          }

        if ( baryCentre ) {
          baryCalculator.getBary( baryCentre );
          _t->reverseApply( baryCentre, baryCentre );
        }
      }
    delete [] cellNodes;
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

  /**
   * Calculates the volume of intersection of this tetrahedron with another one.
   */
  template<class MyMeshType>
  double SplitterTetra<MyMeshType>::intersectTetra(const double** tetraCorners)
  {
    //{ could be done on outside?
    // check if we have planar tetra element
    if(_t->determinant() == 0.0)
    {
      // tetra is planar
      LOG(2, "Planar tetra -- volume 0");
      return 0.0;
    }

    const unsigned nbOfNodes4Type=4;
    // halfspace filtering
    bool isOutside[8] = {true, true, true, true, true, true, true, true};
    bool isTargetOutside = false;

    // calculate the transformed coordinates of the nodes
    double nodes[nbOfNodes4Type][3];
    for(int i = 0;i<(int)nbOfNodes4Type;++i)
    {
      _t->apply(nodes[i], tetraCorners[i]);
      checkIsOutside(nodes[i], isOutside);
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
      const CellModel& cellModelCell=CellModel::getCellModel(NORM_TETRA4);
      int cellNodes[4] = { 0, 1, 2, 3 }, faceNodes[3];

      for(unsigned ii = 0 ; ii < 4 ; ++ii)
      {
        cellModelCell.fillSonCellNodalConnectivity(ii,cellNodes,faceNodes);
        
        TransformedTriangle tri(nodes[faceNodes[0]], nodes[faceNodes[1]], nodes[faceNodes[2]]);
        double vol = tri.calculateIntersectionVolume();
        totalVolume += vol;
      }
      
      // reset if it is very small to keep the matrix sparse
      // is this a good idea?
      if(epsilonEqual(totalVolume, 0.0, SPARSE_TRUNCATION_LIMIT))
      {
        totalVolume = 0.0;
      }
    }
    LOG(2, "Volume = " << totalVolume << ", det= " << _t->determinant());

    // NB : fault in article, Grandy, [8] : it is the determinant of the inverse transformation 
    // that should be used (which is equivalent to dividing by the determinant)
    return std::fabs(1.0 / _t->determinant() * totalVolume) ;
  }

  ////////////////////////////////////////////////////////

  template<class MyMeshTypeT, class MyMeshTypeS>
  SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::SplitterTetra2(const MyMeshTypeT& targetMesh, const MyMeshTypeS& srcMesh, SplittingPolicy policy)
    :_target_mesh(targetMesh),_src_mesh(srcMesh),_splitting_pol(policy)
  {
  }

  template<class MyMeshTypeT, class MyMeshTypeS>
  SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::~SplitterTetra2()
  {
    releaseArrays();
  }

  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::releaseArrays()
  {
    // free potential sub-mesh nodes that have been allocated
    typename MyMeshTypeT::MyConnType nbOfNodesT = _node_ids.size();// Issue 0020634.
    if(_nodes.size()>=/*8*/nbOfNodesT)
      {
        std::vector<const double*>::iterator iter = _nodes.begin() + /*8*/nbOfNodesT;
        while(iter != _nodes.end())
          {
            delete[] *iter;
            ++iter;
          }
      }
    _nodes.clear();
  }

  /*!
   * @param targetCell in C mode.
   * @param tetra is the output result tetra containers.
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::splitTargetCell(typename MyMeshTypeT::MyConnType targetCell,
                                                                 typename MyMeshTypeT::MyConnType nbOfNodesT,
                                                                 typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    typedef typename MyMeshTypeT::MyConnType ConnType;
    const NumberingPolicy numPol=MyMeshTypeT::My_numPol;
    const int numTetra = static_cast<int>(_splitting_pol);
    if(nbOfNodesT==4)
      {
        _nodes.resize(8);
        _node_ids.resize(8);
        tetra.reserve(1);
        const double *nodes[4];
        int conn[4];
        for(int node = 0; node < 4 ; ++node)
          {
            nodes[node]=getCoordsOfNode2(node, OTT<ConnType,numPol>::indFC(targetCell),_target_mesh,conn[node]);
          }
        std::copy(conn,conn+4,_node_ids.begin());
        SplitterTetra<MyMeshTypeS>* t = new SplitterTetra<MyMeshTypeS>(_src_mesh, nodes,conn);
        tetra.push_back(t);
        return ;
      }
    // Issue 0020634. To pass nbOfNodesT to calculateSubNodes (don't want to add an arg)
    _node_ids.resize(nbOfNodesT);

    // pre-calculate nodes
    calculateSubNodes(_target_mesh, OTT<ConnType,numPol>::indFC(targetCell));

    tetra.reserve(numTetra);
    _nodes.reserve(30); // we never have more than this

    switch ( nbOfNodesT )
      {
      case 8:
        {
          switch(_splitting_pol)
            {
            case PLANAR_FACE_5:
              {
                const int subZone[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
                fiveSplit(subZone,tetra);
              }
              break;

            case PLANAR_FACE_6:
              {
                const int subZone[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
                sixSplit(subZone,tetra);
              }
              break;

            case GENERAL_24:
              {
                calculateGeneral24Tetra(tetra);
              }
              break;

            case GENERAL_48:
              {
                calculateGeneral48Tetra(tetra);
              }
              break;
            default:
              assert(false);
            }
          break;
        }
      case 5:
        {
          splitPenta5(tetra);
          break;
        }
      default:
        {
          splitConvex(targetCell, tetra);
        }
      }
  }

  /**
   * Splits the hexahedron into five tetrahedra.
   * This method adds five SplitterTetra objects to the vector tetra. 
   *
   * @param subZone  the local node numbers corresponding to the hexahedron corners - these are mapped onto {0,..,7}. Providing this allows the 
   *                 splitting to be reused on the subzones of the GENERAL_* types of splitting
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::fiveSplit(const int* const subZone, typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    // Schema according to which the splitting is performed.
    // Each line represents one tetrahedron. The numbering is as follows :
    //
    //          7 ------ 6
    //         /|       /|
    //        / |      / |
    //       3 ------ 2  |
    //       |  |     |  |
    //       |  |     |  |
    //       |  4-----|- 5
    //       | /      | /
    //       0 ------ 1


    static const int SPLIT_NODES_5[20] = 
      {
        0, 1, 5, 2,
        0, 4, 5, 7,
        0, 3, 7, 2,
        5, 6, 7, 2,
        0, 2, 5, 7
      };
    
    // create tetrahedra
    for(int i = 0; i < 5; ++i)
      {
        const double* nodes[4];
        int conn[4];
        for(int j = 0; j < 4; ++j)
          {
            nodes[j] = getCoordsOfSubNode2(subZone[ SPLIT_NODES_5[4*i+j] ],conn[j]);
          }
        SplitterTetra<MyMeshTypeS>* t = new SplitterTetra<MyMeshTypeS>(_src_mesh, nodes,conn);
        tetra.push_back(t);
      }
  }

  /**
   * Splits the hexahedron into six tetrahedra.
   * This method adds six SplitterTetra objects to the vector tetra. 
   *
   * @param subZone  the local node numbers corresponding to the hexahedron corners - these are mapped onto {0,..,7}. Providing this allows the 
   *                 splitting to be reused on the subzones of the GENERAL_* types of splitting
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::sixSplit(const int* const subZone, typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    // Schema according to which the splitting is performed.
    // Each line represents one tetrahedron. The numbering is as follows :
    //
    //          7 ------ 6
    //         /|       /|
    //        / |      / |
    //       3 ------ 2  |
    //       |  |     |  |
    //       |  |     |  |
    //       |  4-----|- 5
    //       | /      | /
    //       0 ------ 1

    static const int SPLIT_NODES_6[24] = 
      {
        0, 1, 5, 6,
        0, 2, 1, 6,
        0, 5, 4, 6,
        0, 4, 7, 6,
        0, 3, 2, 6,
        0, 7, 3, 6
      };

    for(int i = 0; i < 6; ++i)
      {
        const double* nodes[4];
        int conn[4];
        for(int j = 0; j < 4; ++j)
          {
            nodes[j] = getCoordsOfSubNode2(subZone[ SPLIT_NODES_6[4*i+j] ],conn[j]);
          }
        SplitterTetra<MyMeshTypeS>* t = new SplitterTetra<MyMeshTypeS>(_src_mesh, nodes,conn);
        tetra.push_back(t);
      }
  }

  /**
   * Splits the hexahedron into 24 tetrahedra.
   * The splitting is done by combining the barycenter of the tetrahedron, the barycenter of each face 
   * and the nodes of each edge of the face. This creates 6 faces * 4 edges / face = 24 tetrahedra.
   * The submesh nodes introduced are the barycenters of the faces and the barycenter of the cell.
   * 
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::calculateGeneral24Tetra(typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    // The two nodes of the original mesh cell used in each tetrahedron.
    // The tetrahedra all have nodes (cellCenter, faceCenter, edgeNode1, edgeNode2)
    // For the correspondance of the nodes, see the GENERAL_48_SUB_NODES table in calculateSubNodes
    static const int TETRA_EDGES[48] = 
      {
        // face with center 9
        0,1,
        1,5,
        5,4,
        4,0,
        // face with center 10
        0,1,
        1,2,
        2,3,
        3,0,
        // face with center 11
        0,4,
        4,7,
        7,3,
        3,0,
        // face with center 12
        1,5,
        5,6,
        6,2,
        2,1,
        // face with center 13
        5,6,
        6,7,
        7,4,
        4,5,
        // face with center 14
        2,6,
        6,7,
        7,3,
        3,2
      };
    
    // nodes to use for tetrahedron
    const double* nodes[4];
    int conn[4];
    // get the cell center
    nodes[0] = getCoordsOfSubNode2(14,conn[0]);
    
    for(int faceCenterNode = 8; faceCenterNode < 14; ++faceCenterNode)
      {
        // get the face center
        nodes[1] = getCoordsOfSubNode2(faceCenterNode,conn[1]);
        for(int j = 0; j < 4; ++j)
          {
            const int row = 4*(faceCenterNode - 9) + j;
            nodes[2] = getCoordsOfSubNode2(TETRA_EDGES[2*row],conn[2]);
            nodes[3] = getCoordsOfSubNode2(TETRA_EDGES[2*row + 1],conn[3]);
           
            SplitterTetra<MyMeshTypeS>* t = new SplitterTetra<MyMeshTypeS>(_src_mesh, nodes, conn);
            tetra.push_back(t);
          }
      }
  }


  /**
   * Splits the hexahedron into 48 tetrahedra.
   * The splitting is done by introducing the midpoints of all the edges 
   * and the barycenter of the element as submesh nodes. The 8 hexahedral subzones thus defined
   * are then split into 6 tetrahedra each, as in Grandy, p. 449. The division of the subzones 
   * is done by calling sixSplit().
   * 
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::calculateGeneral48Tetra(typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    // Define 8 hexahedral subzones as in Grandy, p449
    // the values correspond to the nodes that correspond to nodes 1,2,3,4,5,6,7,8 in the subcell
    // For the correspondance of the nodes, see the GENERAL_48_SUB_NODES table in calculateSubNodes
    static const int subZones[64] = 
      {
        0,8,21,12,9,20,26,22,
        8,1,13,21,20,10,23,26,
        12,21,16,3,22,26,25,17,
        21,13,2,16,26,23,18,25,
        9,20,26,22,4,11,24,14,
        20,10,23,26,11,5,15,24,
        22,26,25,17,14,24,19,7,
        26,23,18,25,24,15,6,19
      };
    
    for(int i = 0; i < 8; ++i)
      {
        sixSplit(&subZones[8*i],tetra);
      }
  }
  
  /**
   * Splits the NORM_PYRA5 into 2 tetrahedra.
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::splitPenta5(typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    static const int SPLIT_PENTA_5[2][4] = 
      {
        {
          0, 1, 2, 4
        },
        {
          0, 2, 3, 4
        }
      };
    
    // create tetrahedra
    const double* nodes[4];
    int conn[4];
    for(int i = 0; i < 2; ++i)
      {
        for(int j = 0; j < 4; ++j)
          nodes[j] = getCoordsOfSubNode2(SPLIT_PENTA_5[i][j],conn[j]);
        SplitterTetra<MyMeshTypeS>* t = new SplitterTetra<MyMeshTypeS>(_src_mesh, nodes,conn);
        tetra.push_back(t);
      }
  }
  
  /**
   * Splits a convex cell into tetrahedra.
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::splitConvex(typename MyMeshTypeT::MyConnType targetCell,
                                                             typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    // Each face of a cell is split into triangles and
    // each of triangles and a cell barycenter form a tetrahedron.

    typedef typename MyMeshTypeT::MyConnType ConnType;
    const NumberingPolicy numPol=MyMeshTypeT::My_numPol;

    // get type of cell and nb of cell nodes
    NormalizedCellType normCellType=_target_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(targetCell));
    const CellModel& cellModelCell=CellModel::getCellModel(normCellType);
    unsigned nbOfCellNodes=cellModelCell.isDynamic() ? _target_mesh.getNumberOfNodesOfElement(OTT<ConnType,numPol>::indFC(targetCell)) : cellModelCell.getNumberOfNodes();

    // get nb of cell sons (faces)
    const ConnType* rawCellConn = _target_mesh.getConnectivityPtr() + OTT<ConnType,numPol>::conn2C( _target_mesh.getConnectivityIndexPtr()[ targetCell ]);
    const int rawNbCellNodes = _target_mesh.getConnectivityIndexPtr()[ targetCell+1 ] - _target_mesh.getConnectivityIndexPtr()[ targetCell ];
    unsigned nbOfSons = cellModelCell.getNumberOfSons2(rawCellConn, rawNbCellNodes);

    // indices of nodes of a son
    static vector<int> allNodeIndices; // == 0,1,2,...,nbOfCellNodes-1
    while ( allNodeIndices.size() < nbOfCellNodes )
      allNodeIndices.push_back( allNodeIndices.size() );
    vector<int> classicFaceNodes(4);
    int* faceNodes = cellModelCell.isDynamic() ? &allNodeIndices[0] : &classicFaceNodes[0];

    // nodes of tetrahedron
    int conn[4];
    const double* nodes[4];
    nodes[3] = getCoordsOfSubNode2( nbOfCellNodes,conn[3]); // barycenter

    for(unsigned ii = 0 ; ii < nbOfSons; ++ii)
      {
        // get indices of son's nodes: it's just next portion of allNodeIndices for polyhedron
        // and some of allNodeIndices accodring to cell model for a classsic cell 
        unsigned nbFaceNodes = cellModelCell.getNumberOfNodesConstituentTheSon2(ii, rawCellConn, rawNbCellNodes);
        if ( normCellType != NORM_POLYHED )
          cellModelCell.fillSonCellNodalConnectivity(ii,&allNodeIndices[0],faceNodes);

        int nbTetra = nbFaceNodes - 2; // split polygon into nbTetra triangles

        // create tetrahedra
        for(int i = 0; i < nbTetra; ++i)
          {
            nodes[0] = getCoordsOfSubNode2( faceNodes[0],  conn[0]);
            nodes[1] = getCoordsOfSubNode2( faceNodes[1+i],conn[1]);
            nodes[2] = getCoordsOfSubNode2( faceNodes[2+i],conn[2]);
            SplitterTetra<MyMeshTypeS>* t = new SplitterTetra<MyMeshTypeS>(_src_mesh, nodes,conn);
            tetra.push_back(t);
          }

        if ( normCellType == NORM_POLYHED )
          faceNodes += nbFaceNodes; // go to the next face
      }
  }
  
  /**
   * Precalculates all the nodes.
   * Retrieves the mesh nodes and allocates the necessary sub-mesh 
   * nodes according to the splitting policy used.
   * This method is meant to be called once by the constructor.
   *
   * @param targetMesh  the target mesh
   * @param targetCell  the global number of the cell that the object represents, in targetMesh mode.
   * @param policy      the splitting policy of the object
   *
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::calculateSubNodes(const MyMeshTypeT& targetMesh, typename MyMeshTypeT::MyConnType targetCell)
  {
    // retrieve real mesh nodes
    
    typename MyMeshTypeT::MyConnType nbOfNodesT = _node_ids.size();// Issue 0020634. _node_ids.resize(8);
    for(int node = 0; node < nbOfNodesT ; ++node)
      {
        // calculate only normal nodes
        _nodes.push_back(getCoordsOfNode2(node, targetCell, targetMesh,_node_ids[node]));
      }

    switch ( nbOfNodesT )
      {
      case 8:

        // create sub-mesh nodes if needed
        switch(_splitting_pol)
          {
          case GENERAL_24:
            {
              // Each sub-node is the barycenter of 4 other nodes.
              // For the faces, these are on the orignal mesh.
              // For the barycenter, the four face sub-nodes are used.
              static const int GENERAL_24_SUB_NODES[28] = 
                {
                  0,1,4,5,// sub-node 9 (face)
                  0,1,2,3,// sub-node 10 (face)
                  0,3,4,7,// sub-node 11 (face)
                  1,2,5,6,// sub-node 12 (face)
                  4,5,6,7,// sub-node 13 (face)
                  2,3,6,7,// sub-node 14 (face)
                  9,10,11,12// sub-node 15 (cell)
                };

              for(int i = 0; i < 7; ++i)
                {
                  double* barycenter = new double[3];
                  calcBarycenter(4, barycenter, &GENERAL_24_SUB_NODES[4*i]);
                  _nodes.push_back(barycenter);
                }
            }
            break;

          case GENERAL_48:
            {
              // Each sub-node is the barycenter of two other nodes.
              // For the edges, these lie on the original mesh.
              // For the faces, these are the edge sub-nodes.
              // For the cell these are two face sub-nodes.
              static const int GENERAL_48_SUB_NODES[38] = 
                {
                  0,1,   // sub-node 9 (edge)
                  0,4,   // sub-node 10 (edge)
                  1,5,   // sub-node 11 (edge)
                  4,5,   // sub-node 12 (edge)
                  0,3,   // sub-node 13 (edge)
                  1,2,   // sub-node 14 (edge)
                  4,7,   // sub-node 15 (edge)
                  5,6,   // sub-node 16 (edge)
                  2,3,   // sub-node 17 (edge)
                  3,7,   // sub-node 18 (edge)
                  2,6,   // sub-node 19 (edge)
                  6,7,   // sub-node 20 (edge)
                  8,11,  // sub-node 21 (face)
                  12,13, // sub-node 22 (face)
                  9,17,  // sub-node 23 (face)
                  10,18, // sub-node 24 (face)
                  14,15, // sub-node 25 (face)
                  16,19, // sub-node 26 (face)
                  20,25  // sub-node 27 (cell)
                };

              for(int i = 0; i < 19; ++i)
                {
                  double* barycenter = new double[3];
                  calcBarycenter(2, barycenter, &GENERAL_48_SUB_NODES[2*i]);
                  _nodes.push_back(barycenter);
                }
            }
            break;

          default:
            break;
          }

      case 5: // NORM_PYRA5
        break;

      default: // convex 3d cell
        {
          // add barycenter of a cell
          vector<int> allIndices(nbOfNodesT);
          for ( int i = 0; i < nbOfNodesT; ++i ) allIndices[i] = i;
          double* barycenter = new double[3];
          calcBarycenter(nbOfNodesT, barycenter, &allIndices[0]);
          _nodes.push_back(barycenter);
        }
      }

  }
}

#endif
