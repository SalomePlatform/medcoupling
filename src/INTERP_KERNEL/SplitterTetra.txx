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
#include "VolSurfFormulae.hxx"

#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <vector>

namespace INTERP_KERNEL
{
  template<class MyMeshType> 
  const double SplitterTetra<MyMeshType>::SPARSE_TRUNCATION_LIMIT=1.0e-14;

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
   * SplitterTetra class computes for a list of cell ids of a given mesh \a srcMesh (badly named) the intersection with a 
   * single TETRA4 cell given by \a tetraCorners (of length 4) and \a nodesId (of length 4 too). \a nodedIds is given only to establish
   * if a partial computation of a triangle has already been performed (to increase performance).
   *
   * The \a srcMesh can contain polyhedron cells.
   * 
   * 
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
    _t=new TetraAffineTransform(_coords);
  }

  /**
   * SplitterTetra class computes for a list of cell ids of a given mesh \a srcMesh (badly named) the intersection with a 
   * single TETRA4 cell given by \a tetraCorners (of length 4) and \a nodesId (of length 4 too). \a nodedIds is given only to establish
   * if a partial computation of a triangle has already been performed (to increase performance).
   *
   * The \a srcMesh can contain polyhedron cells.
   * 
   * 
   * Constructor creating object from the four corners of the tetrahedron.
   *
   * \param [in] srcMesh       mesh containing the source elements
   * \param [in] tetraCorners  array 4*3 doubles containing corners of input tetrahedron (P0X,P0Y,P0Y,P1X,P1Y,P1Z,P2X,P2Y,P2Z,P3X,P3Y,P3Z).
   */
  template<class MyMeshType>
  SplitterTetra<MyMeshType>::SplitterTetra(const MyMeshType& srcMesh, const double tetraCorners[12], const int *conn): _t(0),_src_mesh(srcMesh)
  {
    if(!conn)
      { _conn[0]=0; _conn[1]=1; _conn[2]=2; _conn[3]=3; }
    else
      { _conn[0]=conn[0]; _conn[1]=conn[1]; _conn[2]=conn[2]; _conn[3]=conn[3]; }
    _coords[0]=tetraCorners[0]; _coords[1]=tetraCorners[1]; _coords[2]=tetraCorners[2]; _coords[3]=tetraCorners[3]; _coords[4]=tetraCorners[4]; _coords[5]=tetraCorners[5];
    _coords[6]=tetraCorners[6]; _coords[7]=tetraCorners[7]; _coords[8]=tetraCorners[8]; _coords[9]=tetraCorners[9]; _coords[10]=tetraCorners[10]; _coords[11]=tetraCorners[11];
    // create the affine transform
    _t=new TetraAffineTransform(_coords);
  }
  
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
    for(HashMap< int, double* >::iterator iter = _nodes.begin(); iter != _nodes.end() ; ++iter)
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
    const CellModel& cellModelCell=CellModel::GetCellModel(normCellType);
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
            //for(HashMap< int , double* >::iterator iter3=_nodes.begin();iter3!=_nodes.end();iter3++)
            //  std::cout << (*iter3).first << " ";
            //std::cout << std::endl << "*** " << globalNodeNum << std::endl;
            calculateNode(globalNodeNum);
          }
        CheckIsOutside(_nodes[globalNodeNum], isOutside);       
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
            int *faceNodes, nbFaceNodes=-1;
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
                const CellModel& faceModel=CellModel::GetCellModel(faceType);
                assert(faceModel.getDimension() == 2);
                nbFaceNodes = cellModelCell.getNumberOfNodesConstituentTheSon(ii);
                faceNodes = new int[nbFaceNodes];
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
   * Calculates the intersection surface of two coplanar triangles.
   *
   * @param palneNormal normal of the plane for the first triangle
   * @param planeConstant constant of the equation of the plane for the first triangle
   * @param p1 coordinates of the first  node of the first  triangle
   * @param p2 coordinates of the second node of the first  triangle
   * @param p3 coordinates of the third  node of the first  triangle
   * @param p4 coordinates of the first  node of the second triangle
   * @param p5 coordinates of the second node of the second triangle
   * @param p6 coordinates of the third  node of the second triangle
   * @param dimCaracteristic characteristic size of the meshes containing the triangles
   * @param precision precision for double float data used for comparison
   */
  template<class MyMeshType>
  double SplitterTetra<MyMeshType>::CalculateIntersectionSurfaceOfCoplanarTriangles(const double *const planeNormal,
                                                                                    const double planeConstant,
                                                                                    const double *const p1, const double *const p2, const double *const p3,
                                                                                    const double *const p4, const double *const p5, const double *const p6,
                                                                                    const double dimCaracteristic, const double precision)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    typedef double Vect2[2];
    typedef double Triangle2[3][2];

    const double *const tri0[3] = {p1, p2, p3};
    const double *const tri1[3] = {p4, p5, p6};

    // Plane of the first triangle defined by the normal of the triangle and the constant
    // Project triangles onto coordinate plane most aligned with plane normal
    int maxNormal = 0;
    double fmax = std::abs(planeNormal[0]);
    double absMax = std::abs(planeNormal[1]);
    if (absMax > fmax)
      {
        maxNormal = 1;
        fmax = absMax;
      }
    absMax = std::abs(planeNormal[2]);
    if (absMax > fmax)
      {
        maxNormal = 2;
      }

    Triangle2 projTri0, projTri1;
    int i;

    if (maxNormal == 0)
      {
        // Project onto yz-plane.
        for (i = 0; i < 3; ++i)
          {
            projTri0[i][0] = tri0[i][1];
            projTri0[i][1] = tri0[i][2];
            projTri1[i][0] = tri1[i][1];
            projTri1[i][1] = tri1[i][2];
          }
      }
    else if (maxNormal == 1)
      {
        // Project onto xz-plane.
        for (i = 0; i < 3; ++i)
          {
            projTri0[i][0] = tri0[i][0];
            projTri0[i][1] = tri0[i][2];
            projTri1[i][0] = tri1[i][0];
            projTri1[i][1] = tri1[i][2];
          }
      }
    else
      {
        // Project onto xy-plane.
        for (i = 0; i < 3; ++i)
          {
            projTri0[i][0] = tri0[i][0];
            projTri0[i][1] = tri0[i][1];
            projTri1[i][0] = tri1[i][0];
            projTri1[i][1] = tri1[i][1];
          }
      }

    // 2D triangle intersection routines require counterclockwise ordering.
    Vect2 save;
    Vect2 edge0;
    Vect2 edge1;
    for (int ii = 0; ii < 2; ++ii)
      {
        edge0[ii] = projTri0[1][ii] - projTri0[0][ii];
        edge1[ii] = projTri0[2][ii] - projTri0[0][ii];
      }
    if ((edge0[0] * edge1[1] - edge0[1] * edge1[0]) < (double) 0.)
      {
        // Triangle is clockwise, reorder it.
        for (int ii = 0; ii < 2; ++ii)
          {
            save[ii] = projTri0[1][ii];
            projTri0[1][ii] = projTri0[2][ii];
            projTri0[2][ii] = save[ii];
          }
      }

    for (int ii = 0; ii < 2; ++ii)
      {
        edge0[ii] = projTri1[1][ii] - projTri1[0][ii];
        edge1[ii] = projTri1[2][ii] - projTri1[0][ii];
      }
    if ((edge0[0] * edge1[1] - edge0[1] * edge1[0]) < (double) 0.)
      {
        // Triangle is clockwise, reorder it.
      for (int ii = 0; ii < 2; ++ii)
        {
          save[ii] = projTri1[1][ii];
          projTri1[1][ii] = projTri1[2][ii];
          projTri1[2][ii] = save[ii];
        }
      }

    std::vector<double> inter2;
    intersec_de_triangle(projTri0[0], projTri0[1], projTri0[2],
                         projTri1[0], projTri1[1], projTri1[2],
                         inter2,
                         dimCaracteristic, precision);
    ConnType nb_inter=((ConnType)inter2.size())/2;
    double surface = 0.;
    if(nb_inter >3) inter2=reconstruct_polygon(inter2);
    if (nb_inter > 0)
      {
        std::vector<double> inter3;
        inter3.resize(3 * nb_inter);
        // Map 2D intersections back to the 3D triangle space.
        if (maxNormal == 0)
          {
            double invNX = ((double) 1.) / planeNormal[0];
            for (i = 0; i < nb_inter; i++)
              {
                inter3[3 * i + 1] = inter2[2 * i];
                inter3[3 * i + 2] = inter2[2 * i + 1];
                inter3[3 * i] = invNX * (planeConstant - planeNormal[1] * inter3[3 * i + 1] - planeNormal[2] * inter3[3 * i + 2]);
              }
          }
        else if (maxNormal == 1)
          {
            double invNY = ((double) 1.) / planeNormal[1];
            for (i = 0; i < nb_inter; i++)
              {
                inter3[3 * i] = inter2[2 * i];
                inter3[3 * i + 2] = inter2[2 * i + 1];
                inter3[3 * i + 1] = invNY * (planeConstant - planeNormal[0] * inter3[3 * i] - planeNormal[2] * inter3[3 * i + 2]);
              }
          }
        else
          {
            double invNZ = ((double) 1.) / planeNormal[2];
            for (i = 0; i < nb_inter; i++)
              {
                inter3[3 * i] = inter2[2 * i];
                inter3[3 * i + 1] = inter2[2 * i + 1];
                inter3[3 * i + 2] = invNZ * (planeConstant - planeNormal[0] * inter3[3 * i] - planeNormal[1] * inter3[3 * i + 1]);
              }
          }
        surface = polygon_area<3>(inter3);
      }
    return surface;
  }

  /**
   * Determine if a face is coplanar with a triangle.
   * The first face is characterized by the equation of her plane
   *
   * @param palneNormal normal of the plane for the first triangle
   * @param planeConstant constant of the equation of the plane for the first triangle
   * @param coordsFace coordinates of the triangle face
   * @param precision precision for double float data used for comparison
   */
  template<class MyMeshType>
  bool SplitterTetra<MyMeshType>::IsFacesCoplanar(const double *const planeNormal,
                                                  const double planeConstant,
                                                  const double *const *const coordsFace,
                                                  const double precision)
  {
      // Compute the signed distances of triangle vertices to the plane. Use an epsilon-thick plane test.
      // For faces not left
      int counter = 0;
      for (int i = 0; i < 3; ++i)
        {
          const double distance = dot(planeNormal, coordsFace[i]) - planeConstant;
          if (epsilonEqual(distance, precision))
            {
              counter++;
            }
        }
      return counter == 3;
  }

  /**
   * Calculates the surface of intersection of a polygon face in the source mesh and a cell of the target mesh.
   * It first calculates the transformation that takes the target tetrahedron into the unit tetrahedron. After that, the
   * faces of the source element are triangulated and the calculated transformation is applied
   * to each triangle.
   * The algorithm is based on the algorithm of Grandy used in intersectSourceCell to compute
   * the volume of intersection of two cell elements.
   * The case with a source face colinear to one of the face of tetrahedrons is taking into account:
   * the contribution of the face must not be counted two times.
   *
   * The class will cache the intermediary calculations of transformed nodes of source faces and surfaces associated
   * with triangulated faces to avoid having to recalculate these.
   *
   * @param polyType type of the polygon source face
   * @param polyNodesNbr number of the nodes of the polygon source face
   * @param polyNodes numbers of the nodes of the polygon source face
   * @param polyCoords coordinates of the nodes of the polygon source face
   * @param dimCaracteristic characteristic size of the meshes containing the triangles
   * @param precision precision for double float data used for comparison
   * @param listOfTetraFacesTreated list of tetra faces treated
   * @param listOfTetraFacesColinear list of tetra faces colinear with the polygon source faces
   */
  template<class MyMeshType>
  double SplitterTetra<MyMeshType>::intersectSourceFace(const NormalizedCellType polyType,
                                                        const int polyNodesNbr,
                                                        const int *const polyNodes,
                                                        const double *const *const polyCoords,
                                                        const double dimCaracteristic,
                                                        const double precision,
                                                        std::multiset<TriangleFaceKey>& listOfTetraFacesTreated,
                                                        std::set<TriangleFaceKey>& listOfTetraFacesColinear)
  {
    double totalSurface = 0.0;

    // check if we have planar tetra element
    if(_t->determinant() == 0.0)
      {
        // tetra is planar
        LOG(2, "Planar tetra -- volume 0");
        return 0.0;
      }

    // halfspace filtering
    bool isOutside[8] = {true, true, true, true, true, true, true, true};
    bool isStrictlyOutside[8] = {true, true, true, true, true, true, true, true};
    bool isTargetStrictlyOutside = false;
    bool isTargetOutside = false;

    // calculate the coordinates of the nodes
    for(int i = 0;i<(int)polyNodesNbr;++i)
      {
        const int globalNodeNum = polyNodes[i];
        if(_nodes.find(globalNodeNum) == _nodes.end())
          {
            calculateNode2(globalNodeNum, polyCoords[i]);
          }

        CheckIsStrictlyOutside(_nodes[globalNodeNum], isStrictlyOutside, precision);
        CheckIsOutside(_nodes[globalNodeNum], isOutside, precision);
      }

    // halfspace filtering check
    // NB : might not be beneficial for caching of triangles
    for(int i = 0; i < 8; ++i)
      {
        if(isStrictlyOutside[i])
          {
            isTargetStrictlyOutside = true;
            break;
          }
        else if (isOutside[i])
          {
            isTargetOutside = true;
          }
      }

    if (!isTargetStrictlyOutside)
      {

        if (isTargetOutside)
          {
            // Faces are parallel
            const int tetraFacesNodesConn[4][3] = {
                { 0, 1, 2 },
                { 0, 2, 3 },
                { 0, 3, 1 },
                { 1, 2, 3 } };
            double planeNormal[3];
            for (int iTetraFace = 0; iTetraFace < 4; ++iTetraFace)
              {
                const int * const tetraFaceNodesConn = tetraFacesNodesConn[iTetraFace];
                TriangleFaceKey key = TriangleFaceKey(_conn[tetraFaceNodesConn[0]],
                                                      _conn[tetraFaceNodesConn[1]],
                                                      _conn[tetraFaceNodesConn[2]]);
                if (listOfTetraFacesTreated.find(key) == listOfTetraFacesTreated.end())
                  {
                    const double * const coordsTetraTriNode1 = _coords + tetraFaceNodesConn[0] * MyMeshType::MY_SPACEDIM;
                    const double * const coordsTetraTriNode2 = _coords + tetraFaceNodesConn[1] * MyMeshType::MY_SPACEDIM;
                    const double * const coordsTetraTriNode3 = _coords + tetraFaceNodesConn[2] * MyMeshType::MY_SPACEDIM;
                    calculateNormalForTria(coordsTetraTriNode1, coordsTetraTriNode2, coordsTetraTriNode3, planeNormal);
                    const double normOfTetraTriNormal = norm(planeNormal);
                    if (epsilonEqual(normOfTetraTriNormal, 0.))
                      {
                        for (int i = 0; i < 3; ++i)
                          {
                            planeNormal[i] = 0.;
                          }
                      }
                    else
                      {
                        const double invNormOfTetraTriNormal = 1. / normOfTetraTriNormal;
                        for (int i = 0; i < 3; ++i)
                          {
                            planeNormal[i] *= invNormOfTetraTriNormal;
                          }
                      }
                    double planeConstant = dot(planeNormal, coordsTetraTriNode1);
                    if (IsFacesCoplanar(planeNormal, planeConstant, polyCoords, precision))
                      {
                        int nbrPolyTri = polyNodesNbr - 2; // split polygon into nbrPolyTri triangles
                        for (int iTri = 0; iTri < nbrPolyTri; ++iTri)
                          {
                            double volume = CalculateIntersectionSurfaceOfCoplanarTriangles(planeNormal,
                                                                                            planeConstant,
                                                                                            polyCoords[0],
                                                                                            polyCoords[1 + iTri],
                                                                                            polyCoords[2 + iTri],
                                                                                            coordsTetraTriNode1,
                                                                                            coordsTetraTriNode2,
                                                                                            coordsTetraTriNode3,
                                                                                            dimCaracteristic,
                                                                                            precision);
                            if (!epsilonEqual(volume, 0.))
                              {
                                totalSurface += volume;
                                listOfTetraFacesColinear.insert(key);
                              }
                          }
                      }
                  }
                listOfTetraFacesTreated.insert(key);
              }
          }
        else
          {
              // intersect a son with the unit tetra
              switch (polyType)
                {
                case NORM_TRI3:
                  {
                    // create the face key
                    TriangleFaceKey key = TriangleFaceKey(polyNodes[0], polyNodes[1], polyNodes[2]);

                    // calculate the triangle if needed
                    if (_volumes.find(key) == _volumes.end())
                      {
                        TransformedTriangle tri(_nodes[polyNodes[0]], _nodes[polyNodes[1]], _nodes[polyNodes[2]]);
                        calculateSurface(tri, key);
                        totalSurface += _volumes[key];
                      }
                    else
                      {
                        // count negative as face has reversed orientation
                        totalSurface -= _volumes[key];
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
                    TriangleFaceKey key1 = TriangleFaceKey(polyNodes[0], polyNodes[1], polyNodes[2]);
                    if (_volumes.find(key1) == _volumes.end())
                      {
                        TransformedTriangle tri(_nodes[polyNodes[0]], _nodes[polyNodes[1]], _nodes[polyNodes[2]]);
                        calculateSurface(tri, key1);
                        totalSurface += _volumes[key1];
                      }
                    else
                      {
                        // count negative as face has reversed orientation
                        totalSurface -= _volumes[key1];
                      }

                    // local nodes 1, 3, 4
                    TriangleFaceKey key2 = TriangleFaceKey(polyNodes[0], polyNodes[2], polyNodes[3]);
                    if (_volumes.find(key2) == _volumes.end())
                      {
                        TransformedTriangle tri(_nodes[polyNodes[0]], _nodes[polyNodes[2]], _nodes[polyNodes[3]]);
                        calculateSurface(tri, key2);
                        totalSurface += _volumes[key2];
                      }
                    else
                      {
                        // count negative as face has reversed orientation
                        totalSurface -= _volumes[key2];
                      }
                  }
                  break;

                case NORM_POLYGON:
                  {
                    int nbrPolyTri = polyNodesNbr - 2; // split polygon into nbrPolyTri triangles
                    for (int iTri = 0; iTri < nbrPolyTri; ++iTri)
                      {
                        TriangleFaceKey key = TriangleFaceKey(polyNodes[0], polyNodes[1 + iTri], polyNodes[2 + iTri]);
                        if (_volumes.find(key) == _volumes.end())
                          {
                            TransformedTriangle tri(_nodes[polyNodes[0]], _nodes[polyNodes[1 + iTri]],
                                _nodes[polyNodes[2 + iTri]]);
                            calculateSurface(tri, key);
                            totalSurface += _volumes[key];
                          }
                        else
                          {
                            totalSurface -= _volumes[key];
                          }
                      }
                  }
                  break;

                default:
                  std::cout
                      << "+++ Error : Only elements with triangular and quadratilateral faces are supported at the moment."
                      << std::endl;
                  assert(false);
                }

          }
      }

    // reset if it is very small to keep the matrix sparse
    // is this a good idea?
    if(epsilonEqual(totalSurface, 0.0, SPARSE_TRUNCATION_LIMIT))
      {
        totalSurface = 0.0;
      }
    
    LOG(2, "Volume = " << totalSurface << ", det= " << _t->determinant());

    return totalSurface;
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
      CheckIsOutside(nodes[i], isOutside);
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
      const CellModel& cellModelCell=CellModel::GetCellModel(NORM_TETRA4);
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
    if((int)_nodes.size()>=/*8*/nbOfNodesT)
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
   * \param [in] targetCell in C mode.
   * \param [out] tetra is the output result tetra containers.
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::splitTargetCell2(typename MyMeshTypeT::MyConnType targetCell, typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    const int *refConn(_target_mesh.getConnectivityPtr());
    const int *cellConn(refConn+_target_mesh.getConnectivityIndexPtr()[targetCell]);
    INTERP_KERNEL::NormalizedCellType gt(_target_mesh.getTypeOfElement(targetCell));
    std::vector<int> tetrasNodalConn;
    std::vector<double> addCoords;
    const double *coords(_target_mesh.getCoordinatesPtr());
    SplitIntoTetras(_splitting_pol,gt,cellConn,refConn+_target_mesh.getConnectivityIndexPtr()[targetCell+1],coords,tetrasNodalConn,addCoords);
    std::size_t nbTetras(tetrasNodalConn.size()/4); tetra.resize(nbTetras);
    double tmp[12];
    int tmp2[4];
    for(std::size_t i=0;i<nbTetras;i++)
      {
        for(int j=0;j<4;j++)
          {
            int cellId(tetrasNodalConn[4*i+j]);
            tmp2[j]=cellId;
            if(cellId>=0)
              {
                tmp[j*3+0]=coords[3*cellId+0];
                tmp[j*3+1]=coords[3*cellId+1];
                tmp[j*3+2]=coords[3*cellId+2];
              }
            else
              {
                tmp[j*3+0]=addCoords[3*(-cellId-1)+0];
                tmp[j*3+1]=addCoords[3*(-cellId-1)+1];
                tmp[j*3+2]=addCoords[3*(-cellId-1)+2];
              }
          }
        tetra[i]=new SplitterTetra<MyMeshTypeS>(_src_mesh,tmp,tmp2);
      }
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
          splitPyram5(tetra);
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
    // create tetrahedra
    for(int i = 0; i < 5; ++i)
      {
        const double* nodes[4];
        int conn[4];
        for(int j = 0; j < 4; ++j)
          {
            conn[j] = subZone[ SPLIT_NODES_5[4*i+j] ];
            nodes[j] = getCoordsOfSubNode(conn[j]);
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
    for(int i = 0; i < 6; ++i)
      {
        const double* nodes[4];
        int conn[4];
        for(int j = 0; j < 4; ++j)
          {
            conn[j] = subZone[SPLIT_NODES_6[4*i+j]];
            nodes[j] = getCoordsOfSubNode(conn[j]);
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
    
    // nodes to use for tetrahedron
    const double* nodes[4];
    int conn[4];
    // get the cell center
    conn[0] = 14;
    nodes[0] = getCoordsOfSubNode(conn[0]);

    for(int faceCenterNode = 8; faceCenterNode < 14; ++faceCenterNode)
      {
        // get the face center
        conn[1] = faceCenterNode;
        nodes[1] = getCoordsOfSubNode(conn[1]);
        for(int j = 0; j < 4; ++j)
          {
            const int row = 4*(faceCenterNode - 8) + j;
            conn[2] = TETRA_EDGES_GENERAL_24[2*row];
            conn[3] = TETRA_EDGES_GENERAL_24[2*row + 1];
            nodes[2] = getCoordsOfSubNode(conn[2]);
            nodes[3] = getCoordsOfSubNode(conn[3]);

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
    for(int i = 0; i < 8; ++i)
      {
        sixSplit(GENERAL_48_SUBZONES+8*i,tetra);
      }
  }
  
  /**
   * Splits the NORM_PYRA5 into 2 tetrahedra.
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::splitPyram5(typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra)
  {
    static const int SPLIT_PYPA5[2][4] = 
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
          nodes[j] = getCoordsOfSubNode2(SPLIT_PYPA5[i][j],conn[j]);
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
    const CellModel& cellModelCell=CellModel::GetCellModel(normCellType);
    unsigned nbOfCellNodes=cellModelCell.isDynamic() ? _target_mesh.getNumberOfNodesOfElement(OTT<ConnType,numPol>::indFC(targetCell)) : cellModelCell.getNumberOfNodes();

    // get nb of cell sons (faces)
    const ConnType* rawCellConn = _target_mesh.getConnectivityPtr() + OTT<ConnType,numPol>::conn2C( _target_mesh.getConnectivityIndexPtr()[ targetCell ]);
    const int rawNbCellNodes = _target_mesh.getConnectivityIndexPtr()[ targetCell+1 ] - _target_mesh.getConnectivityIndexPtr()[ targetCell ];
    unsigned nbOfSons = cellModelCell.getNumberOfSons2(rawCellConn, rawNbCellNodes);

    // indices of nodes of a son
    static std::vector<int> allNodeIndices; // == 0,1,2,...,nbOfCellNodes-1
    while ( allNodeIndices.size() < nbOfCellNodes )
      allNodeIndices.push_back( allNodeIndices.size() );
    std::vector<int> classicFaceNodes(4);
    if(cellModelCell.isQuadratic())
      throw INTERP_KERNEL::Exception("SplitterTetra2::splitConvex : quadratic 3D cells are not implemented yet !");
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
          std::vector<int> allIndices(nbOfNodesT);
          for ( int i = 0; i < nbOfNodesT; ++i ) allIndices[i] = i;
          double* barycenter = new double[3];
          calcBarycenter(nbOfNodesT, barycenter, &allIndices[0]);
          _nodes.push_back(barycenter);
        }
      }

  }
}

#endif
