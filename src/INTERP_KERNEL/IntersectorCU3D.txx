// Copyright (C) 2009-2016  OPEN CASCADE
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
// File      : IntersectorCU3D.txx
// Created   : Thu Dec 17 14:17:49 2009
// Author    : Edward AGAPOV (eap)

#ifndef __IntersectorCU3D_TXX__
#define __IntersectorCU3D_TXX__

#include "IntersectorCU3D.hxx"
#include "IntersectorCU.txx"
#include "SplitterTetra.txx"

#define  IntersectorCU3D_TEMPLATE template<class MyCMeshType, class MyUMeshType, class MyMatrix>
#define  INTERSECTOR_CU3D IntersectorCU3D<MyCMeshType,MyUMeshType,MyMatrix >
#define _INTERSECTOR_CU   IntersectorCU  <MyCMeshType,MyUMeshType,MyMatrix,IntersectorCU3D<MyCMeshType,MyUMeshType,MyMatrix> >

namespace INTERP_KERNEL
{
  //================================================================================
  /*!
   * \brief Unstructured hexahedral mesh derived from cartesian 3D mesh.
   *        Mesh contains one HEXA8 element
   */
  //================================================================================

  class _Cartesian3D2UnstructHexMesh
  {
  public:
    static const int MY_SPACEDIM=3;
    static const int MY_MESHDIM=3;
    typedef int MyConnType;
    static const NumberingPolicy My_numPol=ALL_C_MODE;

    _Cartesian3D2UnstructHexMesh(const double * coords[3]): _coordsC(coords) {}
    void setHexa(int I, int J, int K) // indices in C mode
    {
      double* pCoord = _coordsU;
      for ( int k = K; k < K+2; ++k )
        for ( int j = J; j < J+2; ++j )
          for ( int i = I; i < I+2; ++i )
            {
              *pCoord++ = _coordsC[0][i];
              *pCoord++ = _coordsC[1][j];
              *pCoord++ = _coordsC[2][k];
            }
    }
    const int *getConnectivityPtr() const
    {
      static int conn[] = { 1,0,2,3,5,4,6,7 };
      return conn;
    }
    const int *getConnectivityIndexPtr() const
    {
      static int conInd[] = { 0,8 };
      return conInd;
    }
    void getBoundingBox(double *boundingBox) const
    {
      boundingBox[BoundingBox::XMIN] = _coordsU[0];
      boundingBox[BoundingBox::XMAX] = _coordsU[0+1*MY_SPACEDIM];
      boundingBox[BoundingBox::YMIN] = _coordsU[1];
      boundingBox[BoundingBox::YMAX] = _coordsU[1+2*MY_SPACEDIM];
      boundingBox[BoundingBox::ZMIN] = _coordsU[2];
      boundingBox[BoundingBox::ZMAX] = _coordsU[2+4*MY_SPACEDIM];
    }
    NormalizedCellType getTypeOfElement(int eltId) const { return NORM_HEXA8; }
    unsigned char getNumberOfNodesOfElement(int eltId) const { return 8; }
    unsigned long getNumberOfElements() const { return 1; }
    unsigned long getNumberOfNodes()    const { return 8; }
    const double *getCoordinatesPtr()   const { return _coordsU; }
    void releaseTempArrays() {}
  private:
    const double** _coordsC;
    double         _coordsU[3*8];
  };

  //================================================================================
  /*!
   * \brief intersector of the unstructured mesh and the cartesian mesh in 3D
   */
  //================================================================================

  IntersectorCU3D_TEMPLATE
  INTERSECTOR_CU3D::IntersectorCU3D(const MyCMeshType& meshS,
                                    const MyUMeshType& meshT,
                                    SplittingPolicy    splitting_policy):
    _INTERSECTOR_CU( meshS, meshT )
  {
    if ( MyCMeshType::MY_SPACEDIM != 3 || MyCMeshType::MY_MESHDIM != 3 ||
         MyUMeshType::MY_SPACEDIM != 3 || MyUMeshType::MY_MESHDIM != 3 )
      throw Exception("IntersectorCU3D(): Invalid mesh dimension, it must be 3");

    _uHexMesh = new _Cartesian3D2UnstructHexMesh( _INTERSECTOR_CU::_coordsC );
    _split    = new TSplitter( meshT, *_uHexMesh, splitting_policy );
  }

  //================================================================================
  /*!
   * \brief destructor
   */
  //================================================================================

  IntersectorCU3D_TEMPLATE
  INTERSECTOR_CU3D::~IntersectorCU3D()
  {
    delete _uHexMesh; _uHexMesh=0;
    delete _split; _split=0;
  }

  //================================================================================
  /*!
   * \brief Calculate volume of intersection of an unstructured cell and a cartesian one.
   * The cartesian cell is given by its [i,j,k] indices
   */
  //================================================================================

  IntersectorCU3D_TEMPLATE
  double INTERSECTOR_CU3D::intersectGeometry(UConnType                     icellT,
                                             const std::vector<CConnType>& icellS)
  {
    // split an unstructured cell into tetra
    std::vector< TTetra* > tetra;
    UConnType nb_nodes =
      _INTERSECTOR_CU::_connIndexU[icellT+1] - _INTERSECTOR_CU::_connIndexU[icellT];
    _split->releaseArrays();
    _split->splitTargetCell( icellT, nb_nodes, tetra);

    // intersect a cartesian 3d cell with tetra
    _uHexMesh->setHexa( _FMIC(icellS[0]),_FMIC(icellS[1]),_FMIC(icellS[2])); // set cell at i,j,k
    double res = 0;
    for ( unsigned int t = 0; t < tetra.size(); ++t )
      {
        res += tetra[t]->intersectSourceCell( 0 );
        delete tetra[t];
      }
    return res;
  }
}
#endif
