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

#ifndef __POLYHEDRON3D2DINTERSECTORP0P0_HXX__
#define __POLYHEDRON3D2DINTERSECTORP0P0_HXX__

#include "Intersector3DP0P0.hxx"
#include "SplitterTetra.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{


  /** 
   * \brief Class responsible for calculating intersection between a hexahedron target element and  
   * the source elements.
   *
   */
  template<class MyMeshType, class MyMatrixType>
  class Polyhedron3D2DIntersectorP0P0 : public Intersector3DP0P0<MyMeshType,MyMatrixType>
  {
    typedef typename std::map<int,std::set<int> > DuplicateFacesType;

  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;

  public:

    Polyhedron3D2DIntersectorP0P0(const MyMeshType& targetMesh,
                                  const MyMeshType& srcMesh,
                                  const double dimCaracteristic,
                                  const double precision,
                                  DuplicateFacesType& intersectFaces,
                                  SplittingPolicy policy = PLANAR_FACE_5);

    ~Polyhedron3D2DIntersectorP0P0();

    void intersectCells(ConnType targetCell,
                        const std::vector<ConnType>& srcCells,
                        MyMatrixType& matrix);

  private:
    void releaseArrays();
  private:
    /// pointers to the SplitterTetra objects representing the tetrahedra 
    /// that result from the splitting of the hexahedron target cell
    std::vector< SplitterTetra<MyMeshType>* > _tetra;
    
    SplitterTetra2<MyMeshType> _split;

    double _dim_caracteristic;
    double _precision;

    DuplicateFacesType& _intersect_faces;

  };
}

#endif
