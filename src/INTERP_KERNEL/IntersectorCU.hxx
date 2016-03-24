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
// File      : IntersectorCU.hxx
// Created   : Thu Dec 17 12:30:17 2009
// Author    : Edward AGAPOV (eap)
//

#ifndef __IntersectorCU_HXX__
#define __IntersectorCU_HXX__

#include "TargetIntersector.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{
  template<class MyCMeshType, class MyUMeshType, class MyMatrix> class _StabIntersector;

  template<class MyCMeshType, class MyUMeshType, class MyMatrix, class ConcreteIntersector=_StabIntersector<MyCMeshType,MyUMeshType,MyMatrix> >
  class IntersectorCU : public TargetIntersector<MyCMeshType,MyMatrix>
  {
  public:
    static const int SPACEDIM=MyCMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyCMeshType::MY_MESHDIM;
    typedef typename MyUMeshType::MyConnType UConnType;
    typedef typename MyCMeshType::MyConnType CConnType;
  public:
    //! \addtogroup InterpKerGrpIntCU @{
    IntersectorCU(const MyCMeshType& meshS, const MyUMeshType& meshT);
    //! @}
    virtual ~IntersectorCU();
    void getUElemBB(double* bb, UConnType iP);
    void getUCoordinates(UConnType icell, std::vector<double>& coords);

    int getNumberOfRowsOfResMatrix() const;
    int getNumberOfColsOfResMatrix() const;
    void intersectCells(CConnType icellU, const std::vector<CConnType>& icellC, MyMatrix& res);
    double intersectGeometry(CConnType icellT, const std::vector<CConnType>& icellC) { return asLeaf().intersectGeometry(icellT,icellC); }
  protected:
    ConcreteIntersector& asLeaf() { return static_cast<ConcreteIntersector&>(*this); }

  protected:
    const UConnType *_connectU;
    const UConnType *_connIndexU;
    const double *  _coordsU;
    const MyUMeshType& _meshU;

    const double *     _coordsC[SPACEDIM];
    int                _nbCellsC[SPACEDIM];
    const MyCMeshType& _meshC;
  };

  // class to enable usage of IntersectorCU not for intersection but for access to data it encapsulates
  template<class MyCMeshType, class MyUMeshType, class MyMatrix>
  class _StabIntersector: public IntersectorCU<MyCMeshType, MyUMeshType, MyMatrix, _StabIntersector<MyCMeshType, MyUMeshType, MyMatrix> >
  {
  public:
    _StabIntersector(const MyCMeshType& meshS, const MyUMeshType& meshT) : IntersectorCU<MyCMeshType, MyUMeshType, MyMatrix, _StabIntersector<MyCMeshType, MyUMeshType, MyMatrix> >(meshS, meshT) {}
    double intersectGeometry(typename MyUMeshType::MyConnType icellT, const std::vector<typename MyCMeshType::MyConnType>& icellC) { throw Exception("You must provide an intersector as the 4-th template argument of IntersectorCU"); return 0; }
  };
}

#endif
