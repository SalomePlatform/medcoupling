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
// Author : Anthony Geay (CEA/DEN)
#ifndef __INTEGRALUNIFORMINTERSECTOR_TXX__
#define __INTEGRALUNIFORMINTERSECTOR_TXX__

#include "IntegralUniformIntersector.hxx"
#include "VolSurfUser.txx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  IntegralUniformIntersector<MyMeshType,MyMatrix>::IntegralUniformIntersector(const MyMeshType& mesh, bool isAbs):_mesh(mesh),_from_to(false),_is_abs(isAbs)
  {
  }

  template<class MyMeshType, class MyMatrix>
  void IntegralUniformIntersector<MyMeshType,MyMatrix>::putValueIn(ConnType iInCMode, double val1, MyMatrix& res) const
  {
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
    double val=performNormalization(val1);
    if(_from_to)
      {
        typename MyMatrix::value_type& resRow=res[0];
        typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(iInCMode));
        if(iterRes==resRow.end())
          resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(iInCMode),val));
        else
          {
            double val2=(*iterRes).second+val;
            resRow.erase(OTT<ConnType,numPol>::indFC(iInCMode));
            resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(iInCMode),val2));
          }
      }
    else
      {
        typename MyMatrix::value_type& resRow=res[iInCMode];
        typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(0));
        if(iterRes==resRow.end())
          resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(0),val));
        else
          {
            double val2=(*iterRes).second+val;
            resRow.erase(OTT<ConnType,numPol>::indFC(0));
            resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(0),val2));
          }
      }
  }

  template<class MyMeshType, class MyMatrix>
  IntegralUniformIntersectorP0<MyMeshType,MyMatrix>::IntegralUniformIntersectorP0(const MyMeshType& mesh, bool isAbs):IntegralUniformIntersector<MyMeshType,MyMatrix>(mesh,isAbs)
  {
  }

  template<class MyMeshType, class MyMatrix>
  int IntegralUniformIntersectorP0<MyMeshType,MyMatrix>::getNumberOfRowsOfResMatrix() const
  {
    if(IntegralUniformIntersector<MyMeshType,MyMatrix>::_from_to)
      return 1;
    else
      return IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getNumberOfElements();
  }
  
  template<class MyMeshType, class MyMatrix>
  int IntegralUniformIntersectorP0<MyMeshType,MyMatrix>::getNumberOfColsOfResMatrix() const
  {
    if(IntegralUniformIntersector<MyMeshType,MyMatrix>::_from_to)
      return IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getNumberOfElements();
    else
      return 1;
  }
  
  template<class MyMeshType, class MyMatrix>
  void IntegralUniformIntersectorP0<MyMeshType,MyMatrix>::intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res)
  {
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
    res.resize(getNumberOfRowsOfResMatrix());
    unsigned long nbelem=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getNumberOfElements();
    const ConnType *connIndx=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getConnectivityIndexPtr();
    const ConnType *conn=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getConnectivityPtr();
    const double *coords=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getCoordinatesPtr();
    for(unsigned long i=0;i<nbelem;i++)
      {
        INTERP_KERNEL::NormalizedCellType t=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(i));
        double val=computeVolSurfOfCell<ConnType,numPol,MyMeshType::MY_SPACEDIM>(t,conn+OTT<ConnType,numPol>::ind2C(connIndx[i]),connIndx[i+1]-connIndx[i],coords);
        IntegralUniformIntersector<MyMeshType,MyMatrix>::putValueIn(i,val,res);
      }
  }

  template<class MyMeshType, class MyMatrix>
  IntegralUniformIntersectorP1<MyMeshType,MyMatrix>::IntegralUniformIntersectorP1(const MyMeshType& mesh, bool isAbs):IntegralUniformIntersector<MyMeshType,MyMatrix>(mesh,isAbs)
  {
  }

  template<class MyMeshType, class MyMatrix>
  int IntegralUniformIntersectorP1<MyMeshType,MyMatrix>::getNumberOfRowsOfResMatrix() const
  {
    if(IntegralUniformIntersector<MyMeshType,MyMatrix>::_from_to)
      return 1;
    else
      return IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getNumberOfNodes();
  }
  
  template<class MyMeshType, class MyMatrix>
  int IntegralUniformIntersectorP1<MyMeshType,MyMatrix>::getNumberOfColsOfResMatrix() const
  {
    if(IntegralUniformIntersector<MyMeshType,MyMatrix>::_from_to)
      return IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getNumberOfNodes();
    else
      return 1;
  }
  
  template<class MyMeshType, class MyMatrix>
  void IntegralUniformIntersectorP1<MyMeshType,MyMatrix>::intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res)
  {
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
    res.resize(getNumberOfRowsOfResMatrix());
    unsigned long nbelem=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getNumberOfElements();
    const ConnType *connIndx=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getConnectivityIndexPtr();
    const ConnType *conn=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getConnectivityPtr();
    const double *coords=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getCoordinatesPtr();
    for(unsigned long i=0;i<nbelem;i++)
      {
        INTERP_KERNEL::NormalizedCellType t=IntegralUniformIntersector<MyMeshType,MyMatrix>::_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(i));
        int lgth=connIndx[i+1]-connIndx[i];
        const ConnType *locConn=conn+OTT<ConnType,numPol>::ind2C(connIndx[i]);
        double val=computeVolSurfOfCell<ConnType,numPol,MyMeshType::MY_SPACEDIM>(t,locConn,lgth,coords);
        if(t==NORM_TRI3)
          val/=3.;
        else if(t==NORM_TETRA4)
          val/=4.;
        else
          throw INTERP_KERNEL::Exception("Invalid cell type detected : must be TRI3 or TETRA4 ! ");
        for(int j=0;j<lgth;j++)
          IntegralUniformIntersector<MyMeshType,MyMatrix>::putValueIn(OTT<ConnType,numPol>::coo2C(locConn[j]),val,res);
      }
  }
}

#endif
