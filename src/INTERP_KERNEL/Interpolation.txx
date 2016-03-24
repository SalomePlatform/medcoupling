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
#ifndef __INTERPOLATION_TXX__
#define __INTERPOLATION_TXX__

#include "Interpolation.hxx"
#include "IntegralUniformIntersector.hxx"
#include "IntegralUniformIntersector.txx"
#include "VectorUtils.hxx"

namespace INTERP_KERNEL
{ 
  template<class TrueMainInterpolator>
  template<class MyMeshType, class MatrixType>
  int Interpolation<TrueMainInterpolator>::fromToIntegralUniform(bool fromTo, const MyMeshType& mesh, MatrixType& result, const std::string& method)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    int ret=-1;
    if(method=="P0")
      {
        IntegralUniformIntersectorP0<MyMeshType,MatrixType> intersector(mesh,InterpolationOptions::getMeasureAbsStatus());
        intersector.setFromTo(fromTo);
        std::vector<ConnType> tmp;
        intersector.intersectCells(0,tmp,result);
        ret=intersector.getNumberOfColsOfResMatrix();
      }
    else if(method=="P1")
      {
        IntegralUniformIntersectorP1<MyMeshType,MatrixType> intersector(mesh,InterpolationOptions::getMeasureAbsStatus());
        intersector.setFromTo(fromTo);
        std::vector<ConnType> tmp;
        intersector.intersectCells(0,tmp,result);
        ret=intersector.getNumberOfColsOfResMatrix();
      }
    else
      throw INTERP_KERNEL::Exception("Invalid method specified in fromIntegralUniform : must be in { \"P0\", \"P1\"}");
    return ret;
  }

  template<class TrueMainInterpolator>
  template<class MyMeshType>
  double Interpolation<TrueMainInterpolator>::CalculateCharacteristicSizeOfMeshes(const MyMeshType& myMeshS, const MyMeshType& myMeshT, const int printLevel)
  {
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;

    long nbMailleS=myMeshS.getNumberOfElements();
    long nbMailleT=myMeshT.getNumberOfElements();

    /**************************************************/
    /* Search the characteristic size of the meshes   */
    /**************************************************/

    double BoxS[2*SPACEDIM]; myMeshS.getBoundingBox(BoxS);
    double BoxT[2*SPACEDIM]; myMeshT.getBoundingBox(BoxT);
    double diagonalS,dimCaracteristicS=std::numeric_limits<double>::max();
    if(nbMailleS!=0)
      {
        diagonalS=getDistanceBtw2Pts<SPACEDIM>(BoxS+SPACEDIM,BoxS);
        dimCaracteristicS=diagonalS/nbMailleS;
      }
    double diagonalT,dimCaracteristicT=std::numeric_limits<double>::max();
    if(nbMailleT!=0)
      {
        diagonalT=getDistanceBtw2Pts<SPACEDIM>(BoxT+SPACEDIM,BoxT);
        dimCaracteristicT=diagonalT/nbMailleT;
      }
    if (printLevel>=1)
      {
        std::cout << "  - Characteristic size of the source mesh : " << dimCaracteristicS << std::endl;
        std::cout << "  - Characteristic size of the target mesh: " << dimCaracteristicT << std::endl;
      }

    return std::min(dimCaracteristicS, dimCaracteristicT);

  }

}

#endif

