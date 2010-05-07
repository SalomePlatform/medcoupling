//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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
#ifndef __INTERPOLATION_TXX__
#define __INTERPOLATION_TXX__

#include "Interpolation.hxx"
#include "IntegralUniformIntersector.hxx"
#include "IntegralUniformIntersector.txx"

namespace INTERP_KERNEL
{ 
  template<class TrueMainInterpolator>
  template<class MyMeshType, class MatrixType>
  int Interpolation<TrueMainInterpolator>::fromToIntegralUniform(bool fromTo, const MyMeshType& mesh, MatrixType& result, const char *method)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    std::string methodCPP(method);
    int ret=-1;
    if(methodCPP=="P0")
      {
        IntegralUniformIntersectorP0<MyMeshType,MatrixType> intersector(mesh,InterpolationOptions::getMeasureAbsStatus());
        intersector.setFromTo(fromTo);
        std::vector<ConnType> tmp;
        intersector.intersectCells(0,tmp,result);
        ret=intersector.getNumberOfColsOfResMatrix();
      }
    else if(methodCPP=="P1")
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
  void Interpolation<TrueMainInterpolator>::checkAndSplitInterpolationMethod(const char *method, std::string& srcMeth, std::string& trgMeth) throw(INTERP_KERNEL::Exception)
  {
    const int NB_OF_METH_MANAGED=4;
    const char *METH_MANAGED[NB_OF_METH_MANAGED]={"P0P0","P0P1","P1P0","P1P1"};
    std::string methodC(method);
    bool found=false;
    for(int i=0;i<NB_OF_METH_MANAGED && !found;i++)
      found=(methodC==METH_MANAGED[i]);
    if(!found)
      {
        std::string msg("The interpolation method : \'"); msg+=method; msg+="\' not managed !";
        throw INTERP_KERNEL::Exception(msg.c_str());
      }
    srcMeth=methodC.substr(0,2);
    trgMeth=methodC.substr(2);
  }
}

#endif

