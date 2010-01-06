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
#ifndef __PARAMEDMEM_MEDCOUPLINGPOINTSET_TXX__
#define __PARAMEDMEM_MEDCOUPLINGPOINTSET_TXX__

#include "MEDCouplingPointSet.txx"
#include "InterpolationUtils.hxx"
#include "BBTree.txx"

#include <vector>

namespace ParaMEDMEM
{
  template<int SPACEDIM>
  void MEDCouplingPointSet::findCommonNodesAlg(std::vector<double>& bbox,
                                               int nbNodes, double prec,
                                               std::vector<int>& c, std::vector<int>& cI) const
  {
    const double *coordsPtr=_coords->getConstPointer();
    BBTree<SPACEDIM,int> myTree(&bbox[0],0,0,nbNodes,prec);
    double bb[2*SPACEDIM];
    double prec2=prec*prec;
    for(int i=0;i<nbNodes;i++)
      {
        if(std::find(c.begin(),c.end(),i)!=c.end())
          continue;
        for(int j=0;j<SPACEDIM;j++)
          {
            bb[2*j]=coordsPtr[SPACEDIM*i+j];
            bb[2*j+1]=coordsPtr[SPACEDIM*i+j];
          }
        std::vector<int> intersectingElems;
        myTree.getIntersectingElems(bb,intersectingElems);
        if(intersectingElems.size()>1)
          {
            std::vector<int> commonNodes;
            for(std::vector<int>::const_iterator it=intersectingElems.begin();it!=intersectingElems.end();it++)
              if(*it!=i)
                if(INTERP_KERNEL::distance2<SPACEDIM>(coordsPtr+SPACEDIM*i,coordsPtr+SPACEDIM*(*it))<prec2)
                  commonNodes.push_back(*it);
            if(!commonNodes.empty())
              {
                cI.push_back(cI.back()+commonNodes.size()+1);
                c.push_back(i);
                c.insert(c.end(),commonNodes.begin(),commonNodes.end());
              }
          }
      }
  }
}

#endif
