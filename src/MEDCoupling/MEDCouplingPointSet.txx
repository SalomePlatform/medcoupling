// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
#ifndef __PARAMEDMEM_MEDCOUPLINGPOINTSET_TXX__
#define __PARAMEDMEM_MEDCOUPLINGPOINTSET_TXX__

#include "MEDCouplingPointSet.txx"
#include "InterpolationUtils.hxx"
#include "BBTree.txx"

#include <vector>

namespace ParaMEDMEM
{
  template<int SPACEDIM>
  void MEDCouplingPointSet::findNodeIdsNearPointAlg(std::vector<double>& bbox, const double *pos, int nbNodes, double eps,
                                                    std::vector<int>& c, std::vector<int>& cI) const
  {
    const double *coordsPtr=_coords->getConstPointer();
    BBTree<SPACEDIM,int> myTree(&bbox[0],0,0,getNumberOfNodes(),-eps);
    double bb[2*SPACEDIM];
    double eps2=eps*eps;
    for(int i=0;i<nbNodes;i++)
      {
        for(int j=0;j<SPACEDIM;j++)
          {
            bb[2*j]=pos[SPACEDIM*i+j];
            bb[2*j+1]=pos[SPACEDIM*i+j];
          }
        std::vector<int> intersectingElems;
        myTree.getIntersectingElems(bb,intersectingElems);
        std::vector<int> commonNodes;
        for(std::vector<int>::const_iterator it=intersectingElems.begin();it!=intersectingElems.end();it++)
          if(INTERP_KERNEL::distance2<SPACEDIM>(pos+SPACEDIM*i,coordsPtr+SPACEDIM*(*it))<eps2)
            commonNodes.push_back(*it);
        cI.push_back(cI.back()+(int)commonNodes.size());
        c.insert(c.end(),commonNodes.begin(),commonNodes.end());
      }
  }
}

#endif
