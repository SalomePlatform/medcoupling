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

#ifndef __BBTREEDST_TXX__
#define __BBTREEDST_TXX__

#include <vector>
#include <algorithm>

#include <iostream>
#include <limits>
#include <cmath>

template <int dim>
class BBTreeDst
{
private:
  BBTreeDst* _left;
  BBTreeDst* _right;
  int _level;
  double _max_left;
  double _min_right;
  const double *_bb;
  std::vector<int> _elems;
  double  *_terminal;
  int _nbelems;

  static const int MIN_NB_ELEMS=15;
  static const int MAX_LEVEL=20;
public:
  BBTreeDst(const double* bbs, int* elems, int level, int nbelems):
    _left(0),_right(0),_level(level),_bb(bbs),_terminal(0),_nbelems(nbelems)
  {
    if((nbelems < MIN_NB_ELEMS || level> MAX_LEVEL))
      _terminal=new double[2*dim];
    _elems.resize(nbelems);
    for (int i=0; i<nbelems; i++)
      _elems[i]=elems?elems[i]:i;
    if(_terminal)
      {
        fillBBoxTerminal(bbs);
        return ;
      }
    double *nodes=new double[nbelems];
    for (int i=0; i<nbelems; i++)
      nodes[i]=bbs[_elems[i]*dim*2+(level%dim)*2];
    std::nth_element<double*>(nodes, nodes+nbelems/2, nodes+nbelems);
    double median = *(nodes+nbelems/2);
    delete [] nodes;
    std::vector<int> new_elems_left;
    std::vector<int> new_elems_right;
 
    new_elems_left.reserve(nbelems/2+1);
    new_elems_right.reserve(nbelems/2+1);
    double max_left = -std::numeric_limits<double>::max();
    double min_right=  std::numeric_limits<double>::max();
    for(int i=0; i<nbelems;i++)
      {
        int elem;
        if (elems!=0)
          elem= elems[i];
        else
          elem=i;
        double max=bbs[elem*dim*2+(level%dim)*2+1];
        double min = bbs[elem*dim*2+(level%dim)*2];
      
        if (min>median)
          {
            new_elems_right.push_back(elem);
            if (min<min_right) min_right = min;
          }
        else
          {
            new_elems_left.push_back(elem);
            if (max>max_left) max_left = max;
          }
      }
    _max_left=max_left;
    _min_right=min_right;
    int *tmp;
    tmp=0;
    if(!new_elems_left.empty())
      tmp=&(new_elems_left[0]);
    _left=new BBTreeDst(bbs, tmp, level+1, (int)new_elems_left.size());
    tmp=0;
    if(!new_elems_right.empty())
      tmp=&(new_elems_right[0]);
    _right=new BBTreeDst(bbs, tmp, level+1, (int)new_elems_right.size());
  }

  ~BBTreeDst()
  {
    delete _left;
    delete _right;
    delete [] _terminal;
  }

  void getElemsWhoseMinDistanceToPtSmallerThan(const double *pt, double minOfMaxDstsSq, std::vector<int>& elems) const
  {
    if(_terminal)
      {
        for(int i=0; i<_nbelems; i++)
          {
            if(GetMinDistanceFromBBoxToPt(_bb+_elems[i]*2*dim,pt)<minOfMaxDstsSq)
              elems.push_back(_elems[i]);
          }
      }
    else
      {
        double minOfMaxDsts=sqrt(minOfMaxDstsSq);
        if(_min_right-pt[_level%dim]>minOfMaxDsts)
          { _left->getElemsWhoseMinDistanceToPtSmallerThan(pt,minOfMaxDstsSq,elems); return ; }
        if(pt[_level%dim]-_max_left>minOfMaxDsts)
          { _right->getElemsWhoseMinDistanceToPtSmallerThan(pt,minOfMaxDstsSq,elems); return ; }
        _left->getElemsWhoseMinDistanceToPtSmallerThan(pt,minOfMaxDstsSq,elems);
        _right->getElemsWhoseMinDistanceToPtSmallerThan(pt,minOfMaxDstsSq,elems);
      }
  }
  
  /** Get the minimal (square) distance between a point and all the available bounding boxes in the tree.
    The (square) distance to a bbox is the true geometric distance between the point and a face
    (or an edge, or a corner) of the bbox.
  */
  void getMinDistanceOfMax(const double *pt, double& minOfMaxDstsSq) const
  {
    if(_terminal)
      {
        if(GetMinDistanceFromBBoxToPt(_terminal,pt)>minOfMaxDstsSq)//min it is not a bug
          return ;
        for(int i=0; i<_nbelems; i++)
          {
            minOfMaxDstsSq=std::min(minOfMaxDstsSq,GetMaxDistanceFromBBoxToPt(_bb+_elems[i]*2*dim,pt));
          }
      }
    else
      {
        double minOfMaxDsts=sqrt(minOfMaxDstsSq);
        if(_min_right-pt[_level%dim]>minOfMaxDsts)
          { _left->getMinDistanceOfMax(pt,minOfMaxDstsSq); return ; }
        if(pt[_level%dim]-_max_left>minOfMaxDsts)
          { _right->getMinDistanceOfMax(pt,minOfMaxDstsSq); return ; }
        _left->getMinDistanceOfMax(pt,minOfMaxDstsSq);
        _right->getMinDistanceOfMax(pt,minOfMaxDstsSq);
      }
  }

  void fillBBoxTerminal(const double* bbs)
  {
    for(int j=0;j<dim;j++)
      {
        _terminal[2*j]=std::numeric_limits<double>::max();
        _terminal[2*j+1]=-std::numeric_limits<double>::max();
      }
    for(int i=0;i<_nbelems;i++)
      {
        for(int j=0;j<dim;j++)
          {
            _terminal[2*j]=std::min(_terminal[2*j],bbs[2*dim*_elems[i]+2*j]);
            _terminal[2*j+1]=std::max(_terminal[2*j+1],bbs[2*dim*_elems[i]+2*j+1]);
          }
      }
  }

  static double GetMaxDistanceFromBBoxToPt(const double *bbox, const double *pt)
  {
    if(bbox[0]<=bbox[1])
      {
        double zeRes=0.;
        for (int idim=0; idim<dim; idim++)
          {
            double val1=pt[idim]-bbox[idim*2],val2=pt[idim]-bbox[idim*2+1];
            double x=std::max(fabs(val1),fabs(val2));
            zeRes+=x*x;
          }
        return zeRes;
      }
    else//min>max -> no cells in this
      return std::numeric_limits<double>::max();
    
  }
  
  static double GetMinDistanceFromBBoxToPt(const double *bbox, const double *pt)
  {
    if(bbox[0]<=bbox[1])
      {
        double zeRes=0.;
        for (int idim=0; idim<dim; idim++)
          {
            double val1=pt[idim]-bbox[idim*2],val2=pt[idim]-bbox[idim*2+1];
            char pos=(( (0.<val1)-(val1<0.) )+( (0.<val2)-(val2<0.) ))/2;// sign(val) = (0.<val)-(val<0.)
            if(pos!=0)
              {
                double x=pos==1?val2:val1;
                zeRes+=x*x;
              }
          }
        return zeRes;
      }
    else//min>max -> no cells in this
      return std::numeric_limits<double>::max();
  }
};

#endif
