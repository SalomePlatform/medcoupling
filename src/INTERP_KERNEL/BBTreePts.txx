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
#ifndef __BBTREEPTS_TXX__
#define __BBTREEPTS_TXX__

#include <vector>
#include <algorithm>

#include <iostream>
#include <limits>
#include <cmath>

template <int dim, class ConnType = int>
class BBTreePts
{

private:
  BBTreePts* _left;
  BBTreePts* _right;
  int _level;
  double _max_left;
  double _min_right;
  const double *_pts;
  typename std::vector<ConnType> _elems;
  bool _terminal;
  ConnType _nbelems;
  double _epsilon;

  static const int MIN_NB_ELEMS=15;
  static const int MAX_LEVEL=20;
public:

  /*!
    Constructor of the bounding box tree
    \param [in] pts pointer to the array containing the points that are to be indexed.
    \param [in] elems array to the indices of the elements contained in the BBTreePts
    \param [in] level level in the BBTreePts recursive structure
    \param [in] nbelems nb of elements in the BBTreePts
    \param [in] epsilon precision to which points are decided to be coincident. Contrary to BBTree, the absolute epsilon is computed. So the internal epsilon is always positive. 

    Parameters \a elems and \a level are used only by BBTreePts itself for creating trees recursively. A typical use is therefore :
    \code
    int nbelems=...
    double* pts= new double[dim*nbelems];
    // filling pts ...
    ...
    BBTreePts<2> tree = new BBTreePts<2>(elems,0,0,nbelems,1e-12);
    \endcode
  */
  BBTreePts(const double *pts, const ConnType *elems, int level, ConnType nbelems, double epsilon=1e-12):
    _left(0),_right(0),_level(level),_pts(pts),_terminal(nbelems < MIN_NB_ELEMS || level> MAX_LEVEL),_nbelems(nbelems),_epsilon(std::abs(epsilon))
  {
    double *nodes=new double[nbelems];
    _elems.resize(nbelems);
    for (ConnType i=0;i<nbelems;i++)
      {
        ConnType elem;
        if (elems!=0)
          elem= elems[i];
        else
          elem=i;

        _elems[i]=elem;
        nodes[i]=pts[elem*dim+(level%dim)];
      }
    if(_terminal) { delete[] nodes; return; }
    //
    std::nth_element<double*>(nodes, nodes+nbelems/2, nodes+nbelems);
    double median=*(nodes+nbelems/2);
    delete [] nodes;
    std::vector<ConnType> new_elems_left,new_elems_right;
 
    new_elems_left.reserve(nbelems/2+1);
    new_elems_right.reserve(nbelems/2+1);
    double max_left = -std::numeric_limits<double>::max();
    double min_right=  std::numeric_limits<double>::max();
    for(int i=0;i<nbelems;i++)
      {
        int elem;
        if(elems!=0)
          elem= elems[i];
        else
          elem=i;
        double mx=pts[elem*dim+(level%dim)];
        if(mx>median)
          {
            new_elems_right.push_back(elem);
            if(mx<min_right) min_right=mx;
          }
        else
          {
            new_elems_left.push_back(elem);
            if(mx>max_left) max_left=mx;
          }
      }
    _max_left=max_left+_epsilon;
    _min_right=min_right-_epsilon;
    ConnType *tmp;
    tmp=0;
    if(!new_elems_left.empty())
      tmp=&(new_elems_left[0]);
    _left=new BBTreePts(pts, tmp, level+1, (int)new_elems_left.size(),_epsilon);
    tmp=0;
    if(!new_elems_right.empty())
      tmp=&(new_elems_right[0]);
    _right=new BBTreePts(pts, tmp, level+1, (int)new_elems_right.size(),_epsilon);
  }


 
  ~BBTreePts()
  {
    delete _left;
    delete _right;
  }

  /*! returns in \a elems the list of elements potentially containing the point pointed to by \a xx
      Contrary to BBTreePts::getElementsAroundPoint the norm 2 is used here.

    \param [in] xx pointer to query point coords
    \param [in] threshold
    \param elems list of elements (given in 0-indexing) intersecting the bounding box
    \sa BBTreePts::getElementsAroundPoint
  */
  double getElementsAroundPoint2(const double *xx, double threshold, ConnType& elem) const
  {
    //  terminal node : return list of elements intersecting bb
    if(_terminal)
      {
        double ret=std::numeric_limits<double>::max();
        for(ConnType i=0;i<_nbelems;i++)
          {
            const double* const bb_ptr=_pts+_elems[i]*dim;
            double tmp=0.;
            for(int idim=0;idim<dim;idim++)
              tmp+=(bb_ptr[idim]-xx[idim])*(bb_ptr[idim]-xx[idim]);
            if(tmp<threshold)
              {
                if(tmp<ret)
                  { ret=tmp; elem=_elems[i]; }
              }
          }
        return ret;
      }
    //non terminal node
    double s=sqrt(threshold*dim);
    if(xx[_level%dim]+s<_min_right)
      return _left->getElementsAroundPoint2(xx,threshold,elem);
    if(xx[_level%dim]-s>_max_left)
      return _right->getElementsAroundPoint2(xx,threshold,elem);
    int eleml,elemr;
    double retl=_left->getElementsAroundPoint2(xx,threshold,eleml);
    double retr=_right->getElementsAroundPoint2(xx,threshold,elemr);
    if(retl<retr)
      { elem=eleml; return retl; }
    else
      { elem=elemr; return retr; }
  }
 
  /*! returns in \a elems the list of elements potentially containing the point pointed to by \a xx
   * ** Infinite norm is used here not norm 2 ! ***
   * 
   *  \param xx pointer to query point coords
   *  \param elems list of elements (given in 0-indexing) intersecting the bounding box
   * \sa BBTreePts::getElementsAroundPoint2
   */
  void getElementsAroundPoint(const double* xx, std::vector<ConnType>& elems) const
  {
    //  terminal node : return list of elements intersecting bb
    if(_terminal)
      {
        for(ConnType i=0;i<_nbelems;i++)
          {
            const double* const bb_ptr=_pts+_elems[i]*dim;
            bool intersects = true;
            for(int idim=0;idim<dim;idim++)
              intersects=intersects && (std::abs(bb_ptr[idim]-xx[idim])<=_epsilon);
            if(intersects)
              elems.push_back(_elems[i]);
          }
        return;
      }
    //non terminal node 
    if(xx[_level%dim]<_min_right)
      {
        _left->getElementsAroundPoint(xx,elems);
        return;
      }
    if(xx[_level%dim]>_max_left)
      {
        _right->getElementsAroundPoint(xx,elems);
        return;
      }
    _left->getElementsAroundPoint(xx,elems);
    _right->getElementsAroundPoint(xx,elems);
  }
  
  int size() const
  {
    if(_terminal)
      return _nbelems;
    return _left->size()+_right->size();
  }
};

#endif
