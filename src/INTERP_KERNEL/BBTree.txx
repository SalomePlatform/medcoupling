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
#ifndef __BBTREE_TXX__
#define __BBTREE_TXX__

#include <vector>
#include <algorithm>

#include <iostream>
#include <limits>
#include <cmath>

template <int dim, class ConnType = int>
class BBTree
{

private:
  BBTree* _left;
  BBTree* _right;
  int _level;
  double _max_left;
  double _min_right;
  const double *_bb;
  typename std::vector<ConnType> _elems;
  bool _terminal;
  ConnType _nbelems;
  double _epsilon;

  static const int MIN_NB_ELEMS=15;
  static const int MAX_LEVEL=20;
public:

  /*!
    Constructor of the bounding box tree
    \param bbs pointer to the [xmin1 xmax1 ymin1 ymax1 xmin2 xmax2 ...] array containing the bounding boxes that are to be indexed.
    \param elems array to the indices of the elements contained in the BBTree
    \param level level in the BBTree recursive structure
    \param nbelems nb of elements in the BBTree
    \param epsilon precision to which points are decided to be coincident. Epsilon can be positive or negative.
           If \a epsilon is positive the request method will enlarge the computed bounding box (more matching elems return).
           If negative the given bounding box will be tighten (less matching elems return).

    Parameters \a elems and \a level are used only by BBTree itself for creating trees recursively. A typical use is therefore :
    \code
    int nbelems=...
    double* bbs= new double[2*2*nbelems];
    // filling bbs ...
    ...
    BBTree<2> tree = new BBTree<2>(elems,0,0,nbelems,1e-12);
    \endcode
  */

  BBTree(const double* bbs, ConnType* elems, int level, ConnType nbelems, double epsilon=1e-12):
    _left(0), _right(0), _level(level), _bb(bbs), _terminal(false),_nbelems(nbelems),_epsilon(epsilon)
  {
    if (nbelems < MIN_NB_ELEMS || level> MAX_LEVEL)
      {
        _terminal=true;
      
      }
    double* nodes=new double [nbelems];
    _elems.resize(nbelems);
    for (ConnType i=0; i<nbelems; i++)
      {
        ConnType elem;
        if (elems!=0)
          elem= elems[i];
        else
          elem=i;

        _elems[i]=elem;
        nodes[i]=bbs[elem*dim*2+(level%dim)*2];
      }
    if (_terminal) { delete[] nodes; return;}

    std::nth_element<double*>(nodes, nodes+nbelems/2, nodes+nbelems);
    double median = *(nodes+nbelems/2);
    delete[] nodes;
    // std:: cout << *median <<std::endl;

    std::vector<ConnType> new_elems_left;
    std::vector<ConnType> new_elems_right;
 
    new_elems_left.reserve(nbelems/2+1);
    new_elems_right.reserve(nbelems/2+1);
    double max_left = -std::numeric_limits<double>::max();
    double min_right=  std::numeric_limits<double>::max();
    for (int i=0; i<nbelems;i++)
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
    _max_left=max_left+std::abs(_epsilon);
    _min_right=min_right-std::abs(_epsilon);
    ConnType *tmp;
    tmp=0;
    if(!new_elems_left.empty())
      tmp=&(new_elems_left[0]);
    _left=new BBTree(bbs, tmp, level+1, (int)new_elems_left.size(),_epsilon);
    tmp=0;
    if(!new_elems_right.empty())
      tmp=&(new_elems_right[0]);
    _right=new BBTree(bbs, tmp, level+1, (int)new_elems_right.size(),_epsilon);
  
  }


 
  ~BBTree()
  {
    if (_left!=0)  delete _left;
    if (_right!=0) delete _right;

  }

  
  /*! returns in \a elems the list of elements potentially intersecting the bounding box pointed to by \a bb
    
    \param bb pointer to query bounding box
    \param elems list of elements (given in 0-indexing that is to say in \b C \b mode) intersecting the bounding box
  */

  void getIntersectingElems(const double* bb, std::vector<ConnType>& elems) const
  {
    //  terminal node : return list of elements intersecting bb
    if (_terminal)
      {
        for (int i=0; i<_nbelems; i++)
          {
            const double* const  bb_ptr=_bb+_elems[i]*2*dim;
            bool intersects = true;
            for (int idim=0; idim<dim; idim++)
              {
                if (bb_ptr[idim*2]-bb[idim*2+1]>-_epsilon|| bb_ptr[idim*2+1]-bb[idim*2]<_epsilon)
                  intersects=false;
              }
            if (intersects)
              {
                elems.push_back(_elems[i]);
              }
          }
        return;
      }

    //non terminal node 
    double min = bb[(_level%dim)*2];
    double max = bb[(_level%dim)*2+1];
    if (max < _min_right)
      {
        _left->getIntersectingElems(bb, elems);
        return;
      }
    if (min> _max_left)
      {
        _right->getIntersectingElems(bb,elems);
        return;
      }
    _left->getIntersectingElems(bb,elems);
    _right->getIntersectingElems(bb,elems);
  }

  /*!
   * This method is very close to getIntersectingElems except that it returns number of elems instead of elems themselves.
   */
  int getNbOfIntersectingElems(const double* bb)
  {
    //  terminal node : return list of elements intersecting bb
    int ret(0);
    if (_terminal)
      {
        for (int i=0; i<_nbelems; i++)
          {
            const double* const  bb_ptr=_bb+_elems[i]*2*dim;
            bool intersects = true;
            for (int idim=0; idim<dim; idim++)
              {
                if (bb_ptr[idim*2]-bb[idim*2+1]>-_epsilon|| bb_ptr[idim*2+1]-bb[idim*2]<_epsilon)
                  intersects=false;
              }
            if (intersects)
              ret++;
          }
        return ret;
      }
    //non terminal node 
    double min = bb[(_level%dim)*2];
    double max = bb[(_level%dim)*2+1];
    if (max < _min_right)
      return _left->getNbOfIntersectingElems(bb);
    if (min> _max_left)
      return _right->getNbOfIntersectingElems(bb);
    return _left->getNbOfIntersectingElems(bb)+_right->getNbOfIntersectingElems(bb);
  }

 
  /*! returns in \a elems the list of elements potentially containing the point pointed to by \a xx
    \param xx pointer to query point coords
    \param elems list of elements (given in 0-indexing) intersecting the bounding box
  */

  void getElementsAroundPoint(const double* xx, std::vector<ConnType>& elems) const
  {
    //  terminal node : return list of elements intersecting bb
    if (_terminal)
      {
        for (ConnType i=0; i<_nbelems; i++)
          {
            const double* const  bb_ptr=_bb+_elems[i]*2*dim;
            bool intersects = true;
            for (int idim=0; idim<dim; idim++)
              {
                if (bb_ptr[idim*2]-xx[idim]>_epsilon|| bb_ptr[idim*2+1]-xx[idim]<-_epsilon)
                  intersects=false;
              }
            if (intersects)
              {
                elems.push_back(_elems[i]);
              }
          }
        return;
      }

    //non terminal node 
    if (xx[_level%dim] < _min_right)
      {
        _left->getElementsAroundPoint(xx, elems);
        return;
      }
    if (xx[_level%dim]> _max_left)
      {
        _right->getElementsAroundPoint(xx,elems);
        return;
      }
    _left->getElementsAroundPoint(xx,elems);
    _right->getElementsAroundPoint(xx,elems);
  }



  int size()
  {
    if (_terminal) return _nbelems;
    return _left->size()+_right->size();
  }
};
#endif
