#ifndef __BBTREE_H__
#define __BBTREE_H__

#include <vector>
#include <algorithm>

#include <iostream>
#include <math.h> // for HUGE constant
const int MIN_NB_ELEMS =15;
const int MAX_LEVEL=20;

template <int dim>
class BBTree
{

private:
  BBTree* _left;
  BBTree* _right;
  int _level;
  double _max_left;
  double _min_right;
  double* _bb;
  std::vector<int> _elems;
  bool _terminal;
  int _nbelems;

public:


  BBTree(double* bbs, int* elems, int level, int nbelems):
    _left(0), _right(0), _level(level), _bb(bbs), _terminal(false),_nbelems(nbelems)
  {
    if (nbelems < MIN_NB_ELEMS || level> MAX_LEVEL)
      {
        _terminal=true;
      
      }
    double* nodes=new double [nbelems];
    _elems.resize(nbelems);
    for (int i=0; i<nbelems; i++)
      {
        int elem;
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

    std::vector<int> new_elems_left;
    std::vector<int> new_elems_right;
 
    new_elems_left.reserve(nbelems/2+1);
    new_elems_right.reserve(nbelems/2+1);
    double max_left=-HUGE;
    double min_right=HUGE;
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
    _max_left=max_left+1e-12;
    _min_right=min_right-1e-12;
    _left=new BBTree(bbs, &(new_elems_left[0]), level+1, new_elems_left.size());
    _right=new BBTree(bbs, &(new_elems_right[0]), level+1, new_elems_right.size());
  

  }


 
  ~BBTree()
  {
    if (_left!=0)  delete _left;
    if (_right!=0) delete _right;

  }

  
  /*! returns in \a elems the list of elements potentially intersecting the bounding box pointed to by \a bb
    
    \param bb pointer to query bounding box
    \param elems list of elements (given in 0-indexing) intersecting the bounding box
  */

  void getIntersectingElems(const double* bb, std::vector<int>& elems) const
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
                if (bb_ptr[idim*2]-bb[idim*2+1]>-1e-12|| bb_ptr[idim*2+1]-bb[idim*2]<1e-12)
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

 
  /*! returns in \a elems the list of elements potentially containing the point pointed to by \a xx
    \param xx pointer to query point coords
    \param elems list of elements (given in 0-indexing) intersecting the bounding box
  */

  void getElementsAroundPoint(const double* xx, std::vector<int>& elems) const
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
                if (bb_ptr[idim*2]-xx[idim]>1e-12|| bb_ptr[idim*2+1]-xx[idim]<-1e-12)
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
