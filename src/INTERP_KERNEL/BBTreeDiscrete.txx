// Copyright (C) 2024  CEA, EDF
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

#pragma once

#include <vector>
#include <algorithm>

#include <memory>
#include <limits>
#include <cmath>


// DataType is the type of data to locate
// ConnType is the type of IDs returned
template <int dim, class DataType, class ConnType>
class BBTreeDiscrete
{
private:
  std::unique_ptr<BBTreeDiscrete>  _left;
  std::unique_ptr<BBTreeDiscrete> _right;
  int _level;
  DataType _max_left;
  DataType _min_right;
  const DataType *_bb;
  typename std::vector<ConnType> _elems;
  bool _terminal;
  ConnType _nbelems;

  static const int MIN_NB_ELEMS=15;
  static const int MAX_LEVEL=20;
public:
  BBTreeDiscrete() = default;
  /*!
    Constructor of the bounding box tree
    \param bbs pointer to the [x1 y1 x2 y2 ...] array containing the bounding boxes that are to be indexed.
    \param elems array to the indices of the elements contained in the BBTreeDiscrete
    \param level level in the BBTreeDiscrete recursive structure
    \param nbelems nb of elements in the BBTreeDiscrete
  */
  BBTreeDiscrete(const DataType* bbs, ConnType *elems, int level, ConnType nbelems):
    _level(level), _bb(bbs), _terminal(false),_nbelems(nbelems)
  {
    if (nbelems < MIN_NB_ELEMS || level> MAX_LEVEL)
      {
        _terminal=true;
      
      }
    DataType median = std::numeric_limits<DataType>::max();
    {
      std::unique_ptr<DataType[]> nodes( new DataType [nbelems] );
      _elems.resize(nbelems);
      for (ConnType i=0; i<nbelems; i++)
        {
          ConnType elem;
          if (elems)
            elem= elems[i];
          else
            elem=i;

          _elems[i]=elem;
          nodes[i]=bbs[elem*dim+(level%dim)];
        }
      if (_terminal) { return; }

      std::nth_element<DataType*>(nodes.get(), nodes.get()+nbelems/2, nodes.get()+nbelems);
      median = nodes[nbelems/2];
    }

    std::vector<ConnType> new_elems_left;
    std::vector<ConnType> new_elems_right;
 
    new_elems_left.reserve(nbelems/2+1);
    new_elems_right.reserve(nbelems/2+1);
    DataType max_left = -std::numeric_limits<DataType>::max();
    DataType min_right=  std::numeric_limits<DataType>::max();
    for (ConnType i=0; i<nbelems;i++)
      {
        ConnType elem;
        if( elems )
          elem= elems[i];
        else
          elem=i;
        
        DataType value = bbs[elem*dim+(level%dim)];
      
        if (value >= median)
          {
            new_elems_right.push_back(elem);
            if (value<min_right) min_right = value;
          }
        else
          {
            new_elems_left.push_back(elem);
            if (value>max_left) max_left = value;
          }
      }
    _max_left = max_left;
    _min_right = min_right;
    ConnType *tmp( nullptr );
    if(!new_elems_left.empty())
      tmp = new_elems_left.data();
    _left.reset(new BBTreeDiscrete(bbs, tmp, level+1, (ConnType)new_elems_left.size()) );
    tmp = nullptr;
    if(!new_elems_right.empty())
      tmp = new_elems_right.data();
    _right.reset(new BBTreeDiscrete(bbs, tmp, level+1, (ConnType)new_elems_right.size()) );
  
  }
 
  ~BBTreeDiscrete() = default;

  /*! returns in \a elems the list of elements potentially intersecting the bounding box pointed to by \a bb
    
    \param bb pointer to query bounding box
    \param elems list of elements (given in 0-indexing that is to say in \b C \b mode) intersecting the bounding box
  */
  void getIntersectingElems(const DataType *bb, std::vector<ConnType>& elems) const
  {
    //  terminal node : return list of elements intersecting bb
    if (_terminal)
      {
        for (ConnType i=0; i<_nbelems; i++)
          {
            const DataType * const  bb_ptr = _bb + _elems[i]*dim;
            bool intersects = true;
            for (int idim=0; idim<dim; idim++)
              {
                if( bb_ptr[idim] != bb[idim] )
                  intersects=false;
              }
            if(intersects)
              {
                elems.push_back(_elems[i]);
              }
          }
        return ;
      }

    //non terminal node 
    DataType value = bb[_level%dim];
    if (value < _min_right)
      {
        _left->getIntersectingElems(bb, elems);
        return;
      }
    if (value > _max_left)
      {
        _right->getIntersectingElems(bb,elems);
        return;
      }
    _left->getIntersectingElems(bb,elems);
    _right->getIntersectingElems(bb,elems);
  }

  ConnType size()
  {
    if (_terminal) return _nbelems;
    return _left->size()+_right->size();
  }
};
