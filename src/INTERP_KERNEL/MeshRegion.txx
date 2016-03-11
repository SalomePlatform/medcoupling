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
#ifndef __MESHREGION_TXX__
#define __MESHREGION_TXX__

#include "MeshRegion.hxx"

#include "MeshElement.txx"
#include "MeshUtils.hxx"

namespace INTERP_KERNEL
{
    
  /**
   * Default constructor
   * 
   */
  template<class ConnType>
  MeshRegion<ConnType>::MeshRegion():_box(0)
  {
  }
    
  /**
   * Destructor
   *
   */
  template<class ConnType>
  MeshRegion<ConnType>::~MeshRegion()
  {
    delete _box;
  }
  
  /**
   * Adds an element to the region, updating the bounding box. If the bounding box does not yet
   * exist, it is created here. This creation is delayed to make it possible to have empty MeshRegions
   *
   * @param element pointer to element to add to region
   * @param mesh    mesh to which element belongs
   *
   */
  template<class ConnType>
  template<class MyMeshType>
  void MeshRegion<ConnType>::addElement(MeshElement<ConnType>* const element, const MyMeshType& mesh)
  {
    _elements.push_back(element);

    const unsigned char numNodes = element->getNumberOfNodes();
    const ConnType elemIdx = element->getIndex();
       
    if(_box == 0)
      {
        const double** pts = new const double*[numNodes];

        // get coordinates of the nodes of the element
        for(unsigned char i = 0 ; i < numNodes ; ++i)
          {
            pts[i] = getCoordsOfNode(i, OTT<typename MyMeshType::MyConnType,MyMeshType::My_numPol>::indFC(elemIdx), mesh);
          }
           
        _box = new BoundingBox(pts, numNodes);
        delete [] pts;

      } else {

      for(unsigned char i = 0 ; i < numNodes ; ++i)
        {
          const double* pt = getCoordsOfNode(i, OTT<typename MyMeshType::MyConnType,MyMeshType::My_numPol>::indFC(elemIdx), mesh);
          _box->updateWithPoint(pt);
        }
    }
  }

  /**
   * Splits the region in two along the given axis, copying the elements with bounding boxes whose maximum
   * coordinate along the axis are smaller than the middle of the bounding box of this region in region1. The
   * rest of the elements are copied to region2.
   *
   * @param region1 region in which to store one half of this region
   * @param region2 region in which to store the other of this region
   * @param coord   coordinate of BoundingBox to use when splitting the region
   * @param mesh    mesh to which region belongs
   *
   */
  template<class ConnType>
  template<class MyMeshType>
  void MeshRegion<ConnType>::split(MeshRegion<ConnType>& region1, MeshRegion<ConnType>& region2, BoundingBox::BoxCoord coord, const MyMeshType& mesh)
  {
    // create ordering
    ElementBBoxOrder cmp(coord);

    // sort elements by their bounding boxes
    std::sort(_elements.begin(), _elements.end(), cmp);

    // put the first half of the elements in region1 and the 
    // rest in region2
    typename std::vector< MeshElement<ConnType> *>::const_iterator iter = _elements.begin();
    int elemCount = 0;

    while(elemCount < static_cast<int>(_elements.size() / 2))
      {
        region1.addElement(*iter, mesh);
        ++iter;
        ++elemCount;
      }

    while(iter != _elements.end())
      {
        region2.addElement(*iter, mesh);
        ++iter;
      }
  }

  /**
   * Determines if a given element can intersect the elements of this region by 
   * testing whether the bounding box of the region intersects the bounding box of the element.
   * Note that the test is only true in one direction : if the bounding boxes are disjoint, the
   * element cannot intersect any of the elements in the region, but if they are not disjoint, the 
   * element may or may not do so.
   *
   * @param   elem  Element with which to test for disjoint-ness
   * @return  true if the bounding box of the element is disjoint with the bounding box of the region, false otherwise
   */
  template<class ConnType>
  bool MeshRegion<ConnType>::isDisjointWithElementBoundingBox(const MeshElement<ConnType>& elem) const
  {
    const BoundingBox* elemBox = elem.getBoundingBox();

    assert(_box != 0);
    assert(elemBox != 0);

    return _box->isDisjointWith(*elemBox);
  }  

  
}

#endif
