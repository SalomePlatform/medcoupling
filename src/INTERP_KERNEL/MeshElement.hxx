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

#ifndef __MESHELEMENT_HXX__
#define __MESHELEMENT_HXX__

#include "BoundingBox.hxx"

namespace INTERP_KERNEL
{

  /**
   * \brief Class representing a single element of a mesh together with its bounding box.
   * It gives access to the element's global number, type and bounding box and allows
   * easy bounding box intersection tests between MeshElements and collections of MeshElement (MeshRegions)
   */
  template<class ConnType>
  class MeshElement
  {

  public:
    template<class MyMeshType>
    MeshElement(const ConnType index, const MyMeshType& mesh);
    
    ~MeshElement();
    
    ConnType getIndex() const { return _index; }
    
    unsigned char getNumberOfNodes() const { return _number; }
    
    const BoundingBox* getBoundingBox() const { return _box; }

  private:
    /// disallow copying
    MeshElement(const MeshElement& elem);

    /// disallow assignment
    MeshElement& operator=(const MeshElement& elem);

    /// global number of the element
    const ConnType _index;

    const unsigned char _number;
    
    /// bounding box of the element - does not change after having been initialised
    BoundingBox* _box;
  };

  /**
   * \brief Class defining an order for MeshElements based on their bounding boxes.
   * The order defined between two elements is that between a given coordinate of 
   * their bounding boxes. For instance, if the order is based on YMIN, an element whose boxes
   * has a smaller YMIN is sorted before one with a larger YMIN.
   *
   */
  class ElementBBoxOrder
  {
  public : 
    
    ElementBBoxOrder(BoundingBox::BoxCoord coord);
    template<class ConnType>
    bool operator()(MeshElement<ConnType>* elem1, MeshElement<ConnType>* elem2);
    
  private :
    /// BoundingBox coordinate (XMIN, XMAX, etc) on which to base the ordering
    BoundingBox::BoxCoord _coord;  
  };

}

#endif
