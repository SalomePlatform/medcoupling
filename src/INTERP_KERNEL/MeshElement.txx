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
#ifndef __MESHELEMENT_TXX__
#define __MESHELEMENT_TXX__

#include "MeshElement.hxx"

#include "TetraAffineTransform.hxx"
#include "TransformedTriangle.hxx"
#include "MeshUtils.hxx"
#include "BoundingBox.hxx"
#include <assert.h>

namespace INTERP_KERNEL
{

  /**
   * Constructor
   *
   * @param index   global number of element in the mesh in C mode.
   * @param mesh    mesh that the element belongs to
   */
  template<class ConnType>
  template<class MyMeshType>
  MeshElement<ConnType>::MeshElement(const ConnType index, const MyMeshType& mesh)
    : _index(index), _number(mesh.getNumberOfNodesOfElement(OTT<typename MyMeshType::MyConnType,MyMeshType::My_numPol>::indFC(index))), _box(0)
  {
    const double**vertices = new const double*[_number];

    for(unsigned char i = 0 ; i < _number ; ++i)
      vertices[i] = getCoordsOfNode(i , OTT<typename MyMeshType::MyConnType,MyMeshType::My_numPol>::indFC(index), mesh);

    // create bounding box
    _box = new BoundingBox(vertices,_number);
    delete [] vertices;
  }
    
  /**
   * Destructor
   *
   */
  template<class ConnType>
  MeshElement<ConnType>::~MeshElement()
  {
    delete _box;
  }

  

  /////////////////////////////////////////////////////////////////////
  /// ElementBBoxOrder                                    /////////////
  /////////////////////////////////////////////////////////////////////

  /**
   * Comparison operator based on the bounding boxes of the elements
   *
   * @return true if the coordinate _coord of the bounding box of elem1 is 
   *          strictly smaller than that of the bounding box of elem2
   */
  template<class ConnType>
  bool ElementBBoxOrder::operator()( MeshElement<ConnType>* elem1, MeshElement<ConnType>* elem2)
  {
    const BoundingBox* box1 = elem1->getBoundingBox();
    const BoundingBox* box2 = elem2->getBoundingBox();

    assert(elem1 != 0);
    assert(elem2 != 0);
    assert(box1 != 0);
    assert(box2 != 0);
    
    const double coord1 = box1->getCoordinate(_coord);
    const double coord2 = box2->getCoordinate(_coord);
    
    return coord1 < coord2;
  }

}

#endif
