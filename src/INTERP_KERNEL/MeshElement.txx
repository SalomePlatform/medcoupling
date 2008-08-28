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
   * @param index   global number of element in the mesh
   * @param mesh    mesh that the element belongs to
   */
  template<class ConnType>
  template<int SPACEDIM, int MESHDIM, NumberingPolicy numPol, class MyMeshType>
  MeshElement<ConnType>::MeshElement(const ConnType index, const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh)
    : _index(index), _box(0), _number(mesh.getNumberOfNodesOfElement(index))
  {
    const double**vertices = new const double*[_number];

    for(unsigned char i = 0 ; i < _number ; ++i)
      vertices[i] = getCoordsOfNode(i , index, mesh);

    // create bounding box
    _box = new BoundingBox(vertices,_number);
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
   * Constructor
   *
   * @param  coord   BoundingBox coordinate (XMIN, XMAX, etc) on which to base the ordering
   */
  ElementBBoxOrder::ElementBBoxOrder(BoundingBox::BoxCoord coord)
    : _coord(coord)
  {
  }

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
