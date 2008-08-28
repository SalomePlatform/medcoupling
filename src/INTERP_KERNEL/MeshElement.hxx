#ifndef __MESHELEMENT_HXX__
#define __MESHELEMENT_HXX__

#include "BoundingBox.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{

  /**
   * \brief Class representing a single element of a mesh together with its bounding box.
   * It gives access to the element's global number, type and bounding box and allows
   * easy bounding box intersection tests between MeshElements and collections of MeshElement (MeshRegions)
   */
  template<class ConnType>
  class INTERPKERNEL_EXPORT MeshElement
  {

  public:
    template<int SPACEDIM, int MESHDIM, NumberingPolicy numPol, class MyMeshType>
    MeshElement(const ConnType index, const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh);
    
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
  class INTERPKERNEL_EXPORT ElementBBoxOrder
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
