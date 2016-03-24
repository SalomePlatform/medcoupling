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

#ifndef __MESHREGION_HXX__
#define __MESHREGION_HXX__

#include "MeshElement.hxx"
#include "BoundingBox.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  /**
   * \brief Class representing a set of elements in a mesh together with their bounding box.
   * It permits to split itself in two, which is used in the depth-first search filtering process.
   *
   */
  template<class ConnType>
  class MeshRegion
  {
  public:
    
    MeshRegion();

    ~MeshRegion();

    template<class MyMeshType>
    void addElement(MeshElement<ConnType>* const element, const MyMeshType& mesh);

    template<class MyMeshType>
    void split(MeshRegion<ConnType>& region1, MeshRegion<ConnType>& region2, BoundingBox::BoxCoord coord, const MyMeshType& mesh);

    bool isDisjointWithElementBoundingBox(const MeshElement<ConnType>& elem) const;
    /**
     * Accessor to beginning of elements vector
     *
     * @return  constant iterator pointing at the beginning of the vector or elements
     */
    typename std::vector< MeshElement<ConnType>* >::const_iterator getBeginElements() const { return _elements.begin(); }

    /**
     * Accessor to end of elements vector
     *
     * @return  constant iterator pointing at the end of the vector or elements
     */
    typename std::vector< MeshElement<ConnType>* >::const_iterator getEndElements() const { return _elements.end(); }
    
    /**
     * Gives information on how many elements are contained in the region.
     *
     * @return  the number of elements contained in the region
     */
    unsigned getNumberOfElements() const { return _elements.size(); }

  private:

    /// disallow copying
    MeshRegion(const MeshRegion& m);

    /// disallow assignment
    MeshRegion<ConnType>& operator=(const MeshRegion<ConnType>& m);

    /// Vector of pointers to contained MeshElements. 
    /// NB : these pointers are not owned by the region object, and are thus
    /// neither allocated or liberated in this class. The elements must therefore be allocated and liberated outside the class.
    std::vector< MeshElement<ConnType>* > _elements;

    /// BoundingBox containing all the nodes of all the elements in the region.
    BoundingBox* _box;
  
  };

}

#endif
