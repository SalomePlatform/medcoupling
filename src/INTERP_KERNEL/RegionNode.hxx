#ifndef __REGIONNODE_HXX__
#define __REGIONNODE_HXX__

#include "MeshRegion.hxx"

namespace INTERP_KERNEL
{

  /**
   * \brief Class containing a tuplet of a source region and a target region. 
   * This is used as the object to put on the stack in the depth-first search
   * in the bounding-box filtering process.
   */
  template<class ConnType>
  class RegionNode
  {
  public:
    
    RegionNode() { }
    
    ~RegionNode() { }
    
    /**
     *  Accessor to source region
     *
     * @return   reference to source region
     */
    MeshRegion<ConnType>& getSrcRegion() { return _srcRegion; }

    /**
     *  Accessor to target region
     *
     * @return   reference to target region
     */
    MeshRegion<ConnType>& getTargetRegion() { return _targetRegion; }

  private:
    
    /// source region
    MeshRegion<ConnType> _srcRegion;          
    
    /// target region
    MeshRegion<ConnType> _targetRegion;       

  };

}

#endif
