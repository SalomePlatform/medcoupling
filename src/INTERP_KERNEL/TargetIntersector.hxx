#ifndef __TARGETINTERSECTOR__HXX__
#define __TARGETINTERSECTOR__HXX__

namespace INTERP_KERNEL
{
  /**
   * \brief Abstract base class of Intersector classes. 
   * These classes represent a target element and calculate its intersection
   * with the source elements.
   */
  template<class ConnType>
  class TargetIntersector
  {
  public:

    /// Virtual destructor
    virtual ~TargetIntersector() {}
    
    /**
     * Calculate the volume of the intersection between target cell 
     * and the given source cell.
     *
     * @param srcCell     global number of the source cell
     */
    virtual double intersectSourceCell(ConnType srcCell) = 0;
  };
};

#endif
