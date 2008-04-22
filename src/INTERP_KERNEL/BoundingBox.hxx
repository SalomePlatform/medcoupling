#ifndef __BOUNDINGBOX_HXX__
#define __BOUNDINGBOX_HXX__

#include <iostream>

namespace INTERP_KERNEL
{

  /**
   * \brief Class representing the bounding box of a number of points.
   *
   */
  class BoundingBox 
  {
  public:

    /// Enumeration representing the six coordinates that define the bounding box
    enum BoxCoord { XMIN = 0, YMIN = 1, ZMIN = 2, XMAX = 3, YMAX = 4, ZMAX = 5 };
        
    BoundingBox(const double** pts, const unsigned numPts);

    BoundingBox(const BoundingBox& box1, const BoundingBox& box2);

    ~BoundingBox();

    bool isDisjointWith(const BoundingBox& box) const;
    
    inline void setCoordinate(const BoxCoord coord, double value);

    inline double getCoordinate(const BoxCoord coord) const;

    void updateWithPoint(const double* pt);

    inline void dumpCoords() const;

  private:
    
    bool isValid() const;

    /// disallow copying
    BoundingBox(const BoundingBox& box);
    
    /// disallow assignment
    BoundingBox& operator=(const BoundingBox& box);
    
    /// Vector containing the coordinates of the box
    /// interlaced in the order XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX
    double* _coords;

  };

  /**
   * Sets a coordinate of the box to a given value.
   * 
   * @param coord coordinate to set
   * @param value new value for coordinate
   *
   */
  inline void BoundingBox::setCoordinate(const BoxCoord coord, double value)
  {
    _coords[coord] = value;
  }

  /**
   * Gets a coordinate of the box
   * 
   * @param coord coordinate to get
   * @return value of coordinate
   *
   */
  inline double BoundingBox::getCoordinate(const BoxCoord coord) const
  {
    return _coords[coord];
  }

  /**
   * Prints the coordinates of the box to std::cout
   *
   */
  inline void BoundingBox::dumpCoords() const
  {
    std::cout << "[xmin, xmax] = [" << _coords[XMIN] << ", " << _coords[XMAX] << "]" << " | ";
    std::cout << "[ymin, ymax] = [" << _coords[YMIN] << ", " << _coords[YMAX] << "]" << " | ";
    std::cout << "[zmin, zmax] = [" << _coords[ZMIN] << ", " << _coords[ZMAX] << "]";
    std::cout << std::endl;
  }

};

#endif
