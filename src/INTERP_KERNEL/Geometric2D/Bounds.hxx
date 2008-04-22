#ifndef __BOUNDS_HXX__
#define __BOUNDS_HXX__

namespace INTERP_KERNEL
{
  /*!
   * Relative LOC
   */
  typedef enum
    {
      IN              = 0,
      OUT             = 1,
      ON_BOUNDARY_POS = 2,
      ON_BOUNDARY_NEG = 3
    } Position;
  
  class Bounds
  {
  public:
    Bounds():_xMin(0),_xMax(0.),_yMin(0.),_yMax(0.) { }
    double &operator[](int i);
    const double& operator[](int i) const;
    Bounds& operator=(const Bounds& other) { _xMin=other._xMin; _xMax=other._xMax; _yMin=other._yMin; _yMax=other._yMax; return *this; }
    Bounds(double xMin, double xMax, double yMin, double yMax):_xMin(xMin),_xMax(xMax),_yMin(yMin),_yMax(yMax) { }
    void setValues(double xMin, double xMax, double yMin, double yMax) { _xMin=xMin; _xMax=xMax; _yMin=yMin; _yMax=yMax; }
    void prepareForAggregation();
    int fitXForXFig(double val, int res) const { return (int)fitXForXFigD(val,res); }
    int fitYForXFig(double val, int res) const { return (int)fitYForXFigD(val,res); }
    double fitXForXFigD(double val, int res) const;
    double fitYForXFigD(double val, int res) const;
    Bounds *nearlyAmIIntersectingWith(const Bounds& other) const;
    Bounds *amIIntersectingWith(const Bounds& other) const;
    //! No approximations.
    Position where(double x, double y) const;
    //! Idem where method but with approximations.
    Position nearlyWhere(double x, double y) const;
    void aggregate(const Bounds& other);
    double getCaracteristicDim() const;
  protected:
    double _xMin;
    double _xMax;
    double _yMin;
    double _yMax;
  };
}

#endif
