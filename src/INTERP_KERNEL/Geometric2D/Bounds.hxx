#ifndef __BOUNDS_HXX__
#define __BOUNDS_HXX__

#include "Geometric2D_defines.hxx"
#include <cmath>

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
  
  class GEOMETRIC2D_EXPORT Bounds
  {
  public:
    Bounds():_xMin(0.),_xMax(0.),_yMin(0.),_yMax(0.) { }
    double &operator[](int i);
    const double& operator[](int i) const;
    double getDiagonal() const;
    void getBarycenter(double& xBary, double& yBary) const;
    void applySimilarity(double xBary, double yBary, double dimChar);
    Bounds& operator=(const Bounds& other) { _xMin=other._xMin; _xMax=other._xMax; _yMin=other._yMin; _yMax=other._yMax; return *this; }
    Bounds(double xMin, double xMax, double yMin, double yMax):_xMin(xMin),_xMax(xMax),_yMin(yMin),_yMax(yMax) { }
    void setValues(double xMin, double xMax, double yMin, double yMax) { _xMin=xMin; _xMax=xMax; _yMin=yMin; _yMax=yMax; }
    void prepareForAggregation();
    void getInterceptedArc(const double *center, double radius, double& intrcptArcAngle0, double& intrcptArcDelta) const;
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
    double getCaracteristicDim() const { return fmax(_xMax-_xMin,_yMax-_yMin); }
  protected:
    double _xMin;
    double _xMax;
    double _yMin;
    double _yMax;
  };
}

#endif
