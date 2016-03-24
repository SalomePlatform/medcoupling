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
// Author : Anthony Geay (CEA/DEN)

#ifndef __INTERPKERNELGEO2DBOUNDS_HXX__
#define __INTERPKERNELGEO2DBOUNDS_HXX__

#include "INTERPKERNELDefines.hxx"

#include <algorithm>

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

  class INTERPKERNEL_EXPORT Bounds
  {
  public:
    Bounds():_x_min(0.),_x_max(0.),_y_min(0.),_y_max(0.) { }
    double &operator[](int i);
    const double& operator[](int i) const;
    double getXMin() const { return _x_min; }
    double getXMax() const { return _x_max; }
    double getYMin() const { return _y_min; }
    double getYMax() const { return _y_max; }
    double getDiagonal() const;
    void getBarycenter(double& xBary, double& yBary) const;
    void applySimilarity(double xBary, double yBary, double dimChar);
    void unApplySimilarity(double xBary, double yBary, double dimChar);
    Bounds& operator=(const Bounds& other) { _x_min=other._x_min; _x_max=other._x_max; _y_min=other._y_min; _y_max=other._y_max; return *this; }
    Bounds(double xMin, double xMax, double yMin, double yMax):_x_min(xMin),_x_max(xMax),_y_min(yMin),_y_max(yMax) { }
    void setValues(double xMin, double xMax, double yMin, double yMax) { _x_min=xMin; _x_max=xMax; _y_min=yMin; _y_max=yMax; }
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
    double getCaracteristicDim() const { return std::max(_x_max-_x_min,_y_max-_y_min); }
  protected:
    double _x_min;
    double _x_max;
    double _y_min;
    double _y_max;
  };
}

#endif
