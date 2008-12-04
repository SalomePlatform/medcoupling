//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
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
