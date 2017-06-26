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

#include "InterpKernelGeo2DBounds.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DNode.hxx"

using namespace INTERP_KERNEL;

const double& Bounds::operator[](int i) const
{
  switch(i)
  {
    case 0:
      return _x_min;
    case 1:
      return _x_max;
    case 2:
      return _y_min;
    case 3:
      return _y_max;
  }
  throw Exception("internal error occurs !");
}

double &Bounds::operator[](int i)
{
  switch(i)
  {
    case 0:
      return _x_min;
    case 1:
      return _x_max;
    case 2:
      return _y_min;
    case 3:
      return _y_max;
  }
  throw Exception("internal error occurs !");
}

double Bounds::getDiagonal() const
{
  double a=_x_max-_x_min;
  double b=_y_max-_y_min;
  return sqrt(a*a+b*b);
}

/*!
 * See Node::applySimilarity to see signification of params.
 */
void Bounds::applySimilarity(double xBary, double yBary, double dimChar)
{
  _x_min=(_x_min-xBary)/dimChar;
  _x_max=(_x_max-xBary)/dimChar;
  _y_min=(_y_min-yBary)/dimChar;
  _y_max=(_y_max-yBary)/dimChar;
}

/*!
 * See Node::unApplySimilarity to see signification of params.
 */
void Bounds::unApplySimilarity(double xBary, double yBary, double dimChar)
{
  _x_min=_x_min*dimChar+xBary;
  _x_max=_x_max*dimChar+xBary;
  _y_min=_y_min*dimChar+yBary;
  _y_max=_y_max*dimChar+yBary;
}

void Bounds::getBarycenter(double& xBary, double& yBary) const
{
  xBary=(_x_min+_x_max)/2.;
  yBary=(_y_max+_y_min)/2.;
}

void Bounds::prepareForAggregation()
{
  _x_min=1e200; _x_max=-1e200; _y_min=1e200; _y_max=-1e200;
}

/*! 
 * Given an arc defined by 'center', 'radius' and 'intrcptArcDelta' in radian, returns (by outputs intrcptArcAngle0 and intrcptArcDelta)
 * the intercepted angle of 'this' from 'center' point of view.
 * If diagonal of 'this' is the same order of 2*radius, intrcptArcAngle0 and intrcptArcDelta remains unchanged.
 * @param center IN parameter.
 * @param radius IN parameter.
 * @param [out] intrcptArcAngle0 OUT parameter.
 * @param [out] intrcptArcDelta OUT parameter.
 */
void Bounds::getInterceptedArc(const double *center, double radius, double& intrcptArcAngle0, double& intrcptArcDelta) const
{
  double diag=getDiagonal();
  if(diag<2.*radius)
    {
      double v1[2],v2[2],w1[2],w2[2];
      v1[0]=_x_min-center[0]; v1[1]=_y_max-center[1]; v2[0]=_x_max-center[0]; v2[1]=_y_min-center[1];
      w1[0]=v1[0]; w1[1]=_y_min-center[1];           w2[0]=v2[0]; w2[1]=_y_max-center[1];
      double delta1=EdgeArcCircle::SafeAsin(v1[0]*v2[1]-v1[1]*v2[0]);
      double delta2=EdgeArcCircle::SafeAsin(w1[0]*w2[1]-w1[1]*w2[0]);
      double tmp;
      if(fabs(delta1)>fabs(delta2))
        {
          intrcptArcDelta=delta1;
          intrcptArcAngle0=EdgeArcCircle::GetAbsoluteAngle(v1,tmp);
        }
      else
        {
          intrcptArcDelta=delta2;
          intrcptArcAngle0=EdgeArcCircle::GetAbsoluteAngle(w1,tmp);
        }
    }
}

double Bounds::fitXForXFigD(double val, int res) const
{
  double delta=std::max(_x_max-_x_min,_y_max-_y_min)/2.;
  double ret=val-(_x_max+_x_min)/2.+delta;
  delta=11.1375*res/(2.*delta);
  return ret*delta;
}

double Bounds::fitYForXFigD(double val, int res) const
{
  double delta=std::max(_x_max-_x_min,_y_max-_y_min)/2.;
  double ret=(_y_max+_y_min)/2.-val+delta;
  delta=11.1375*res/(2.*delta);
  return ret*delta;
}

Bounds *Bounds::nearlyAmIIntersectingWith(const Bounds& other) const
{
  double eps = QuadraticPlanarPrecision::getPrecision();
  if( (other._x_min > _x_max+eps) || (other._x_max < _x_min-eps) || (other._y_min > _y_max+eps)
      || (other._y_max < _y_min-eps) )
    return 0;
  if( (other._x_min >= _x_max ) || (other._x_max <= _x_min) || (other._y_min >= _y_max) || (other._y_max <= _y_min) )
    {
      return new Bounds(std::max(_x_min-eps,other._x_min),
          std::min(_x_max+eps,other._x_max),
          std::max(_y_min-eps,other._y_min),
          std::min(_y_max+eps,other._y_max));//In approx cases.
    }
  else
    return new Bounds(std::max(_x_min,other._x_min),std::min(_x_max,other._x_max),std::max(_y_min,other._y_min),std::min(_y_max,other._y_max));
}

Bounds *Bounds::amIIntersectingWith(const Bounds& other) const
{
  if( (other._x_min > _x_max) || (other._x_max < _x_min) || (other._y_min > _y_max) || (other._y_max < _y_min) )
    return 0;
  return new Bounds(std::max(_x_min,other._x_min),std::min(_x_max,other._x_max),std::max(_y_min,other._y_min),std::min(_y_max,other._y_max));
}

Position Bounds::where(double x, double y) const
{
  if((x>=_x_min && x<=_x_max) && (y>=_y_min && y<=_y_max))
    return IN;
  else
    return OUT;
}

Position Bounds::nearlyWhere(double x, double y) const
{
  bool thinX=Node::areDoubleEquals(_x_min,_x_max);
  bool thinY=Node::areDoubleEquals(_y_min,_y_max);
  if(!thinX)
    {
      if((Node::areDoubleEquals(x,_x_min) || Node::areDoubleEquals(x,_x_max)) && ((y<_y_max+QuadraticPlanarPrecision::getPrecision()) && (y>_y_min-QuadraticPlanarPrecision::getPrecision())))
        return ON_BOUNDARY_POS;
    }
  else
    if(!Node::areDoubleEquals(_x_min,x) && !Node::areDoubleEquals(_x_max,x))
      return OUT;
  if(!thinY)
    {
      if((Node::areDoubleEquals(y,_y_min) || Node::areDoubleEquals(y,_y_max)) && ((x<_x_max+QuadraticPlanarPrecision::getPrecision()) && (x>_x_min-QuadraticPlanarPrecision::getPrecision())))
        return ON_BOUNDARY_POS;
    }
  else
    if(!Node::areDoubleEquals(_y_min,y) && !Node::areDoubleEquals(_y_max,y))
      return OUT;
  if(thinX && thinY)
    return ON_BOUNDARY_POS;
  if((x>=_x_min && x<=_x_max) && (y>=_y_min && y<=_y_max))
    return IN;
  else
    return OUT;
}

void Bounds::aggregate(const Bounds& other)
{
  _x_min=std::min(_x_min,other._x_min); _x_max=std::max(_x_max,other._x_max);
  _y_min=std::min(_y_min,other._y_min); _y_max=std::max(_y_max,other._y_max);
}
