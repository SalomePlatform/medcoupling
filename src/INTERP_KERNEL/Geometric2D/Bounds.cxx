#include "Bounds.hxx"
#include "InterpolationUtils.hxx"
#include "Node.hxx"

#include <cmath>

using namespace INTERP_KERNEL;

const double& Bounds::operator[](int i) const
{
  switch(i)
    {
    case 0:
      return _xMin;
    case 1:
      return _xMax;
    case 2:
      return _yMin;
    case 3:
      return _yMax;
    }
  throw Exception("internal error occurs !");
}

double &Bounds::operator[](int i)
{
  switch(i)
    {
    case 0:
      return _xMin;
    case 1:
      return _xMax;
    case 2:
      return _yMin;
    case 3:
      return _yMax;
    }
  throw Exception("internal error occurs !");
}

void Bounds::prepareForAggregation()
{
  _xMin=1e200; _xMax=-1e200; _yMin=1e200; _yMax=-1e200;
}

double Bounds::fitXForXFigD(double val, int res) const
{
  double delta=fmax(_xMax-_xMin,_yMax-_yMin)/2.;
  double ret=val-(_xMax+_xMin)/2.+delta;
  delta=11.1375*res/(2.*delta);
  return ret*delta;
}

double Bounds::fitYForXFigD(double val, int res) const
{
  double delta=fmax(_xMax-_xMin,_yMax-_yMin)/2.;
  double ret=val-(_yMax+_yMin)/2.+delta;
  delta=11.1375*res/(2.*delta);
  return ret*delta;
}

Bounds *Bounds::nearlyAmIIntersectingWith(const Bounds& other) const
{
  if( (other._xMin > _xMax+QUADRATIC_PLANAR::_precision) || (other._xMax < _xMin-QUADRATIC_PLANAR::_precision) || (other._yMin > _yMax+QUADRATIC_PLANAR::_precision) 
      || (other._yMax < _yMin-QUADRATIC_PLANAR::_precision) )
    return 0;
  if( (other._xMin >= _xMax ) || (other._xMax <= _xMin) || (other._yMin >= _yMax) || (other._yMax <= _yMin) )
    return new Bounds(fmax(_xMin-QUADRATIC_PLANAR::_precision,other._xMin),
                      fmin(_xMax+QUADRATIC_PLANAR::_precision,other._xMax),
                      fmax(_yMin-QUADRATIC_PLANAR::_precision,other._yMin),
                      fmin(_yMax+QUADRATIC_PLANAR::_precision,other._yMax));//In approx cases.
  else
    return new Bounds(fmax(_xMin,other._xMin),fmin(_xMax,other._xMax),fmax(_yMin,other._yMin),fmin(_yMax,other._yMax));
}

Bounds *Bounds::amIIntersectingWith(const Bounds& other) const
{
  if( (other._xMin > _xMax) || (other._xMax < _xMin) || (other._yMin > _yMax) || (other._yMax < _yMin) )
    return 0;
  return new Bounds(fmax(_xMin,other._xMin),fmin(_xMax,other._xMax),fmax(_yMin,other._yMin),fmin(_yMax,other._yMax));
}

Position Bounds::where(double x, double y) const
{
  if((x>=_xMin && x<=_xMax) && (y>=_yMin && y<=_yMax))
    return IN;
  else
    return OUT;
}

Position Bounds::nearlyWhere(double x, double y) const
{
  bool thinX=Node::areDoubleEquals(_xMin,_xMax);
  bool thinY=Node::areDoubleEquals(_yMin,_yMax);
  if(!thinX)
    {
      if(Node::areDoubleEquals(x,_xMin) || Node::areDoubleEquals(x,_xMax) && (y<_yMax+QUADRATIC_PLANAR::_precision) && (y>_yMin-QUADRATIC_PLANAR::_precision))
        return ON_BOUNDARY_POS;
    }
  else
    if(!Node::areDoubleEquals(_xMin,x) && !Node::areDoubleEquals(_xMax,x))
      return OUT;
  if(!thinY)
    {
      if(Node::areDoubleEquals(y,_yMin) || Node::areDoubleEquals(y,_yMax) && (x<_xMax+QUADRATIC_PLANAR::_precision) && (x>_xMin-QUADRATIC_PLANAR::_precision))
        return ON_BOUNDARY_POS;
    }
  else
    if(!Node::areDoubleEquals(_yMin,y) && !Node::areDoubleEquals(_yMax,y))
      return OUT;
  if(thinX && thinY)
    return ON_BOUNDARY_POS;
  if((x>=_xMin && x<=_xMax) && (y>=_yMin && y<=_yMax))
    return IN;
  else
    return OUT;
}

void Bounds::aggregate(const Bounds& other)
{
  _xMin=fmin(_xMin,other._xMin); _xMax=fmax(_xMax,other._xMax);
  _yMin=fmin(_yMin,other._yMin); _yMax=fmax(_yMax,other._yMax);
}

double Bounds::getCaracteristicDim() const
{
  return fmax(_xMax-_xMin,_yMax-_yMin);
}
