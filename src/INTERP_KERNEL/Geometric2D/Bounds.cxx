#include "Bounds.hxx"
#include "InterpolationUtils.hxx"
#include "EdgeArcCircle.hxx"
#include "Node.hxx"

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

double Bounds::getDiagonal() const
{
  double a=_xMax-_xMin;
  double b=_yMax-_yMin;
  return sqrt(a*a+b*b);
}

/*!
 * See Node::applySimilarity to see signification of params.
 */
void Bounds::applySimilarity(double xBary, double yBary, double dimChar)
{
  _xMin=(_xMin-xBary)/dimChar;
  _xMax=(_xMax-xBary)/dimChar;
  _yMin=(_yMin-yBary)/dimChar;
  _yMax=(_yMax-yBary)/dimChar;
}

void Bounds::getBarycenter(double& xBary, double& yBary) const
{
  xBary=(_xMin+_xMax)/2.;
  yBary=(_yMax+_yMin)/2.;
}

void Bounds::prepareForAggregation()
{
  _xMin=1e200; _xMax=-1e200; _yMin=1e200; _yMax=-1e200;
}

/*! 
 * Given an arc defined by 'center', 'radius' and 'intrcptArcDelta' in radian, returns (by outputs intrcptArcAngle0 and intrcptArcDelta)
 * the intercepted angle of 'this' from 'center' point of view.
 * If diagonal of 'this' is the same order of 2*radius, intrcptArcAngle0 and intrcptArcDelta remains unchanged.
 * @param center IN parameter.
 * @param radius IN parameter.
 * @param intrcptArcAngle0 OUT parameter.
 * @param intrcptArcDelta IN/OUT parameter.
 */
void Bounds::getInterceptedArc(const double *center, double radius, double& intrcptArcAngle0, double& intrcptArcDelta) const
{
  double diag=getDiagonal();
  if(diag<2.*radius)
    {
      double v1[2],v2[2],w1[2],w2[2];
      v1[0]=_xMin-center[0]; v1[1]=_yMax-center[1]; v2[0]=_xMax-center[0]; v2[1]=_yMin-center[1];
      w1[0]=v1[0]; w1[1]=_yMin-center[1];           w2[0]=v2[0]; w2[1]=_yMax-center[1];
      double delta1=EdgeArcCircle::safeAsin(v1[0]*v2[1]-v1[1]*v2[0]);
      double delta2=EdgeArcCircle::safeAsin(w1[0]*w2[1]-w1[1]*w2[0]);
      double tmp;
      if(fabs(delta1)>fabs(delta2))
        {
          intrcptArcDelta=delta1;
          intrcptArcAngle0=EdgeArcCircle::getAbsoluteAngle(v1,tmp);
        }
      else
        {
          intrcptArcDelta=delta2;
          intrcptArcAngle0=EdgeArcCircle::getAbsoluteAngle(w1,tmp);
        }
    }
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
