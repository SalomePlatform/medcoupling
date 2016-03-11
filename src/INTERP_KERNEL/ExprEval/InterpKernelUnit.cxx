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

#include "InterpKernelUnit.hxx"
#include "InterpKernelExprParser.hxx"

#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>

using namespace INTERP_KERNEL;

UnitDataBase UnitDataBase::_uniqueMapForExpr;

static const char InterpKernelMuAscii[2]={-0x4B,0x0};

static const char InterpKernelMuUnicode[3]={-0x3E,-0x4B,0x0};

const char *UnitDataBase::PREF_POW10[NB_OF_PREF_POW10]={"y","z","a","f","p","n",InterpKernelMuAscii,InterpKernelMuUnicode,"u","m","c","d",
                                                        "da","h","k","M","G","T","P","E","Z","Y"};

const double UnitDataBase::POW10[NB_OF_PREF_POW10]={1e-24,1e-21,1e-18,1e-15,1e-12,1e-9,1e-6,1e-6,1e-6,1e-3,1e-2,1e-1,
                                                  1e1,1e2,1e3,1e6,1e9,1e12,1e15,1e18,1e21,1e24};

static const char InterpKernelDegreeCAscii[3]={-0x50,0x43,0x0};

static const char InterpKernelDegreeCUnicode[4]={-0x3E,-0x50,0x43,0x0};

static const char InterpKernelDegreeCUnicodeWin[3]={-0x08,0x43,0x0};

const char *UnitDataBase::UNITS_RECOGN[NB_OF_UNITS_RECOGN]={"g","m","s","A","K",
                                                            "W","J","Hz","V","h","min","t","N","dyn",
                                                            "eV","Pa","atm","bar",InterpKernelDegreeCAscii,"C","ohm","F","S",
                                                            "T","H","P","St",InterpKernelDegreeCUnicode,InterpKernelDegreeCUnicodeWin};

const short UnitDataBase::PROJ_IN_BASE[NB_OF_UNITS_RECOGN][SIZE_OF_UNIT_BASE]=
  {
    {1,0,0,0,0},//g
    {0,1,0,0,0},//m
    {0,0,1,0,0},//s
    {0,0,0,1,0},//A
    {0,0,0,0,1},//K
    {1,2,-3,0,0},//W
    {1,2,-2,0,0},//J
    {0,0,-1,0,0},//Hz
    {1,2,-3,-1,0},//V
    {0,0,1,0,0},//h
    {0,0,1,0,0},//min
    {1,0,0,0,0},//t
    {1,1,-2,0,0},//N
    {1,1,-2,0,0},//dyn
    {1,2,-2,0,0},//eV
    {1,-1,-2,0,0},//Pa
    {1,-1,-2,0,0},//atm
    {1,-1,-2,0,0},//bar
    {0,0,0,0,1},//degree C
    {0,0,1,1,0},//C
    {1,2,-3,-2,0},//ohm
    {-1,-2,4,2,0},//F
    {-1,-2,3,2,0},//S
    {1,0,-2,-1,0},//T
    {1,2,-2,-2,0},//H
    {1,-1,-1,0,0},//P
    {0,2,-1,0,0},//St
    {0,0,0,0,1},//degree C
    {0,0,0,0,1}//degree C
  };

const double UnitDataBase::MUL_COEFF[NB_OF_UNITS_RECOGN]=
  { 1.,1.,1.,1.,1.,
    1000.,1000.,1.,1000.,3600.,3600.,1e6,1000.,1e-2,
    1.60217733e-16,1000.,1.01325e8,1e8,1.,1.,1000.,1e-3,
    1000.,1000.,100.,1.,1.,1.,1.};

const double UnitDataBase::ADD_COEFF[NB_OF_UNITS_RECOGN]=
  { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 273.15, 0., 0., 0., 0., 0., 0., 0., 0., 273.15 ,273.15};

UnitDataBase::UnitDataBase()
{
  for(int i=0;i<NB_OF_PREF_POW10;i++)
    _prefix_pow_10[PREF_POW10[i]]=POW10[i];
  for(int i=0;i<NB_OF_UNITS_RECOGN;i++)
    {
      _units_semantic[UNITS_RECOGN[i]]=PROJ_IN_BASE[i];
      _units_mul[UNITS_RECOGN[i]]=MUL_COEFF[i];
      _units_add[UNITS_RECOGN[i]]=ADD_COEFF[i];
    }
}

const short *UnitDataBase::getInfoForUnit(const std::string& unit, double& addFact, double& mFact) const
{
  std::size_t lgth=unit.length();
  std::string work,work2;
  const short *ret=0;
  for(std::size_t i=0;i<lgth && !ret;i++)
    {
      work=unit.substr(i);
      std::map<std::string,const short *>::const_iterator iter=_units_semantic.find(work);
      if(iter!=_units_semantic.end())
        {
          ret=(*iter).second;
          std::map<std::string,double>::const_iterator iter2=_units_add.find(work);
          addFact=(*iter2).second;
          std::map<std::string,double>::const_iterator iter3=_units_mul.find(work);
          mFact=(*iter3).second;
          work2=unit.substr(0,i);
        }
    }
  if(!ret)
    {
      std::ostringstream os;
      os << "Unit : " << unit << " not recognized !";
      throw INTERP_KERNEL::Exception(os.str().c_str());
    }
  if(!work2.empty())
    {
      std::map<std::string,double>::const_iterator iter4=_prefix_pow_10.find(work2);
      if(iter4==_prefix_pow_10.end())
        {
          std::ostringstream os;
          os << "Unit : " << unit << " not fully recognized : \"" << work << "\" detected as core unit and \"";
          os << work2 << "\" not recognized prefix !";
          throw INTERP_KERNEL::Exception(os.str().c_str());
        }
      addFact=0.;
      mFact*=(*iter4).second;
    }
  return ret;
}

DecompositionInUnitBase::DecompositionInUnitBase():_add_to_base(0.),_mult_fact_to_base(1.)
{
  _value[0]=0;
  _value[1]=0;
  _value[2]=0;
  _value[3]=0;
  _value[4]=0;
}

void DecompositionInUnitBase::setInfo(const short *vals, double addFact, double mFact)
{
  _add_to_base=addFact;
  _mult_fact_to_base=mFact;
  _value[0]=vals[0];
  _value[1]=vals[1];
  _value[2]=vals[2];
  _value[3]=vals[3];
  _value[4]=vals[4];
}

bool DecompositionInUnitBase::operator==(const DecompositionInUnitBase& other) const
{
  return _value[0]==other._value[0] && _value[1]==other._value[1] && _value[2]==other._value[2] && _value[3]==other._value[3] && _value[4]==other._value[4];
}

void DecompositionInUnitBase::getTranslationParams(const DecompositionInUnitBase& other, double& mul, double& add) const
{
  if((*this)==other)
    {
      mul=_mult_fact_to_base/other._mult_fact_to_base;
      add=_add_to_base/other._mult_fact_to_base-other._add_to_base;
    }
  else
    {
      mul=std::numeric_limits<double>::max();
      add=std::numeric_limits<double>::max();
    }
}

bool DecompositionInUnitBase::isEqual(short mass, short lgth, short time, short intensity, short temp, double add, double mult)
{
  bool ret1=mass==_value[0];
  bool ret2=lgth==_value[1];
  bool ret3=time==_value[2];
  bool ret4=intensity==_value[3];
  bool ret5=temp==_value[4];
  bool ret6=areDoubleEquals(add,_add_to_base);
  bool ret7=areDoubleEquals(mult,_mult_fact_to_base);
  return ret1 && ret2 && ret3 && ret4 && ret5 && ret6 && ret7;
}

void DecompositionInUnitBase::negate()
{
  _mult_fact_to_base=-_mult_fact_to_base;
}

bool DecompositionInUnitBase::isAdimensional() const
{
  return _value[0]==0 && _value[1]==0 && _value[2]==0 && _value[3]==0 && _value[4]==0;
}

bool DecompositionInUnitBase::isUnitary() const
{
  return areDoubleEquals(_add_to_base,0.) && areDoubleEquals(_mult_fact_to_base,1.);
}

void DecompositionInUnitBase::tryToConvertInUnit(double val)
{
  int valI=(int)val;
  if((val-(double)valI)!=0.)
    {
      std::ostringstream os;
      os << "Double value " << val << " can't be considered as integer. Not admitable for units !";
      throw INTERP_KERNEL::Exception(os.str().c_str());
    }
  _value[0]=0;
  _value[1]=0;
  _value[2]=0;
  _value[3]=0;
  _value[4]=0;
  _add_to_base=0;
  _mult_fact_to_base=valI;
}

DecompositionInUnitBase &DecompositionInUnitBase::operator*(const DecompositionInUnitBase& other)
{
  _value[0]+=other._value[0]; _value[1]+=other._value[1]; _value[2]+=other._value[2]; _value[3]+=other._value[3]; _value[4]+=other._value[4];
  _mult_fact_to_base*=other._mult_fact_to_base;
  _add_to_base=0.;
  return *this;
}

DecompositionInUnitBase &DecompositionInUnitBase::operator/(const DecompositionInUnitBase& other)
{
  _value[0]-=other._value[0]; _value[1]-=other._value[1]; _value[2]-=other._value[2]; _value[3]-=other._value[3]; _value[4]-=other._value[4];
  _mult_fact_to_base/=other._mult_fact_to_base;
 _add_to_base=0.;
 return *this;
}

DecompositionInUnitBase &DecompositionInUnitBase::operator^(const DecompositionInUnitBase& other)
{
  if(!other.isAdimensional())
    throw INTERP_KERNEL::Exception("Trying to execute operator ^ with a second member not adimensionnal");
  int exp=couldItBeConsideredAsInt(other._mult_fact_to_base);
  _value[0]*=exp; _value[1]*=exp; _value[2]*=exp; _value[3]*=exp; _value[4]*=exp;
  _mult_fact_to_base=powInt(_mult_fact_to_base,exp);
  _add_to_base=0.;
  return *this;
}

void DecompositionInUnitBase::dealWithAddFactor(const DecompositionInUnitBase& other)
{
  if(!areDoubleEquals(_add_to_base,0.))
    if(other.isAdimensional())
      if(areDoubleEquals(other._mult_fact_to_base,1.))
        return ;
  if(!other.areDoubleEquals(_add_to_base,0.))
    if(isAdimensional())
      if(areDoubleEquals(_mult_fact_to_base,1.))
        return ;
  _add_to_base=0.;
}

double DecompositionInUnitBase::powInt(double val, int exp)
{
  double work=1.;
  if(exp==0)
    return 1.;
  if(exp>0)
    for(int i=0;i<exp;i++)
      work*=val;
  else
    {
      int tmp=-exp;
      for(int i=0;i<tmp;i++)
        work*=1/val;
    }
  return work;
}

bool DecompositionInUnitBase::areDoubleEquals(double a, double b)
{
  if(a==0. || b==0.)
    return a==b;
  double ref=std::max(a,b);
  return fabs((a-b)/ref)<1e-7;
}

int DecompositionInUnitBase::couldItBeConsideredAsInt(double val)
{
  int ret=(int)val;
  double valT=(double) ret;
  if(valT==val)
    return ret;
  else
    {
      std::ostringstream stream; stream << "Invalid double number " << std::setprecision(16) << val << " can's be considered for ^ operation on unit.";
      throw INTERP_KERNEL::Exception(stream.str().c_str());
    }
}

Unit::Unit(const char *reprC, bool tryToInterp):_coarse_repr(reprC),
                                                _is_interpreted(false),
                                                _is_interpretation_ok(false)
{
  if(tryToInterp)
    tryToInterprate();
}

Unit::Unit(const char *reprFortran, int sizeOfRepr, bool tryToInterp):_coarse_repr(ExprParser::buildStringFromFortran(reprFortran,sizeOfRepr)),
                                                                      _is_interpreted(false),
                                                                      _is_interpretation_ok(false)
{
}

void Unit::tryToInterprate() const
{
  if(!_is_interpreted)
    {
      _is_interpreted=true;
      _is_interpretation_ok=false;
      try
        {
          ExprParser expr(_coarse_repr.c_str());
          expr.parse();
          _decomp_in_base=expr.evaluateUnit();
          _is_interpretation_ok=true;
        }
      catch(INTERP_KERNEL::Exception&) { }
    }
}

bool Unit::isInterpretationOK() const
{
  return _is_interpretation_ok;
}

bool Unit::isCompatibleWith(const Unit& other) const
{
  tryToInterprate();
  other.tryToInterprate();
  if(_is_interpretation_ok && other._is_interpretation_ok)
    return _decomp_in_base==other._decomp_in_base;
  else
    return false;
}

double Unit::convert(const Unit& target, double sourceVal) const
{
  if(isCompatibleWith(target))
    {
      double mult,add;
      _decomp_in_base.getTranslationParams(target._decomp_in_base,mult,add);
      return mult*sourceVal+add;
    }
  else
    return std::numeric_limits<double>::max();
}

std::string Unit::getCoarseRepr() const
{
  return _coarse_repr;
}
