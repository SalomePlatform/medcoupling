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

#include "InterpKernelValue.hxx"
#include "InterpKernelFunction.hxx"

#include <cmath>
#include <limits>
#include <algorithm>

using namespace INTERP_KERNEL;

ValueDouble::ValueDouble():_data(std::numeric_limits<double>::max())
{
}

Value *ValueDouble::newInstance() const
{
  return new ValueDouble;
}

ValueDouble::ValueDouble(double val):_data(val)
{
}

void ValueDouble::setDouble(double val)
{
  _data=val;
}

void ValueDouble::setVarname(int fastPos, const std::string& var)
{
  std::string msg("Error var : "); msg+=var; msg+=" not numeric : use another expression evaluator !";
  throw INTERP_KERNEL::Exception(msg.c_str());
}

void ValueDouble::positive()
{
}

void ValueDouble::negate()
{
  _data=-_data;
}

void ValueDouble::sqrt()
{
  _data=std::sqrt(_data);
}

void ValueDouble::cos()
{
  _data=std::cos(_data);
}

void ValueDouble::sin()
{
  _data=std::sin(_data);
}

void ValueDouble::tan()
{
  _data=std::tan(_data);
}

void ValueDouble::acos()
{
  _data=std::acos(_data);
}

void ValueDouble::asin()
{
  _data=std::asin(_data);
}

void ValueDouble::atan()
{
  _data=std::atan(_data);
}

void ValueDouble::cosh()
{
  _data=std::cosh(_data);
}

void ValueDouble::sinh()
{
  _data=std::sinh(_data);
}

void ValueDouble::tanh()
{
  _data=std::tanh(_data);
}

void ValueDouble::abs()
{
  if(_data<0.)
    _data=-_data;
}

void ValueDouble::exp()
{
  _data=std::exp(_data);
}

void ValueDouble::ln()
{
  _data=std::log(_data);
}

void ValueDouble::log10()
{
  _data=std::log10(_data);
}

Value *ValueDouble::plus(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(_data+valC->_data);
}

Value *ValueDouble::minus(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(_data-valC->_data);
}

Value *ValueDouble::mult(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(_data*valC->_data);
}

Value *ValueDouble::div(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(_data/valC->_data);
}

Value *ValueDouble::pow(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(std::pow(_data,valC->_data));
}

Value *ValueDouble::max(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(std::max(_data,valC->_data));
}

Value *ValueDouble::min(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(std::min(_data,valC->_data));
}

Value *ValueDouble::greaterThan(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(_data>valC->_data?std::numeric_limits<double>::max():-std::numeric_limits<double>::max());
}

Value *ValueDouble::lowerThan(const Value *other) const
{
  const ValueDouble *valC=checkSameType(other);
  return new ValueDouble(_data<valC->_data?std::numeric_limits<double>::max():-std::numeric_limits<double>::max());
}

Value *ValueDouble::ifFunc(const Value *the, const Value *els) const
{
  const ValueDouble *theC=checkSameType(the);
  const ValueDouble *elsC=checkSameType(els);
  if(_data==std::numeric_limits<double>::max())
    return new ValueDouble(theC->_data);
  if(_data==-std::numeric_limits<double>::max())
    return new ValueDouble(elsC->_data);
  throw INTERP_KERNEL::Exception("ValueDouble::ifFunc : The fist element of ternary function if is not a binary op !");
}

const ValueDouble *ValueDouble::checkSameType(const Value *val)
{
  const ValueDouble *valC=dynamic_cast<const ValueDouble *>(val);
  if(!valC)
    throw INTERP_KERNEL::Exception("Trying to operate on non homogeneous Values (double with other type) !");
  return valC;
}

ValueUnit::ValueUnit()
{
}

Value *ValueUnit::newInstance() const
{
  return new ValueUnit;
}

ValueUnit::ValueUnit(const DecompositionInUnitBase& unit):_data(unit)
{
}

void ValueUnit::setDouble(double val)
{
  _data.tryToConvertInUnit(val);
}

void ValueUnit::setVarname(int fastPos, const std::string& var)
{
  double add,mul;
  const short *projInBase=UnitDataBase::_uniqueMapForExpr.getInfoForUnit(var,add,mul);
  _data.setInfo(projInBase,add,mul);
}

void ValueUnit::positive()
{
  unsupportedOp(PositiveFunction::REPR);
}

void ValueUnit::negate()
{
  _data.negate();
}

void ValueUnit::sqrt()
{
  unsupportedOp(SqrtFunction::REPR);
}

void ValueUnit::cos()
{
  unsupportedOp(CosFunction::REPR);
}

void ValueUnit::sin()
{
  unsupportedOp(SinFunction::REPR);
}

void ValueUnit::tan()
{
  unsupportedOp(TanFunction::REPR);
}

void ValueUnit::acos()
{
  unsupportedOp(ACosFunction::REPR);
}

void ValueUnit::asin()
{
  unsupportedOp(ASinFunction::REPR);
}

void ValueUnit::atan()
{
  unsupportedOp(ATanFunction::REPR);
}

void ValueUnit::cosh()
{
  unsupportedOp(CoshFunction::REPR);
}

void ValueUnit::sinh()
{
  unsupportedOp(SinhFunction::REPR);
}

void ValueUnit::tanh()
{
  unsupportedOp(TanhFunction::REPR);
}

void ValueUnit::abs()
{
  unsupportedOp(AbsFunction::REPR);
}

void ValueUnit::exp()
{
  unsupportedOp(ExpFunction::REPR);
}

void ValueUnit::ln()
{
  unsupportedOp(LnFunction::REPR);
}

void ValueUnit::log10()
{
  unsupportedOp(Log10Function::REPR);
}

Value *ValueUnit::plus(const Value *other) const
{
  unsupportedOp(PlusFunction::REPR);
  return 0;
}

Value *ValueUnit::minus(const Value *other) const
{
  unsupportedOp(MinusFunction::REPR);
  return 0;
}

Value *ValueUnit::greaterThan(const Value *other) const
{
  unsupportedOp(GreaterThanFunction::REPR);
  return 0;
}

Value *ValueUnit::lowerThan(const Value *other) const
{
  unsupportedOp(LowerThanFunction::REPR);
  return 0;
}

Value *ValueUnit::ifFunc(const Value *the, const Value *els) const
{
  unsupportedOp(IfFunction::REPR);
  return 0;
}

Value *ValueUnit::mult(const Value *other) const
{
  const ValueUnit *valC=checkSameType(other);
  DecompositionInUnitBase tmp=_data;
  tmp*valC->getData();
  return new ValueUnit(tmp);
}

Value *ValueUnit::div(const Value *other) const
{
  const ValueUnit *valC=checkSameType(other);
  DecompositionInUnitBase tmp=_data;
  tmp/valC->getData();
  return new ValueUnit(tmp);
}

Value *ValueUnit::pow(const Value *other) const
{
  const ValueUnit *valC=checkSameType(other);
  DecompositionInUnitBase tmp=_data;
  tmp^valC->getData();
  return new ValueUnit(tmp);
}

Value *ValueUnit::max(const Value *other) const
{
  unsupportedOp(MaxFunction::REPR);
  return 0;
}

Value *ValueUnit::min(const Value *other) const
{
  unsupportedOp(MinFunction::REPR);
  return 0;
}

const ValueUnit *ValueUnit::checkSameType(const Value *val)
{
  const ValueUnit *valC=dynamic_cast<const ValueUnit *>(val);
  if(!valC)
    throw INTERP_KERNEL::Exception("Trying to operate on non homogeneous Values (Units with other type) !");
  return valC;
}

void ValueUnit::unsupportedOp(const char *type)
{
  const char msg[]="Unsupported operation for units :";
  std::string msgStr(msg);
  msgStr+=type;
  throw INTERP_KERNEL::Exception(msgStr.c_str());
}

ValueDoubleExpr::ValueDoubleExpr(int szDestData, const double *srcData):_sz_dest_data(szDestData),_dest_data(new double[_sz_dest_data]),_src_data(srcData)
{
}

ValueDoubleExpr::~ValueDoubleExpr()
{
  delete [] _dest_data;
}

Value *ValueDoubleExpr::newInstance() const
{
  return new ValueDoubleExpr(_sz_dest_data,_src_data);
}

void ValueDoubleExpr::setDouble(double val)
{
  std::fill(_dest_data,_dest_data+_sz_dest_data,val);
}

void ValueDoubleExpr::setVarname(int fastPos, const std::string& var)
{
  if(fastPos==-2)
    std::copy(_src_data,_src_data+_sz_dest_data,_dest_data);
  else if(fastPos>-2)
    std::fill(_dest_data,_dest_data+_sz_dest_data,_src_data[fastPos]);
  else
    {
      std::fill(_dest_data,_dest_data+_sz_dest_data,0.);
      _dest_data[-7-fastPos]=1.;
    }
}

void ValueDoubleExpr::positive()
{
}

void ValueDoubleExpr::negate()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::negate<double>());
}

void ValueDoubleExpr::sqrt()
{
  double *it=std::find_if(_dest_data,_dest_data+_sz_dest_data,std::bind2nd(std::less<double>(),0.));
  if(it!=_dest_data+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to apply sqrt on < 0. value !");
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::sqrt));
}

void ValueDoubleExpr::cos()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::cos));
}

void ValueDoubleExpr::sin()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::sin));
}

void ValueDoubleExpr::tan()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::tan));
}

void ValueDoubleExpr::acos()
{
  double *it=std::find_if(_dest_data,_dest_data+_sz_dest_data,std::bind2nd(std::less<double>(),-1.));
  if(it!=_dest_data+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to apply acos on < 1. value !");
  it=std::find_if(_dest_data,_dest_data+_sz_dest_data,std::bind2nd(std::greater<double>(),1.));
  if(it!=_dest_data+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to apply acos on > 1. value !");
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::acos));
}

void ValueDoubleExpr::asin()
{
   double *it=std::find_if(_dest_data,_dest_data+_sz_dest_data,std::bind2nd(std::less<double>(),-1.));
   if(it!=_dest_data+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to apply asin on < 1. value !");
  it=std::find_if(_dest_data,_dest_data+_sz_dest_data,std::bind2nd(std::greater<double>(),1.));
  if(it!=_dest_data+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to apply asin on > 1. value !");
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::asin));
}

void ValueDoubleExpr::atan()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::atan));
}

void ValueDoubleExpr::cosh()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::cosh));
}

void ValueDoubleExpr::sinh()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::sinh));
}

void ValueDoubleExpr::tanh()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::tanh));
}

void ValueDoubleExpr::abs()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(fabs));
}

void ValueDoubleExpr::exp()
{
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::exp));
}

void ValueDoubleExpr::ln()
{
  double *it=std::find_if(_dest_data,_dest_data+_sz_dest_data,std::bind2nd(std::less_equal<double>(),0.));
  if(it!=_dest_data+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to apply neperian/natural log on <= 0. value !");
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::log));
}

void ValueDoubleExpr::log10()
{
  double *it=std::find_if(_dest_data,_dest_data+_sz_dest_data,std::bind2nd(std::less_equal<double>(),0.));
  if(it!=_dest_data+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to apply log10 on <= 0. value !");
  std::transform(_dest_data,_dest_data+_sz_dest_data,_dest_data,std::ptr_fun<double,double>(std::log10));
}

Value *ValueDoubleExpr::plus(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  std::transform(_dest_data,_dest_data+_sz_dest_data,otherC->getData(),ret->getData(),std::plus<double>());
  return ret;
}

Value *ValueDoubleExpr::minus(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  std::transform(_dest_data,_dest_data+_sz_dest_data,otherC->getData(),ret->getData(),std::minus<double>());
  return ret;
}

Value *ValueDoubleExpr::mult(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  std::transform(_dest_data,_dest_data+_sz_dest_data,otherC->getData(),ret->getData(),std::multiplies<double>());
  return ret;
}

Value *ValueDoubleExpr::div(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  double *it=std::find(otherC->getData(),otherC->getData()+_sz_dest_data,0.);
  if(it!=otherC->getData()+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to operate division by 0. !");
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  std::transform(_dest_data,_dest_data+_sz_dest_data,otherC->getData(),ret->getData(),std::divides<double>());
  return ret;
}

Value *ValueDoubleExpr::pow(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  double p=otherC->getData()[0];
  double *it=std::find_if(_dest_data,_dest_data+_sz_dest_data,std::bind2nd(std::less<double>(),0.));
  if(it!=_dest_data+_sz_dest_data)
    throw INTERP_KERNEL::Exception("Trying to operate pow(a,b) with a<0. !");
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  std::transform(_dest_data,_dest_data+_sz_dest_data,ret->getData(),std::bind2nd(std::ptr_fun<double,double,double>(std::pow),p));
  return ret;
}

Value *ValueDoubleExpr::max(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  std::transform(_dest_data,_dest_data+_sz_dest_data,otherC->getData(),ret->getData(),std::ptr_fun<const double&, const double&, const double& >(std::max));
  return ret;
}

Value *ValueDoubleExpr::min(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  std::transform(_dest_data,_dest_data+_sz_dest_data,otherC->getData(),ret->getData(),std::ptr_fun<const double&, const double&, const double& >(std::min));
  return ret;
}

Value *ValueDoubleExpr::greaterThan(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  for(int i=0;i<_sz_dest_data;i++)
    if(_dest_data[i]<=otherC->getData()[i])
      {
        std::fill(ret->getData(),ret->getData()+_sz_dest_data,-std::numeric_limits<double>::max());
        return ret;
      }
  std::fill(ret->getData(),ret->getData()+_sz_dest_data,std::numeric_limits<double>::max());
  return ret;
}

Value *ValueDoubleExpr::lowerThan(const Value *other) const
{
  const ValueDoubleExpr *otherC=static_cast<const ValueDoubleExpr *>(other);
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  for(int i=0;i<_sz_dest_data;i++)
    if(_dest_data[i]>=otherC->getData()[i])
      {
        std::fill(ret->getData(),ret->getData()+_sz_dest_data,-std::numeric_limits<double>::max());
        return ret;
      }
  std::fill(ret->getData(),ret->getData()+_sz_dest_data,std::numeric_limits<double>::max());
  return ret;
}

Value *ValueDoubleExpr::ifFunc(const Value *the, const Value *els) const
{
  const ValueDoubleExpr *theC=static_cast<const ValueDoubleExpr *>(the);
  const ValueDoubleExpr *elsC=static_cast<const ValueDoubleExpr *>(els);
  ValueDoubleExpr *ret=new ValueDoubleExpr(_sz_dest_data,_src_data);
  bool okmax=true;
  bool okmin=true;
  for(int i=0;i<_sz_dest_data && (okmax || okmin);i++)
    {
      okmax=_dest_data[i]==std::numeric_limits<double>::max();
      okmin=_dest_data[i]==-std::numeric_limits<double>::max();
    }
  if(okmax || okmin)
    {
      if(okmax)
        std::copy(theC->getData(),theC->getData()+_sz_dest_data,ret->getData());
      else
        std::copy(elsC->getData(),elsC->getData()+_sz_dest_data,ret->getData());
      return ret;
    }
  else
    {
      throw INTERP_KERNEL::Exception("ValueDoubleExpr::ifFunc : first parameter of ternary func is NOT a consequence of a boolean op !");
    }
}
