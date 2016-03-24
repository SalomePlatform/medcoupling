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

#include "InterpKernelFunction.hxx"
#include "InterpKernelValue.hxx"

#include <cmath>
#include <limits>

using namespace INTERP_KERNEL;

const char IdentityFunction::REPR[]="Id";

const char PositiveFunction::REPR[]="+";

const char NegateFunction::REPR[]="-";

const char CosFunction::REPR[]="cos";

const char SinFunction::REPR[]="sin";

const char TanFunction::REPR[]="tan";

const char ACosFunction::REPR[]="acos";

const char ASinFunction::REPR[]="asin";

const char ATanFunction::REPR[]="atan";

const char CoshFunction::REPR[]="cosh";

const char SinhFunction::REPR[]="sinh";

const char TanhFunction::REPR[]="tanh";

const char SqrtFunction::REPR[]="sqrt";

const char AbsFunction::REPR[]="abs";

const char PlusFunction::REPR[]="+";

const char MinusFunction::REPR[]="-";

const char MultFunction::REPR[]="*";

const char DivFunction::REPR[]="/";

const char PowFunction::REPR[]="^";

const char ExpFunction::REPR[]="exp";

const char LnFunction::REPR[]="ln";

const char LogFunction::REPR[]="log";

const char Log10Function::REPR[]="log10";

const char MaxFunction::REPR[]="max";

const char MinFunction::REPR[]="min";

const char GreaterThanFunction::REPR[]=">";

const char LowerThanFunction::REPR[]="<";

const char IfFunction::REPR[]="if";

Function *FunctionsFactory::buildFuncFromString(const char *type, int nbOfParams)
{
  switch(nbOfParams)
    {
    case 1:
      return buildUnaryFuncFromString(type);
    case 2:
      return buildBinaryFuncFromString(type);
    case 3:
      return buildTernaryFuncFromString(type);
    default:
      throw INTERP_KERNEL::Exception("Invalid number of params detected : limited to 2 !");
    }
}

Function *FunctionsFactory::buildUnaryFuncFromString(const char *type)
{
  std::string tmp(type);
  if(tmp.empty())
    return new IdentityFunction;
  if(tmp==CosFunction::REPR)
    return new CosFunction;
  if(tmp==SinFunction::REPR)
    return new SinFunction;
  if(tmp==TanFunction::REPR)
    return new TanFunction;
  if(tmp==ACosFunction::REPR)
    return new ACosFunction;
  if(tmp==ASinFunction::REPR)
    return new ASinFunction;
  if(tmp==ATanFunction::REPR)
    return new ATanFunction;
  if(tmp==CoshFunction::REPR)
    return new CoshFunction;
  if(tmp==SinhFunction::REPR)
    return new SinhFunction;
  if(tmp==TanhFunction::REPR)
    return new TanhFunction;
  if(tmp==SqrtFunction::REPR)
    return new SqrtFunction;
  if(tmp==AbsFunction::REPR)
    return new AbsFunction;
  if(tmp==PositiveFunction::REPR)
    return new PositiveFunction;
  if(tmp==NegateFunction::REPR)
    return new NegateFunction;
  if(tmp==ExpFunction::REPR)
    return new ExpFunction;
  if(tmp==LnFunction::REPR)
    return new LnFunction;
  if(tmp==LogFunction::REPR)
    return new LogFunction;
  if(tmp==Log10Function::REPR)
    return new Log10Function;
  //
  std::string msg("Invalid unary function detected : \"");
  msg+=type; msg+="\"";
  throw INTERP_KERNEL::Exception(msg.c_str());
}

Function *FunctionsFactory::buildBinaryFuncFromString(const char *type)
{
  std::string tmp(type);
  if(tmp==PositiveFunction::REPR)
    return new PlusFunction;
  if(tmp==NegateFunction::REPR)
    return new MinusFunction;
  if(tmp==MultFunction::REPR)
    return new MultFunction;
  if(tmp==DivFunction::REPR)
    return new DivFunction;
  if(tmp==PowFunction::REPR)
    return new PowFunction;
  if(tmp==MaxFunction::REPR)
    return new MaxFunction;
  if(tmp==MinFunction::REPR)
    return new MinFunction;
  if(tmp==GreaterThanFunction::REPR)
    return new GreaterThanFunction;
  if(tmp==LowerThanFunction::REPR)
    return new LowerThanFunction;
  std::string msg("Invalid binary function detected : \"");
  msg+=type; msg+="\"";
  throw INTERP_KERNEL::Exception(msg.c_str());
}

Function *FunctionsFactory::buildTernaryFuncFromString(const char *type)
{
  std::string tmp(type);
  if(tmp==IfFunction::REPR)
    return new IfFunction();
  std::string msg("Invalid ternary function detected : \"");
  msg+=type; msg+="\"";
  throw INTERP_KERNEL::Exception(msg.c_str());
}

Function *FunctionsFactory::buildBinaryFuncFromString(char type)
{
  char tmp[2]; tmp[0]=type; tmp[1]='\0';
  return buildBinaryFuncFromString(tmp);
}

Function::~Function()
{
}

IdentityFunction::~IdentityFunction()
{
}

void IdentityFunction::operate(std::vector<Value *>& stck) const
{
}

void IdentityFunction::operateX86(std::vector<std::string>& asmb) const
{
}

void IdentityFunction::operateStackOfDouble(std::vector<double>& stck) const
{
}

const char *IdentityFunction::getRepr() const
{
  return REPR;
}

bool IdentityFunction::isACall() const
{
  return false;
}

PositiveFunction::~PositiveFunction()
{
}

int UnaryFunction::getNbInputParams() const
{
  return 1;
}

void PositiveFunction::operate(std::vector<Value *>& stck) const
{
}

void PositiveFunction::operateX86(std::vector<std::string>& asmb) const
{
}

void PositiveFunction::operateStackOfDouble(std::vector<double>& stck) const
{
}

const char *PositiveFunction::getRepr() const
{
  return REPR;
}

bool PositiveFunction::isACall() const
{
  return false;
}

NegateFunction::~NegateFunction()
{
}

void NegateFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->negate();
}

void NegateFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("fchs");
}

void NegateFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=-v;
}

const char *NegateFunction::getRepr() const
{
  return REPR;
}

bool NegateFunction::isACall() const
{
  return false;
}

CosFunction::~CosFunction()
{
}

void CosFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->cos();
}

void CosFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("fcos");
}

void CosFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=cos(v);
}

const char *CosFunction::getRepr() const
{
  return REPR;
}

bool CosFunction::isACall() const
{
  return true;
}

SinFunction::~SinFunction()
{
}

void SinFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->sin();
}

void SinFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("fsin");
}

void SinFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=sin(v);
}

const char *SinFunction::getRepr() const
{
  return REPR;
}

bool SinFunction::isACall() const
{
  return true;
}

TanFunction::~TanFunction()
{
}

void TanFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->tan();
}

void TanFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void TanFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=tan(v);
}

const char *TanFunction::getRepr() const
{
  return REPR;
}

bool TanFunction::isACall() const
{
  return true;
}

ACosFunction::~ACosFunction()
{
}

void ACosFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->acos();
}

void ACosFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void ACosFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=acos(v);
}

void ACosFunction::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double v(stck.back());
  if(fabs(v)>1.)
    throw INTERP_KERNEL::Exception("acos on a value which absolute is > 1 !");
  stck.back()=acos(v);
}

const char *ACosFunction::getRepr() const
{
  return REPR;
}

bool ACosFunction::isACall() const
{
  return true;
}

ASinFunction::~ASinFunction()
{
}

void ASinFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->asin();
}

void ASinFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void ASinFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=asin(v);
}

void ASinFunction::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double v(stck.back());
  if(fabs(v)>1.)
    throw INTERP_KERNEL::Exception("asin on a value which absolute is > 1 !");
  stck.back()=asin(v);
}

const char *ASinFunction::getRepr() const
{
  return REPR;
}

bool ASinFunction::isACall() const
{
  return true;
}

ATanFunction::~ATanFunction()
{
}

void ATanFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->atan();
}

void ATanFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void ATanFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=atan(v);
}

const char *ATanFunction::getRepr() const
{
  return REPR;
}

bool ATanFunction::isACall() const
{
  return true;
}

CoshFunction::~CoshFunction()
{
}

void CoshFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->cosh();
}

void CoshFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void CoshFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=cosh(v);
}

const char *CoshFunction::getRepr() const
{
  return REPR;
}

bool CoshFunction::isACall() const
{
  return true;
}

SinhFunction::~SinhFunction()
{
}

void SinhFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->sinh();
}

void SinhFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void SinhFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=sinh(v);
}

const char *SinhFunction::getRepr() const
{
  return REPR;
}

bool SinhFunction::isACall() const
{
  return true;
}

TanhFunction::~TanhFunction()
{
}

void TanhFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->tanh();
}

void TanhFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void TanhFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=tanh(v);
}

const char *TanhFunction::getRepr() const
{
  return REPR;
}

bool TanhFunction::isACall() const
{
  return true;
}

SqrtFunction::~SqrtFunction()
{
}

void SqrtFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->sqrt();
}

void SqrtFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("fsqrt");
}

void SqrtFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=sqrt(v);
}

void SqrtFunction::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double v(stck.back());
  if(v<0.)
    throw INTERP_KERNEL::Exception("sqrt on a value < 0. !");
  stck.back()=sqrt(v);
}

const char *SqrtFunction::getRepr() const
{
  return REPR;
}

bool SqrtFunction::isACall() const
{
  return true;
}

AbsFunction::~AbsFunction()
{
}

void AbsFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->abs();
}

void AbsFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("fabs");
}

void AbsFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=fabs(v);
}

const char *AbsFunction::getRepr() const
{
  return REPR;
}

bool AbsFunction::isACall() const
{
  return false;
}

void ExpFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->exp();
}

void ExpFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void ExpFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=std::exp(v);
}

const char *ExpFunction::getRepr() const
{
  return REPR;
}

bool ExpFunction::isACall() const
{
  return true;
}

LnFunction::~LnFunction()
{
}

void LnFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->ln();
}

void LnFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void LnFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=std::log(v);
}

void LnFunction::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double v(stck.back());
  if(v<0.)
    throw INTERP_KERNEL::Exception("ln on a value < 0. !");
  stck.back()=std::log(v);
}

const char *LnFunction::getRepr() const
{
  return REPR;
}

bool LnFunction::isACall() const
{
  return true;
}

LogFunction::~LogFunction()
{
}

void LogFunction::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->ln();
}

void LogFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly for log Not implemented yet !");
}

void LogFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=std::log(v);
}

void LogFunction::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double v(stck.back());
  if(v<0.)
    throw INTERP_KERNEL::Exception("log on a value < 0. !");
  stck.back()=std::log(v);
}

const char *LogFunction::getRepr() const
{
  return REPR;
}

bool LogFunction::isACall() const
{
  return true;
}

Log10Function::~Log10Function()
{
}

void Log10Function::operate(std::vector<Value *>& stck) const
{
  Value *val=stck.back();
  val->log10();
}

void Log10Function::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly for log Not implemented yet !");
}

void Log10Function::operateStackOfDouble(std::vector<double>& stck) const
{
  double v(stck.back());
  stck.back()=std::log10(v);
}

void Log10Function::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double v(stck.back());
  if(v<0.)
    throw INTERP_KERNEL::Exception("log10 on a value < 0. !");
  stck.back()=std::log10(v);
}

const char *Log10Function::getRepr() const
{
  return REPR;
}

bool Log10Function::isACall() const
{
  return true;
}

int BinaryFunction::getNbInputParams() const
{
  return 2;
}

PlusFunction::~PlusFunction()
{
}

void PlusFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3;
  try
    {
      val3=val1->plus(val2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      throw e;
    }
  delete val1;
  delete val2;
  val2=val3;
}

void PlusFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("faddp st1");
}

void PlusFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  stck.back()=a+stck.back();
}

const char *PlusFunction::getRepr() const
{
  return REPR;
}

bool PlusFunction::isACall() const
{
  return false;
}

MinusFunction::~MinusFunction()
{
}

void MinusFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3;
  try
    {
      val3=val1->minus(val2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      throw e;
    }
  delete val1;
  delete val2;
  val2=val3;
}

void MinusFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("fsubp st1");
}

void MinusFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  stck.back()=a-stck.back();
}

const char *MinusFunction::getRepr() const
{
  return REPR;
}

bool MinusFunction::isACall() const
{
  return false;
}

MultFunction::~MultFunction()
{
}

void MultFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3=val1->mult(val2);
  delete val1;
  delete val2;
  val2=val3;
}

void MultFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("fmulp st1");
}

void MultFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  stck.back()=a*stck.back();
}

const char *MultFunction::getRepr() const
{
  return REPR;
}

bool MultFunction::isACall() const
{
  return false;
}

DivFunction::~DivFunction()
{
}

void DivFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3;
  try
    {
      val3=val1->div(val2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      throw e;
    }
  delete val1;
  delete val2;
  val2=val3;
}

void DivFunction::operateX86(std::vector<std::string>& asmb) const
{
  asmb.push_back("fdivp st1");
}

void DivFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  stck.back()=a/stck.back();
}

void DivFunction::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  if(stck.back()==0.)
    throw INTERP_KERNEL::Exception("division by 0. !");
  stck.back()=a/stck.back();
}

const char *DivFunction::getRepr() const
{
  return REPR;
}

bool DivFunction::isACall() const
{
  return false;
}

PowFunction::~PowFunction()
{
}

void PowFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3;
  try
    {
      val3=val1->pow(val2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      throw e;
    }
  delete val1;
  delete val2;
  val2=val3;
}

void PowFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void PowFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  stck.back()=std::pow(a,stck.back());
}

void PowFunction::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  double b(stck.back());
  if(a<0.)
    throw INTERP_KERNEL::Exception("pow with val < 0. !");
  stck.back()=std::pow(a,b);
}

const char *PowFunction::getRepr() const
{
  return REPR;
}

bool PowFunction::isACall() const
{
  return true;
}

ExpFunction::~ExpFunction()
{
}

MaxFunction::~MaxFunction()
{
}

void MaxFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3;
  try
    {
      val3=val1->max(val2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      throw e;
    }
  delete val1;
  delete val2;
  val2=val3;
}

void MaxFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void MaxFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  stck.back()=std::max(stck.back(),a);
}

const char *MaxFunction::getRepr() const
{
  return REPR;
}

bool MaxFunction::isACall() const
{
  return false;
}

MinFunction::~MinFunction()
{
}

void MinFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3;
  try
    {
      val3=val1->min(val2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      throw e;
    }
  delete val1;
  delete val2;
  val2=val3;
}

void MinFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void MinFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  stck.back()=std::min(stck.back(),a);
}

const char *MinFunction::getRepr() const
{
  return REPR;
}

bool MinFunction::isACall() const
{
  return false;
}

GreaterThanFunction::~GreaterThanFunction()
{
}

void GreaterThanFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3;
  try
    {
      val3=val1->greaterThan(val2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      throw e;
    }
  delete val1;
  delete val2;
  val2=val3;
}

void GreaterThanFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void GreaterThanFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  double b(stck.back());
  stck.back()=a>b?std::numeric_limits<double>::max():-std::numeric_limits<double>::max();
}

const char *GreaterThanFunction::getRepr() const
{
  return REPR;
}

bool GreaterThanFunction::isACall() const
{
  return false;
}

LowerThanFunction::~LowerThanFunction()
{
}

void LowerThanFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *& val2=stck.back();
  Value *val3;
  try
    {
      val3=val1->lowerThan(val2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      throw e;
    }
  delete val1;
  delete val2;
  val2=val3;
}

void LowerThanFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void LowerThanFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double a(stck.back());
  stck.pop_back();
  double b(stck.back());
  stck.back()=a<b?std::numeric_limits<double>::max():-std::numeric_limits<double>::max();
}

const char *LowerThanFunction::getRepr() const
{
  return REPR;
}

bool LowerThanFunction::isACall() const
{
  return false;
}

int TernaryFunction::getNbInputParams() const
{
  return 3;
}

IfFunction::~IfFunction()
{
}

void IfFunction::operate(std::vector<Value *>& stck) const
{
  Value *val1=stck.back();
  stck.pop_back();
  Value *val2=stck.back();
  stck.pop_back();
  Value *&val3=stck.back();
  Value *val4;
  try
    {
      val4=val1->ifFunc(val2,val3);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete val1;
      delete val2;
      throw e;
    }
  delete val1;
  delete val2;
  delete val3;
  val3=val4;
}

void IfFunction::operateX86(std::vector<std::string>& asmb) const
{
  throw INTERP_KERNEL::Exception("Assembly Not implemented yet !");
}

void IfFunction::operateStackOfDouble(std::vector<double>& stck) const
{
  double cond(stck.back());
  stck.pop_back();
  double the(stck.back());
  stck.pop_back();
  if(cond==std::numeric_limits<double>::max())
    stck.back()=the;
}

void IfFunction::operateStackOfDoubleSafe(std::vector<double>& stck) const
{
  double cond(stck.back());
  stck.pop_back();
  double the(stck.back());
  stck.pop_back();
  if(cond!=std::numeric_limits<double>::max() && cond!=-std::numeric_limits<double>::max())
    throw INTERP_KERNEL::Exception("ifFunc : first parameter of ternary func is NOT a consequence of a boolean op !");
  if(cond==std::numeric_limits<double>::max())
    stck.back()=the;
}

const char *IfFunction::getRepr() const
{
  return REPR;
}

bool IfFunction::isACall() const
{
  return false;
}

