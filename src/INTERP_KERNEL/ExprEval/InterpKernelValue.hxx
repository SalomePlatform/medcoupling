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

#ifndef __INTERPKERNELVALUE_HXX__
#define __INTERPKERNELVALUE_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelUnit.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT Value
  {
  public:
    virtual Value *newInstance() const = 0;
    virtual ~Value() { }
    virtual void setDouble(double val) = 0;
    virtual void setVarname(int fastPos, const std::string& var) = 0;
    //unary
    virtual void positive() = 0;
    virtual void negate() = 0;
    virtual void sqrt() = 0;
    virtual void cos() = 0;
    virtual void sin() = 0;
    virtual void tan() = 0;
    virtual void acos() = 0;
    virtual void asin() = 0;
    virtual void atan() = 0;
    virtual void cosh() = 0;
    virtual void sinh() = 0;
    virtual void tanh() = 0;
    virtual void abs() = 0;
    virtual void exp() = 0;
    virtual void ln() = 0;
    virtual void log10() = 0;
    //binary
    virtual Value *plus(const Value *other) const = 0;
    virtual Value *minus(const Value *other) const = 0;
    virtual Value *mult(const Value *other) const = 0;
    virtual Value *div(const Value *other) const = 0;
    virtual Value *pow(const Value *other) const = 0;
    virtual Value *max(const Value *other) const = 0;
    virtual Value *min(const Value *other) const = 0;
    virtual Value *greaterThan(const Value *other) const = 0;
    virtual Value *lowerThan(const Value *other) const = 0;
    //ternary
    virtual Value *ifFunc(const Value *the, const Value *els) const = 0;
  };

  class INTERPKERNEL_EXPORT ValueDouble : public Value
  {
  public:
    ValueDouble();
    Value *newInstance() const;
    void setDouble(double val);
    void setVarname(int fastPos, const std::string& var);
    //
    double getData() const { return _data; }
    void positive();
    void negate();
    void sqrt();
    void cos();
    void sin();
    void tan();
    void acos();
    void asin();
    void atan();
    void cosh();
    void sinh();
    void tanh();
    void abs();
    void exp();
    void ln();
    void log10();
    //
    Value *plus(const Value *other) const;
    Value *minus(const Value *other) const;
    Value *mult(const Value *other) const;
    Value *div(const Value *other) const;
    Value *pow(const Value *other) const;
    Value *max(const Value *other) const;
    Value *min(const Value *other) const;
    Value *greaterThan(const Value *other) const;
    Value *lowerThan(const Value *other) const;
    //
    Value *ifFunc(const Value *the, const Value *els) const;
  private:
    ValueDouble(double val);
    static const ValueDouble *checkSameType(const Value *val);
  private:
    double _data;
  };

  class ValueUnit : public Value
  {
  public:
    INTERPKERNEL_EXPORT ValueUnit();
    INTERPKERNEL_EXPORT Value *newInstance() const;
    INTERPKERNEL_EXPORT void setDouble(double val);
    INTERPKERNEL_EXPORT void setVarname(int fastPos, const std::string& var);
    //
    INTERPKERNEL_EXPORT DecompositionInUnitBase getData() const { return _data; }
    INTERPKERNEL_EXPORT void positive();
    INTERPKERNEL_EXPORT void negate();
    INTERPKERNEL_EXPORT void sqrt();
    INTERPKERNEL_EXPORT void cos();
    INTERPKERNEL_EXPORT void sin();
    INTERPKERNEL_EXPORT void tan();
    INTERPKERNEL_EXPORT void acos();
    INTERPKERNEL_EXPORT void asin();
    INTERPKERNEL_EXPORT void atan();
    INTERPKERNEL_EXPORT void cosh();
    INTERPKERNEL_EXPORT void sinh();
    INTERPKERNEL_EXPORT void tanh();
    INTERPKERNEL_EXPORT void abs();
    INTERPKERNEL_EXPORT void exp();
    INTERPKERNEL_EXPORT void ln();
    INTERPKERNEL_EXPORT void log10();
    //
    INTERPKERNEL_EXPORT Value *plus(const Value *other) const;
    INTERPKERNEL_EXPORT Value *minus(const Value *other) const;
    INTERPKERNEL_EXPORT Value *mult(const Value *other) const;
    INTERPKERNEL_EXPORT Value *div(const Value *other) const;
    INTERPKERNEL_EXPORT Value *pow(const Value *other) const;
    INTERPKERNEL_EXPORT Value *max(const Value *other) const;
    INTERPKERNEL_EXPORT Value *min(const Value *other) const;
    INTERPKERNEL_EXPORT Value *greaterThan(const Value *other) const;
    INTERPKERNEL_EXPORT Value *lowerThan(const Value *other) const;
    //
    INTERPKERNEL_EXPORT Value *ifFunc(const Value *the, const Value *els) const;
  private:
    ValueUnit(const DecompositionInUnitBase& unit);
    static void unsupportedOp(const char *type);
    static const ValueUnit *checkSameType(const Value *val);
  private:
    DecompositionInUnitBase _data;
  };

  class INTERPKERNEL_EXPORT ValueDoubleExpr : public Value
  {
  public:
    ValueDoubleExpr(int szDestData, const double *srcData);
    ~ValueDoubleExpr();
    double *getData() const { return _dest_data; }
    Value *newInstance() const;
    void setDouble(double val);
    void setVarname(int fastPos, const std::string& var);
    //
    void positive();
    void negate();
    void sqrt();
    void cos();
    void sin();
    void tan();
    void acos();
    void asin();
    void atan();
    void cosh();
    void sinh();
    void tanh();
    void abs();
    void exp();
    void ln();
    void log10();
    //
    Value *plus(const Value *other) const;
    Value *minus(const Value *other) const;
    Value *mult(const Value *other) const;
    Value *div(const Value *other) const;
    Value *pow(const Value *other) const;
    Value *max(const Value *other) const;
    Value *min(const Value *other) const;
    Value *greaterThan(const Value *other) const;
    Value *lowerThan(const Value *other) const;
    //
    Value *ifFunc(const Value *the, const Value *els) const;
  private:
    int _sz_dest_data;
    double *_dest_data;
    const double *_src_data;
  };
}

#endif
