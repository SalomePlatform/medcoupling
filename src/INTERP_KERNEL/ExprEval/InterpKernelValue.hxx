// Copyright (C) 2007-2024  CEA, EDF
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
#include "InterpKernelUnit.hxx"
#include <string>

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT Value
  {
  public:
    virtual Value *newInstance() const = 0;
    virtual ~Value() = default;
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
    Value *newInstance() const override;
    void setDouble(double val) override;
    void setVarname(int fastPos, const std::string& var) override;
    //
    double getData() const { return _data; }
    void positive() override;
    void negate() override;
    void sqrt() override;
    void cos() override;
    void sin() override;
    void tan() override;
    void acos() override;
    void asin() override;
    void atan() override;
    void cosh() override;
    void sinh() override;
    void tanh() override;
    void abs() override;
    void exp() override;
    void ln() override;
    void log10() override;
    //
    Value *plus(const Value *other) const override;
    Value *minus(const Value *other) const override;
    Value *mult(const Value *other) const override;
    Value *div(const Value *other) const override;
    Value *pow(const Value *other) const override;
    Value *max(const Value *other) const override;
    Value *min(const Value *other) const override;
    Value *greaterThan(const Value *other) const override;
    Value *lowerThan(const Value *other) const override;
    //
    Value *ifFunc(const Value *the, const Value *els) const override;
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
    INTERPKERNEL_EXPORT Value *newInstance() const override;
    INTERPKERNEL_EXPORT void setDouble(double val) override;
    INTERPKERNEL_EXPORT void setVarname(int fastPos, const std::string& var) override;
    //
    INTERPKERNEL_EXPORT DecompositionInUnitBase getData() const { return _data; }
    INTERPKERNEL_EXPORT void positive() override;
    INTERPKERNEL_EXPORT void negate() override;
    INTERPKERNEL_EXPORT void sqrt() override;
    INTERPKERNEL_EXPORT void cos() override;
    INTERPKERNEL_EXPORT void sin() override;
    INTERPKERNEL_EXPORT void tan() override;
    INTERPKERNEL_EXPORT void acos() override;
    INTERPKERNEL_EXPORT void asin() override;
    INTERPKERNEL_EXPORT void atan() override;
    INTERPKERNEL_EXPORT void cosh() override;
    INTERPKERNEL_EXPORT void sinh() override;
    INTERPKERNEL_EXPORT void tanh() override;
    INTERPKERNEL_EXPORT void abs() override;
    INTERPKERNEL_EXPORT void exp() override;
    INTERPKERNEL_EXPORT void ln() override;
    INTERPKERNEL_EXPORT void log10() override;
    //
    INTERPKERNEL_EXPORT Value *plus(const Value *other) const override;
    INTERPKERNEL_EXPORT Value *minus(const Value *other) const override;
    INTERPKERNEL_EXPORT Value *mult(const Value *other) const override;
    INTERPKERNEL_EXPORT Value *div(const Value *other) const override;
    INTERPKERNEL_EXPORT Value *pow(const Value *other) const override;
    INTERPKERNEL_EXPORT Value *max(const Value *other) const override;
    INTERPKERNEL_EXPORT Value *min(const Value *other) const override;
    INTERPKERNEL_EXPORT Value *greaterThan(const Value *other) const override;
    INTERPKERNEL_EXPORT Value *lowerThan(const Value *other) const override;
    //
    INTERPKERNEL_EXPORT Value *ifFunc(const Value *the, const Value *els) const override;
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
    ~ValueDoubleExpr() override;
    double *getData() const { return _dest_data; }
    Value *newInstance() const override;
    void setDouble(double val) override;
    void setVarname(int fastPos, const std::string& var) override;
    //
    void positive() override;
    void negate() override;
    void sqrt() override;
    void cos() override;
    void sin() override;
    void tan() override;
    void acos() override;
    void asin() override;
    void atan() override;
    void cosh() override;
    void sinh() override;
    void tanh() override;
    void abs() override;
    void exp() override;
    void ln() override;
    void log10() override;
    //
    Value *plus(const Value *other) const override;
    Value *minus(const Value *other) const override;
    Value *mult(const Value *other) const override;
    Value *div(const Value *other) const override;
    Value *pow(const Value *other) const override;
    Value *max(const Value *other) const override;
    Value *min(const Value *other) const override;
    Value *greaterThan(const Value *other) const override;
    Value *lowerThan(const Value *other) const override;
    //
    Value *ifFunc(const Value *the, const Value *els) const override;
  private:
    int _sz_dest_data;
    double *_dest_data;
    const double *_src_data;
  };
}

#endif
