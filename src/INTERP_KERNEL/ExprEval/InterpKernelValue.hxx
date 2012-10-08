// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
    virtual void setDouble(double val) throw(INTERP_KERNEL::Exception) = 0;
    virtual void setVarname(int fastPos, const std::string& var) throw(INTERP_KERNEL::Exception) = 0;
    //unary
    virtual void positive() throw(INTERP_KERNEL::Exception) = 0;
    virtual void negate() throw(INTERP_KERNEL::Exception) = 0;
    virtual void sqrt() throw(INTERP_KERNEL::Exception) = 0;
    virtual void cos() throw(INTERP_KERNEL::Exception) = 0;
    virtual void sin() throw(INTERP_KERNEL::Exception) = 0;
    virtual void tan() throw(INTERP_KERNEL::Exception) = 0;
    virtual void abs() throw(INTERP_KERNEL::Exception) = 0;
    virtual void exp() throw(INTERP_KERNEL::Exception) = 0;
    virtual void ln() throw(INTERP_KERNEL::Exception) = 0;
    virtual void log10() throw(INTERP_KERNEL::Exception) = 0;
    //binary
    virtual Value *plus(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *minus(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *mult(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *div(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *pow(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *max(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *min(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *greaterThan(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *lowerThan(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    //ternary
    virtual Value *ifFunc(const Value *the, const Value *els) const throw(INTERP_KERNEL::Exception) = 0;
  };

  class INTERPKERNEL_EXPORT ValueDouble : public Value
  {
  public:
    ValueDouble();
    Value *newInstance() const;
    void setDouble(double val) throw(INTERP_KERNEL::Exception);
    void setVarname(int fastPos, const std::string& var) throw(INTERP_KERNEL::Exception);
    //
    double getData() const { return _data; }
    void positive() throw(INTERP_KERNEL::Exception);
    void negate() throw(INTERP_KERNEL::Exception);
    void sqrt() throw(INTERP_KERNEL::Exception);
    void cos() throw(INTERP_KERNEL::Exception);
    void sin() throw(INTERP_KERNEL::Exception);
    void tan() throw(INTERP_KERNEL::Exception);
    void abs() throw(INTERP_KERNEL::Exception);
    void exp() throw(INTERP_KERNEL::Exception);
    void ln() throw(INTERP_KERNEL::Exception);
    void log10() throw(INTERP_KERNEL::Exception);
    //
    Value *plus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *minus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *mult(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *div(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *pow(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *max(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *min(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *greaterThan(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *lowerThan(const Value *other) const throw(INTERP_KERNEL::Exception);
    //
    Value *ifFunc(const Value *the, const Value *els) const throw(INTERP_KERNEL::Exception);
  private:
    ValueDouble(double val);
    static const ValueDouble *checkSameType(const Value *val) throw(INTERP_KERNEL::Exception);
  private:
    double _data;
  };

  class INTERPKERNEL_EXPORT ValueUnit : public Value
  {
  public:
    ValueUnit();
    Value *newInstance() const;
    void setDouble(double val) throw(INTERP_KERNEL::Exception);
    void setVarname(int fastPos, const std::string& var) throw(INTERP_KERNEL::Exception);
    //
    DecompositionInUnitBase getData() const { return _data; }
    void positive() throw(INTERP_KERNEL::Exception);
    void negate() throw(INTERP_KERNEL::Exception);
    void sqrt() throw(INTERP_KERNEL::Exception);
    void cos() throw(INTERP_KERNEL::Exception);
    void sin() throw(INTERP_KERNEL::Exception);
    void tan() throw(INTERP_KERNEL::Exception);
    void abs() throw(INTERP_KERNEL::Exception);
    void exp() throw(INTERP_KERNEL::Exception);
    void ln() throw(INTERP_KERNEL::Exception);
    void log10() throw(INTERP_KERNEL::Exception);
    //
    Value *plus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *minus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *mult(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *div(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *pow(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *max(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *min(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *greaterThan(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *lowerThan(const Value *other) const throw(INTERP_KERNEL::Exception);
    //
    Value *ifFunc(const Value *the, const Value *els) const throw(INTERP_KERNEL::Exception);
  private:
    ValueUnit(const DecompositionInUnitBase& unit);
    static void unsupportedOp(const char *type) throw(INTERP_KERNEL::Exception);
    static const ValueUnit *checkSameType(const Value *val) throw(INTERP_KERNEL::Exception);
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
    void setDouble(double val) throw(INTERP_KERNEL::Exception);
    void setVarname(int fastPos, const std::string& var) throw(INTERP_KERNEL::Exception);
    //
    void positive() throw(INTERP_KERNEL::Exception);
    void negate() throw(INTERP_KERNEL::Exception);
    void sqrt() throw(INTERP_KERNEL::Exception);
    void cos() throw(INTERP_KERNEL::Exception);
    void sin() throw(INTERP_KERNEL::Exception);
    void tan() throw(INTERP_KERNEL::Exception);
    void abs() throw(INTERP_KERNEL::Exception);
    void exp() throw(INTERP_KERNEL::Exception);
    void ln() throw(INTERP_KERNEL::Exception);
    void log10() throw(INTERP_KERNEL::Exception);
    //
    Value *plus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *minus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *mult(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *div(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *pow(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *max(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *min(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *greaterThan(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *lowerThan(const Value *other) const throw(INTERP_KERNEL::Exception);
    //
    Value *ifFunc(const Value *the, const Value *els) const throw(INTERP_KERNEL::Exception);
  private:
    int _sz_dest_data;
    double *_dest_data;
    const double *_src_data;
  };
}

#endif
