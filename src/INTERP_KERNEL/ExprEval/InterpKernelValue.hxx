#ifndef __INTERPKERNELVALUE_HXX__
#define __INTERPKERNELVALUE_HXX__

#include "INTERPKERNELEXPREVALDefines.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelUnit.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNELEXPREVAL_EXPORT Value
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
    //binary
    virtual Value *plus(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *minus(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *mult(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *div(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *pow(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *max(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual Value *min(const Value *other) const throw(INTERP_KERNEL::Exception) = 0;
  };

  class INTERPKERNELEXPREVAL_EXPORT ValueDouble : public Value
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
    //
    Value *plus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *minus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *mult(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *div(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *pow(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *max(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *min(const Value *other) const throw(INTERP_KERNEL::Exception);
  private:
    ValueDouble(double val);
    static const ValueDouble *checkSameType(const Value *val) throw(INTERP_KERNEL::Exception);
  private:
    double _data;
  };

  class INTERPKERNELEXPREVAL_EXPORT ValueUnit : public Value
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
    //
    Value *plus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *minus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *mult(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *div(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *pow(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *max(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *min(const Value *other) const throw(INTERP_KERNEL::Exception);
  private:
    ValueUnit(const DecompositionInUnitBase& unit);
    static void unsupportedOp(const char *type) throw(INTERP_KERNEL::Exception);
    static const ValueUnit *checkSameType(const Value *val) throw(INTERP_KERNEL::Exception);
  private:
    DecompositionInUnitBase _data;
  };

  class INTERPKERNELEXPREVAL_EXPORT ValueDoubleExpr : public Value
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
    //
    Value *plus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *minus(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *mult(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *div(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *pow(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *max(const Value *other) const throw(INTERP_KERNEL::Exception);
    Value *min(const Value *other) const throw(INTERP_KERNEL::Exception);
  private:
    int _sz_dest_data;
    double *_dest_data;
    const double *_src_data;
  };
}

#endif
