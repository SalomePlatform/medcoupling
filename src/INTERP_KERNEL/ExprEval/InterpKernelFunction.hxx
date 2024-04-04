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

#ifndef __INTERPKERNELFUNCTION_HXX__
#define __INTERPKERNELFUNCTION_HXX__

#include "INTERPKERNELDefines.hxx"

#include <string>
#include <vector>

namespace INTERP_KERNEL
{
  class Value;
  class Function;

  class INTERPKERNEL_EXPORT FunctionsFactory
  {
  public:
    static Function *buildFuncFromString(const char *type, int nbOfParams);
    static Function *buildUnaryFuncFromString(const char *type);
    //static Function *buildUnaryFuncFromString(char type);
    static Function *buildBinaryFuncFromString(const char *type);
    static Function *buildBinaryFuncFromString(char type);
    static Function *buildTernaryFuncFromString(const char *type);
  };

  class INTERPKERNEL_EXPORT Function
  {
  public:
    virtual ~Function();
    virtual int getNbInputParams() const = 0;
    virtual void operate(std::vector<Value *>& stck) const = 0;
    virtual void operateX86(std::vector<std::string>& asmb) const = 0;
    virtual void operateStackOfDouble(std::vector<double>& stck) const = 0;
    virtual void operateStackOfDoubleSafe(std::vector<double>& stck) const { operateStackOfDouble(stck); }
    virtual const char *getRepr() const = 0;
    virtual bool isACall() const = 0;
    virtual Function *deepCopy() const = 0;
  };

  class INTERPKERNEL_EXPORT UnaryFunction : public Function
  { 
  public:
    int getNbInputParams() const override;
  };

  class INTERPKERNEL_EXPORT IdentityFunction : public UnaryFunction
  {
  public:
    ~IdentityFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    IdentityFunction *deepCopy() const override { return new IdentityFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT PositiveFunction : public UnaryFunction
  {
  public:
    ~PositiveFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    PositiveFunction *deepCopy() const override { return new PositiveFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT NegateFunction : public UnaryFunction
  {
  public:
    ~NegateFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    NegateFunction *deepCopy() const override { return new NegateFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT CosFunction : public UnaryFunction
  {
  public:
    ~CosFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    CosFunction *deepCopy() const override { return new CosFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT SinFunction : public UnaryFunction
  {
  public:
    ~SinFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    SinFunction *deepCopy() const override { return new SinFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT TanFunction : public UnaryFunction
  {
  public:
    ~TanFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    TanFunction *deepCopy() const override { return new TanFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT ACosFunction : public UnaryFunction
  {
  public:
    ~ACosFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    ACosFunction *deepCopy() const override { return new ACosFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT ASinFunction : public UnaryFunction
  {
  public:
    ~ASinFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    ASinFunction *deepCopy() const override { return new ASinFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT ATanFunction : public UnaryFunction
  {
  public:
    ~ATanFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    ATanFunction *deepCopy() const override { return new ATanFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT CoshFunction : public UnaryFunction
  {
  public:
    ~CoshFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    CoshFunction *deepCopy() const override { return new CoshFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT SinhFunction : public UnaryFunction
  {
  public:
    ~SinhFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    SinhFunction *deepCopy() const override { return new SinhFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT TanhFunction : public UnaryFunction
  {
  public:
    ~TanhFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    TanhFunction *deepCopy() const override { return new TanhFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT SqrtFunction : public UnaryFunction
  {
  public:
    ~SqrtFunction() override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operate(std::vector<Value *>& stck) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    SqrtFunction *deepCopy() const override { return new SqrtFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT AbsFunction : public UnaryFunction
  {
  public:
    ~AbsFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    AbsFunction *deepCopy() const override { return new AbsFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT ExpFunction : public UnaryFunction
  {
  public:
    ~ExpFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    ExpFunction *deepCopy() const override { return new ExpFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT LnFunction : public UnaryFunction
  {
  public:
    ~LnFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    LnFunction *deepCopy() const override { return new LnFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT LogFunction : public UnaryFunction
  {
  public:
    ~LogFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    LogFunction *deepCopy() const override { return new LogFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT Log10Function : public UnaryFunction
  {
  public:
    ~Log10Function() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    Log10Function *deepCopy() const override { return new Log10Function; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT BinaryFunction : public Function
  {
  public:
    int getNbInputParams() const override;
  };

  class PlusFunction : public BinaryFunction
  {
  public:
    ~PlusFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    PlusFunction *deepCopy() const override { return new PlusFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT MinusFunction : public BinaryFunction
  {
  public:
    ~MinusFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    MinusFunction *deepCopy() const override { return new MinusFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT MultFunction : public BinaryFunction
  {
  public:
    ~MultFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    MultFunction *deepCopy() const override { return new MultFunction; }
  public:
    static const char REPR[];
  };
  
  class INTERPKERNEL_EXPORT DivFunction : public BinaryFunction
  {
  public:
    ~DivFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    DivFunction *deepCopy() const override { return new DivFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT PowFunction : public BinaryFunction
  {
  public:
    ~PowFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    PowFunction *deepCopy() const override { return new PowFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT MaxFunction : public BinaryFunction
  {
  public:
    ~MaxFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    MaxFunction *deepCopy() const override { return new MaxFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT MinFunction : public BinaryFunction
  {
  public:
    ~MinFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    MinFunction *deepCopy() const override { return new MinFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT GreaterThanFunction : public BinaryFunction
  {
  public:
    ~GreaterThanFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    GreaterThanFunction *deepCopy() const override { return new GreaterThanFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT LowerThanFunction : public BinaryFunction
  {
  public:
    ~LowerThanFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    LowerThanFunction *deepCopy() const override { return new LowerThanFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT TernaryFunction : public Function
  {
  public:
    int getNbInputParams() const override;
  };

  class INTERPKERNEL_EXPORT IfFunction : public TernaryFunction
  {
  public:
    ~IfFunction() override;
    void operate(std::vector<Value *>& stck) const override;
    void operateX86(std::vector<std::string>& asmb) const override;
    void operateStackOfDouble(std::vector<double>& stck) const override;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const override;
    const char *getRepr() const override;
    bool isACall() const override;
    IfFunction *deepCopy() const override { return new IfFunction; }
  public:
    static const char REPR[];
  };
}

#endif
