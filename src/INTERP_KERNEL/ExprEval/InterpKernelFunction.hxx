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

#ifndef __INTERPKERNELFUNCTION_HXX__
#define __INTERPKERNELFUNCTION_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"

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
    int getNbInputParams() const;
  };

  class INTERPKERNEL_EXPORT IdentityFunction : public UnaryFunction
  {
  public:
    ~IdentityFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    IdentityFunction *deepCopy() const { return new IdentityFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT PositiveFunction : public UnaryFunction
  {
  public:
    ~PositiveFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    PositiveFunction *deepCopy() const { return new PositiveFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT NegateFunction : public UnaryFunction
  {
  public:
    ~NegateFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    NegateFunction *deepCopy() const { return new NegateFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT CosFunction : public UnaryFunction
  {
  public:
    ~CosFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    CosFunction *deepCopy() const { return new CosFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT SinFunction : public UnaryFunction
  {
  public:
    ~SinFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    SinFunction *deepCopy() const { return new SinFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT TanFunction : public UnaryFunction
  {
  public:
    ~TanFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    TanFunction *deepCopy() const { return new TanFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT ACosFunction : public UnaryFunction
  {
  public:
    ~ACosFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    ACosFunction *deepCopy() const { return new ACosFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT ASinFunction : public UnaryFunction
  {
  public:
    ~ASinFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    ASinFunction *deepCopy() const { return new ASinFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT ATanFunction : public UnaryFunction
  {
  public:
    ~ATanFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    ATanFunction *deepCopy() const { return new ATanFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT CoshFunction : public UnaryFunction
  {
  public:
    ~CoshFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    CoshFunction *deepCopy() const { return new CoshFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT SinhFunction : public UnaryFunction
  {
  public:
    ~SinhFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    SinhFunction *deepCopy() const { return new SinhFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT TanhFunction : public UnaryFunction
  {
  public:
    ~TanhFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    TanhFunction *deepCopy() const { return new TanhFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT SqrtFunction : public UnaryFunction
  {
  public:
    ~SqrtFunction();
    void operateX86(std::vector<std::string>& asmb) const;
    void operate(std::vector<Value *>& stck) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    SqrtFunction *deepCopy() const { return new SqrtFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT AbsFunction : public UnaryFunction
  {
  public:
    ~AbsFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    AbsFunction *deepCopy() const { return new AbsFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT ExpFunction : public UnaryFunction
  {
  public:
    ~ExpFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    ExpFunction *deepCopy() const { return new ExpFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT LnFunction : public UnaryFunction
  {
  public:
    ~LnFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    LnFunction *deepCopy() const { return new LnFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT LogFunction : public UnaryFunction
  {
  public:
    ~LogFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    LogFunction *deepCopy() const { return new LogFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT Log10Function : public UnaryFunction
  {
  public:
    ~Log10Function();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    Log10Function *deepCopy() const { return new Log10Function; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT BinaryFunction : public Function
  {
  public:
    int getNbInputParams() const;
  };

  class PlusFunction : public BinaryFunction
  {
  public:
    ~PlusFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    PlusFunction *deepCopy() const { return new PlusFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT MinusFunction : public BinaryFunction
  {
  public:
    ~MinusFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    MinusFunction *deepCopy() const { return new MinusFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT MultFunction : public BinaryFunction
  {
  public:
    ~MultFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    MultFunction *deepCopy() const { return new MultFunction; }
  public:
    static const char REPR[];
  };
  
  class INTERPKERNEL_EXPORT DivFunction : public BinaryFunction
  {
  public:
    ~DivFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    DivFunction *deepCopy() const { return new DivFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT PowFunction : public BinaryFunction
  {
  public:
    ~PowFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    PowFunction *deepCopy() const { return new PowFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT MaxFunction : public BinaryFunction
  {
  public:
    ~MaxFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    MaxFunction *deepCopy() const { return new MaxFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT MinFunction : public BinaryFunction
  {
  public:
    ~MinFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    MinFunction *deepCopy() const { return new MinFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT GreaterThanFunction : public BinaryFunction
  {
  public:
    ~GreaterThanFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    GreaterThanFunction *deepCopy() const { return new GreaterThanFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT LowerThanFunction : public BinaryFunction
  {
  public:
    ~LowerThanFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    LowerThanFunction *deepCopy() const { return new LowerThanFunction; }
  public:
    static const char REPR[];
  };

  class INTERPKERNEL_EXPORT TernaryFunction : public Function
  {
  public:
    int getNbInputParams() const;
  };

  class INTERPKERNEL_EXPORT IfFunction : public TernaryFunction
  {
  public:
    ~IfFunction();
    void operate(std::vector<Value *>& stck) const;
    void operateX86(std::vector<std::string>& asmb) const;
    void operateStackOfDouble(std::vector<double>& stck) const;
    void operateStackOfDoubleSafe(std::vector<double>& stck) const;
    const char *getRepr() const;
    bool isACall() const;
    IfFunction *deepCopy() const { return new IfFunction; }
  public:
    static const char REPR[];
  };
}

#endif
