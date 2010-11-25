//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __INTERPKERNELFUNCTION_HXX__
#define __INTERPKERNELFUNCTION_HXX__

#include "INTERPKERNELEXPREVALDefines.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  class Value;
  class Function;

  class INTERPKERNELEXPREVAL_EXPORT FunctionsFactory
  {
  public:
    static Function *buildFuncFromString(const char *type, int nbOfParams) throw(INTERP_KERNEL::Exception);
    static Function *buildUnaryFuncFromString(const char *type) throw(INTERP_KERNEL::Exception);
    //static Function *buildUnaryFuncFromString(char type) throw(INTERP_KERNEL::Exception);
    static Function *buildBinaryFuncFromString(const char *type) throw(INTERP_KERNEL::Exception);
    static Function *buildBinaryFuncFromString(char type) throw(INTERP_KERNEL::Exception);
  };

  class INTERPKERNELEXPREVAL_EXPORT Function
  {
  public:
    virtual ~Function();
    virtual int getNbInputParams() const = 0;
    virtual void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception) = 0;
    virtual const char *getRepr() const = 0;
    virtual bool isACall() const = 0;
  };

  class INTERPKERNELEXPREVAL_EXPORT UnaryFunction : public Function
  { 
  public:
    int getNbInputParams() const;
  };

  class INTERPKERNELEXPREVAL_EXPORT IdentityFunction : public UnaryFunction
  {
  public:
    ~IdentityFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT PositiveFunction : public UnaryFunction
  {
  public:
    ~PositiveFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT NegateFunction : public UnaryFunction
  {
  public:
    ~NegateFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT CosFunction : public UnaryFunction
  {
  public:
    ~CosFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT SinFunction : public UnaryFunction
  {
  public:
    ~SinFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT TanFunction : public UnaryFunction
  {
  public:
    ~TanFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT SqrtFunction : public UnaryFunction
  {
  public:
    ~SqrtFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT AbsFunction : public UnaryFunction
  {
  public:
    ~AbsFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT ExpFunction : public UnaryFunction
  {
  public:
    ~ExpFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT LnFunction : public UnaryFunction
  {
  public:
    ~LnFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT BinaryFunction : public Function
  {
  public:
    int getNbInputParams() const;
  };

  class PlusFunction : public BinaryFunction
  {
  public:
    ~PlusFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT MinusFunction : public BinaryFunction
  {
  public:
    ~MinusFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT MultFunction : public BinaryFunction
  {
  public:
    ~MultFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };
  
  class INTERPKERNELEXPREVAL_EXPORT DivFunction : public BinaryFunction
  {
  public:
    ~DivFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT PowFunction : public BinaryFunction
  {
  public:
    ~PowFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT MaxFunction : public BinaryFunction
  {
  public:
    ~MaxFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };

  class INTERPKERNELEXPREVAL_EXPORT MinFunction : public BinaryFunction
  {
  public:
    ~MinFunction();
    void operate(std::vector<Value *>& stack) const throw(INTERP_KERNEL::Exception);
    void operateX86(std::vector<std::string>& asmb) const throw(INTERP_KERNEL::Exception);
    const char *getRepr() const;
    bool isACall() const;
  public:
    static const char REPR[];
  };
}

#endif
