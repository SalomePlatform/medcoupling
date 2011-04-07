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

#ifndef __INTERPKERNELEXPRPARSER_HXX__
#define __INTERPKERNELEXPRPARSER_HXX__

#include "INTERPKERNELEXPREVALDefines.hxx"
#include "InterpKernelUnit.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelFunction.hxx"

#include <string>
#include <list>
#include <map>
#include <set>

namespace INTERP_KERNEL
{
  class ValueDouble;

  class INTERPKERNELEXPREVAL_EXPORT LeafExpr
  {
  public:
    virtual ~LeafExpr();
    virtual void fillValue(Value *val) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void compileX86(std::vector<std::string>& ass) const = 0;
    virtual void compileX86_64(std::vector<std::string>& ass) const = 0;
    virtual void replaceValues(const std::vector<double>& valuesInExpr) throw(INTERP_KERNEL::Exception) = 0;
    static LeafExpr *buildInstanceFrom(const std::string& expr) throw(INTERP_KERNEL::Exception);
  };

  class INTERPKERNELEXPREVAL_EXPORT LeafExprVal : public LeafExpr
  {
  public:
    LeafExprVal(double value);
    ~LeafExprVal();
    void compileX86(std::vector<std::string>& ass) const;
    void compileX86_64(std::vector<std::string>& ass) const;
    void fillValue(Value *val) const throw(INTERP_KERNEL::Exception);
    void replaceValues(const std::vector<double>& valuesInExpr) throw(INTERP_KERNEL::Exception);
  private:
    double _value;
  };

  class INTERPKERNELEXPREVAL_EXPORT LeafExprVar : public LeafExpr
  {
  public:
    LeafExprVar(const std::string& var);
    ~LeafExprVar();
    void compileX86(std::vector<std::string>& ass) const;
    void compileX86_64(std::vector<std::string>& ass) const;
    void fillValue(Value *val) const throw(INTERP_KERNEL::Exception);
    std::string getVar() const { return _var_name; }
    void prepareExprEvaluation(const std::vector<std::string>& vars) const throw(INTERP_KERNEL::Exception);
    void prepareExprEvaluationVec() const throw(INTERP_KERNEL::Exception);
    void replaceValues(const std::vector<double>& valuesInExpr) throw(INTERP_KERNEL::Exception);
    static bool isRecognizedKeyVar(const std::string& var, int& pos);
  public:
    static const char END_OF_RECOGNIZED_VAR[];
  private:
    mutable int _fast_pos;
    std::string _var_name;
  };

  class INTERPKERNELEXPREVAL_EXPORT ExprParser
  {
  public:
    ExprParser(const char *expr, ExprParser *father=0);
    ExprParser(const char *expr, int lgth, ExprParser *father=0);
    ~ExprParser();
    void parse() throw(INTERP_KERNEL::Exception);
    bool isParsingSuccessfull() const { return _is_parsing_ok; }
    double evaluate() const throw(INTERP_KERNEL::Exception);
    DecompositionInUnitBase evaluateUnit() const throw(INTERP_KERNEL::Exception);
    void prepareExprEvaluation(const std::vector<std::string>& vars) const throw(INTERP_KERNEL::Exception);
    void evaluateExpr(int szOfOutParam, const double *inParam, double *outParam) const throw(INTERP_KERNEL::Exception);
    void prepareExprEvaluationVec() const throw(INTERP_KERNEL::Exception);
    void getSetOfVars(std::set<std::string>& vars) const;
    void getTrueSetOfVars(std::set<std::string>& vars) const;
    //
    char *compileX86() const;
    char *compileX86_64() const;
    void compileX86LowLev(std::vector<std::string>& ass) const;
    void compileX86_64LowLev(std::vector<std::string>& ass) const;
    int getStackSizeToPlayX86(const ExprParser *asker) const;
    //
    static std::string buildStringFromFortran(const char *expr, int lgth);
    static std::string deleteWhiteSpaces(const std::string& expr);
  private:
    Value *evaluateLowLev(Value *valGen) const throw(INTERP_KERNEL::Exception);
  private:
    void prepareExprEvaluationVecLowLev() const throw(INTERP_KERNEL::Exception);
    bool tryToInterpALeaf() throw(INTERP_KERNEL::Exception);
    void parseUnaryFunc() throw(INTERP_KERNEL::Exception);
    void parseForCmp() throw(INTERP_KERNEL::Exception);
    void parseForAddMin() throw(INTERP_KERNEL::Exception);
    void parseForMulDiv() throw(INTERP_KERNEL::Exception);
    void parseForPow() throw(INTERP_KERNEL::Exception);
    void parseDeeper() throw(INTERP_KERNEL::Exception);
    bool simplify() throw(INTERP_KERNEL::Exception);
    void releaseFunctions();
    void checkBracketsParity() const throw(INTERP_KERNEL::Exception);
    void fillValuesInExpr(std::vector<double>& valuesInExpr) throw(INTERP_KERNEL::Exception);
    void replaceValues(const std::vector<double>& valuesInExpr) throw(INTERP_KERNEL::Exception);
    static double ReplaceAndTraduce(std::string& expr, int id, std::size_t bg, std::size_t end, int& delta) throw(INTERP_KERNEL::Exception);
    static std::size_t FindCorrespondingOpenBracket(const std::string& expr, std::size_t posOfCloseBracket);
    static void LocateError(std::ostream& stringToDisp, const std::string& srcOfErr, int posOfErr);
  private:
    ExprParser *_father;
    bool _is_parsed;
    LeafExpr *_leaf;
    bool _is_parsing_ok;
    std::string _expr;
    std::list<ExprParser> _sub_expr;
    std::list<Function *> _func_btw_sub_expr;
  private:
    static const int MAX_X86_FP_ST=8;
    static const char WHITE_SPACES[];
    static const char EXPR_PARSE_ERR_MSG[];
  };
}

#endif
