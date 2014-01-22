// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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

#ifndef __INTERPKERNELEXPRPARSER_HXX__
#define __INTERPKERNELEXPRPARSER_HXX__

#include "INTERPKERNELDefines.hxx"
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

  class LeafExpr
  {
  public:
    INTERPKERNEL_EXPORT virtual ~LeafExpr();
    INTERPKERNEL_EXPORT virtual void fillValue(Value *val) const = 0;
    INTERPKERNEL_EXPORT virtual void compileX86(std::vector<std::string>& ass) const = 0;
    INTERPKERNEL_EXPORT virtual void compileX86_64(std::vector<std::string>& ass) const = 0;
    INTERPKERNEL_EXPORT virtual void replaceValues(const std::vector<double>& valuesInExpr) = 0;
    INTERPKERNEL_EXPORT static LeafExpr *buildInstanceFrom(const std::string& expr);
  };

  class LeafExprVal : public LeafExpr
  {
  public:
    INTERPKERNEL_EXPORT LeafExprVal(double value);
    INTERPKERNEL_EXPORT ~LeafExprVal();
    INTERPKERNEL_EXPORT void compileX86(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void compileX86_64(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void fillValue(Value *val) const;
    INTERPKERNEL_EXPORT void replaceValues(const std::vector<double>& valuesInExpr);
  private:
    double _value;
  };

  class LeafExprVar : public LeafExpr
  {
  public:
    INTERPKERNEL_EXPORT LeafExprVar(const std::string& var);
    INTERPKERNEL_EXPORT ~LeafExprVar();
    INTERPKERNEL_EXPORT void compileX86(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void compileX86_64(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void fillValue(Value *val) const;
    INTERPKERNEL_EXPORT std::string getVar() const { return _var_name; }
    INTERPKERNEL_EXPORT void prepareExprEvaluation(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo) const;
    INTERPKERNEL_EXPORT void prepareExprEvaluationVec() const;
    INTERPKERNEL_EXPORT void replaceValues(const std::vector<double>& valuesInExpr);
    INTERPKERNEL_EXPORT static bool isRecognizedKeyVar(const std::string& var, int& pos);
  public:
    static const char END_OF_RECOGNIZED_VAR[];
  private:
    mutable int _fast_pos;
    std::string _var_name;
  };

  class ExprParser
  {
  public:
    INTERPKERNEL_EXPORT ExprParser(const std::string& expr, ExprParser *father=0);
    INTERPKERNEL_EXPORT ExprParser(const char *expr, int lgth, ExprParser *father=0);
    INTERPKERNEL_EXPORT ~ExprParser();
    INTERPKERNEL_EXPORT void parse();
    INTERPKERNEL_EXPORT bool isParsingSuccessfull() const { return _is_parsing_ok; }
    INTERPKERNEL_EXPORT double evaluate() const;
    INTERPKERNEL_EXPORT DecompositionInUnitBase evaluateUnit() const;
    INTERPKERNEL_EXPORT void prepareExprEvaluation(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo) const;
    INTERPKERNEL_EXPORT void evaluateExpr(int szOfOutParam, const double *inParam, double *outParam) const;
    INTERPKERNEL_EXPORT void prepareExprEvaluationVec() const;
    INTERPKERNEL_EXPORT void getSetOfVars(std::set<std::string>& vars) const;
    INTERPKERNEL_EXPORT void getTrueSetOfVars(std::set<std::string>& vars) const;
    //
    INTERPKERNEL_EXPORT char *compileX86() const;
    INTERPKERNEL_EXPORT char *compileX86_64() const;
    INTERPKERNEL_EXPORT void compileX86LowLev(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void compileX86_64LowLev(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT int getStackSizeToPlayX86(const ExprParser *asker) const;
    //
    INTERPKERNEL_EXPORT static std::string buildStringFromFortran(const char *expr, int lgth);
    INTERPKERNEL_EXPORT static std::string deleteWhiteSpaces(const std::string& expr);
  private:
    Value *evaluateLowLev(Value *valGen) const;
  private:
    void prepareExprEvaluationVecLowLev() const;
    bool tryToInterpALeaf();
    void parseUnaryFunc();
    void parseForCmp();
    void parseForAddMin();
    void parseForMulDiv();
    void parseForPow();
    void parseDeeper();
    bool simplify();
    void releaseFunctions();
    void checkBracketsParity() const;
    void fillValuesInExpr(std::vector<double>& valuesInExpr);
    void replaceValues(const std::vector<double>& valuesInExpr);
    static double ReplaceAndTraduce(std::string& expr, int id, std::size_t bg, std::size_t end, int& delta);
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
