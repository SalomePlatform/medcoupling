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
    INTERPKERNEL_EXPORT virtual double getDoubleValue() const = 0;
    INTERPKERNEL_EXPORT virtual void fillValue(Value *val) const = 0;
    INTERPKERNEL_EXPORT virtual void compileX86(std::vector<std::string>& ass) const = 0;
    INTERPKERNEL_EXPORT virtual void compileX86_64(std::vector<std::string>& ass) const = 0;
    INTERPKERNEL_EXPORT virtual void replaceValues(const std::vector<double>& valuesInExpr) = 0;
    INTERPKERNEL_EXPORT virtual LeafExpr *deepCopy() const = 0;
    INTERPKERNEL_EXPORT static LeafExpr *buildInstanceFrom(const std::string& expr);
  };

  class LeafExprVal : public LeafExpr
  {
  public:
    INTERPKERNEL_EXPORT LeafExprVal(double value);
    INTERPKERNEL_EXPORT ~LeafExprVal();
    INTERPKERNEL_EXPORT double getDoubleValue() const;
    INTERPKERNEL_EXPORT void compileX86(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void compileX86_64(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void fillValue(Value *val) const;
    INTERPKERNEL_EXPORT void replaceValues(const std::vector<double>& valuesInExpr);
    INTERPKERNEL_EXPORT LeafExprVal *deepCopy() const;
  private:
    double _value;
  };

  class LeafExprVar : public LeafExpr
  {
  public:
    INTERPKERNEL_EXPORT LeafExprVar(const LeafExprVar& other):_fast_pos(other._fast_pos),_ref_pos(other._ref_pos),_var_name(other._var_name),_val(other._val) { }
    INTERPKERNEL_EXPORT LeafExprVar(const std::string& var);
    INTERPKERNEL_EXPORT ~LeafExprVar();
    INTERPKERNEL_EXPORT double getDoubleValue() const;
    INTERPKERNEL_EXPORT void compileX86(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void compileX86_64(std::vector<std::string>& ass) const;
    INTERPKERNEL_EXPORT void fillValue(Value *val) const;
    INTERPKERNEL_EXPORT std::string getVar() const { return _var_name; }
    INTERPKERNEL_EXPORT void prepareExprEvaluation(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo) const;
    INTERPKERNEL_EXPORT void prepareExprEvaluationDouble(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo, int refPos, const double *ptOfInputStart, const double *ptOfInputEnd) const;
    INTERPKERNEL_EXPORT void prepareExprEvaluationVec() const;
    INTERPKERNEL_EXPORT void replaceValues(const std::vector<double>& valuesInExpr);
    INTERPKERNEL_EXPORT static bool isRecognizedKeyVar(const std::string& var, int& pos);
    INTERPKERNEL_EXPORT LeafExprVar *deepCopy() const;
  public:
    static const char END_OF_RECOGNIZED_VAR[];
  private:
    mutable int _fast_pos;
    mutable int _ref_pos;
    std::string _var_name;
    mutable const double *_val;
  };

  class ExprParserOfEval
  {
  public:
    ExprParserOfEval():_leaf(0) { }
    ExprParserOfEval(LeafExpr *leaf, const std::vector<ExprParserOfEval>& subParts, const std::vector<Function *>& funcs):_leaf(leaf),_sub_parts(subParts),_funcs(funcs) { }
    void evaluateDoubleInternal(std::vector<double>& stck) const
    {
      if(_leaf)
        stck.push_back(_leaf->getDoubleValue());
      else
        for(std::vector<ExprParserOfEval>::const_iterator iter=_sub_parts.begin();iter!=_sub_parts.end();iter++)
          (*iter).evaluateDoubleInternal(stck);
      for(std::vector<Function *>::const_iterator iter3=_funcs.begin();iter3!=_funcs.end();iter3++)
        (*iter3)->operateStackOfDouble(stck);
    }
    void evaluateDoubleInternalSafe(std::vector<double>& stck) const
    {
      if(_leaf)
        stck.push_back(_leaf->getDoubleValue());
      else
        for(std::vector<ExprParserOfEval>::const_iterator iter=_sub_parts.begin();iter!=_sub_parts.end();iter++)
          (*iter).evaluateDoubleInternalSafe(stck);
      for(std::vector<Function *>::const_iterator iter3=_funcs.begin();iter3!=_funcs.end();iter3++)
        (*iter3)->operateStackOfDoubleSafe(stck);
    }
    void clearSortedMemory();
    void sortMemory();
  private:
    LeafExpr *_leaf;
    std::vector<ExprParserOfEval> _sub_parts;
    std::vector<Function *> _funcs;
  };

  class ExprParser
  {
  public:
#if __cplusplus >= 201103L
    INTERPKERNEL_EXPORT ExprParser(ExprParser&& other);
    INTERPKERNEL_EXPORT ExprParser& operator=(ExprParser&& other);
#endif
    INTERPKERNEL_EXPORT ExprParser(const std::string& expr, ExprParser *father=0);
    INTERPKERNEL_EXPORT ExprParser(const char *expr, int lgth, ExprParser *father=0);
    INTERPKERNEL_EXPORT ~ExprParser();
    INTERPKERNEL_EXPORT void parse();
    INTERPKERNEL_EXPORT bool isParsingSuccessfull() const { return _is_parsing_ok; }
    INTERPKERNEL_EXPORT double evaluate() const;
    INTERPKERNEL_EXPORT DecompositionInUnitBase evaluateUnit() const;
    INTERPKERNEL_EXPORT void prepareExprEvaluation(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo) const;
    INTERPKERNEL_EXPORT void prepareExprEvaluationDouble(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo, int refPos, const double *ptOfInputStart, const double *ptOfInputEnd) const;
    INTERPKERNEL_EXPORT void prepareFastEvaluator() const;
    INTERPKERNEL_EXPORT void prepareExprEvaluationVec() const;
    INTERPKERNEL_EXPORT double evaluateDouble() const;
    INTERPKERNEL_EXPORT void evaluateDoubleInternal(std::vector<double>& stck) const { _for_eval.evaluateDoubleInternal(stck); }
    INTERPKERNEL_EXPORT void evaluateDoubleInternalSafe(std::vector<double>& stck) const { _for_eval.evaluateDoubleInternalSafe(stck); }
    INTERPKERNEL_EXPORT void checkForEvaluation() const;
    INTERPKERNEL_EXPORT void evaluateExpr(int szOfOutParam, const double *inParam, double *outParam) const;
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
    void reverseThis();
    ExprParserOfEval convertMeTo() const;
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
    mutable ExprParserOfEval _for_eval;
    std::vector<ExprParser> _sub_expr;
    std::vector<Function *> _func_btw_sub_expr;
  private:
    static const int MAX_X86_FP_ST=8;
    static const char WHITE_SPACES[];
    static const char EXPR_PARSE_ERR_MSG[];
  };
}

#endif
