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

#include "InterpKernelExprParser.hxx"
#include "InterpKernelValue.hxx"
#include "InterpKernelAsmX86.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <cctype>
#include <sstream>
#include <limits>
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>

using namespace INTERP_KERNEL;

const char LeafExprVar::END_OF_RECOGNIZED_VAR[]="Vec";

const char ExprParser::WHITE_SPACES[]=" \n";

const char ExprParser::EXPR_PARSE_ERR_MSG[]="Invalid expression detected : ";

LeafExpr *LeafExpr::buildInstanceFrom(const std::string& expr)
{
  std::istringstream stream;
  stream.str(expr);
  double val;
  stream >> val;
  if(!stream.fail())
    if(stream.eof())
      return new LeafExprVal(val);
    else
      {
        std::ostringstream errMsg;
        char MSGTYP6[]="Error following expression is not consedered as a double value : ";
        errMsg << MSGTYP6 << expr;
        throw INTERP_KERNEL::Exception(errMsg.str().c_str());
      }
  else
    return new LeafExprVar(expr);
}

LeafExpr::~LeafExpr()
{
}

LeafExprVal::LeafExprVal(double value):_value(value)
{
}

LeafExprVal::~LeafExprVal()
{
}

void LeafExprVal::fillValue(Value *val) const
{
  val->setDouble(_value);
}

void LeafExprVal::replaceValues(const std::vector<double>& valuesInExpr)
{
  int pos=(int)_value;
  int lgth=(int)valuesInExpr.size();
  if(pos>=lgth || pos<0)
    throw INTERP_KERNEL::Exception("LeafExprVal::replaceValues : Big Problem detected! Send a mail to Salome support with expression.");
  _value=valuesInExpr[pos];
}

LeafExprVal *LeafExprVal::deepCopy() const
{
  return new LeafExprVal(*this);
}

LeafExprVar::LeafExprVar(const std::string& var):_fast_pos(-1),_var_name(var),_val(0)
{
}

void LeafExprVar::fillValue(Value *val) const
{
  if(_val)
    val->setDouble(_val[_fast_pos]);
  else
    val->setVarname(_fast_pos,_var_name);

}

void LeafExprVar::prepareExprEvaluation(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo) const
{
  std::vector<std::string>::const_iterator iter=std::find(vars.begin(),vars.end(),_var_name);
  if(iter==vars.end())
    {
      if(!isRecognizedKeyVar(_var_name,_fast_pos))
        {
          std::ostringstream oss; oss << "Var : " << _var_name << " not in : ";
          std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss,", "));
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      else
        {
          int relPos=-7-_fast_pos;
          if(relPos>=targetNbOfCompo)
            {
              std::ostringstream oss; oss << "LeafExprVar::prepareExprEvaluation : Found recognized unitary vector \"" << _var_name << "\" which implies that component #" << relPos;
              oss << " exists, but it is not the case component id should be in [0," << targetNbOfCompo << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          else
            return;
        }
    }
  _fast_pos=(int)std::distance(vars.begin(),iter);
  if(_fast_pos>=nbOfCompo)
    {
      std::ostringstream oss; oss << "LeafExprVar::prepareExprEvaluation : Found var \"" << _var_name << "\" on place " << _fast_pos << " whereas only must be in [0," << nbOfCompo << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * \param [in] vars - the sorted list of vars
 * \param [in] nbOfCompo - the size of the input tuples (it is used to scan if no problem occurs)
 * \param [in] targetNbOfCompo - the size of the output tuple (it is used to check that no problem occurs)
 * \param [in] refPos - is an integer in [0,targetNbOfCompo), that tell the id of \a this. It is for multi interpreters.
 * \sa evaluateDouble
 */
void LeafExprVar::prepareExprEvaluationDouble(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo, int refPos, const double *ptOfInputStart, const double *ptOfInputEnd) const
{
  if((int)vars.size()!=std::distance(ptOfInputStart,ptOfInputEnd))
    throw INTERP_KERNEL::Exception("LeafExprVar::prepareExprEvaluationDouble : size of input vector must be equal to the input vector !");
  prepareExprEvaluation(vars,nbOfCompo,targetNbOfCompo);
  _ref_pos=refPos;
  _val=ptOfInputStart;
}

void LeafExprVar::prepareExprEvaluationVec() const
{
  if(!isRecognizedKeyVar(_var_name,_fast_pos))
    _fast_pos=-2;
}

bool LeafExprVar::isRecognizedKeyVar(const std::string& var, int& pos)
{
  if(var.length()!=sizeof(END_OF_RECOGNIZED_VAR))
    return false;
  std::string end=var.substr(1);
  if(end!=END_OF_RECOGNIZED_VAR)
    return false;
  char first=var[0];
  if(first<'I' || first>'Z')
    return false;
  pos=-7-(first-'I');
  return true;
}

LeafExprVar *LeafExprVar::deepCopy() const
{
  return new LeafExprVar(*this);
}

/*!
 * Nothing to do it is not a bug.
 */
void LeafExprVar::replaceValues(const std::vector<double>& valuesInExpr)
{
}

LeafExprVar::~LeafExprVar()
{
}

void ExprParserOfEval::clearSortedMemory()
{
  delete _leaf;
  for(std::vector<ExprParserOfEval>::iterator it=_sub_parts.begin();it!=_sub_parts.end();it++)
    (*it).clearSortedMemory();
  for(std::vector<Function *>::iterator it=_funcs.begin();it!=_funcs.end();it++)
    delete *it;
}

void ExprParserOfEval::sortMemory()
{
  for(std::vector<ExprParserOfEval>::iterator it=_sub_parts.begin();it!=_sub_parts.end();it++)
    (*it).sortMemory();
  if(_leaf)
    _leaf=_leaf->deepCopy();
  for(std::vector<Function *>::iterator it=_funcs.begin();it!=_funcs.end();it++)
    if(*it)
      *it=(*it)->deepCopy();
}

ExprParser::ExprParser(const std::string& expr, ExprParser *father):_father(father),_is_parsed(false),_leaf(0),_is_parsing_ok(false),_expr(expr)
{
  _expr=deleteWhiteSpaces(_expr);
}

//! For \b NOT null terminated strings coming from FORTRAN.
ExprParser::ExprParser(const char *expr, int lgth, ExprParser *father):_father(father),_is_parsed(false),_leaf(0),_is_parsing_ok(false)
{
  _expr=buildStringFromFortran(expr,lgth);
  _expr=deleteWhiteSpaces(_expr);
}

ExprParser::~ExprParser()
{
  delete _leaf;
  _for_eval.clearSortedMemory();
  releaseFunctions();
}

std::size_t ExprParser::FindCorrespondingOpenBracket(const std::string& expr, std::size_t posOfCloseBracket)
{
  int level=0;
  for(std::size_t iter=0;iter<posOfCloseBracket;iter++)
    {
      std::size_t iter2=posOfCloseBracket-1-iter;
      if(expr[iter2]==')')
        level++;
      else if(expr[iter2]=='(')
        {
          if(level==0)
            return iter2;
          else
            level--;
        }
    }
  return std::string::npos;
}

std::string ExprParser::buildStringFromFortran(const char *expr, int lgth)
{
  std::string ret(expr,lgth);
  std::string whiteSpaces(WHITE_SPACES);
  std::size_t found=ret.find_last_not_of(whiteSpaces);
  if (found!=std::string::npos)
    ret.erase(found+1);
  else
    ret.clear();//ret is all whitespace
  return ret;
}

std::string ExprParser::deleteWhiteSpaces(const std::string& expr)
{
  std::string ret(expr);
  std::string whiteSpaces(WHITE_SPACES);
  std::size_t where1=0,where2=0;
  while(where2!=std::string::npos && where1!=std::string::npos)
    {
      where1=ret.find_first_of(whiteSpaces.c_str(),where1,whiteSpaces.length());
      if(where1!=std::string::npos)
        {
          where2=ret.find_first_not_of(whiteSpaces,where1);
          if(where2!=std::string::npos)
            ret.erase(ret.begin()+where1,ret.begin()+where2);
          else
            ret.erase(ret.begin()+where1,ret.end());
        }
    }
  return ret;
}

void ExprParser::parse()
{
  _is_parsed=true;
  _is_parsing_ok=false;
  _sub_expr.clear();
  releaseFunctions();
  if(!_expr.empty())
    {
      std::string tmp(_expr);
      std::vector<double> valuesInExpr;
      fillValuesInExpr(valuesInExpr);
      checkBracketsParity();
      if(!simplify())
        parseDeeper();
      replaceValues(valuesInExpr);
      _expr=tmp;
    }
  reverseThis();
  _is_parsing_ok=true;
}

double ExprParser::evaluate() const
{
  AutoCppPtr<Value> gen(new ValueDouble);
  AutoCppPtr<ValueDouble> res(static_cast<ValueDouble *>(evaluateLowLev(gen)));
  return res->getData();
}

DecompositionInUnitBase ExprParser::evaluateUnit() const
{
  Value *gen=new ValueUnit;
  ValueUnit *res=0;
  try
    {
      res=(ValueUnit *)evaluateLowLev(gen);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete gen;
      throw e;
    }
  delete gen;
  DecompositionInUnitBase ret=res->getData();
  delete res;
  return ret;
}

void ExprParser::evaluateExpr(int szOfOutParam, const double *inParam, double *outParam) const
{
  AutoCppPtr<Value> gen(new ValueDoubleExpr(szOfOutParam,inParam));
  AutoCppPtr<ValueDoubleExpr> res(static_cast<ValueDoubleExpr *>(evaluateLowLev(gen)));
  std::copy(res->getData(),res->getData()+szOfOutParam,outParam);
}

void ExprParser::prepareExprEvaluation(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo) const
{
  if(_leaf)
    {
      LeafExprVar *leafC=dynamic_cast<LeafExprVar *>(_leaf);
      if(leafC)
        leafC->prepareExprEvaluation(vars,nbOfCompo,targetNbOfCompo);
    }
  else
    for(std::vector<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
      (*iter).prepareExprEvaluation(vars,nbOfCompo,targetNbOfCompo);
}

/*!
 * \param [in] vars - the sorted list of vars
 * \param [in] nbOfCompo - the size of the input tuples (it is used to scan if no problem occurs)
 * \param [in] targetNbOfCompo - the size of the output tuple (it is used to check that no problem occurs)
 * \param [in] refPos - is an integer in [0,targetNbOfCompo), that tell the id of \a this. It is for multi interpreters.
 * \sa evaluateDouble
 */
void ExprParser::prepareExprEvaluationDouble(const std::vector<std::string>& vars, int nbOfCompo, int targetNbOfCompo, int refPos, const double *ptOfInputStart, const double *ptOfInputEnd) const
{
  if((int)vars.size()!=std::distance(ptOfInputStart,ptOfInputEnd))
    throw INTERP_KERNEL::Exception("ExprParser::prepareExprEvaluationDouble : size of input vector must be equal to the input vector !");
  if(_leaf)
    {
      LeafExprVar *leafC=dynamic_cast<LeafExprVar *>(_leaf);
      if(leafC)
        leafC->prepareExprEvaluationDouble(vars,nbOfCompo,targetNbOfCompo,refPos,ptOfInputStart,ptOfInputEnd);
    }
  else
    for(std::vector<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
      (*iter).prepareExprEvaluationDouble(vars,nbOfCompo,targetNbOfCompo,refPos,ptOfInputStart,ptOfInputEnd);
}

void ExprParser::prepareFastEvaluator() const
{
  _for_eval.clearSortedMemory();
  _for_eval=convertMeTo();
  _for_eval.sortMemory();
}

/*!
 * \sa prepareExprEvaluationDouble
 */
double ExprParser::evaluateDouble() const
{
  checkForEvaluation();
  std::vector<double> stackOfVal;
  evaluateDoubleInternal(stackOfVal);
  return stackOfVal.back();
}

void ExprParser::checkForEvaluation() const
{
  if(!_is_parsing_ok)
    throw INTERP_KERNEL::Exception("checkForEvaluation : Parsing fails ! Invalid expression !");
  if(_sub_expr.empty() && !_leaf)
    throw INTERP_KERNEL::Exception("checkForEvaluation : Empty expression !");
}

void ExprParser::prepareExprEvaluationVec() const
{
  std::set<std::string> trueVars;
  getTrueSetOfVars(trueVars);
  if(trueVars.size()>1)
    {
      std::ostringstream oss; oss << "For this type of evaluation only one not keyword variable authorized : ";
      oss << "having " << trueVars.size() << " : ";
      std::copy(trueVars.begin(),trueVars.end(),std::ostream_iterator<std::string>(oss," ")); oss << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  prepareExprEvaluationVecLowLev();
}

void ExprParser::prepareExprEvaluationVecLowLev() const
{
  if(_leaf)
    {
      LeafExprVar *leafC=dynamic_cast<LeafExprVar *>(_leaf);
      if(leafC)
        leafC->prepareExprEvaluationVec();
    }
  else
    for(std::vector<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
      (*iter).prepareExprEvaluationVecLowLev();
}

Value *ExprParser::evaluateLowLev(Value *valGen) const
{
  checkForEvaluation();
  std::vector<Value *> stackOfVal;
  try
    {
      if(_leaf)
        {
          Value *ret=valGen->newInstance();
          try
            {
              _leaf->fillValue(ret);
            }
          catch(INTERP_KERNEL::Exception& e)
            {
              delete ret;
              throw e;
            }
          stackOfVal.resize(1);
          stackOfVal[0]=ret;
        }
      else
        {
          stackOfVal.resize(_sub_expr.size());
          std::vector<Value *>::iterator iter2=stackOfVal.begin();
          for(std::vector<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++,iter2++)
            *iter2=(*iter).evaluateLowLev(valGen);
        }
      for(std::vector<Function *>::const_iterator iter3=_func_btw_sub_expr.begin();iter3!=_func_btw_sub_expr.end();iter3++)
        (*iter3)->operate(stackOfVal);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      for(std::vector<Value *>::iterator iter4=stackOfVal.begin();iter4!=stackOfVal.end();iter4++)
        delete *iter4;
      throw e;
    }
  return stackOfVal.back();
}

#if __cplusplus >= 201103L

ExprParser::ExprParser(ExprParser&& other):_father(other._father),_leaf(other._leaf),_is_parsing_ok(std::move(other._is_parsing_ok)),_expr(std::move(other._expr)),_sub_expr(std::move(other._sub_expr)),_func_btw_sub_expr(std::move(other._func_btw_sub_expr))
{
  other._leaf=0;
}

ExprParser& ExprParser::operator=(ExprParser&& other)
{
  _father=other._father;
  _is_parsing_ok=std::move(other._is_parsing_ok);
  _leaf=other._leaf;
  _expr=std::move(other._expr);
  _sub_expr=std::move(other._sub_expr);
  _func_btw_sub_expr=std::move(other._func_btw_sub_expr);
  other._leaf=other._leaf;
  other._leaf=0;
  return *this;
}

#endif

void ExprParser::reverseThis()
{
  if(_leaf)
    return ;
  for(std::vector<ExprParser>::iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
    (*iter).reverseThis();
  std::size_t sz(_sub_expr.size());
  std::size_t nbOfTurn(sz/2);
#if __cplusplus >= 201103L
  for(std::size_t i=0;i<nbOfTurn;i++)
    std::swap(_sub_expr[i],_sub_expr[sz-i-1]);
#else
  AutoPtr<char> buf(new char[sizeof(ExprParser)]);
  char *loc(reinterpret_cast<char *>(&_sub_expr[0])),*bufPtr(buf);
  for(std::size_t i=0;i<nbOfTurn;i++)
    {
      std::copy(loc+i*sizeof(ExprParser),loc+(i+1)*sizeof(ExprParser),bufPtr);
      std::copy(loc+(sz-i-1)*sizeof(ExprParser),loc+(sz-i)*sizeof(ExprParser),loc+i*sizeof(ExprParser));
      std::copy(bufPtr,bufPtr+sizeof(ExprParser),loc+(sz-i-1)*sizeof(ExprParser));
    }
#endif
}

ExprParserOfEval ExprParser::convertMeTo() const
{
  std::size_t sz(_sub_expr.size());
  std::vector<ExprParserOfEval> subExpr(sz);
  for(std::size_t i=0;i<sz;i++)
    subExpr[i]=_sub_expr[i].convertMeTo();
  return ExprParserOfEval(_leaf,subExpr,_func_btw_sub_expr);
}

void ExprParser::getSetOfVars(std::set<std::string>& vars) const
{
  if(_leaf)
    {
      LeafExprVar *leafC=dynamic_cast<LeafExprVar *>(_leaf);
      if(leafC)
        vars.insert(leafC->getVar());
    }
  else
    for(std::vector<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
      (*iter).getSetOfVars(vars);
}

void ExprParser::getTrueSetOfVars(std::set<std::string>& trueVars) const
{
  std::set<std::string> vars;
  getSetOfVars(vars);
  trueVars.clear();
  for(std::set<std::string>::const_iterator iter=vars.begin();iter!=vars.end();iter++)
    {
      int tmp;
      if(!LeafExprVar::isRecognizedKeyVar(*iter,tmp))
        trueVars.insert(*iter);
    }
}

void ExprParser::parseDeeper()
{
  for(std::vector<ExprParser>::iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
    if(!(*iter).simplify())
      (*iter).parseDeeper();
}

/*!
 * This method has the responsability to see if this->_expr can be seen as a unary function of something.
 * Something defined as the contain of highest level barckets.
 * Typically '(3*x+2)' and 'cos(4*l+p*n)' will be intercepted by this method whereas '3*x+2' not...etc..
 */
void ExprParser::parseUnaryFunc()
{
  if(_expr[_expr.length()-1]!=')')
    return ;
  //at this level of code _expr
  std::size_t pos1=_expr.find_first_of('(');
  std::size_t pos4=FindCorrespondingOpenBracket(_expr,_expr.length()-1);
  if(pos4!=pos1)
    return ;
  std::string funcName=_expr.substr(0,pos1);
  std::size_t pos2=funcName.find_first_of("+-*/^><",0,7);
  std::size_t pos3=funcName.find_first_not_of("+-*/^><",0,7);
  if(pos2!=std::string::npos && pos3!=std::string::npos)
    return ;//Bracket group is not alone, can't conclude not recursively.
  std::string newExp2=_expr.substr(pos1+1,_expr.length()-pos1-2);
  std::size_t nbOfParamsInFunc=std::count(newExp2.begin(),newExp2.end(),',')+1;
  if(pos3!=std::string::npos)
    _func_btw_sub_expr.push_back(FunctionsFactory::buildFuncFromString(funcName.c_str(),(int)nbOfParamsInFunc));
  else
    {
      std::size_t lgth=funcName.length();
      char tmp[2]; tmp[1]='\0';
      for(std::size_t i=0;i<lgth;i++)
        {
          tmp[0]=funcName[i];
          _func_btw_sub_expr.push_back(FunctionsFactory::buildFuncFromString(tmp,(int)nbOfParamsInFunc));
        }
    }
  std::size_t pos6=0;
  for(std::size_t i=0;i<nbOfParamsInFunc;i++)
    {
      std::size_t pos5=newExp2.find_first_of(',',pos6);
      std::size_t len=std::string::npos;
      if(pos5!=std::string::npos)
        len=pos5-pos6;
      std::string newExp3=newExp2.substr(pos6,len);
      _sub_expr.push_back(ExprParser(newExp3.c_str(),this));
      pos6=pos5+1;
    }
  _is_parsing_ok=true;
}

/*!
 *  This method has the responsability to see if this->_expr is interpretable without any recursion.
 * \return true if no recursion needed, false if this->_expr is too complex to be interpreted at this level.
 * \throw exception if this->_expr is simple enough to try to interprate this and this expression contains an error.
 */
bool ExprParser::tryToInterpALeaf()
{
  std::size_t pos=_expr.find_first_not_of("+-",0,2);
  std::string minimizedExpr=_expr.substr(pos);
  std::size_t pos2=minimizedExpr.find_first_of("+-*/^()<>",0,9);
  if(pos2!=std::string::npos)
    return false;
  delete _leaf;
  _leaf=LeafExpr::buildInstanceFrom(minimizedExpr);
  int nbOfNegs=0;
  for(std::size_t i=0;i<pos;i++)
    if(_expr[i]=='-')
      nbOfNegs++;
  if(nbOfNegs%2)
    _func_btw_sub_expr.push_back(FunctionsFactory::buildUnaryFuncFromString("-"));
  _is_parsing_ok=true;
  return true;
}

void ExprParser::parseForCmp()
{
  std::string::const_iterator iter;
  int curLevel=0;
  std::string curPart;
  bool isParsingSucceed=false;
  for(iter=_expr.begin();iter!=_expr.end();iter++)
    {
      switch(*iter)
        {
        case '>':
        case '<':
          {
            isParsingSucceed=true;
            if(!curPart.empty())
              {
                _sub_expr.push_back(ExprParser(curPart.c_str(),this));
                curPart.clear();
                _func_btw_sub_expr.push_back(FunctionsFactory::buildBinaryFuncFromString(*iter));
              }
            else
              {
                std::ostringstream errMsg;
                char MSGTYP1[]="Error non unary function for '";
                errMsg << EXPR_PARSE_ERR_MSG << MSGTYP1 << *iter << "'";
                std::string tmp=_expr.substr(iter-_expr.begin());
                LocateError(errMsg,tmp,0);
                throw INTERP_KERNEL::Exception(errMsg.str().c_str());
              }
            break;
          }
        case '(':
          curLevel++;
          curPart+=*iter;
          break;
        case ')':
          curLevel--;
          curPart+=*iter;
          break;
        default:
          curPart+=*iter;
        }
    }
  if(isParsingSucceed)
    {
      if(!curPart.empty())
        {
          _sub_expr.push_back(ExprParser(curPart.c_str(),this));
          _is_parsing_ok=true;
        }
      else
        {
          std::ostringstream errMsg;
          char MSGTYP4[]="Error following expression finished by > / < without right part.";
          errMsg << EXPR_PARSE_ERR_MSG << MSGTYP4 << _expr;
          throw INTERP_KERNEL::Exception(errMsg.str().c_str());
        }
    }
}

void ExprParser::parseForAddMin()
{
  std::string::const_iterator iter;
  int curLevel=0;
  std::string curPart;
  bool isParsingSucceed=false;
  for(iter=_expr.begin();iter!=_expr.end();iter++)
    {
      switch(*iter)
        {
        case '+':
        case '-':
          if(curLevel!=0)
            curPart+=*iter;
          else
            {
              if(!curPart.empty())
                {
                  std::string::reverse_iterator accessor=curPart.rbegin();
                  if(*accessor!='*' && *accessor!='/' && *accessor!='^')
                    {
                      isParsingSucceed=true;
                      _sub_expr.push_back(ExprParser(curPart.c_str(),this));
                      curPart.clear();
                      _func_btw_sub_expr.push_back(FunctionsFactory::buildBinaryFuncFromString(*iter));
                    }
                  else
                    curPart+=*iter;
                }
              else
                curPart+=*iter;
            }
        break;
        case '(':
          curLevel++;
          curPart+=*iter;
          break;
        case ')':
          curLevel--;
          curPart+=*iter;
          break;
        default:
          curPart+=*iter;
        }
    }
  if(isParsingSucceed)
    {
      if(!curPart.empty())
        {
          _sub_expr.push_back(ExprParser(curPart.c_str(),this));
          _is_parsing_ok=true;
        }
      else
        {
          std::ostringstream errMsg;
          char MSGTYP4[]="Error following expression finished by +/- without right part.";
          errMsg << EXPR_PARSE_ERR_MSG << MSGTYP4 << _expr;
          throw INTERP_KERNEL::Exception(errMsg.str().c_str());
        }
    }
}

void ExprParser::parseForMulDiv()
{
  std::string::const_iterator iter;
  int curLevel=0;
  std::string curPart;
  bool isParsingSucceed=false;
  for(iter=_expr.begin();iter!=_expr.end();iter++)
    {
      switch(*iter)
        {
        case '/':
        case '*':
          if(curLevel!=0)
            curPart+=*iter;
          else
            {
              isParsingSucceed=true;
              if(!curPart.empty())
                {
                  _sub_expr.push_back(ExprParser(curPart.c_str(),this));
                  curPart.clear();
                  _func_btw_sub_expr.push_back(FunctionsFactory::buildBinaryFuncFromString(*iter));
                }
              else
                {
                  std::ostringstream errMsg;
                  char MSGTYP1[]="Error non unary function for '";
                  errMsg << EXPR_PARSE_ERR_MSG << MSGTYP1 << *iter << "'";
                  std::string tmp=_expr.substr(iter-_expr.begin());
                  LocateError(errMsg,tmp,0);
                  throw INTERP_KERNEL::Exception(errMsg.str().c_str());
                }
            }
        break;
        case '(':
          curLevel++;
          curPart+=*iter;
          break;
        case ')':
          curLevel--;
          curPart+=*iter;
          break;
        default:
          curPart+=*iter;
        }
    }
  if(isParsingSucceed)
    {
      if(!curPart.empty())
        {
          _sub_expr.push_back(ExprParser(curPart.c_str(),this));
          _is_parsing_ok=true;
        }
      else
        {
          std::ostringstream errMsg;
          char MSGTYP5[]="Error following expression finished by *// without right part.";
          errMsg << EXPR_PARSE_ERR_MSG << MSGTYP5 << _expr;
          throw INTERP_KERNEL::Exception(errMsg.str().c_str());
        }
    }
}

void ExprParser::parseForPow()
{
  std::string::const_iterator iter;
  int curLevel=0;
  std::string curPart;
  bool isParsingSucceed=false;
  for(iter=_expr.begin();iter!=_expr.end();iter++)
    {
      switch(*iter)
        {
        case '^':
          if(curLevel!=0)
            curPart+=*iter;
          else
            if(!curPart.empty())
              {
                isParsingSucceed=true;
                _sub_expr.push_back(ExprParser(curPart.c_str(),this));
                curPart.clear();
                _func_btw_sub_expr.push_back(FunctionsFactory::buildBinaryFuncFromString(*iter));
              }
            else
              {
                std::ostringstream errMsg;
                char MSGTYP1[]="Error non unary function for '";
                errMsg << EXPR_PARSE_ERR_MSG << MSGTYP1 << *iter << "'";
                std::string tmp=_expr.substr(iter-_expr.begin());
                LocateError(errMsg,tmp,0);curPart+=*iter;
                throw INTERP_KERNEL::Exception(errMsg.str().c_str());
              }
          break;
        case '(':
          curLevel++;
          curPart+=*iter;
          break;
        case ')':
          curLevel--;
          curPart+=*iter;
          break;
        default:
          curPart+=*iter;
        }
    }
  if(isParsingSucceed)
    {
      if(!curPart.empty())
        {
          _sub_expr.push_back(ExprParser(curPart.c_str(),this));
          _is_parsing_ok=true;
        }
      else
        {
          std::ostringstream errMsg;
          char MSGTYP6[]="Error following expression finished by ^ without right part.";
          errMsg << EXPR_PARSE_ERR_MSG << MSGTYP6 << _expr;
          throw INTERP_KERNEL::Exception(errMsg.str().c_str());
        }
    }
}

void ExprParser::releaseFunctions()
{
  for(std::vector<Function *>::iterator iter=_func_btw_sub_expr.begin();iter!=_func_btw_sub_expr.end();iter++)
    delete *iter;
  _func_btw_sub_expr.clear();
}

/*!
 * This method parse this->_expr at the current level.
 * This method first try to see if this->_expr is a leaf, if not it try a unary function of something (see INTERP_KERNEL::ExprParser::parseUnaryFunc method)
 * If true is returned, no deeper parsing needed, if false is returned for a full parsing of this->_expr INTERP_KERNEL::ExprParser::parseDeeper call needed.
 */
bool ExprParser::simplify()
{
  if(tryToInterpALeaf())
    return true;
  parseUnaryFunc();
  if(!_is_parsing_ok)
    {
      parseForCmp();
      if(!_is_parsing_ok)
        {
          parseForAddMin();
          if(!_is_parsing_ok)
            {
              parseForMulDiv();
              if(!_is_parsing_ok)
                parseForPow();
            }
        }
    }
  if(!_is_parsing_ok)
    {
      std::ostringstream errMsg;
      char MSGTYP3[]="Error in interpreting : ";
      errMsg << EXPR_PARSE_ERR_MSG << MSGTYP3 << _expr;
      LocateError(errMsg,_expr,0);
      throw INTERP_KERNEL::Exception(errMsg.str().c_str());
    }
  return false;
}

void ExprParser::checkBracketsParity() const
{
  std::string::const_iterator iter;
  int curLevel=0;
  for(iter=_expr.begin();iter!=_expr.end();iter++)
    {
      if(*iter=='(')
        curLevel++;
      else if(*iter==')')
        {
          if(curLevel==0)
            {
              std::ostringstream errMsg;
              char MSGTYP1[]="Error in brackets : closing brackets ')' before openning '('";
              errMsg << EXPR_PARSE_ERR_MSG << MSGTYP1;
              LocateError(errMsg,_expr,(int)std::distance(_expr.begin(),iter));
              throw INTERP_KERNEL::Exception(errMsg.str().c_str());
            }
          curLevel--;
        }
    }
  if(curLevel!=0)
    {
      std::ostringstream errMsg;
      char MSGTYP2[]="Error in brackets : not finally closed expr.";
      errMsg << EXPR_PARSE_ERR_MSG << MSGTYP2;
      throw INTERP_KERNEL::Exception(errMsg.str().c_str());
    }
}

/*!
 * This method substitutes part in [bg,end) in expr by the content of (str(id)) and returns the double value representation in expr[bg,end).
 * If double representation is invalid an exception is thrown.
 * This method returns a delta that is the delta to operate to pos in expr after substitution.
 */
double ExprParser::ReplaceAndTraduce(std::string& expr, int id, std::size_t bg, std::size_t end, int& delta)
{
  static const char MSG[]="Interal error : A string expected to be a float is not one ! Bug to signal !";
  std::istringstream stream;
  std::ostringstream oss;
  std::size_t end2=end!=std::string::npos?end-bg:end;
  std::string tmp=expr.substr(bg,end2);
  stream.str(tmp);
  double ret=std::numeric_limits<double>::max();
  stream >> ret;
  if(stream.fail())
    throw INTERP_KERNEL::Exception(MSG);
  if(!stream.eof())
    throw INTERP_KERNEL::Exception(MSG);
  oss << id;
  std::string tmp2(oss.str());
  std::size_t l1=tmp.length();
  delta=(int)tmp2.length()-(int)l1;
  expr.replace(bg,l1,tmp2);
  return ret;
}

/*!
 * This method makes the assumption that _expr has no white space.
 * This method scans _expr finding in greedy mode the following pattern :
 * {0..9}+{.}?{0..9}*{{eE}{-}?{0..9}+}?
 */
void ExprParser::fillValuesInExpr(std::vector<double>& valuesInExpr)
{
  const char FIGURES[]="0123456789";
  const std::string other("+-*^/(<>,");
  std::size_t lgth=_expr.length();
  int id=0,delta;
  for(std::size_t pos=0;pos!=std::string::npos;id++)
    {
      std::size_t pos2=_expr.find_first_of(FIGURES,pos,10);
      if(pos2==std::string::npos)
        break;
      if(pos2>0)
        {//treat case of "x*log10(x)" -> "10" should NOT be intercepted by this
          if(other.find_first_of(_expr[pos2-1])==std::string::npos)
            {
              pos=_expr.find_first_not_of(FIGURES,pos2,10);
              id--;
              continue;
            }
          if(_expr[pos2-1]==')')
            {
              pos=_expr.find_first_not_of(FIGURES,pos2,10);
              std::ostringstream oss; oss << "Problem on parsing : Number \"" << _expr.substr(pos2,pos!=std::string::npos?pos2-pos:std::string::npos);
              oss << "\" is right after close parenthesis... ')'";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      std::size_t pos3=_expr.find_first_not_of(FIGURES,pos2,10);
      if(pos3==std::string::npos)
        {//"x+1223442320"
          valuesInExpr.push_back(ReplaceAndTraduce(_expr,id,pos2,std::string::npos,delta));
          break;
        }
      if(_expr[pos3]=='.')
        pos3++;
      if(pos3<lgth)
        {
          std::size_t pos4=_expr.find_first_not_of(FIGURES,pos3,10);
          if(pos4==std::string::npos)
            {//"x+1223334.223"
              valuesInExpr.push_back(ReplaceAndTraduce(_expr,id,pos2,std::string::npos,delta));
              break;
            }
          else
            {
              if(_expr[pos4]!='e' && _expr[pos4]!='E')
                {//"x+1223334.223+x"
                  valuesInExpr.push_back(ReplaceAndTraduce(_expr,id,pos2,pos4,delta));
                  pos=pos4+delta;
                  continue;
                }
              else
                {
                  if(++pos4<lgth)
                    {
                      if(_expr[pos4]=='+' || _expr[pos4]=='-')
                        pos4++;
                      if(pos4>=lgth)
                        {//"x+1223334.223e+" or "1223334.223E-"
                          std::ostringstream oss; oss << "Invalid expr : float number at the end of expr is invalid lacking number after exponential and sign ! -> \"" << _expr.substr(pos2) << "\"";
                          throw INTERP_KERNEL::Exception(oss.str().c_str());
                        }
                      std::size_t pos5=_expr.find_first_not_of(FIGURES,pos4,10);
                      if(pos4==pos5)
                        {//"x+1223334.223e+x" or "1223334.223E-y"
                          std::ostringstream oss; oss << "Invalid expr : float number in expr is invalid lacking number after exponential ! -> \"" << _expr.substr(pos2,pos4-pos2) << "\"";
                          throw INTERP_KERNEL::Exception(oss.str().c_str());
                        }
                      //OK, normal case
                      valuesInExpr.push_back(ReplaceAndTraduce(_expr,id,pos2,pos5,delta));
                      pos=pos5+delta;
                      continue;
                    }
                  else//"x+1223334.223e"
                    {
                      std::ostringstream oss; oss << "Invalid expr : float number at the end of expr is invalid lacking number after exponential ! " << _expr.substr(pos2);
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                }
            }
        }
      else
        {//"x+1223334."
          valuesInExpr.push_back(ReplaceAndTraduce(_expr,id,pos2,std::string::npos,delta));
          break;
        }
    }
}

void ExprParser::replaceValues(const std::vector<double>& valuesInExpr)
{
  if(_leaf)
    _leaf->replaceValues(valuesInExpr);
  else
    {
      for(std::vector<ExprParser>::iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
        (*iter).replaceValues(valuesInExpr);
    }
}

void ExprParser::LocateError(std::ostream& stringToDisp, const std::string& srcOfErr, int posOfErr)
{
  stringToDisp << "Position is " << posOfErr << " of string : \"" <<  srcOfErr << "\"" << std::endl;
}

char *ExprParser::compileX86() const
{
  std::vector<std::string> ass;
  //need in stack
  ass.push_back("push ebp");
  ass.push_back("mov ebp,esp");
  compileX86LowLev(ass);
  ass.push_back("pop ebp");
  ass.push_back("ret");
  std::cout << std::endl;
  for(std::vector<std::string>::const_iterator iter=ass.begin();iter!=ass.end();iter++)
    std::cout << "        " << *iter << std::endl;
  AsmX86 asmb;
  std::vector<char> output=asmb.convertIntoMachineLangage(ass);
  for(std::vector<char>::const_iterator iter=output.begin();iter!=output.end();iter++)
    std::cout << std::hex << (int)((unsigned char)(*iter)) << " ";
  std::cout << std::endl;
  unsigned offset;
  return asmb.copyToExecMemZone(output,offset);
}

char *ExprParser::compileX86_64() const
{
  std::vector<std::string> ass;
  //need in stack
  ass.push_back("push rbp");
  ass.push_back("mov rbp,rsp");
  compileX86_64LowLev(ass);
  ass.push_back("sub rsp,8");
  ass.push_back("fst qword [rsp]");
  ass.push_back("movsd xmm0,[rsp]");
  ass.push_back("add rsp,8");
  ass.push_back("leave");
  ass.push_back("ret");
  std::cout << std::endl;
  for(std::vector<std::string>::const_iterator iter=ass.begin();iter!=ass.end();iter++)
    std::cout << "        " << *iter << std::endl;
  AsmX86 asmb;
  std::vector<char> output=asmb.convertIntoMachineLangage(ass);
  for(std::vector<char>::const_iterator iter=output.begin();iter!=output.end();iter++)
    std::cout << std::hex << (int)((unsigned char)(*iter)) << " ";
  std::cout << std::endl;
  unsigned offset;
  return asmb.copyToExecMemZone(output,offset);
}

void ExprParser::compileX86LowLev(std::vector<std::string>& ass) const
{
  if(_leaf)
    _leaf->compileX86(ass);
  else
    {
      for(std::vector<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
        (*iter).compileX86LowLev(ass);
    }
  for(std::vector<Function *>::const_iterator iter2=_func_btw_sub_expr.begin();iter2!=_func_btw_sub_expr.end();iter2++)
    (*iter2)->operateX86(ass);
}

void ExprParser::compileX86_64LowLev(std::vector<std::string>& ass) const
{
  if(_leaf)
    _leaf->compileX86_64(ass);
  else
    {
      for(std::vector<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
        (*iter).compileX86_64LowLev(ass);
    }
  for(std::vector<Function *>::const_iterator iter2=_func_btw_sub_expr.begin();iter2!=_func_btw_sub_expr.end();iter2++)
    (*iter2)->operateX86(ass);
}

double LeafExprVal::getDoubleValue() const
{
  return _value;
}

void LeafExprVal::compileX86(std::vector<std::string>& ass) const
{
  ass.push_back("sub esp,8");
  const int *b=reinterpret_cast<const int *>(&_value),*c=reinterpret_cast<const int *>(&_value);
  c++;
  std::ostringstream oss;
  oss << std::hex;
  oss << "mov dword [esp+4],0x" << *c;
  ass.push_back(oss.str());
  oss.str("");
  oss << "mov dword [esp],0x" << *b;
  ass.push_back(oss.str());
  ass.push_back("fld qword [esp]");
  ass.push_back("add esp,8");
}

void LeafExprVal::compileX86_64(std::vector<std::string>& ass) const
{
  ass.push_back("sub rsp,8");
  const int *b=reinterpret_cast<const int *>(&_value),*c=reinterpret_cast<const int *>(&_value);
  c++;
  std::ostringstream oss;
  oss << std::hex;
  oss << "mov dword [rsp+4],0x" << *c;
  ass.push_back(oss.str());
  oss.str("");
  oss << "mov dword [rsp],0x" << *b;
  ass.push_back(oss.str());
  ass.push_back("fld qword [rsp]");
  ass.push_back("add rsp,8");
}

double LeafExprVar::getDoubleValue() const
{
  if(_fast_pos>=0)
    return _val[_fast_pos];
  else
    {
      int pos(-7-_fast_pos);
      return pos==_ref_pos?1.:0.;
    }
}

void LeafExprVar::compileX86(std::vector<std::string>& ass) const
{
  ass.push_back("fld qword [ebp+8]");
}

void LeafExprVar::compileX86_64(std::vector<std::string>& ass) const
{
  ass.push_back("sub rsp,8");
  ass.push_back("movsd [rsp],xmm0");
  ass.push_back("fld qword [rsp]");
  ass.push_back("add rsp,8");
}

int ExprParser::getStackSizeToPlayX86(const ExprParser *asker) const
{
  if(asker)
    {
      int sz=_father->getStackSizeToPlayX86(this);
      int i=0;
      for(std::vector<ExprParser>::const_reverse_iterator iter=_sub_expr.rbegin();iter!=_sub_expr.rend();iter++,i++)
        {
          const ExprParser& obj=(*iter);
          const ExprParser *pt=&obj;
          if(pt==asker)
            return sz-i;
        }
      throw INTERP_KERNEL::Exception("error getStackSizeToPlayX86 an object ExprParser called as father, whereas it is not one !");
    }
  else
    {
      if(!_father)
        return MAX_X86_FP_ST;
      return _father->getStackSizeToPlayX86(this);
    }
}
