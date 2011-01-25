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

#include "InterpKernelExprParser.hxx"
#include "InterpKernelValue.hxx"
#include "InterpKernelAsmX86.hxx"

#include <cctype>
#include <sstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>

using namespace INTERP_KERNEL;

const char LeafExprVar::END_OF_RECOGNIZED_VAR[]="Vec";

const char ExprParser::WHITE_SPACES[]=" \n";

const char ExprParser::EXPR_PARSE_ERR_MSG[]="Invalid expression detected : ";

LeafExpr *LeafExpr::buildInstanceFrom(const std::string& expr) throw(INTERP_KERNEL::Exception)
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

void LeafExprVal::fillValue(Value *val) const throw(INTERP_KERNEL::Exception)
{
  val->setDouble(_value);
}

LeafExprVar::LeafExprVar(const std::string& var):_fast_pos(-1),_var_name(var)
{
}

void LeafExprVar::fillValue(Value *val) const throw(INTERP_KERNEL::Exception)
{
  val->setVarname(_fast_pos,_var_name);
}

void LeafExprVar::prepareExprEvaluation(const std::vector<std::string>& vars) const throw(INTERP_KERNEL::Exception)
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
      return;
    }
  _fast_pos=iter-vars.begin();
}

void LeafExprVar::prepareExprEvaluationVec() const throw(INTERP_KERNEL::Exception)
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

LeafExprVar::~LeafExprVar()
{
}

ExprParser::ExprParser(const char *expr, ExprParser *father):_father(father),_is_parsed(false),_leaf(0),_is_parsing_ok(false),_expr(expr)
{
}

//! For \b NOT null terminated strings coming from FORTRAN.
ExprParser::ExprParser(const char *expr, int lgth, ExprParser *father):_father(father),_is_parsed(false),_leaf(0),_is_parsing_ok(false)
{
  _expr=buildStringFromFortran(expr,lgth);
}

ExprParser::~ExprParser()
{
  delete _leaf;
  releaseFunctions();
}

std::size_t ExprParser::findCorrespondingOpenBracket(const std::string& expr, std::size_t posOfCloseBracket)
{
  int level=0;
  for(std::size_t iter=posOfCloseBracket-1;iter>=0;iter--)
    if(expr[iter]==')')
      level++;
    else if(expr[iter]=='(')
      {
        if(level==0)
          return iter;
        else
          level--;
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

void ExprParser::parse() throw(INTERP_KERNEL::Exception)
{
  _is_parsed=true;
  _is_parsing_ok=false;
  _sub_expr.clear();
  releaseFunctions();
  if(!_expr.empty())
    {
      checkBracketsParity();
      if(!simplify())
        parseDeeper();
    }
  _is_parsing_ok=true;
}

double ExprParser::evaluate() const throw(INTERP_KERNEL::Exception)
{
  Value *gen=new ValueDouble;
  ValueDouble *res=(ValueDouble *)evaluateLowLev(gen);
  delete gen;
  double ret=res->getData();
  delete res;
  return ret;
}

DecompositionInUnitBase ExprParser::evaluateUnit() const throw(INTERP_KERNEL::Exception)
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

void ExprParser::evaluateExpr(int szOfOutParam, const double *inParam, double *outParam) const throw(INTERP_KERNEL::Exception)
{
  Value *gen=new ValueDoubleExpr(szOfOutParam,inParam);
  ValueDoubleExpr *res=0;
  try
    {
      res=(ValueDoubleExpr *)evaluateLowLev(gen);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete gen;
      throw e;
    }
  delete gen;
  std::copy(res->getData(),res->getData()+szOfOutParam,outParam);
  delete res;
}

void ExprParser::prepareExprEvaluation(const std::vector<std::string>& vars) const throw(INTERP_KERNEL::Exception)
{
  if(_leaf)
    {
      LeafExprVar *leafC=dynamic_cast<LeafExprVar *>(_leaf);
      if(leafC)
        leafC->prepareExprEvaluation(vars);
    }
  else
    for(std::list<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
      (*iter).prepareExprEvaluation(vars);
}

void ExprParser::prepareExprEvaluationVec() const throw(INTERP_KERNEL::Exception)
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

void ExprParser::prepareExprEvaluationVecLowLev() const throw(INTERP_KERNEL::Exception)
{
  if(_leaf)
    {
      LeafExprVar *leafC=dynamic_cast<LeafExprVar *>(_leaf);
      if(leafC)
        leafC->prepareExprEvaluationVec();
    }
  else
    for(std::list<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
      (*iter).prepareExprEvaluationVecLowLev();
}

Value *ExprParser::evaluateLowLev(Value *valGen) const throw(INTERP_KERNEL::Exception)
{
  if(!_is_parsing_ok)
    throw INTERP_KERNEL::Exception("Parsing fails ! Invalid expression !");
  if(_sub_expr.empty() && !_leaf)
    throw INTERP_KERNEL::Exception("Empty expression !");
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
          std::vector<Value *>::reverse_iterator iter2=stackOfVal.rbegin();
          for(std::list<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++,iter2++)
            *iter2=(*iter).evaluateLowLev(valGen);
        }
      std::list<Function *>::const_iterator iter3;
      for(iter3=_func_btw_sub_expr.begin();iter3!=_func_btw_sub_expr.end();iter3++)
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

void ExprParser::getSetOfVars(std::set<std::string>& vars) const
{
  if(_leaf)
    {
      LeafExprVar *leafC=dynamic_cast<LeafExprVar *>(_leaf);
      if(leafC)
        vars.insert(leafC->getVar());
    }
  else
    for(std::list<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
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

void ExprParser::parseDeeper() throw(INTERP_KERNEL::Exception)
{
  for(std::list<ExprParser>::iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
    if(!(*iter).simplify())
      (*iter).parseDeeper();
}

/*!
 * This method has the responsability to see if this->_expr can be seen as a unary function of something.
 * Something defined as the contain of highest level barckets.
 * Typically '(3*x+2)' and 'cos(4*l+p*n)' will be intercepted by this method whereas '3*x+2' not...etc..
 */
void ExprParser::parseUnaryFunc() throw(INTERP_KERNEL::Exception)
{
  if(_expr[_expr.length()-1]!=')')
    return ;
  //at this level of code _expr 
  std::size_t pos1=_expr.find_first_of('(');
  std::size_t pos4=findCorrespondingOpenBracket(_expr,_expr.length()-1);
  if(pos4!=pos1)
    return ;
  std::string funcName=_expr.substr(0,pos1);
  std::size_t pos2=funcName.find_first_of("+-*/^><",0,7);
  std::size_t pos3=funcName.find_first_not_of("+-*/^><",0,7);
  if(pos2!=std::string::npos && pos3!=std::string::npos)
    return ;//Bracket group is not alone, can't conclude not recursively.
  std::string newExp2=_expr.substr(pos1+1,_expr.length()-pos1-2);
  int nbOfParamsInFunc=std::count(newExp2.begin(),newExp2.end(),',')+1;
  if(pos3!=std::string::npos)
    _func_btw_sub_expr.push_back(FunctionsFactory::buildFuncFromString(funcName.c_str(),nbOfParamsInFunc));
  else
    {
      int lgth=funcName.length();
      char tmp[2]; tmp[1]='\0';
      for(int i=0;i<lgth;i++)
        {
          tmp[0]=funcName[i];
          _func_btw_sub_expr.push_back(FunctionsFactory::buildFuncFromString(tmp,nbOfParamsInFunc));
        }
    }
  std::size_t pos6=0;
  for(int i=0;i<nbOfParamsInFunc;i++)
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
bool ExprParser::tryToInterpALeaf() throw(INTERP_KERNEL::Exception)
{
  std::size_t pos=_expr.find_first_not_of("+-",0,2);
  std::string minimizedExpr=_expr.substr(pos);
  std::size_t pos2=minimizedExpr.find_first_of("+-*/^()",0,7);
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

void ExprParser::parseForCmp() throw(INTERP_KERNEL::Exception)
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
                locateError(errMsg,tmp,0);
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

void ExprParser::parseForAddMin() throw(INTERP_KERNEL::Exception)
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

void ExprParser::parseForMulDiv() throw(INTERP_KERNEL::Exception)
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
                  locateError(errMsg,tmp,0);
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

void ExprParser::parseForPow() throw(INTERP_KERNEL::Exception)
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
                locateError(errMsg,tmp,0);curPart+=*iter;
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
  for(std::list<Function *>::iterator iter=_func_btw_sub_expr.begin();iter!=_func_btw_sub_expr.end();iter++)
    delete *iter;
  _func_btw_sub_expr.clear();
}

/*!
 * This method parse this->_expr at the current level.
 * This method first try to see if this->_expr is a leaf, if not it try a unary function of something (see INTERP_KERNEL::ExprParser::parseUnaryFunc method)
 * If true is returned, no deeper parsing needed, if false is returned for a full parsing of this->_expr INTERP_KERNEL::ExprParser::parseDeeper call needed.
 */
bool ExprParser::simplify() throw(INTERP_KERNEL::Exception)
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
      locateError(errMsg,_expr,0);
      throw INTERP_KERNEL::Exception(errMsg.str().c_str());
    }
  return false;
}

void ExprParser::checkBracketsParity() const throw(INTERP_KERNEL::Exception)
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
              locateError(errMsg,_expr,iter-_expr.begin());
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

void ExprParser::locateError(std::ostream& stringToDisp, const std::string& srcOfErr, int posOfErr)
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
      for(std::list<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
        (*iter).compileX86LowLev(ass);
    }
  for(std::list<Function *>::const_iterator iter2=_func_btw_sub_expr.begin();iter2!=_func_btw_sub_expr.end();iter2++)
    (*iter2)->operateX86(ass);
}

void ExprParser::compileX86_64LowLev(std::vector<std::string>& ass) const
{
  if(_leaf)
    _leaf->compileX86_64(ass);
  else
    {
      for(std::list<ExprParser>::const_iterator iter=_sub_expr.begin();iter!=_sub_expr.end();iter++)
        (*iter).compileX86_64LowLev(ass);
    }
  for(std::list<Function *>::const_iterator iter2=_func_btw_sub_expr.begin();iter2!=_func_btw_sub_expr.end();iter2++)
    (*iter2)->operateX86(ass);
}

void LeafExprVal::compileX86(std::vector<std::string>& ass) const
{
  ass.push_back("sub esp,8");
  int *b=(int *)&_value,*c=(int *)&_value;
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
  int *b=(int *)&_value,*c=(int *)&_value;
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
      for(std::list<ExprParser>::const_reverse_iterator iter=_sub_expr.rbegin();iter!=_sub_expr.rend();iter++,i++)
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
