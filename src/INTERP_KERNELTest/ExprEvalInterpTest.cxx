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

#include "ExprEvalInterpTest.hxx"
#include "InterpKernelExprParser.hxx"

#include <limits>
#include <iterator>

using namespace INTERP_TEST;

void ExprEvalInterpTest::testBuildStringFromFortran()
{
  char toto1[]="123456  ";
  char result[]="123456";
  std::string titi;
  titi=INTERP_KERNEL::ExprParser::buildStringFromFortran(toto1,8);
  CPPUNIT_ASSERT_EQUAL(6,(int)titi.length());
  CPPUNIT_ASSERT(titi==result);
  //
  char toto2[]=" 23456  ";
  char result2[]=" 23456";
  titi=INTERP_KERNEL::ExprParser::buildStringFromFortran(toto2,8);
  CPPUNIT_ASSERT(titi==result2);
  CPPUNIT_ASSERT_EQUAL(6,(int)titi.length());
  //
  char toto3[]="  3456  ";
  char result3[]="  3456";
  titi=INTERP_KERNEL::ExprParser::buildStringFromFortran(toto3,8);
  CPPUNIT_ASSERT(titi==result3);
  CPPUNIT_ASSERT_EQUAL(6,(int)titi.length());
  //
  char toto4[]="        ";
  titi=INTERP_KERNEL::ExprParser::buildStringFromFortran(toto4,8);
  CPPUNIT_ASSERT_EQUAL(0,(int)titi.length());
  //
  char toto5[]="  345677";
  titi=INTERP_KERNEL::ExprParser::buildStringFromFortran(toto5,8);
  CPPUNIT_ASSERT(titi==toto5);
  CPPUNIT_ASSERT_EQUAL(8,(int)titi.length());
}

void ExprEvalInterpTest::testDeleteWhiteSpaces()
{
  char toto[]=" jkhjkh ooooppp l ";
  char result[]="jkhjkhoooopppl";
  std::string totoS(toto);
  std::string totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
  //
  char toto2[]=" jkhjkh     ooooppp    l ";
  totoS=toto2;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
  //
  char toto3[]=" jkhjkh     oooo pppl ";
  totoS=toto3;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
  //
  char toto4[]=" jkhjkh     oooo pppl";
  totoS=toto4;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
  //
  char toto5[]="jkhjkh     oooo pppl";
  totoS=toto5;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
  //
  totoS=result;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
  //
  char toto6[]="j k h j k h o o o o p p p l";
  totoS=toto6;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
  //
  char toto7[]="j  k  h j    k h   o  o  o  o  p  pp    l";
  totoS=toto7;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
  //
  char toto8[]="           ";
  totoS=toto8;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR.empty());
  //
  char toto9[]="";
  totoS=toto9;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR.empty());
  //
  char toto10[]="j\n k \nh\nj \n\n  k\nh \n o \no\n o\n o \np\n\npp \n\n l";
  totoS=toto10;
  totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
  CPPUNIT_ASSERT(totoR==result);
  CPPUNIT_ASSERT_EQUAL(14,(int)totoR.length());
}

void ExprEvalInterpTest::testInterpreter0()
{
  INTERP_KERNEL::ExprParser expr1("3*-2");
  expr1.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-6.,expr1.evaluate(),1e-15);
  INTERP_KERNEL::ExprParser expr2("-2.3");
  expr2.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-2.3,expr2.evaluate(),1e-15);
  INTERP_KERNEL::ExprParser expr3("--2.3");
  expr3.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.3,expr3.evaluate(),1e-15);
  INTERP_KERNEL::ExprParser expr4("-++2.3");
  expr4.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-2.3,expr4.evaluate(),1e-15);
  INTERP_KERNEL::ExprParser expr5("+2.3");
  expr5.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.3,expr5.evaluate(),1e-15);
  INTERP_KERNEL::ExprParser expr6("3^-1");
  expr6.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333333,expr6.evaluate(),1e-15);
}

void ExprEvalInterpTest::testInterpreter1()
{
  INTERP_KERNEL::ExprParser expr1("3+2*5");
  expr1.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(13.,expr1.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr2("3+2^3*5");
  expr2.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(43.,expr2.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr3("3+2^(2*5)");
  expr3.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1027.,expr3.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr4("(3.2+4.3)*(1.3+2.3*7.8)");
  expr4.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(144.3,expr4.evaluate(),1e-10);
  INTERP_KERNEL::ExprParser expr5("(3.2+4.3)*cos(1.3+2.3*7.8)");
  expr5.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.9355510138337619,expr5.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr6("3+2-4-7+4.3");
  expr6.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.7,expr6.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr7("3.2*4.5/3.3/2.2");
  expr7.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.9834710743801653,expr7.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr8("3.2*4.5/3.3/2.2");
  expr8.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.9834710743801653,expr8.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr9("(((1.23456789)))");
  expr9.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.23456789,expr9.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr10("3.2*((2*5.2+6.)+(1.2*1.2+3.))");
  expr10.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(66.688,expr10.evaluate(),1e-13);
  INTERP_KERNEL::ExprParser expr11("((3.2*(((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr11.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(66.688,expr11.evaluate(),1e-13);
  INTERP_KERNEL::ExprParser expr12("((3.2*(cos((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr12.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.3038041398761016,expr12.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr13("((3.2*(sin((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr13.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.9223440531261784,expr13.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr14("((3.2*(tan((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr14.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-7.1724737512280257,expr14.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr15("((3.2*(sqrt((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr15.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.608271629457059,expr15.evaluate(),1e-13);
  INTERP_KERNEL::ExprParser expr16("-((3.2*(sqrt((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr16.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-14.608271629457059,expr16.evaluate(),1e-13);
  INTERP_KERNEL::ExprParser expr17("(-(3.2*(sqrt((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr17.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-14.608271629457059,expr17.evaluate(),1e-13);
  INTERP_KERNEL::ExprParser expr18("((-3.2*(sqrt((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr18.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-14.608271629457059,expr18.evaluate(),1e-13);
  INTERP_KERNEL::ExprParser expr19("((3.2*(exp((6.+2*5.2)+(1.2*1.2+3.)))))");
  expr19.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3596226038.1784945,expr19.evaluate(),1e-6);
  INTERP_KERNEL::ExprParser expr20("((3.2*(ln((2*5.2+6.)+(1.2*1.2+3.)))))");
  expr20.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9.7179974940325309,expr20.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr21("max(3.2,4.5)");
  expr21.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.5,expr21.evaluate(),1e-14);
  INTERP_KERNEL::ExprParser expr22("3.*max(((3.2*(ln((2*5.2+6.)+(1.2*1.2+3.))))),((3.2*(exp((6.+2*5.2)+(1.2*1.2+3.))))))");
  expr22.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10788678114.535484,expr22.evaluate(),1e-5);
  INTERP_KERNEL::ExprParser expr23("min(3.2,4.5)");
  expr23.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2,expr23.evaluate(),1e-14);
}

void ExprEvalInterpTest::testInterpreter2()
{
  INTERP_KERNEL::ExprParser expr1("3.5*x+x*x*x/(2+x)+2*5*y");
  expr1.parse();
  std::set<std::string> res,expected;
  expr1.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
  expected.insert("x"); expected.insert("y");
  CPPUNIT_ASSERT(std::equal(res.begin(),res.end(),expected.begin()));
  double xyValue[2]={1.,3.};
  double res1;
  std::vector<std::string> vars; vars.push_back("x"); vars.push_back("y");
  expr1.prepareExprEvaluation(vars,2,1);
  expr1.evaluateExpr(1,xyValue,&res1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(33.833333333333336,res1,1e-13);
  xyValue[0]=-2.;
  CPPUNIT_ASSERT_THROW(expr1.evaluateExpr(1,xyValue,&res1),INTERP_KERNEL::Exception);
  double res2[2];
  xyValue[0]=1.;
  expr1.evaluateExpr(2,xyValue,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(33.833333333333336,res2[0],1e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(33.833333333333336,res2[1],1e-13);
  INTERP_KERNEL::ExprParser expr2("3.5*tan(2.3*x)*IVec+(cos(1.2+y/x)*JVec)");
  expr2.parse();
  res.clear(); expected.clear();
  expr2.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(4,(int)res.size());
  expected.insert("x"); expected.insert("y"); expected.insert("IVec"); expected.insert("JVec");
  CPPUNIT_ASSERT(std::equal(res.begin(),res.end(),expected.begin()));
  expr2.prepareExprEvaluation(vars,2,2);
  expr2.evaluateExpr(2,xyValue,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.9172477460694637,res2[0],1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.49026082134069943,res2[1],1e-14);
  INTERP_KERNEL::ExprParser expr3("3.5*u+u^2.4+2.");
  expr3.parse();
  expr3.prepareExprEvaluationVec();
  expr3.evaluateExpr(2,xyValue,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.5,res2[0],1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(26.466610165238237,res2[1],1e-14);
  INTERP_KERNEL::ExprParser expr4("3.5*v+u^2.4+2.");
  expr4.parse();
  CPPUNIT_ASSERT_THROW(expr4.prepareExprEvaluationVec(),INTERP_KERNEL::Exception);
}

void ExprEvalInterpTest::testInterpreterUnit0()
{
  INTERP_KERNEL::ExprParser expr1("kg");
  expr1.parse();
  INTERP_KERNEL::DecompositionInUnitBase unit=expr1.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(1,0,0,0,0,0.,1000.));
  INTERP_KERNEL::ExprParser expr2("kgT");
  expr2.parse();
  CPPUNIT_ASSERT_THROW(expr2.evaluateUnit(),INTERP_KERNEL::Exception);
  INTERP_KERNEL::ExprParser expr3("g");
  expr3.parse();
  unit=expr3.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(1,0,0,0,0,0.,1.));
  INTERP_KERNEL::ExprParser expr4("g*m");
  expr4.parse();
  unit=expr4.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(1,1,0,0,0,0.,1.));
  INTERP_KERNEL::ExprParser expr5("g*m/K");
  expr5.parse();
  unit=expr5.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(1,1,0,0,-1,0.,1.));
  INTERP_KERNEL::ExprParser expr6("g*m/K^2");
  expr6.parse();
  unit=expr6.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(1,1,0,0,-2,0.,1.));
  INTERP_KERNEL::ExprParser expr7("g/K^2*m");
  expr7.parse();
  unit=expr7.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(1,1,0,0,-2,0.,1.));
  INTERP_KERNEL::ExprParser expr8("g/(K^2*m)");
  expr8.parse();
  unit=expr8.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(1,-1,0,0,-2,0.,1.));
  INTERP_KERNEL::ExprParser expr9("km/h");
  expr9.parse();
  unit=expr9.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,1,-1,0,0,0.,0.27777777777777779));
  INTERP_KERNEL::ExprParser expr10("m/s");
  expr10.parse();
  unit=expr10.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,1,-1,0,0,0.,1.));
  INTERP_KERNEL::ExprParser expr11("m+s");
  expr11.parse();
  CPPUNIT_ASSERT_THROW(expr11.evaluateUnit(),INTERP_KERNEL::Exception);
  INTERP_KERNEL::ExprParser expr12("m-m");
  expr12.parse();
  CPPUNIT_ASSERT_THROW(expr12.evaluateUnit(),INTERP_KERNEL::Exception);
  const char expr13C[3]={-0x50,0x43,0x0};
  INTERP_KERNEL::ExprParser expr13(expr13C);
  expr13.parse();
  unit=expr13.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,0,0,0,1,273.15,1.));
  const char expr14C[4]={-0x3E,-0x50,0x43,0x0};
  INTERP_KERNEL::ExprParser expr14(expr14C);
  expr14.parse();
  unit=expr14.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,0,0,0,1,273.15,1.));
  INTERP_KERNEL::ExprParser expr15("kN/kg");
  expr15.parse();
  unit=expr15.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,1,-2,0,0,0.,1000.));
  INTERP_KERNEL::ExprParser expr16("cm");
  expr16.parse();
  unit=expr16.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,1,0,0,0,0.,0.01));
  INTERP_KERNEL::ExprParser expr17("m");
  expr17.parse();
  unit=expr17.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,1,0,0,0,0.,1));
  const char expr18C[3]={-0x08,0x43,0x0};
  INTERP_KERNEL::ExprParser expr18(expr18C);
  expr18.parse();
  unit=expr18.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,0,0,0,1,273.15,1.));
  const char expr19C[6]={-0x50,0x43,0x2F,-0x50,0x43,0x0};
  INTERP_KERNEL::ExprParser expr19(expr19C);
  expr19.parse();
  unit=expr19.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,0,0,0,0,0.,1.));
  const char expr20C[9]={-0x50,0x43,0x2A,-0x50,0x43,0x2F,-0x50,0x43,0x0};
  INTERP_KERNEL::ExprParser expr20(expr20C);
  expr20.parse();
  unit=expr20.evaluateUnit();
  CPPUNIT_ASSERT(unit.isEqual(0,0,0,0,1,0.,1.));
}

void ExprEvalInterpTest::testInterpreterUnit1()
{
  INTERP_KERNEL::Unit unit1("m/s");
  INTERP_KERNEL::Unit unit2("km/h");
  CPPUNIT_ASSERT(unit1.isCompatibleWith(unit2) && unit2.isCompatibleWith(unit1));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(360,unit1.convert(unit2,100.),1e-10);
  INTERP_KERNEL::Unit unit3("J/s");
  INTERP_KERNEL::Unit unit4("kW");
  CPPUNIT_ASSERT(unit3.isCompatibleWith(unit4) && unit4.isCompatibleWith(unit3));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,unit3.convert(unit4,1000.),1e-10);
  CPPUNIT_ASSERT(unit4.getCoarseRepr()=="kW");
  INTERP_KERNEL::Unit unit5("kpT");
  CPPUNIT_ASSERT(!unit5.isInterpretationOK());
  CPPUNIT_ASSERT(unit5.getCoarseRepr()=="kpT");
  INTERP_KERNEL::Unit unit6("m*kpT");
  CPPUNIT_ASSERT(!unit6.isInterpretationOK());
  INTERP_KERNEL::Unit unit7("m*s^-1");
  CPPUNIT_ASSERT(unit7.isCompatibleWith(unit2) && unit2.isCompatibleWith(unit7));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(360,unit7.convert(unit2,100.),1e-10);
  const char unit8C[3]={-0x50,0x43,0x0};
  INTERP_KERNEL::Unit unit8(unit8C);
  INTERP_KERNEL::Unit unit9("K");
  CPPUNIT_ASSERT(unit9.isCompatibleWith(unit8) && unit8.isCompatibleWith(unit9));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(335.15,unit8.convert(unit9,62.),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-16.37,unit9.convert(unit8,256.78),1e-10);
  INTERP_KERNEL::Unit unit10("m");
  INTERP_KERNEL::Unit unit11("cm");
  CPPUNIT_ASSERT(unit10.isCompatibleWith(unit11) && unit11.isCompatibleWith(unit10));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6200.,unit10.convert(unit11,62.),1e-8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.62,unit11.convert(unit10,62.),1e-15);
  INTERP_KERNEL::Unit unit12("m-m");
  CPPUNIT_ASSERT(!unit12.isInterpretationOK());
}

void ExprEvalInterpTest::testInterpreter3()
{
  std::set<std::string> res;
  double input[3];
  double res2[3];
  INTERP_KERNEL::ExprParser expr1("2.3+x>5.");
  expr1.parse();
  expr1.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT(*(res.begin())=="x");
  expr1.prepareExprEvaluationVec();
  input[0]=0.;
  expr1.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT(-std::numeric_limits<double>::max()==res2[0]);
  input[0]=2.8;
  expr1.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT(std::numeric_limits<double>::max()==res2[0]);
  input[0]=2.6;
  expr1.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT(-std::numeric_limits<double>::max()==res2[0]);
  //
  INTERP_KERNEL::ExprParser expr2("2.3+x<5.");
  expr2.parse();
  res.clear();
  expr2.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT(*(res.begin())=="x");
  expr2.prepareExprEvaluationVec();
  input[0]=0.;
  expr2.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT(std::numeric_limits<double>::max()==res2[0]);
  input[0]=2.8;
  expr2.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT(-std::numeric_limits<double>::max()==res2[0]);
  input[0]=2.6;
  expr2.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT(std::numeric_limits<double>::max()==res2[0]);
  //
  INTERP_KERNEL::ExprParser expr3("if(2.3+x<5.,2+3*x,3+x/2)");
  expr3.parse();
  res.clear();
  expr3.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT(*(res.begin())=="x");
  expr3.prepareExprEvaluationVec();
  input[0]=0.;
  expr3.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res2[0],1e-12);
  input[0]=2.8;
  expr3.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.4,res2[0],1e-12);
  input[0]=2.6;
  expr3.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9.8,res2[0],1e-12);
  //
  INTERP_KERNEL::ExprParser expr4("if(x>1000,2*x,x/3)");
  expr4.parse();
  res.clear();
  expr4.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT(*(res.begin())=="x");
  expr4.prepareExprEvaluationVec();
  input[0]=2.7;
  expr4.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.9,res2[0],1e-12);
  input[0]=999.;
  expr4.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(333.,res2[0],1e-12);
  input[0]=1002.;
  expr4.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2004.,res2[0],1e-12);
  //
  INTERP_KERNEL::ExprParser expr5("4.4*x*log10(x)*10");
  expr5.parse();
  res.clear();
  expr5.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT(*(res.begin())=="x");
  expr5.prepareExprEvaluationVec();
  input[0]=273.15;
  expr5.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(29282.131520617437,res2[0],1e-9);
  input[0]=0.;
  CPPUNIT_ASSERT_THROW(expr5.evaluateExpr(1,input,res2),INTERP_KERNEL::Exception);
}

/*!
 * Bug detected by Marc concerning expressions with blanks.
 */
void ExprEvalInterpTest::testInterpreter4()
{
  INTERP_KERNEL::ExprParser expr("2*x + 1");
  double vals[3]={0.1,0.2,0.3};
  std::vector<std::string> varsV(3);
  varsV[0] = "x";
  varsV[1] = "y";
  varsV[2] = "z";
  expr.parse();
  expr.prepareExprEvaluation(varsV,3,1);
  double result;
  expr.evaluateExpr(1,vals, &result);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.2,result,1e-12);
}

/*!
 * Allowing scientific format for floats.
 */
void ExprEvalInterpTest::testInterpreter5()
{
  std::set<std::string> res;
  double input[3];
  double res2[3];
  INTERP_KERNEL::ExprParser expr1("1.85e-3*x");
  expr1.parse();
  expr1.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT(*(res.begin())=="x");
  input[0]=56.7;
  expr1.prepareExprEvaluationVec();
  expr1.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.104895,res2[0],1e-12);
  input[0]=-65.7;
  expr1.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.121545,res2[0],1e-12);
  //
  INTERP_KERNEL::ExprParser expr2("x*1.85e-3");
  expr2.parse();
  expr2.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT(*(res.begin())=="x");
  input[0]=56.7;
  expr2.prepareExprEvaluationVec();
  expr2.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.104895,res2[0],1e-12);
  input[0]=-65.7;
  expr2.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.121545,res2[0],1e-12);
  //
  INTERP_KERNEL::ExprParser expr3("2.6E+1+x*1.85e-3");
  expr3.parse();
  expr3.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT(*(res.begin())=="x");
  input[0]=56.7;
  expr3.prepareExprEvaluationVec();
  expr3.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(26.104895,res2[0],1e-12);
  input[0]=-65.7;
  expr3.evaluateExpr(1,input,res2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(25.878455,res2[0],1e-12);
  //
  INTERP_KERNEL::ExprParser expr4("3.*max(((3.2e+1*(ln((2*5.2E-02+6.)+(1.2E-001*1.2E+2+3e-4))))),((3.2E-2*(exp((6e-1+2*5.2e-2)+(1.2E001*1.2+3.))))))");
  expr4.parse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6994207.8359543988,expr4.evaluate(),1e-5);
}

/*!
 * Test focusing on fast evaluator.
 */
void ExprEvalInterpTest::testInterpreter6()
{
  std::vector<double> stackOfVal;
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr1("1.-2./3.");
  expr1.parse();
  expr1.prepareFastEvaluator();
  expr1.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333333,stackOfVal.back(),1e-14);
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr2("1.-2.^3.");
  expr2.parse();
  expr2.prepareFastEvaluator();
  expr2.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-7.,stackOfVal.back(),1e-14);
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr3("(7.-2.)^3.");
  expr3.parse();
  expr3.prepareFastEvaluator();
  expr3.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125.,stackOfVal.back(),1e-12);
  // now playing with one parameter - calling several times
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr4("1.2/(7.-2.*cos(x/3))");
  expr4.parse();
  expr4.prepareFastEvaluator();
  double z;
  expr4.prepareExprEvaluationDouble(std::vector<std::string>(1,"x"),1,1,0,&z,&z+1);
  expr4.prepareFastEvaluator();
  z=0.77;
  expr4.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.23689586281558844,stackOfVal.back(),1e-12);
  stackOfVal.pop_back();
  z=0.55;
  expr4.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2384018932069258,stackOfVal.back(),1e-12);
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr5("x-2*cos(y/3.)");
  expr5.parse();
  expr5.prepareFastEvaluator();
  double *aa(new double[2]);
  std::vector<std::string> vv(2); vv[0]="x"; vv[1]="y";
  expr5.prepareExprEvaluationDouble(vv,2,1,0,aa,aa+2);
  expr5.prepareFastEvaluator();
  aa[0]=0.3; aa[1]=0.5;
  expr5.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.67228646312585,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  aa[0]=0.5; aa[1]=0.3;
  expr5.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.4900083305560516,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  delete [] aa;
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr6("x*IVec-2*cos(y/3.)*JVec");
  expr6.parse();
  aa=new double[2];
  vv.resize(2); vv[0]="x"; vv[1]="y";
  expr6.prepareExprEvaluationDouble(vv,2,2,0,aa,aa+2);//emulate 1st interpreter
  expr6.prepareFastEvaluator();
  aa[0]=0.3; aa[1]=0.5;
  expr6.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  expr6.prepareExprEvaluationDouble(vv,2,2,1,aa,aa+2);//emulate 2nd interpreter
  expr6.prepareFastEvaluator();
  expr6.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.97228646312585,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  delete [] aa;
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr7("if(x>3.,-6,7.)");
  expr7.parse();
  expr7.prepareExprEvaluationDouble(std::vector<std::string>(1,"x"),1,1,0,&z,&z+1);
  expr7.prepareFastEvaluator();
  z=3.1;
  expr7.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-6.,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  z=2.8;
  expr7.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr8("if(x<3.,-6,7.)");
  expr8.parse();
  expr8.prepareExprEvaluationDouble(std::vector<std::string>(1,"x"),1,1,0,&z,&z+1);
  expr8.prepareFastEvaluator();
  z=3.1;
  expr8.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  z=2.8;
  expr8.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-6.,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr9("x*x/2");//this test is very important to test for iso priority with more than one
  expr9.parse();
  expr9.prepareExprEvaluationDouble(std::vector<std::string>(1,"x"),1,1,0,&z,&z+1);
  expr9.prepareFastEvaluator();
  z=3.;
  expr9.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.5,stackOfVal.back(),1e-14);
  stackOfVal.pop_back();
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr10("2./3./4./5.");
  expr10.parse();
  expr10.prepareFastEvaluator();
  expr10.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.03333333333333333,stackOfVal.back(),1e-13);
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr11("2./3./4.*5.");
  expr11.parse();
  expr11.prepareFastEvaluator();
  expr11.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.83333333333333333,stackOfVal.back(),1e-14);
  //
  stackOfVal.clear();
  INTERP_KERNEL::ExprParser expr12("2./3.*4.*5.");
  expr12.parse();
  expr12.prepareFastEvaluator();
  expr12.evaluateDoubleInternal(stackOfVal);
  CPPUNIT_ASSERT_EQUAL(1,(int)stackOfVal.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(13.333333333333333,stackOfVal.back(),1e-14);
}
