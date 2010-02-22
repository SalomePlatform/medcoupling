//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#include "ExprEvalInterpTest.hxx"
#include "InterpKernelExprParser.hxx"

using namespace std;
using namespace INTERP_TEST;

void ExprEvalInterpTest::testBuildStringFromFortran()
{
  char toto1[]="123456  ";
  char result[]="123456";
  string titi;
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
  string totoS(toto);
  string totoR=INTERP_KERNEL::ExprParser::deleteWhiteSpaces(totoS);
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
  expr1.prepareExprEvaluation(vars);
  expr1.evaluateExpr(1,&res1,xyValue);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(33.833333333333336,res1,1e-13);
  xyValue[0]=-2.;
  CPPUNIT_ASSERT_THROW(expr1.evaluateExpr(1,&res1,xyValue),INTERP_KERNEL::Exception);
  double res2[2];
  xyValue[0]=1.;
  expr1.evaluateExpr(2,res2,xyValue);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(33.833333333333336,res2[0],1e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(33.833333333333336,res2[1],1e-13);
  INTERP_KERNEL::ExprParser expr2("3.5*tan(2.3*x)*IVec+(cos(1.2+y/x)*JVec)");
  expr2.parse();
  res.clear(); expected.clear();
  expr2.getSetOfVars(res);
  CPPUNIT_ASSERT_EQUAL(4,(int)res.size());
  expected.insert("x"); expected.insert("y"); expected.insert("IVec"); expected.insert("JVec");
  CPPUNIT_ASSERT(std::equal(res.begin(),res.end(),expected.begin()));
  expr2.prepareExprEvaluation(vars);
  expr2.evaluateExpr(2,res2,xyValue);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.9172477460694637,res2[0],1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.49026082134069943,res2[1],1e-14);
  INTERP_KERNEL::ExprParser expr3("3.5*u+u^2.4+2.");
  expr3.parse();
  expr3.prepareExprEvaluationVec();
  expr3.evaluateExpr(2,res2,xyValue);
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
