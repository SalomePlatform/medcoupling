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

#include "QuadraticPlanarInterpTest.hxx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "InterpKernelGeo2DElementaryEdge.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"

#include <cmath>
#include <sstream>
#include <iostream>
#include <iterator>

using namespace INTERP_KERNEL;

namespace INTERP_TEST
{

class DoubleEqual
{
public:
  DoubleEqual(double eps):_eps(eps) { }
  bool operator()(double x, double y) { return fabs(x-y)<_eps; }
private:
  double _eps;
};

void QuadraticPlanarInterpTest::checkNonRegressionOmar0000()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.383022221559489, 0.3213938048432697, -0.5745333323392334, 0.4820907072649046, 0.5745333323392335, 0.4820907072649044, 0.383022221559489, 0.3213938048432696,
    -0.4787777769493612, 0.4017422560540872, 4.592273826833915e-17, 0.75, 0.4787777769493612, 0.401742256054087, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -0.383022221559489, -0.1786061951567303, -0.5745333323392334, -0.01790929273509539, 0.5745333323392335, -0.01790929273509556, 0.383022221559489, -0.1786061951567304,
    -0.4787777769493612, -0.0982577439459128, 4.592273826833915e-17, 0.25, 0.4787777769493612, -0.09825774394591297, 3.061515884555943e-17, 0 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0001()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.383022221559489, 0.3213938048432697, -0.5745333323392334, 0.4820907072649046, 0.5745333323392335, 0.4820907072649044, 0.383022221559489, 0.3213938048432696,
    -0.4787777769493612, 0.4017422560540872, 4.592273826833915e-17, 0.75, 0.4787777769493612, 0.401742256054087, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -0.383022221559489, 0.3213938048432697, -0.5745333323392334, 0.4820907072649046, 0.5745333323392335, 0.4820907072649044, 0.383022221559489, 0.3213938048432696,
    -0.4787777769493612, 0.4017422560540872, 4.592273826833915e-17, 0.75, 0.4787777769493612, 0.401742256054087, 3.061515884555943e-17, 0.5 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.272708,pol1->intersectWith(*pol2),1.e-6);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.272708,pol2->intersectWith(*pol1),1.e-6);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0002()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.383022221559489, 0.3213938048432697, -0.5745333323392334, 0.4820907072649046, 0.5745333323392335, 0.4820907072649044, 0.383022221559489, 0.3213938048432696,
    -0.4787777769493612, 0.4017422560540872, 4.592273826833915e-17, 0.75, 0.4787777769493612, 0.401742256054087, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -0.4979288880273356, 0.4178119462962507, -0.6128355544951823, 0.5142300877492316, 0.6128355544951825, 0.5142300877492314, 0.4979288880273357, 0.4178119462962505,
    -0.555382221261259, 0.4660210170227412, 4.898425415289509e-17, 0.8, 0.5553822212612591, 0.466021017022741, 3.979970649922726e-17, 0.65 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.122173,pol1->intersectWith(*pol2),1.e-6);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.122173,pol2->intersectWith(*pol1),1.e-6);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0003()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.3535533905932737, 0.3535533905932738, -0.5303300858899106, 0.5303300858899107, 0.5303300858899107, 0.5303300858899106, 0.3535533905932738, 0.3535533905932737,
    -0.4419417382415922, 0.4419417382415922, 4.592273826833915e-17, 0.75, 0.4419417382415922, 0.4419417382415922, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -0.4979288880273356, 0.4178119462962507, -0.6128355544951823, 0.5142300877492316, 0.6128355544951825, 0.5142300877492314, 0.4979288880273357, 0.4178119462962505,
    -0.555382221261259, 0.4660210170227412, 4.898425415289509e-17, 0.8, 0.5553822212612591, 0.466021017022741, 3.979970649922726e-17, 0.65 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.109956,pol1->intersectWith(*pol2),1.e-6);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.109956,pol2->intersectWith(*pol1),1.e-6);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0004()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.4596194077712559, 0.4596194077712559, -0.5303300858899106, 0.5303300858899107, 0.5303300858899107, 0.5303300858899106, 0.4596194077712559, 0.4596194077712559,
    -0.4949747468305832, 0.4949747468305833, 4.592273826833915e-17, 0.75, 0.4949747468305833, 0.4949747468305832, 3.979970649922726e-17, 0.65 };

  double coords2[16]={
    -0.383022221559489, 0.3213938048432697, -0.6128355544951823, 0.5142300877492316, 0.6128355544951825, 0.5142300877492314, 0.383022221559489, 0.3213938048432696,
    -0.4979288880273356, 0.4178119462962507, 4.898425415289509e-17, 0.8, 0.4979288880273357, 0.4178119462962505, 3.061515884555943e-17, 0.5 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.109956,pol1->intersectWith(*pol2),1.e-6);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.109956,pol2->intersectWith(*pol1),1.e-6);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0005()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.383022221559489, 0.3213938048432697, -0.6128355544951823, 0.5142300877492316, 0.6128355544951825, 0.5142300877492314, 0.383022221559489, 0.3213938048432696,
    -0.4979288880273356, 0.4178119462962507, 4.898425415289509e-17, 0.8, 0.4979288880273357, 0.4178119462962505, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -0.4596194077712559, 0.4596194077712559, -0.5303300858899106, 0.5303300858899107, 0.5303300858899107, 0.5303300858899106, 0.4596194077712559, 0.4596194077712559,
    -0.4949747468305832, 0.4949747468305833, 4.592273826833915e-17, 0.75, 0.4949747468305833, 0.4949747468305832, 3.979970649922726e-17, 0.65 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.109956,pol1->intersectWith(*pol2),1.e-6);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.109956,pol2->intersectWith(*pol1),1.e-6);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0006()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.383022221559489, 0.3213938048432697, -0.5362311101832845, 0.4499513267805776, 0.5362311101832846, 0.4499513267805774, 0.383022221559489, 0.3213938048432696,
    -0.4596266658713867, 0.3856725658119237, 4.28612223837832e-17, 0.7, 0.4596266658713868, 0.3856725658119236, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -0.1811733315717646, 0.6761480784023478, -0.2070552360820167, 0.7727406610312547, 0.2070552360820166, 0.7727406610312547, 0.1811733315717645, 0.6761480784023478,
    -0.1941142838268906, 0.7244443697168013, 4.898425415289509e-17, 0.8, 0.1941142838268906, 0.7244443697168013, 4.28612223837832e-17, 0.7 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.366519,0.,0.};
  double test2_res[4]={0.,0.,0.,0.366519};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-6)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-6)));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0007()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.383022221559489, 0.3213938048432697, -0.5362311101832845, 0.4499513267805776, 0.5362311101832846, 0.4499513267805774, 0.383022221559489, 0.3213938048432696,
    -0.4596266658713867, 0.3856725658119237, 4.28612223837832e-17, 0.7, 0.4596266658713868, 0.3856725658119236, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -0.4499513267805775, 0.5362311101832846, -0.5142300877492315, 0.6128355544951825, -0.1389185421335442, 0.7878462024097664, -0.1215537243668512, 0.6893654271085455,
    -0.4820907072649045, 0.5745333323392335, -0.3380946093925595, 0.7250462296293201, -0.1302361332501977, 0.738605814759156, -0.2958327832184895, 0.634415450925655 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.366519,0.,0.};
  double test2_res[4]={0.,0.,0.,0.366519};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-6)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-6)));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0008()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.383022221559489, 0.3213938048432697, -0.5362311101832845, 0.4499513267805776, 0.5362311101832846, 0.4499513267805774, 0.383022221559489, 0.3213938048432696,
    -0.4596266658713867, 0.3856725658119237, 4.28612223837832e-17, 0.7, 0.4596266658713868, 0.3856725658119236, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -0.6344154509256549, 0.2958327832184896, -0.72504622962932, 0.3380946093925596, -0.4588611490808367, 0.6553216354311937, -0.401503505445732, 0.5734064310022944,
    -0.6797308402774874, 0.3169636963055246, -0.6128355544951823, 0.5142300877492316, -0.4301823272632844, 0.614364033216744, -0.5362311101832845, 0.4499513267805776 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.18326,0.,0.};
  double test2_res[4]={0.,0.,0.,0.18326};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-5)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-5)));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0009()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.3863703305156274, -0.1035276180410081, -0.4829629131445342, -0.1294095225512602, 0.4829629131445342, -0.1294095225512604, 0.3863703305156274, -0.1035276180410083,
    -0.4346666218300808, -0.1164685702961342, 1.416374613080751e-16, 0.5, 0.4346666218300808, -0.1164685702961343, 1.133099690464601e-16, 0.4 };
  double coords2[16]={
    0.5, -1.224606353822377e-16, 0.6, -1.469527624586853e-16, -0.6, 7.347638122934263e-17, -0.5, 6.123031769111886e-17,
    0.55, -1.347066989204615e-16, -1.102145718440139e-16, -0.6, -0.55, 6.735334946023075e-17, -9.184547653667829e-17, -0.5 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0010()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.3863703305156274, -0.1035276180410081, -0.4829629131445342, -0.1294095225512602, 0.4829629131445342, -0.1294095225512604, 0.3863703305156274, -0.1035276180410083,
-0.4346666218300808, -0.1164685702961342, 1.416374613080751e-16, 0.5, 0.4346666218300808, -0.1164685702961343, 1.133099690464601e-16, 0.4 };
  double coords2[16]={
    0.4346666218300808, -0.1164685702961343, 0.579555495773441, -0.1552914270615124, -0.579555495773441, -0.1552914270615122, -0.4346666218300808, -0.1164685702961342,
0.5071110588017609, -0.1358799986788234, -1.102145718440139e-16, -0.6, -0.507111058801761, -0.1358799986788232, -8.266092888301047e-17, -0.45 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0011()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.3863703305156274, -0.1035276180410081, -0.4829629131445342, -0.1294095225512602, 0.4829629131445342, -0.1294095225512604, 0.3863703305156274, -0.1035276180410083,
-0.4346666218300808, -0.1164685702961342, 1.416374613080751e-16, 0.5, 0.4346666218300808, -0.1164685702961343, 1.133099690464601e-16, 0.4 };
  double coords2[16]={
    0.4829629131445342, -0.1294095225512603, 0.579555495773441, -0.1552914270615124, -0.579555495773441, -0.1552914270615122, -0.4829629131445342, -0.1294095225512602,
0.5312592044589877, -0.1423504748063864, -1.102145718440139e-16, -0.6, -0.5312592044589877, -0.1423504748063862, -9.184547653667829e-17, -0.5 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  double val1,val2,val3;
  pol1->intersectForPerimeter(*pol2,val1,val2,val3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val1,1.e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val2,1.e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val3,1.e-13);
  std::vector<double> val4,val5;
  pol1->intersectForPerimeterAdvanced(*pol2,val4,val5);
  double test1_res[4]={0.,0.,0.,0.};
  CPPUNIT_ASSERT(std::equal(val4.begin(),val4.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val5.begin(),val5.end(),test1_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  pol1->intersectForPerimeter(*pol2,val1,val2,val3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val1,1.e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val2,1.e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val3,1.e-13);
  val4.clear(); val5.clear();
  pol1->intersectForPerimeterAdvanced(*pol2,val4,val5);
  CPPUNIT_ASSERT(std::equal(val4.begin(),val4.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val5.begin(),val5.end(),test1_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar2511()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.3863703305156274, -0.1035276180410081, -0.4829629131445342, -0.1294095225512602, 0.4829629131445342, -0.1294095225512604, 0.3863703305156274, -0.1035276180410083,
    -0.4346666218300808, -0.1164685702961342, 1.416374613080751e-16, 0.5, 0.4346666218300808, -0.1164685702961343, 1.133099690464601e-16, 0.4, };
  
  double coords2[16]={
    0.579555495773441, -0.1552914270615124, -0.579555495773441, -0.1552914270615122, -0.4829629131445342, -0.1294095225512602, 0.4829629131445342, -0.1294095225512603,
    -1.102145718440139e-16, -0.6, -0.5312592044589877, -0.1423504748063862, -9.184547653667829e-17, -0.5, 0.5312592044589877, -0.1423504748063864, };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  double val1,val2,val3;
  pol1->intersectForPerimeter(*pol2,val1,val2,val3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val1,1.e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val2,1.e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val3,1.e-13);
  std::vector<double> val4,val5;
  pol1->intersectForPerimeterAdvanced(*pol2,val4,val5);
  double test1_res[4]={0.,0.,0.,0.};
  CPPUNIT_ASSERT(std::equal(val4.begin(),val4.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val5.begin(),val5.end(),test1_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  pol1->intersectForPerimeter(*pol2,val1,val2,val3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val1,1.e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val2,1.e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,val3,1.e-13);
  val4.clear(); val5.clear();
  pol1->intersectForPerimeterAdvanced(*pol2,val4,val5);
  CPPUNIT_ASSERT(std::equal(val4.begin(),val4.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val5.begin(),val5.end(),test1_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0012()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -1, 1.224606353822377e-16, -1.6, 1.959370166115804e-16, 9.796850830579018e-17, 1.6, 6.123031769111886e-17, 1,
    -1.3, 1.591988259969091e-16, -1.131370849898476, 1.131370849898476, 7.959941299845453e-17, 1.3, -0.7071067811865475, 0.7071067811865476 };
  
  double coords2[16]={
    6.123031769111886e-18, 1.85, 1.224606353822377e-17, 1.95, 1.224606353822377e-17, 1.55, 6.123031769111886e-18, 1.65,
    9.18454765366783e-18, 1.9, 0.2, 1.75, 9.18454765366783e-18, 1.6, 0.1, 1.75 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.05,0.};
  double test2_res[4]={0.,0.,0.05,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,1,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0013()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -1, 1.224606353822377e-16, -1.6, 1.959370166115804e-16, 9.796850830579018e-17, 1.6, 6.123031769111886e-17, 1,
    -1.3, 1.591988259969091e-16, -1.131370849898476, 1.131370849898476, 7.959941299845453e-17, 1.3, -0.7071067811865475, 0.7071067811865476 };
  
  double coords2[16]={
    6.123031769111886e-18, 1.7, 1.224606353822377e-17, 1.8, 1.224606353822377e-17, 1.4, 6.123031769111886e-18, 1.5,
    9.18454765366783e-18, 1.75, 0.2, 1.6, 9.18454765366783e-18, 1.45, 0.1, 1.6 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.1,0.};
  double test2_res[4]={0.,0.,0.1,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,2,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0014()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -1, 1.224606353822377e-16, -1.6, 1.959370166115804e-16, 9.796850830579018e-17, 1.6, 6.123031769111886e-17, 1,
-1.3, 1.591988259969091e-16, -1.131370849898476, 1.131370849898476, 7.959941299845453e-17, 1.3, -0.7071067811865475, 0.7071067811865476 };
  double coords2[16]={
    6.123031769111886e-18, 1.55, 1.224606353822377e-17, 1.65, 1.224606353822377e-17, 1.25, 6.123031769111886e-18, 1.35,
9.18454765366783e-18, 1.6, 0.2, 1.45, 9.18454765366783e-18, 1.3, 0.1, 1.45 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.15,0.};
  double test2_res[4]={0.05,0.,0.1,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,3,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0015()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -1, 1.224606353822377e-16, -1.6, 1.959370166115804e-16, 9.796850830579018e-17, 1.6, 6.123031769111886e-17, 1,
-1.3, 1.591988259969091e-16, -1.131370849898476, 1.131370849898476, 7.959941299845453e-17, 1.3, -0.7071067811865475, 0.7071067811865476 };
  double coords2[16]={
    6.123031769111886e-18, 1.4, 1.224606353822377e-17, 1.5, 1.224606353822377e-17, 1.1, 6.123031769111886e-18, 1.2,
9.18454765366783e-18, 1.45, 0.2, 1.3, 9.18454765366783e-18, 1.15, 0.1, 1.3 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.2,0.};
  double test2_res[4]={0.1,0.,0.1,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,4,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0016()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -1, 1.224606353822377e-16, -1.6, 1.959370166115804e-16, 9.796850830579018e-17, 1.6, 6.123031769111886e-17, 1,
-1.3, 1.591988259969091e-16, -1.131370849898476, 1.131370849898476, 7.959941299845453e-17, 1.3, -0.7071067811865475, 0.7071067811865476 };
  double coords2[16]={
    6.123031769111886e-18, 1.25, 1.224606353822377e-17, 1.35, 1.224606353822377e-17, 0.95, 6.123031769111886e-18, 1.05,
9.18454765366783e-18, 1.3, 0.2, 1.15, 9.18454765366783e-18, 0.9999999999999999, 0.1, 1.15 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.15,0.};
  double test2_res[4]={0.1,0.,0.05,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,3,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0017()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -1, 1.224606353822377e-16, -1.6, 1.959370166115804e-16, 9.796850830579018e-17, 1.6, 6.123031769111886e-17, 1,
    -1.3, 1.591988259969091e-16, -1.131370849898476, 1.131370849898476, 7.959941299845453e-17, 1.3, -0.7071067811865475, 0.7071067811865476 };
  
  double coords2[16]={
    6.123031769111886e-18, 1.1, 1.224606353822377e-17, 1.2, 1.224606353822377e-17, 0.8, 6.123031769111886e-18, 0.9,
    9.18454765366783e-18, 1.15, 0.2, 1, 9.18454765366783e-18, 0.85, 0.1, 1 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.1,0.};
  double test2_res[4]={0.1,0.,0.,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,2,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0018()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -1, 1.224606353822377e-16, -1.6, 1.959370166115804e-16, 9.796850830579018e-17, 1.6, 6.123031769111886e-17, 1,
    -1.3, 1.591988259969091e-16, -1.131370849898476, 1.131370849898476, 7.959941299845453e-17, 1.3, -0.7071067811865475, 0.7071067811865476 };
  
  double coords2[16]={
    6.123031769111886e-18, 0.95, 1.224606353822377e-17, 1.05, 1.224606353822377e-17, 0.6499999999999999, 6.123031769111886e-18, 0.75,
    9.18454765366783e-18, 1, 0.2, 0.85, 9.18454765366783e-18, 0.7, 0.1, 0.85 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.05,0.};
  double test2_res[4]={0.05,0.,0.,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,1,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0019()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.5, 6.123031769111886e-17, -0.8, 9.796850830579018e-17, 0.8, 0, 0.5, 0,
    -0.65, 7.959941299845453e-17, 4.898425415289509e-17, 0.8, 0.65, 0, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    0.9500000000000001, 1.836909530733566e-17, 0.8, 3.673819061467131e-17, 1.4, 0, 1.25, 0,
    0.8750000000000001, 2.755364296100349e-17, 1.1, 0.3, 1.325, 0, 1.1, 0.15 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0020()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.5, 6.123031769111886e-17, -0.8, 9.796850830579018e-17, 0.8, 0, 0.5, 0,
    -0.65, 7.959941299845453e-17, 4.898425415289509e-17, 0.8, 0.65, 0, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    0.05000000000000002, 1.836909530733566e-17, -0.09999999999999998, 3.673819061467131e-17, 0.5, 0, 0.35, 0,
    -0.02499999999999997, 2.755364296100349e-17, 0.2, 0.3, 0.425, 0, 0.2, 0.15 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.,0.};
  double test2_res[4]={0.,0.,0.,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-6)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-6)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,0,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0021()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.5, 6.123031769111886e-17, -0.8, 9.796850830579018e-17, 0.8, 0, 0.5, 0,
    -0.65, 7.959941299845453e-17, 4.898425415289509e-17, 0.8, 0.65, 0, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    -1, -0.07999999999999999, -1.15, -0.07999999999999996, -0.55, -0.08, -0.7, -0.08,
    -1.075, -0.07999999999999997, -0.85, 0.22, -0.625, -0.08, -0.85, 0.06999999999999999 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0311485,pol1->intersectWith(*pol2),1.e-7);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0311485,pol2->intersectWith(*pol1),1.e-7);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.162251,0.151523,0.,0.};
  double test2_res[4]={0.,0.311383,0.,0.0978193};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-6)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-6)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={2,2,0,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}
void QuadraticPlanarInterpTest::checkNonRegressionOmar0022()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.5, 6.123031769111886e-17, -0.8, 9.796850830579018e-17, 0.8, 0, 0.5, 0,
    -0.65, 7.959941299845453e-17, 4.898425415289509e-17, 0.8, 0.65, 0, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    0.15, -0.07999999999999999, 0, -0.07999999999999996, 0.6, -0.08, 0.45, -0.08,
    0.07500000000000001, -0.07999999999999997, 0.3, 0.22, 0.5249999999999999, -0.08, 0.3, 0.06999999999999999 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00902229,pol1->intersectWith(*pol2),1.e-8);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00902229,pol2->intersectWith(*pol1),1.e-8);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0023()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.5, 6.123031769111886e-17, -0.8, 9.796850830579018e-17, 0.8, 0, 0.5, 0,
    -0.65, 7.959941299845453e-17, 4.898425415289509e-17, 0.8, 0.65, 0, 3.061515884555943e-17, 0.5, };
  
  double coords2[16]={
    0.4156854249492381, 0.5656854249492381, 0.2656854249492381, 0.5656854249492381, 0.8656854249492381, 0.5656854249492381, 0.7156854249492381, 0.5656854249492381,
    0.3406854249492381, 0.5656854249492381, 0.5656854249492381, 0.8656854249492381, 0.7906854249492381, 0.5656854249492381, 0.5656854249492381, 0.7156854249492381 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0215659,pol1->intersectWith(*pol2),1.e-7);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0215659,pol2->intersectWith(*pol1),1.e-7);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0024()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.5, 6.123031769111886e-17, -0.8, 9.796850830579018e-17, 0.8, 0, 0.5, 0,
-0.65, 7.959941299845453e-17, 4.898425415289509e-17, 0.8, 0.65, 0, 3.061515884555943e-17, 0.5 };
  double coords2[16]={
    0.5656854249492381, 0.5656854249492381, 0.4156854249492382, 0.5656854249492381, 1.015685424949238, 0.5656854249492381, 0.8656854249492382, 0.5656854249492381,
0.4906854249492382, 0.5656854249492381, 0.7156854249492381, 0.8656854249492381, 0.9406854249492381, 0.5656854249492381, 0.7156854249492381, 0.7156854249492381 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00877657,pol1->intersectWith(*pol2),1.e-8);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00877657,pol2->intersectWith(*pol1),1.e-8);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar2524()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.5, 6.123031769111886e-17, -0.8, 9.796850830579018e-17, 0.8, 0, 0.5, 0,
-0.65, 7.959941299845453e-17, 4.898425415289509e-17, 0.8, 0.65, 0, 3.061515884555943e-17, 0.5 };
  double coords2[16]={
    0.4156854249492382, 0.5656854249492381, 1.015685424949238, 0.5656854249492381, 0.8656854249492382, 0.5656854249492381, 0.5656854249492381, 0.5656854249492381,
0.7156854249492381, 0.8656854249492381, 0.9406854249492381, 0.5656854249492381, 0.7156854249492381, 0.7156854249492381, 0.4906854249492382, 0.5656854249492381 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00877657,pol1->intersectWith(*pol2),1.e-8);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00877657,pol2->intersectWith(*pol1),1.e-8);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0025()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.5, 6.123031769111886e-17, -0.8, 9.796850830579018e-17, 0.8, 0, 0.5, 0,
    -0.65, 7.959941299845453e-17, 4.898425415289509e-17, 0.8, 0.65, 0, 3.061515884555943e-17, 0.5 };
  
  double coords2[16]={
    0.715685424949238, 0.5656854249492381, 0.565685424949238, 0.5656854249492381, 1.165685424949238, 0.5656854249492381, 1.015685424949238, 0.5656854249492381,
    0.6406854249492381, 0.5656854249492381, 0.8656854249492381, 0.8656854249492381, 1.090685424949238, 0.5656854249492381, 0.8656854249492381, 0.7156854249492381 };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,1,0,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0026()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.4, 4.898425415289509e-17, -0.75, 9.184547653667829e-17, 0.75, 0, 0.4, 0,
    -0.575, 7.041486534478669e-17, 4.592273826833915e-17, 0.75, 0.575, 0, 2.449212707644755e-17, 0.4 };
  
  double coords2[16]={
    0.1, 0.95, 0.2, 0.95, -0.2, 0.95, -0.1, 0.95,
    0.15, 0.95, 1.224606353822377e-17, 0.75, -0.15, 0.95, 6.123031769111886e-18, 0.85 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,1,0,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0027()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.4, 4.898425415289509e-17, -0.75, 9.184547653667829e-17, 0.75, 0, 0.4, 0,
    -0.575, 7.041486534478669e-17, 4.592273826833915e-17, 0.75, 0.575, 0, 2.449212707644755e-17, 0.4 };
  
  double coords2[16]={
    -0.1, 0.7, -0.2, 0.7, 0.2, 0.7, 0.1, 0.7,
    -0.15, 0.7, 1.224606353822377e-17, 0.8999999999999999, 0.15, 0.7, 6.123031769111886e-18, 0.7999999999999999 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00712309,pol1->intersectWith(*pol2),1.e-8);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00712309,pol2->intersectWith(*pol1),1.e-8);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.222704,0.,0.};
  double test2_res[4]={0.1,0.0465335,0.1,0.092554};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-6)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-6)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,4,0,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0028()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.4, 4.898425415289509e-17, -0.75, 9.184547653667829e-17, 0.75, 0, 0.4, 0,
    -0.575, 7.041486534478669e-17, 4.592273826833915e-17, 0.75, 0.575, 0, 2.449212707644755e-17, 0.4 };
  
  double coords2[16]={
    -0.07071067811865477, 0.4792893218813453, -0.1414213562373095, 0.4085786437626905, 0.1414213562373095, 0.6914213562373095, 0.07071067811865477, 0.6207106781186548,
    -0.1060660171779822, 0.4439339828220179, -0.1414213562373095, 0.6914213562373096, 0.1060660171779822, 0.6560660171779822, -0.07071067811865475, 0.6207106781186548 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0471239,pol1->intersectWith(*pol2),1.e-7);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0471239,pol2->intersectWith(*pol1),1.e-7);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.,0.};
  double test2_res[4]={0.1,0.628319,0.1,0.314159};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-6)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-6)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,1,0,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0029()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.4, 4.898425415289509e-17, -0.75, 9.184547653667829e-17, 0.75, 0, 0.4, 0,
    -0.575, 7.041486534478669e-17, 4.592273826833915e-17, 0.75, 0.575, 0, 2.449212707644755e-17, 0.4 };
  
  double coords2[16]={
    -0.07071067811865477, 0.1292893218813453, -0.1414213562373095, 0.05857864376269051, 0.1414213562373095, 0.3414213562373095, 0.07071067811865477, 0.2707106781186548,
    -0.1060660171779822, 0.09393398282201787, -0.1414213562373095, 0.3414213562373095, 0.1060660171779822, 0.3060660171779822, -0.07071067811865475, 0.2707106781186548 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.,0.};
  double test2_res[4]={0.,0.,0.,0.};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-13)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-13)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,0,0,1};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegressionOmar0030()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -0.4, 4.898425415289509e-17, -0.75, 9.184547653667829e-17, 0.75, 0, 0.4, 0,
    -0.575, 7.041486534478669e-17, 4.592273826833915e-17, 0.75, 0.575, 0, 2.449212707644755e-17, 0.4 };
  
  double coords2[16]={
    -0.4889087296526012, 0.3889087296526012, -0.5889087296526012, 0.3889087296526012, -0.1889087296526012, 0.3889087296526012, -0.2889087296526012, 0.3889087296526012,
    -0.5389087296526012, 0.3889087296526012, -0.3889087296526012, 0.5889087296526012, -0.2389087296526012, 0.3889087296526012, -0.3889087296526012, 0.4889087296526012 };
  
  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0471239,pol1->intersectWith(*pol2),1.e-7);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0471239,pol2->intersectWith(*pol1),1.e-7);
  delete pol1;
  delete pol2;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  std::vector<double> val1,val2;
  pol1->intersectForPerimeterAdvanced(*pol2,val1,val2);
  double test1_res[4]={0.,0.,0.,0.};
  double test2_res[4]={0.1,0.628319,0.1,0.314159};
  CPPUNIT_ASSERT(std::equal(val1.begin(),val1.end(),test1_res,DoubleEqual(1e-6)));
  CPPUNIT_ASSERT(std::equal(val2.begin(),val2.end(),test2_res,DoubleEqual(1e-6)));
  delete pol1;
  delete pol2;
  std::vector<int> val3;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  pol1->intersectForPoint(*pol2,val3);
  int test3_res[4]={0,1,0,0};
  CPPUNIT_ASSERT(std::equal(val3.begin(),val3.end(),test3_res));
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkIsInOrOut()
{
  double coords[8]={   0.30662641093707971,  -0.47819928619088981,
                      -0.47819928619088964,  0.30662641093707987,
                       0.0, 0.0,
                       0.4, 0.4
  };
  coords[4] = (coords[0] + coords[2]) / 2.0;
  coords[5] = (coords[1] + coords[3]) / 2.0;

  int tab4[4]={ 0, 1, 2, 3};
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab4,4);
  Node * n = new Node(0.3175267678416348, -0.4890996430954449);

  CPPUNIT_ASSERT(! pol1->isInOrOut(n)); // node should be out
  n->decrRef();
  delete pol1;
}

void QuadraticPlanarInterpTest::checkGetMiddleOfPoints()
{
  { // from testIntersect2DMeshWith1DLine6()
    double p1[] = {0.51641754716735844, 2.0};
    double p2[] = {0.0, 1.0};
    double e_center[] = {-0.71, 2.0};
    double mid[] = {0.0,0.0}; // out
    double mide[] = {0.0,0.0}; // expected

    Node * start = new Node(0.,0.); Node * end = new Node(0.,0.); // unused
    // start, end, center_x, center_y, radius, angle0, angle
    EdgeArcCircle e(start, end, e_center, 1.2264175471673588, -0.9533904350433241, 0.95339043504332388);

    e.getMiddleOfPoints(p1, p2, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.37969180470645592, mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4372640310451197, mid[1], 1.e-7);

    e.getMiddleOfPoints(p2, p1, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.37969180470645592, mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4372640310451197, mid[1], 1.e-7);

    start->decrRef(); end->decrRef();
  }
  { // from testSwig2Intersect2DMeshWith1DLine11()
    double p1[] = {-1., 0.23453685964236054};
    double p2[] = {-0.23453685964235979, 1.0};
    double e_center[] = {-4.85, 4.85};
    double mid[] = {0.0,0.0}; // out

    Node * start = new Node(0.,0.); Node * end = new Node(0.,0.); // unused
    // start, end, center_x, center_y, radius, angle0, angle
    EdgeArcCircle e(start, end, e_center, 6.0104076400856474, -0.69522150912422953, -0.18035330854643861);

    e.getMiddleOfPoints(p1, p2, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.6, mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.6, mid[1], 1.e-7);

    e.getMiddleOfPoints(p2, p1, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.6, mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.6, mid[1], 1.e-7);

    start->decrRef(); end->decrRef();
  }
  { // from testSwig2Intersect2DMeshWith1DLine11()
    double p1[] = {-0.1303327636866019, -1.0};
    double p2[] = {-1.0, -0.1303327636866019};
    double e_center[] = {-1.9833333333333298, -1.9833333333333298};
    double mid[] = {0.0,0.0}; // out

    Node * start = new Node(0.,0.); Node * end = new Node(0.,0.); // unused
    // start, end, center_x, center_y, radius, angle0, angle
    EdgeArcCircle e(start, end, e_center, 2.0977501175200861, 1.0829141821052615, -0.59503203741562627);

    e.getMiddleOfPoints(p1, p2, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, mid[1], 1.e-7);

    e.getMiddleOfPoints(p2, p1, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, mid[1], 1.e-7);

    start->decrRef(); end->decrRef();
  }
}

void QuadraticPlanarInterpTest::checkGetMiddleOfPointsOriented()
{
  { // from testSwig2Colinearize2D3()
    double p1[] = {-0.70710678118654746, 0.70710678118654757};
    double p2[] = {-0.70710678118654768, -0.70710678118654746};
    double e_center[] = {0., 0.};
    double mid[] = {0.0,0.0}; // out

    Node * start = new Node(0.,0.); Node * end = new Node(0.,0.); // unused
    // start, end, center_x, center_y, radius, angle0, angle
    EdgeArcCircle e(start, end, e_center, 1.0, -0.7853981633974485, -1.5707963267948966);

    e.getMiddleOfPointsOriented(p1, p2, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1., mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., mid[1], 1.e-7);

    e.getMiddleOfPoints(p1, p2, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1., mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., mid[1], 1.e-7);

    e.getMiddleOfPointsOriented(p2, p1, mid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1., mid[0], 1.e-7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., mid[1], 1.e-7);

    start->decrRef(); end->decrRef();
  }
}

}
