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

#ifndef __SINGLE_ELEMENT_PLANAR_TESTS_HXX_
#define __SINGLE_ELEMENT_PLANAR_TESTS_HXX_ 

#include "InterpKernelTestExport.hxx"
#include "InterpolationPlanarTestSuite.hxx"

namespace INTERP_TEST 
{
  /**
   * \brief Class testing algorithm by intersecting simple meshes having only one planar element each. 
   * This serves mainly to verify that the volume calculations between elements is correct.
   *
   */
  class INTERPKERNELTEST_EXPORT SingleElementPlanarTests : public InterpolationPlanarTestSuite
  {
    CPPUNIT_TEST_SUITE( SingleElementPlanarTests );
    
    CPPUNIT_TEST( diamondsBasic ); 
    CPPUNIT_TEST( tangentDiamonds );
    CPPUNIT_TEST( tangentSquares );
    CPPUNIT_TEST( diamondsSharingVertex1 );
    CPPUNIT_TEST( identicalSquares );
    CPPUNIT_TEST( squareAndDiamondBasic );
    CPPUNIT_TEST( squareAndDiamondCritical );
    CPPUNIT_TEST( diamondsCritical );
    CPPUNIT_TEST( quadranglesCritical );
    CPPUNIT_TEST( quadrangleAndDiamondCritical );
    CPPUNIT_TEST( diamondsCritical2 );
    CPPUNIT_TEST( hexagonsCritical1 );
    CPPUNIT_TEST( hexagonsCritical2 );
    CPPUNIT_TEST( squareAndQuadrangleCritical );
    CPPUNIT_TEST( diamondsSharingVertex2 );
    CPPUNIT_TEST( triangleAndDiamondCritical );
    CPPUNIT_TEST( triangleAndSquareBasic );
    CPPUNIT_TEST( trianglesCritical );
    CPPUNIT_TEST( paralellogramsCritical1 );
    CPPUNIT_TEST( paralellogramsCritical2 );
    CPPUNIT_TEST( trianglesTangencyCritical );
    CPPUNIT_TEST( trianglesTangencyCritical2 );
    CPPUNIT_TEST( trianglesTangencyCritical3 );
    CPPUNIT_TEST( trianglesTangencyCritical4 );
    CPPUNIT_TEST( diamondsBasic_Triangulation ); 
    CPPUNIT_TEST( tangentDiamonds_Triangulation );
    CPPUNIT_TEST( tangentSquares_Triangulation );
    CPPUNIT_TEST( diamondsSharingVertex1_Triangulation );
    CPPUNIT_TEST( identicalSquares_Triangulation );
    //CPPUNIT_TEST( squareAndDiamondBasic_Triangulation );
    //CPPUNIT_TEST( squareAndDiamondCritical_Triangulation );
    CPPUNIT_TEST( diamondsCritical_Triangulation );
    CPPUNIT_TEST( quadranglesCritical_Triangulation );
    CPPUNIT_TEST( quadrangleAndDiamondCritical_Triangulation );
    CPPUNIT_TEST( diamondsCritical2_Triangulation );
    CPPUNIT_TEST( hexagonsCritical1_Triangulation );
    CPPUNIT_TEST( hexagonsCritical2_Triangulation );
    CPPUNIT_TEST( squareAndQuadrangleCritical_Triangulation );
    CPPUNIT_TEST( diamondsSharingVertex2_Triangulation );
    CPPUNIT_TEST( triangleAndDiamondCritical_Triangulation );
    CPPUNIT_TEST( triangleAndSquareBasic_Triangulation );
    CPPUNIT_TEST( trianglesCritical_Triangulation );
    CPPUNIT_TEST( paralellogramsCritical1_Triangulation );
    CPPUNIT_TEST( paralellogramsCritical2_Triangulation );
    CPPUNIT_TEST( trianglesTangencyCritical_Triangulation );
    CPPUNIT_TEST( trianglesTangencyCritical2_Triangulation );
    CPPUNIT_TEST( trianglesTangencyCritical3_Triangulation );
    CPPUNIT_TEST( trianglesTangencyCritical4_Triangulation );

    CPPUNIT_TEST_SUITE_END();
    
  public:

    void diamondsBasic();
    void tangentDiamonds();
    void tangentSquares();
    void diamondsSharingVertex1();
    void identicalSquares();
    void squareAndDiamondBasic();
    void squareAndDiamondCritical();
    void diamondsCritical();
    void quadranglesCritical();  
    void quadrangleAndDiamondCritical();
    void diamondsCritical2();
    void hexagonsCritical1();
    void hexagonsCritical2();
    void squareAndQuadrangleCritical();
    void diamondsSharingVertex2();
    void triangleAndDiamondCritical();
    void triangleAndSquareBasic();
    void trianglesCritical();
    void paralellogramsCritical1();
    void paralellogramsCritical2();
    void trianglesTangencyCritical();
    void trianglesTangencyCritical2();
    void trianglesTangencyCritical3();
    void trianglesTangencyCritical4();
    void diamondsBasic_Triangulation();
    void tangentDiamonds_Triangulation(); 
    void tangentSquares_Triangulation();
    void diamondsSharingVertex1_Triangulation();
    void identicalSquares_Triangulation();
    void squareAndDiamondBasic_Triangulation();
    void squareAndDiamondCritical_Triangulation();
    void diamondsCritical_Triangulation();
    void quadranglesCritical_Triangulation();  
    void quadrangleAndDiamondCritical_Triangulation();
    void diamondsCritical2_Triangulation();
    void hexagonsCritical1_Triangulation();
    void hexagonsCritical2_Triangulation();
    void squareAndQuadrangleCritical_Triangulation();
    void diamondsSharingVertex2_Triangulation();
    void triangleAndDiamondCritical_Triangulation();
    void triangleAndSquareBasic_Triangulation();
    void trianglesCritical_Triangulation();
    void paralellogramsCritical1_Triangulation();
    void paralellogramsCritical2_Triangulation();
    void trianglesTangencyCritical_Triangulation();
    void trianglesTangencyCritical2_Triangulation();
    void trianglesTangencyCritical3_Triangulation();
    void trianglesTangencyCritical4_Triangulation();
  };
}
#endif
