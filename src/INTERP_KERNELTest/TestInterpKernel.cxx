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

#include "CppUnitTest.hxx"
#include "BBTreeTest.hxx"
#include "ExprEvalInterpTest.hxx"
#include "QuadraticPlanarInterpTest.hxx"
#include "SingleElementPlanarTests.hxx"
#include "TransformedTriangleIntersectTest.hxx"
#include "TransformedTriangleTest.hxx"
#include "UnitTetraIntersectionBaryTest.hxx"
#include "UnitTetra3D2DIntersectionTest.hxx"

#ifndef MEDCOUPLING_MICROMED
#include "HexaTests.hxx"
#include "InterpolationOptionsTest.hxx"
#include "MultiElement2DTests.hxx"
#include "MultiElementTetraTests.hxx"
#include "SingleElementTetraTests.hxx"
#include "ThreeDSurfProjectionTest.hxx"
#endif

using namespace INTERP_TEST;

//--- Registers the fixture into the 'registry'

CPPUNIT_TEST_SUITE_REGISTRATION( BBTreeTest );
CPPUNIT_TEST_SUITE_REGISTRATION( ExprEvalInterpTest );
CPPUNIT_TEST_SUITE_REGISTRATION( QuadraticPlanarInterpTest );
CPPUNIT_TEST_SUITE_REGISTRATION( SingleElementPlanarTests );
CPPUNIT_TEST_SUITE_REGISTRATION( TransformedTriangleIntersectTest );
CPPUNIT_TEST_SUITE_REGISTRATION( TransformedTriangleTest );
CPPUNIT_TEST_SUITE_REGISTRATION( UnitTetraIntersectionBaryTest );
CPPUNIT_TEST_SUITE_REGISTRATION( UnitTetra3D2DIntersectionTest );

#ifndef MEDCOUPLING_MICROMED
// These test suites need MEDLoader to load some test files:
CPPUNIT_TEST_SUITE_REGISTRATION( InterpolationOptionsTest );
CPPUNIT_TEST_SUITE_REGISTRATION( HexaTests );
CPPUNIT_TEST_SUITE_REGISTRATION( MultiElement2DTests );
CPPUNIT_TEST_SUITE_REGISTRATION( MultiElementTetraTests );
CPPUNIT_TEST_SUITE_REGISTRATION( SingleElementTetraTests );
CPPUNIT_TEST_SUITE_REGISTRATION( ThreeDSurfProjectionTest );
#endif
// --- generic Main program from KERNEL_SRC/src/Basics/Test

#include "BasicMainTest.hxx"
