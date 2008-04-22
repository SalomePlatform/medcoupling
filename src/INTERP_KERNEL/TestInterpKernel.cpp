// Copyright (C) 2005  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
//
// This library is distributed in the hope that it will be useful
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
#include "TransformedTriangleTest.hxx"
#include "TransformedTriangleIntersectTest.hxx"
#include "MultiElementTetraTests.hxx"
#include "SingleElementTetraTests.hxx"
#include "HexaTests.hxx"

using namespace INTERP_TEST;

// --- Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( HexaTests );
CPPUNIT_TEST_SUITE_REGISTRATION( MultiElementTetraTests );
CPPUNIT_TEST_SUITE_REGISTRATION( SingleElementTetraTests );
CPPUNIT_TEST_SUITE_REGISTRATION( INTERP_TEST::TransformedTriangleIntersectTest );
CPPUNIT_TEST_SUITE_REGISTRATION( INTERP_TEST::TransformedTriangleTest );

// --- generic Main program from KERNEL_SRC/src/Basics/Test

#include "BasicMainTest.hxx"
