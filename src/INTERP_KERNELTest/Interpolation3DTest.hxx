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

#ifndef __TU_INTERPOLATION_3D_TEST_HXX__
#define __TU_INTERPOLATION_3D_TEST_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include "Interpolation3D.hxx"

#define ERR_TOL 1.0e-8

/// \brief OBSOLETE - renamed Interpolation3DTestSuite
class Interpolation3DTest : public CppUnit::TestFixture
{

public:
  void setUp()
  {
    _testTools = new MeshTestToolkit();
  }

  void tearDown()
  {
    delete _testTools;
  }

protected:

  MeshToolkit* _testTools; 

};

#endif
