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

#ifndef __TU_INTERPOLATION_TEST_SUITE_HXX__
#define __TU_INTERPOLATION_TEST_SUITE_HXX__

#include "MeshTestToolkit.txx"

#include <cppunit/extensions/HelperMacros.h>

namespace INTERP_TEST
{

  /**
   * \brief Base class for mesh intersection test suites.
   * 
   */
  template<int SPACEDIM, int MESHDIM>
  class InterpolationTestSuite : public CppUnit::TestFixture
  {

  public:
    /**
     * Sets up the test suite.
     * Creates the MeshTestToolkit object used by the tests.
     *
     */
    void setUp()
    {
      _testTools = new MeshTestToolkit<SPACEDIM,MESHDIM>();
    }

    /**
     * Cleans up after the test suite.
     * Liberates the MeshTestToolkit object used by the tests.
     */
    void tearDown()
    {
      delete _testTools;
    }

    

  protected:
    /// MeshTestToolkit object to which the tests are delegated
    MeshTestToolkit<SPACEDIM,MESHDIM>* _testTools; 
  
  };
}
#endif
