#ifndef __TU_INTERPOLATION_3D_TEST_HXX__
#define __TU_INTERPOLATION_3D_TEST_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include "../Interpolation3D.hxx"

#define ERR_TOL 1.0e-8

using MEDMEM::Interpolation3D;
class MEDMEM::MESH;

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
