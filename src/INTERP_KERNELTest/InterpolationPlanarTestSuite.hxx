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

#ifndef __TU_INTERPOLATION_PLANAR_TEST_SUITE_HXX__
#define __TU_INTERPOLATION_PLANAR_TEST_SUITE_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include <deque>
#include <cmath>
#include <iostream>

namespace INTERP_TEST
{

  /**
   * \brief Base class for planar mesh intersection test suites.
   * 
   */
  class InterpolationPlanarTestSuite : public CppUnit::TestFixture
  {

  public:
    double _Epsilon;
    double _Precision;

    /**
     * Sets up the test suite.
     *
     */
    void setUp()
    {
      _Epsilon = 1.e-6;
      _Precision = 1.e-6;
    }
    void tearDown()  {}

    //     bool checkDequesEqual(std::deque< double > deque1, std::deque< double > deque2, double epsilon);
    //     bool checkVectorsEqual(std::vector< double > Vect1, std::vector< double > Vect2, double epsilon);
    //     void dequePrintOut(std::deque< double > deque1);
    //     void vectPrintOut(std::vector< double > vect);
    //     void tabPrintOut( const double * tab, int size);

    bool checkDequesEqual(std::deque< double > deque1,  
                          std::deque< double > deque2, double epsilon)
    {
      int size1 = deque1.size();
      int size2 = deque2.size();
      bool are_equal = size1 == size2;
    
      if(are_equal)
        for(int i = 0; i < size1 && are_equal; i++)
          are_equal = fabs(deque1[i] - deque2[i]) < epsilon;
      
      return are_equal; 
    }
    bool checkVectorsEqual(std::vector< double > vect1,  
                           std::vector< double > vect2, double epsilon)
    {
      int size1 = vect1.size();
      int size2 = vect2.size();
      bool are_equal = size1 == size2;
      
      if(are_equal)
        for(int i = 0; i < size1 && are_equal; i++)
          are_equal = fabs(vect1[i] - vect2[i]) < epsilon;
      
      return are_equal; 
    }
    void dequePrintOut(std::deque< double > deque1)
    {
      for(int i = 0; i< (int)deque1.size(); i++)
        {
          std::cerr << deque1[i] << " ";
        }
      std::cerr<< std::endl;
    }
    void vectPrintOut(std::vector< double > vect)
    {
      for(int i = 0; i< (int)vect.size(); i++)
        {
          std::cerr << vect[i] << " ";
        }
      std::cerr<< std::endl;
    }  
    void tabPrintOut( const double * tab,int size)
    {
      for(int i = 0; i< size; i++)
        {
          std::cerr << tab[i] << " ";
        }
      std::cerr<< std::endl;
    }  

    /**
     * Cleans up after the test suite.
     * Liberates the MeshTestToolkit object used by the tests.
     */
    //     void tearDown()
    //     {
    //       delete _testTools;
    //     }

    

    //   protected:
    //     /// MeshTestToolkit object to which the tests are delegated
    //     MeshTestToolkit* _testTools; 

  };
}
#endif
