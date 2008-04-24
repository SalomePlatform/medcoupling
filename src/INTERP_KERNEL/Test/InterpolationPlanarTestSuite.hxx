#ifndef __TU_INTERPOLATION_PLANAR_TEST_SUITE_HXX__
#define __TU_INTERPOLATION_PLANAR_TEST_SUITE_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include <deque>
#include <cmath>
#include <iostream>

using namespace std;

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
		void tearDown()	{}

// 		bool checkDequesEqual(std::deque< double > deque1, std::deque< double > deque2, double epsilon);
// 		bool checkVectorsEqual(std::vector< double > Vect1, std::vector< double > Vect2, double epsilon);
// 		void dequePrintOut(std::deque< double > deque1);
// 		void vectPrintOut(std::vector< double > vect);
// 		void tabPrintOut( const double * tab, int size);

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
			std::cerr<< endl;
		}
		void vectPrintOut(std::vector< double > vect)
		{
			for(int i = 0; i< (int)vect.size(); i++)
				{
					std::cerr << vect[i] << " ";
				}
			std::cerr<< endl;
		}	
		void tabPrintOut( const double * tab,int size)
		{
			for(int i = 0; i< size; i++)
				{
					std::cerr << tab[i] << " ";
				}
			std::cerr<< endl;
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
