
#include "TimeInterpolator.hxx"

using namespace std;

namespace ParaMEDMEM {    

TimeInterpolator::TimeInterpolator( double InterpPrecision, int nStepBefore,
                                    int nStepAfter ){
  _InterpPrecision = InterpPrecision ;
  _nStepBefore = nStepBefore ;
  _nStepAfter = nStepAfter ;
}

TimeInterpolator::~TimeInterpolator() {
} 

}
