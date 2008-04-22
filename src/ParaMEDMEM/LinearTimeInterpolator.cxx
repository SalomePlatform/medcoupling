
#include "LinearTimeInterpolator.hxx"

using namespace std;

namespace ParaMEDMEM {    

LinearTimeInterpolator::LinearTimeInterpolator( double InterpPrecision, int nStepBefore,
                                                int nStepAfter ):
  TimeInterpolator( InterpPrecision, nStepBefore, nStepAfter ) {
}

LinearTimeInterpolator::~LinearTimeInterpolator() {
} 

void LinearTimeInterpolator::DoInterp( double time0, double time1, double time,
                                       int recvcount , int nbuff0, int nbuff1,
                                       int **recvbuff0, int **recvbuff1, int *result ) {
  int i ;
  for ( i = 0 ; i < recvcount ; i++ ) {
     result[i] = (int ) ((recvbuff0[0][i]*(time1 - time) + recvbuff1[0][i]*(time - time0))/(time1 - time0) + _InterpPrecision) ;
     //cout << "DoInterpint time " << time << " time0 " << time0 << " time1 " << time1
     //     << " recvbuff0 " << recvbuff0[0][i] << " recvbuff1 " << recvbuff1[0][i]
     //     << " --> " << result[i] << endl ;
  }
}

void LinearTimeInterpolator::DoInterp( double time0, double time1, double time,
                                       int recvcount , int nbuff0, int nbuff1,
                                       double **recvbuff0, double **recvbuff1,
                                       double *result ) {
  int i ;
  for ( i = 0 ; i < recvcount ; i++ ) {
     result[i] = (recvbuff0[0][i]*(time1 - time) + recvbuff1[0][i]*(time - time0))/(time1 - time0) ;
     //cout << "DoInterpdouble time " << time << " time0 " << time0 << " time1 " << time1
     //     << " recvbuff0 " << recvbuff0[0][i] << " recvbuff1 " << recvbuff1[0][i]
     //     << " --> " << result[i] << endl ;
  }
}

}
