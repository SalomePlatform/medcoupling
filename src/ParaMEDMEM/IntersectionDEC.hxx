#ifndef INTERSECTIONDEC_HXX_
#define INTERSECTIONDEC_HXX_

#include "DEC.hxx"
#include "MPI_AccessDEC.hxx"
#include "MxN_Mapping.hxx" 
#include "InterpolationPlanar.hxx"
#include "IntersectorHexa.hxx"
#include "InterpolationOptions.hxx"

namespace ParaMEDMEM
{
  class InterpolationMatrix;

  class IntersectionDEC:public DEC, public INTERP_KERNEL::InterpolationOptions
  {
    public:  
    IntersectionDEC();
    IntersectionDEC(ProcessorGroup& local_group, ProcessorGroup& distant_group);
    
    virtual ~IntersectionDEC();

    void synchronize();

    void recvData();
    void recvData( double time );

    void sendData();

    void sendData( double time , double deltatime );

    void prepareSourceDE(){};
    void prepareTargetDE(){};
    
    private :
    
    //Number of distant points to be located locally 
    int _nb_distant_points;
    
    //coordinates of distant points 
    const double* _distant_coords;
    
    //local element number containing the distant points  
    const int* _distant_locations; 
     
   InterpolationMatrix* _interpolation_matrix;

  };
}


#endif
