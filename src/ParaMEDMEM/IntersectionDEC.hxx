#ifndef INTERSECTIONDEC_HXX_
#define INTERSECTIONDEC_HXX_

#include "MEDMEM_OptionManager.hxx"
#include "MPI_AccessDEC.hxx"
 
namespace ParaMEDMEM
{
  class DEC;
  class InterpolationMatrix;

  class IntersectionDEC:public DEC
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
   
   //inerpolation method
   string _method;
   
   InterpolationMatrix* _interpolation_matrix;

		double _bb_adjustment;

		bool _asynchronous;
    TimeInterpolationMethod _timeinterpolationmethod ;
		AllToAllMethod _allToAllMethod;
  };
}


#endif
