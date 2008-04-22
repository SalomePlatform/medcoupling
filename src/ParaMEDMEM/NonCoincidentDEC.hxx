#ifndef NONCOINCIDENTDEC_HXX_
#define NONCOINCIDENTDEC_HXX_

struct _fvm_locator_t;

typedef enum {NN} InterpolationMethod;

namespace ParaMEDMEM
{
  class DEC;
    
  class NonCoincidentDEC:public DEC
  {
    public:  
    NonCoincidentDEC();
    NonCoincidentDEC(ProcessorGroup& , ProcessorGroup&);

    virtual ~NonCoincidentDEC();

    void synchronize();

    void recvData();

    void sendData();
    
    void prepareSourceDE(){};
    void prepareTargetDE(){};
    
    void setInterpolationMethod(InterpolationMethod method)
    {_method=method;}
    
    private :
    // Structure for computing the localization
    // of remote nodes on local mesh
    _fvm_locator_t* _locator;
    
    //Number of distant points to be located locally 
    int _nb_distant_points;
    
    //coordinates of distant points 
    const double* _distant_coords;
    
    //local element number containing the distant points  
    const int* _distant_locations; 
   
   //inerpolation method
   InterpolationMethod _method;
  };
}


#endif
