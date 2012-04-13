// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __NONCOINCIDENTDEC_HXX__
#define __NONCOINCIDENTDEC_HXX__

#include "DEC.hxx"

struct _fvm_locator_t;

typedef enum {NN} InterpolationMethod;

namespace ParaMEDMEM
{   
  class NonCoincidentDEC : public DEC
  {
    public:  
    NonCoincidentDEC();
    NonCoincidentDEC(ProcessorGroup& , ProcessorGroup&);

    virtual ~NonCoincidentDEC();

    void synchronize();

    void recvData();

    void sendData();
    
    void prepareSourceDE() { }
    void prepareTargetDE() { }
    
    void setInterpolationMethod(InterpolationMethod method) { _method=method; }
    
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
