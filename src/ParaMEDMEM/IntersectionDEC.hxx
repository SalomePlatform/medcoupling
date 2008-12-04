//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
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
