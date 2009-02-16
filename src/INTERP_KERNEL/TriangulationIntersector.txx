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
#ifndef __TRIANGULATIONINTERSECTOR_TXX__
#define __TRIANGULATIONINTERSECTOR_TXX__

#include "TriangulationIntersector.hxx"
#include "PlanarIntersectorP0P0.txx"
#include "PlanarIntersectorP0P1.txx"
#include "PlanarIntersectorP1P0.txx"

#include "InterpolationUtils.hxx"
#include "PlanarIntersector.hxx"

#include <iostream>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  TriangulationIntersector<MyMeshType,MyMatrix,InterpType>::TriangulationIntersector(const MyMeshType& meshT, const MyMeshType& meshS, 
                                                                              double DimCaracteristic, double Precision,
                                                                              double MedianPlane, int orientation, int PrintLevel)
    :InterpType<MyMeshType,MyMatrix,TriangulationIntersector<MyMeshType,MyMatrix,InterpType> >(meshT,meshS,DimCaracteristic, Precision, MedianPlane, true, orientation, PrintLevel)
  {
    if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 1)
      {
        std::cout << "  - intersection type = triangles " << std::endl;
        if(SPACEDIM==3) std::cout << "_do_rotate = true"<< std::endl;
      }
  }
  
  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  double TriangulationIntersector<MyMeshType,MyMatrix,InterpType>::intersectGeometry(ConnType icellT, ConnType icellS, ConnType nbNodesT, ConnType nbNodesS)
  {
    double result = 0.;
    int orientation = 1;
                    
    //Obtain the coordinates of T and S
    std::vector<double> CoordsT;
    std::vector<double> CoordsS;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealCoordinates(icellT,icellS,nbNodesT,nbNodesS,CoordsT,CoordsS,orientation);
    //Compute the intersection area
    double area[SPACEDIM];
    for(ConnType iT = 1; iT<nbNodesT-1; iT++)
      {
        for(ConnType iS = 1; iS<nbNodesS-1; iS++)
          {
            std::vector<double> inter;
            INTERP_KERNEL::intersec_de_triangle(&CoordsT[0],&CoordsT[SPACEDIM*iT],&CoordsT[SPACEDIM*(iT+1)],
                                                &CoordsS[0],&CoordsS[SPACEDIM*iS],&CoordsS[SPACEDIM*(iS+1)],
                                                inter, PlanarIntersector<MyMeshType,MyMatrix>::_dim_caracteristic,
                                                PlanarIntersector<MyMeshType,MyMatrix>::_precision);
            ConnType nb_inter=((ConnType)inter.size())/2;
            if(nb_inter >3) inter=reconstruct_polygon(inter);
            for(ConnType i = 1; i<nb_inter-1; i++)
              {
                INTERP_KERNEL::crossprod<2>(&inter[0],&inter[2*i],&inter[2*(i+1)],area);
                result +=0.5*fabs(area[0]);
              }
            //DEBUG prints
            if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 3)
              {
                std::cout << std::endl << "Number of nodes of the intersection = "<< nb_inter << std::endl;
                for(ConnType i=0; i< nb_inter; i++)
                  {for (int idim=0; idim<2; idim++) std::cout << inter[2*i+idim] << " "; std::cout << std::endl; }
              }
          }
      }
    
    //DEBUG PRINTS    
    if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 3) 
      std::cout << std::endl <<"Intersection area = " << result << std::endl;
    
    return orientation*result;
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  double TriangulationIntersector<MyMeshType,MyMatrix,InterpType>::intersectGeometryWithQuadrangle(const double *quadrangle, const std::vector<double>& sourceCoords, bool isSourceQuad)
  {
    double result = 0.;
    ConnType nbNodesS=sourceCoords.size()/SPACEDIM;
    //Compute the intersection area
    double area[SPACEDIM];
    for(ConnType iT = 1; iT<3; iT++)
      {
        for(ConnType iS = 1; iS<nbNodesS-1; iS++)
          {
            std::vector<double> inter;
            INTERP_KERNEL::intersec_de_triangle(quadrangle,&quadrangle[SPACEDIM*iT],&quadrangle[SPACEDIM*(iT+1)],
                                                &sourceCoords[0],&sourceCoords[SPACEDIM*iS],&sourceCoords[SPACEDIM*(iS+1)],
                                                inter, PlanarIntersector<MyMeshType,MyMatrix>::_dim_caracteristic,
                                                PlanarIntersector<MyMeshType,MyMatrix>::_precision);
            ConnType nb_inter=((ConnType)inter.size())/2;
            if(nb_inter >3) inter=reconstruct_polygon(inter);
            for(ConnType i = 1; i<nb_inter-1; i++)
              {
                INTERP_KERNEL::crossprod<2>(&inter[0],&inter[2*i],&inter[2*(i+1)],area);
                result +=0.5*fabs(area[0]);
              }
            //DEBUG prints
            if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 3)
              {
                std::cout << std::endl << "Number of nodes of the intersection = "<< nb_inter << std::endl;
                for(ConnType i=0; i< nb_inter; i++)
                  {for (int idim=0; idim<2; idim++) std::cout << inter[2*i+idim] << " "; std::cout << std::endl; }
              }
          }
      }
    
    //DEBUG PRINTS    
    if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 3) 
      std::cout << std::endl <<"Intersection area = " << result << std::endl;
    
    return result;
  }
}

#endif
