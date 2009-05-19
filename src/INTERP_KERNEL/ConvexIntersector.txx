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
#ifndef __CONVEXINTERSECTOR_TXX__
#define __CONVEXINTERSECTOR_TXX__

#include "ConvexIntersector.hxx"
#include "PlanarIntersectorP0P0.txx"
#include "PlanarIntersectorP0P1.txx"
#include "PlanarIntersectorP1P0.txx"
#include "PlanarIntersectorP1P1.txx"

#include "PolygonAlgorithms.txx"

#include <iostream>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  ConvexIntersector<MyMeshType,MyMatrix,InterpType>::ConvexIntersector(const MyMeshType& meshT, const MyMeshType& meshS, 
                                                                       double dimCaracteristic, double precision, double md3DSurf,
                                                                       double medianPlane, bool doRotate , int oriantation, int printLevel)
    :InterpType<MyMeshType,MyMatrix,ConvexIntersector<MyMeshType,MyMatrix,InterpType> >(meshT,meshS,dimCaracteristic, precision, md3DSurf, medianPlane, doRotate, oriantation, printLevel),
     _epsilon(precision*dimCaracteristic)
  {
    if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 1)
      {        
        std::cout << " - intersection type = convex " << std::endl;
        if(SPACEDIM==3){
          if(PlanarIntersector<MyMeshType,MyMatrix>::_do_rotate) std::cout << "  _do_rotate = true" << std::endl;
          else std::cout << "  _do_rotate = false" << std::endl;
        }
      }
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  double ConvexIntersector<MyMeshType,MyMatrix,InterpType>::intersectGeometry(ConnType icellT, ConnType icellS, ConnType nbNodesT, ConnType nbNodesS)
  {
    double result = 0;
    int orientation = 1;
    
    /*** Obtain the coordinates of T and S ***/
    std::vector<double> CoordsT;
    std::vector<double> CoordsS;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealCoordinates(icellT,icellS,nbNodesT,nbNodesS,CoordsT,CoordsS,orientation);
    /*** Compute the intersection area ***/
    INTERP_KERNEL::PolygonAlgorithms<SPACEDIM> P(_epsilon, PlanarIntersector<MyMeshType,MyMatrix>::_precision);
    std::deque<double> inter =  P.intersectConvexPolygons(&CoordsT[0], &CoordsS[0],
                                                            CoordsT.size()/SPACEDIM, CoordsS.size()/SPACEDIM);
    double area[SPACEDIM];
    int nb_inter =((int)inter.size())/SPACEDIM;
    for(int i = 1; i<nb_inter-1; i++)
      {
        INTERP_KERNEL::crossprod<SPACEDIM>(&inter[0],&inter[SPACEDIM*i],&inter[SPACEDIM*(i+1)],area);
        result +=0.5*norm<SPACEDIM>(area);
      }

    //DEBUG prints
    if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 3)
      {       
        std::cout << std::endl << "Number of nodes of the intersection = "<<  nb_inter << std::endl;
        for(int i=0; i<  nb_inter; i++)
          {for (int idim=0; idim<SPACEDIM; idim++) std::cout << inter[SPACEDIM*i+idim]<< " "; std::cout << std::endl;}
        std::cout << std::endl <<"Intersection area = " << result << std::endl;
      }
    
    return orientation*result;
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  double ConvexIntersector<MyMeshType,MyMatrix,InterpType>::intersectGeometryWithQuadrangle(const double *quadrangle, const std::vector<double>& sourceCoords, bool isSourceQuad)
  {
    double result = 0;
    int nbOfNodesS=sourceCoords.size()/SPACEDIM;

    /*** Compute the intersection area ***/
    INTERP_KERNEL::PolygonAlgorithms<SPACEDIM> P(_epsilon, PlanarIntersector<MyMeshType,MyMatrix>::_precision);
    std::deque<double> inter =  P.intersectConvexPolygons(quadrangle, &sourceCoords[0],
                                                            4, nbOfNodesS);
    double area[SPACEDIM];
    int nb_inter =((int)inter.size())/SPACEDIM;
    for(int i = 1; i<nb_inter-1; i++)
      {
        INTERP_KERNEL::crossprod<SPACEDIM>(&inter[0],&inter[SPACEDIM*i],&inter[SPACEDIM*(i+1)],area);
        result +=0.5*norm<SPACEDIM>(area);
      }

    //DEBUG prints
    if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 3)
      {       
        std::cout << std::endl << "Number of nodes of the intersection = "<<  nb_inter << std::endl;
        for(int i=0; i<  nb_inter; i++)
          {for (int idim=0; idim<SPACEDIM; idim++) std::cout << inter[SPACEDIM*i+idim]<< " "; std::cout << std::endl;}
        std::cout << std::endl <<"Intersection area = " << result << std::endl;
      }
    
    return result;
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  double ConvexIntersector<MyMeshType,MyMatrix,InterpType>::intersectGeometryGeneral(const std::vector<double>& targetCoords, const std::vector<double>& sourceCoords)
  {
    double result = 0;
    int nbOfNodesS=sourceCoords.size()/SPACEDIM;
    int nbOfNodesT=targetCoords.size()/SPACEDIM;
    /*** Compute the intersection area ***/
    INTERP_KERNEL::PolygonAlgorithms<SPACEDIM> P(_epsilon, PlanarIntersector<MyMeshType,MyMatrix>::_precision);
    std::deque<double> inter =  P.intersectConvexPolygons(&targetCoords[0], &sourceCoords[0],
                                                          nbOfNodesT, nbOfNodesS);
    double area[SPACEDIM];
    int nb_inter =((int)inter.size())/SPACEDIM;
    for(int i = 1; i<nb_inter-1; i++)
      {
        INTERP_KERNEL::crossprod<SPACEDIM>(&inter[0],&inter[SPACEDIM*i],&inter[SPACEDIM*(i+1)],area);
        result +=0.5*norm<SPACEDIM>(area);
      }
    return result;
  }
}

#endif
