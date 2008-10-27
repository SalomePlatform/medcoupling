#ifndef __CONVEXINTERSECTOR_TXX__
#define __CONVEXINTERSECTOR_TXX__

#include "ConvexIntersector.hxx"

#include "PolygonAlgorithms.txx"

#include <iostream>

namespace INTERP_KERNEL
{
  template<class MyMeshType>
  ConvexIntersector<MyMeshType>::ConvexIntersector(const MyMeshType& mesh_A, const MyMeshType& mesh_B, 
                                                   double dimCaracteristic, double precision,
                                                   double medianPlane, bool doRotate , int printLevel)
    :PlanarIntersector<MyMeshType>(dimCaracteristic, precision, medianPlane, doRotate, printLevel),
     _epsilon(precision*dimCaracteristic)
  {
    _connectA= mesh_A.getConnectivityPtr();
    _connectB= mesh_B.getConnectivityPtr();
    _connIndexA= mesh_A.getConnectivityIndexPtr();
    _connIndexB= mesh_B.getConnectivityIndexPtr();
    _coordsA = mesh_A.getCoordinatesPtr();
    _coordsB = mesh_B.getCoordinatesPtr();
    if(PlanarIntersector<MyMeshType>::_printLevel >= 1)
      {        
        std::cout << " - intersection type = convex " << std::endl;
        if(SPACEDIM==3){
          if(PlanarIntersector<MyMeshType>::_doRotate) std::cout << "  _do_rotate = true" << std::endl;
          else std::cout << "  _do_rotate = false" << std::endl;
        }
      }
  }

  template<class MyMeshType>
  double ConvexIntersector<MyMeshType>::intersectCells(ConnType icell_A, ConnType icell_B, 
                                                       int nb_NodesA, int nb_NodesB)
  {
    double result = 0;
		int orientation = 1;

    /*** Obtain the coordinates of A and B ***/
    std::vector<double> Coords_A (SPACEDIM*nb_NodesA);
    std::vector<double> Coords_B (SPACEDIM*nb_NodesB);
    int nb_dist_NodesA=nb_NodesA;
    int nb_dist_NodesB=nb_NodesB;
    int i_last = nb_NodesA - 1;
    const double * Pi_last= _coordsA + SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectA[OTT<ConnType,numPol>::conn2C(_connIndexA[OTT<ConnType,numPol>::ind2C(icell_A)]+i_last)]);

    for (int i_A=0; i_A<nb_NodesA; i_A++)
      {
        const double * Pi_A = _coordsA + SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectA[OTT<ConnType,numPol>::conn2C(_connIndexA[OTT<ConnType,numPol>::ind2C(icell_A)]+i_A)]);
        if(distance2<SPACEDIM>(Pi_last, Pi_A)>_epsilon)
          {
            for (int idim=0; idim<SPACEDIM; idim++)
              Coords_A[SPACEDIM*i_A+idim]=Pi_A[idim];
            i_last=i_A; Pi_last = Pi_A;
          }
        else 
          nb_dist_NodesA--;
      }
    i_last = nb_NodesB - 1;
    Pi_last=_coordsB + SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectB[OTT<ConnType,numPol>::conn2C(_connIndexB[OTT<ConnType,numPol>::ind2C(icell_B)]+i_last)]);
    for (int i_B=0; i_B<nb_NodesB; i_B++)
      {
        const double * Pi_B=_coordsB+SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectB[OTT<ConnType,numPol>::conn2C(_connIndexB[OTT<ConnType,numPol>::ind2C(icell_B)]+i_B)]);
        if(distance2<SPACEDIM>(Pi_last, Pi_B)>_epsilon)
          {
            for (int idim=0; idim<SPACEDIM; idim++)
              Coords_B[SPACEDIM*i_B+idim]=Pi_B[idim];
            i_last=i_B; Pi_last = Pi_B;
          }
        else
          nb_dist_NodesB--;
      }
      
    /*** project cells A and B on the median plane ***/
    /***  and rotate the median plane ***/
    if(SPACEDIM==3) 
			orientation = projection(Coords_A, Coords_B, nb_dist_NodesA, nb_dist_NodesB,_epsilon,
                               PlanarIntersector<MyMeshType>::_medianPlane,
                               PlanarIntersector<MyMeshType>::_doRotate);

    //DEBUG PRINTS
    if(PlanarIntersector<MyMeshType>::_printLevel >= 3) 
      {
        std::cout << std::endl << "Cell coordinates (possibly after projection)" << std::endl;
        std::cout << std::endl << "icell_A= " << icell_A << ", nb nodes A= " <<  nb_NodesA << std::endl;
        for(int i_A =0; i_A< nb_NodesA; i_A++)
          {for (int idim=0; idim<SPACEDIM; idim++) std::cout << Coords_A[SPACEDIM*i_A+idim] << " "; std::cout << std::endl;}
        std::cout << std::endl << "icell_B= " << icell_B << ", nb nodes B= " <<  nb_NodesB << std::endl;
        for(int i_B =0; i_B< nb_NodesB; i_B++)
          {for (int idim=0; idim<SPACEDIM; idim++) std::cout << Coords_B[SPACEDIM*i_B+idim]<< " "; std::cout << std::endl; }
      }

    /*** Compute the intersection area ***/
    INTERP_KERNEL::PolygonAlgorithms<SPACEDIM> P(_epsilon, PlanarIntersector<MyMeshType>::_precision);
    std::deque<double> inter =  P.intersect_convex_polygons(&Coords_A[0], &Coords_B[0],
                                                            nb_dist_NodesA, nb_dist_NodesB);
    double area[SPACEDIM];
    int nb_inter =((int)inter.size())/SPACEDIM;
    for(int i = 1; i<nb_inter-1; i++)
      {
        INTERP_KERNEL::crossprod<SPACEDIM>(&inter[0],&inter[SPACEDIM*i],&inter[SPACEDIM*(i+1)],area);
        result +=0.5*norm<SPACEDIM>(area);
      }

    //DEBUG prints
    if(PlanarIntersector<MyMeshType>::_printLevel >= 3)
      {       
        std::cout << std::endl << "Number of nodes of the intersection = "<<  nb_inter << std::endl;
        for(int i=0; i<  nb_inter; i++)
          {for (int idim=0; idim<SPACEDIM; idim++) std::cout << inter[SPACEDIM*i+idim]<< " "; std::cout << std::endl;}
        std::cout << std::endl <<"Intersection area = " << result << std::endl;
      }
    
    return orientation*result;
  }
}

#endif
