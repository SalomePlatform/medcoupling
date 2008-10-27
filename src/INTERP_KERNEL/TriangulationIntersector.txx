#ifndef __TRIANGULATIONINTERSECTOR_TXX__
#define __TRIANGULATIONINTERSECTOR_TXX__

#include "TriangulationIntersector.hxx"

#include "InterpolationUtils.hxx"
#include "PlanarIntersector.hxx"

#include <iostream>

namespace INTERP_KERNEL
{
  template<class MyMeshType>
  TriangulationIntersector<MyMeshType>::TriangulationIntersector(const MyMeshType& mesh_A, const MyMeshType& mesh_B, 
                                                                 double DimCaracteristic, double Precision,
                                                                 double MedianPlane, int PrintLevel)
    :PlanarIntersector<MyMeshType>(DimCaracteristic, Precision, MedianPlane, true, PrintLevel)
  {
    _connectA= mesh_A.getConnectivityPtr();
    _connectB= mesh_B.getConnectivityPtr();
    _connIndexA= mesh_A.getConnectivityIndexPtr();
    _connIndexB= mesh_B.getConnectivityIndexPtr();
    _coordsA = mesh_A.getCoordinatesPtr();
    _coordsB = mesh_B.getCoordinatesPtr();
    if(PlanarIntersector<MyMeshType>::_printLevel >= 1)
      {
        std::cout << "  - intersection type = triangles " << std::endl;
        if(SPACEDIM==3) std::cout << "_do_rotate = true"<< std::endl;
      }
  }

  template<class MyMeshType>
  double TriangulationIntersector<MyMeshType>::intersectCells(ConnType icell_A, ConnType icell_B,
                                                              int nb_NodesA, int nb_NodesB)
  {
    double result = 0.;
		int orientation = 1;
                    
    //Obtain the coordinates of A and B
    std::vector<double> Coords_A (SPACEDIM*nb_NodesA);
    std::vector<double> Coords_B (SPACEDIM*nb_NodesB);
    for (int idim=0; idim<SPACEDIM; idim++)
      {
        for (ConnType i_A=0; i_A<nb_NodesA; i_A++)
          Coords_A[SPACEDIM*i_A+idim] = _coordsA[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectA[OTT<ConnType,numPol>::conn2C(_connIndexA[OTT<ConnType,numPol>::ind2C(icell_A)]+i_A)])+idim];
        for (ConnType i_B=0; i_B<nb_NodesB; i_B++)
          Coords_B[SPACEDIM*i_B+idim] = _coordsB[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectB[OTT<ConnType,numPol>::conn2C(_connIndexB[OTT<ConnType,numPol>::ind2C(icell_B)]+i_B)])+idim];
      }
    
    //project cells A and B on the median plane and rotate the median plane
    if(SPACEDIM==3)
      orientation = projection(Coords_A, Coords_B, nb_NodesA, nb_NodesB,
															 PlanarIntersector<MyMeshType>::_dimCaracteristic * PlanarIntersector<MyMeshType>::_precision,
															 PlanarIntersector<MyMeshType>::_medianPlane, PlanarIntersector<MyMeshType>::_doRotate);
    
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
                
    //Compute the intersection area
    double area[SPACEDIM];
    for(ConnType i_A = 1; i_A<nb_NodesA-1; i_A++)
      {
        for(ConnType i_B = 1; i_B<nb_NodesB-1; i_B++)
          {
            std::vector<double> inter;
            INTERP_KERNEL::intersec_de_triangle(&Coords_A[0],&Coords_A[SPACEDIM*i_A],&Coords_A[SPACEDIM*(i_A+1)],
                                                &Coords_B[0],&Coords_B[SPACEDIM*i_B],&Coords_B[SPACEDIM*(i_B+1)],
                                                inter, PlanarIntersector<MyMeshType>::_dimCaracteristic,
                                                PlanarIntersector<MyMeshType>::_precision);
            ConnType nb_inter=((ConnType)inter.size())/2;
            if(nb_inter >3) inter=reconstruct_polygon(inter);
            for(ConnType i = 1; i<nb_inter-1; i++)
              {
                INTERP_KERNEL::crossprod<2>(&inter[0],&inter[2*i],&inter[2*(i+1)],area);
                result +=0.5*fabs(area[0]);
              }
            //DEBUG prints
            if(PlanarIntersector<MyMeshType>::_printLevel >= 3)
              {
                std::cout << std::endl << "Number of nodes of the intersection = "<< nb_inter << std::endl;
                for(ConnType i=0; i< nb_inter; i++)
                  {for (int idim=0; idim<2; idim++) std::cout << inter[2*i+idim] << " "; std::cout << std::endl; }
              }
          }
      }
    
    //DEBUG PRINTS    
    if(PlanarIntersector<MyMeshType>::_printLevel >= 3) 
      std::cout << std::endl <<"Intersection area = " << result << std::endl;
    
    return orientation*result;
  }
}

#endif
