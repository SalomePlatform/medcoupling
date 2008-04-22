#ifndef __PLANARINTERSECTOR_TXX__
#define __PLANARINTERSECTOR_TXX__

#include "PlanarIntersector.hxx"
#include "InterpolationUtils.hxx"
#include "TranslationRotationMatrix.hxx"

#include <iostream>

namespace INTERP_KERNEL
{
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::PlanarIntersector(double dimCaracteristic, double precision, double medianPlane, bool doRotate, int printLevel):
    _dimCaracteristic(dimCaracteristic),_precision(precision),_medianPlane(medianPlane),
    _doRotate(doRotate),_printLevel(printLevel)
  {
  }

  /*!
    \brief creates the bounding boxes for all the cells of mesh \a mesh
  
    The method accepts mixed meshes (containing triangles and quadrangles).
    The vector returned is of dimension 6*nb_elems with bounding boxes stored as xmin1, xmax1, ymin1, ymax1, zmin1, zmax1, xmin2, xmax2, ymin2,... 
    The returned pointer must be deleted by the calling code.
  
    \param mesh structure pointing to the mesh
    \param bbox vector containing the bounding boxes
  */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  void PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::createBoundingBoxes(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh,
                                                                                           std::vector<double>& bbox)
  {
    /* We build the segment tree for locating possible matching intersections*/
    long nbelems = mesh.getNumberOfElements();
    bbox.resize(2*SPACEDIM* nbelems);
    const double* coords = mesh.getCoordinatesPtr();
    const ConnType* conn = mesh.getConnectivityPtr();
    const ConnType* conn_index = mesh.getConnectivityIndexPtr();  
    int ibox=0;
    for(long icell=0; icell<nbelems; icell++)
      {
        int nb_nodes_per_elem =conn_index[icell+1]-conn_index[icell];
        //initializing bounding box limits
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            bbox[2*SPACEDIM*ibox+2*idim]   =  HUGE;
            bbox[2*SPACEDIM*ibox+2*idim+1] = -HUGE;
          }
        //updating the bounding box with each node of the element
        for (int j=0; j<nb_nodes_per_elem; j++)
          {
            const double* coord_node=coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(conn[OTT<ConnType,numPol>::conn2C(conn_index[icell]+j)]);
            for(int idim=0; idim<SPACEDIM; idim++)
              {            
                double x=*(coord_node+idim);
                bbox[ibox*2*SPACEDIM + 2*idim]   = (bbox[ibox*2*SPACEDIM + 2*idim]  <x)?bbox[ibox*2*SPACEDIM + 2*idim  ]:x;
                bbox[ibox*2*SPACEDIM + 2*idim+1] = (bbox[ibox*2*SPACEDIM + 2*idim+1]>x)?bbox[ibox*2*SPACEDIM + 2*idim+1]:x;
              }
          }
        ibox++;
      }                        
  }

  /*!
    Computes the bouding box of a given element. iP in numPol mode.
  */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  void PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::getElemBB(double* bb,
                                                                                 const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh,
                                                                                 ConnType iP, ConnType nb_nodes)
  {
    const double* coords = mesh.getCoordinatesPtr();
    const ConnType* conn_index = mesh.getConnectivityIndexPtr();
    const ConnType* conn = mesh.getConnectivityPtr();
    //initializing bounding box limits
    for(int idim=0; idim<SPACEDIM; idim++)
      {
      bb[2*idim  ] =  HUGE;
      bb[2*idim+1] = -HUGE;
      }
  
    for (ConnType i=0; i<nb_nodes; i++)
      {
        //MN: iP= cell index, not node index, use of connectivity array ?
        const double* coord_node=coords+SPACEDIM*(OTT<ConnType,numPol>::coo2C(conn[OTT<ConnType,numPol>::conn2C(conn_index[OTT<ConnType,numPol>::ind2C(iP)]+i)]));
        for(int idim=0; idim<SPACEDIM; idim++)
          {            
            double x = *(coord_node+idim);
            //double y = *(mesh.getCoordinates(MED_FULL_INTERLACE)+3*(iP+i)+1);
            bb[2*idim  ] = (x<bb[2*idim  ])?x:bb[2*idim  ];
            bb[2*idim+1] = (x>bb[2*idim+1])?x:bb[2*idim+1];
          }            
      }
  }

  /*! Readjusts a set of bounding boxes so that they are extended
    in all dimensions for avoiding missing interesting intersections
  
    \param bbox vector containing the bounding boxes
  */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  void PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::adjustBoundingBoxes(std::vector<double>& bbox, double Surf3DAdjustmentEps)
  {
    /* We build the segment tree for locating possible matching intersections*/
  
    long size = bbox.size()/(2*SPACEDIM);
    for (int i=0; i<size; i++)
      {
        double max=-HUGE;
        for(int idim=0; idim<SPACEDIM; idim++)
          {            
            double Dx=bbox[i*2*SPACEDIM+1+2*idim]-bbox[i*2*SPACEDIM+2*idim];
            max=(max<Dx)?Dx:max;
          }
        for(int idim=0; idim<SPACEDIM; idim++)
          {            
            bbox[i*2*SPACEDIM+2*idim  ] -= Surf3DAdjustmentEps*max;
            bbox[i*2*SPACEDIM+2*idim+1] += Surf3DAdjustmentEps*max;
          }
      }
  }
  
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  void PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::projection(std::vector< double>& Coords_A, std::vector< double>& Coords_B, 
                                                                                  int nb_NodesA, int nb_NodesB, double epsilon, double median_plane, bool do_rotate)
  {
    double normal_A[3]={0,0,0};
    double normal_B[3]={0,0,0};
    double linear_comb[3];
    double proj;
    //Find the normal to cells A and B
    int i_A1=1;
    while(i_A1<nb_NodesA && distance2<SPACEDIM>(&Coords_A[0],&Coords_A[SPACEDIM*i_A1])< epsilon) i_A1++;
    int i_A2=i_A1+1;
    crossprod<SPACEDIM>(&Coords_A[0], &Coords_A[SPACEDIM*i_A1], &Coords_A[SPACEDIM*i_A2],normal_A);
    while(i_A2<nb_NodesA && dotprod<SPACEDIM>(normal_A,normal_A)<epsilon)
      {
      crossprod<SPACEDIM>(&Coords_A[0], &Coords_A[SPACEDIM*i_A1], &Coords_A[SPACEDIM*i_A2],normal_A);
      i_A2++;
      }
    int i_B1=1;
    while(i_B1<nb_NodesB && distance2<SPACEDIM>(&Coords_B[0],&Coords_B[SPACEDIM*i_B1])< epsilon) i_B1++;
    int i_B2=i_B1+1;
    crossprod<SPACEDIM>(&Coords_B[0], &Coords_B[SPACEDIM*i_B1], &Coords_B[SPACEDIM*i_B2],normal_B);
    while(i_B2<nb_NodesB && dotprod<SPACEDIM>(normal_B,normal_B)< epsilon)
      {
      crossprod<SPACEDIM>(&Coords_B[0], &Coords_B[SPACEDIM*i_B1], &Coords_B[SPACEDIM*i_B2],normal_B);
      i_B2++;
      }
    if(i_A2<nb_NodesA && i_B2<nb_NodesB)
      {
      //Build the normal of the median plane
      if(dotprod<SPACEDIM>(normal_A,normal_B)<0)
        for(int idim =0; idim< SPACEDIM; idim++) normal_A[idim] *=-1;
      for(int idim =0; idim< SPACEDIM; idim++)
        linear_comb[idim] = median_plane*normal_A[idim] + (1-median_plane)*normal_B[idim];
      double norm= sqrt(dotprod<SPACEDIM>(linear_comb,linear_comb));

      if(norm>epsilon)
        {
          for(int idim =0; idim< SPACEDIM; idim++) linear_comb[idim]/=norm;
          
          //Project the nodes of A and B on the median plane
          for(int i_A=0; i_A<nb_NodesA; i_A++)
            {
            proj = dotprod<SPACEDIM>(&Coords_A[SPACEDIM*i_A],linear_comb);
            for(int idim =0; idim< SPACEDIM; idim++)
              Coords_A[SPACEDIM*i_A+idim] -=  proj*linear_comb[idim];
            }
          for(int i_B=0; i_B<nb_NodesB; i_B++)
            {
            proj = dotprod<SPACEDIM>(&Coords_B[SPACEDIM*i_B],linear_comb);
            for(int idim =0; idim< SPACEDIM; idim++)
              Coords_B[SPACEDIM*i_B+idim] -=  proj*linear_comb[idim];
            }
      
          //Buid the matrix sending  A into the Oxy plane and apply it to A and B  
          if(do_rotate)
            {
            TranslationRotationMatrix rotation;
            //rotate3DTriangle(&Coords_A[0], &Coords_A[SPACEDIM*i_A1], &Coords_A[SPACEDIM*i_A2], rotation);
            rotate3DTriangle(&Coords_B[0], &Coords_B[SPACEDIM*i_B1], &Coords_B[SPACEDIM*i_B2], rotation);
            for (int i=0; i<nb_NodesA; i++)    rotation.transform_vector(&Coords_A[SPACEDIM*i]);
            for (int i=0; i<nb_NodesB; i++)    rotation.transform_vector(&Coords_B[SPACEDIM*i]);
            }
        }
      }
    else
      {
        std::cout << " Maille d�g�n�r�e " << "epsilon = " << epsilon << std::endl;
				std::cout << " i_A1= " << i_A1 << " i_A2= " << i_A2 << std::endl;
				std::cout << " distance2<SPACEDIM>(&Coords_A[0],&Coords_A[i_A1])= " <<  distance2<SPACEDIM>(&Coords_A[0],&Coords_A[i_A1]) << std::endl;
				std::cout << "abs(normal_A) = " << fabs(normal_A[0]) << " ; " <<fabs( normal_A[1]) << " ; " << fabs(normal_A[2]) << std::endl;
				std::cout << " i_B1= " << i_B1 << " i_B2= " << i_B2 << std::endl; 
				std::cout << " distance2<SPACEDIM>(&Coords_B[0],&Coords_B[i_B1])= " <<  distance2<SPACEDIM>(&Coords_B[0],&Coords_B[i_B1]) << std::endl;
				std::cout << "normal_B = " << normal_B[0] << " ; " << normal_B[1] << " ; " << normal_B[2] << std::endl;
      }
  }
  
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  void PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::rotate3DTriangle(double* PP1, double*PP2, double*PP3,
                                                                                        TranslationRotationMatrix& rotation_matrix)
  {
    //initializes
    rotation_matrix.translate(PP1);
    
    double P2w[3];
    double P3w[3];
    P2w[0]=PP2[0]; P2w[1]=PP2[1];P2w[2]=PP2[2];
    P3w[0]=PP3[0]; P3w[1]=PP3[1];P3w[2]=PP3[2];
    
    // translating to set P1 at the origin
    for (int i=0; i<3; i++)
      {
      P2w[i]-=PP1[i];
      P3w[i]-=PP1[i];
      }
   
    // rotating to set P2 on the Oxy plane
    TranslationRotationMatrix A;
    A.rotate_x(P2w);
    A.rotate_vector(P3w);
    rotation_matrix.multiply(A);

    //rotating to set P2 on the Ox axis
    TranslationRotationMatrix B;
    B.rotate_z(P2w);
    B.rotate_vector(P3w);
    rotation_matrix.multiply(B);
  
    //rotating to set P3 on the Oxy plane
    TranslationRotationMatrix C;
    C.rotate_x(P3w);
    rotation_matrix.multiply(C);
  }
}

#endif
