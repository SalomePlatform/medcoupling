// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
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
// Author : Anthony Geay (CEA/DEN)
#ifndef __PLANARINTERSECTOR_TXX__
#define __PLANARINTERSECTOR_TXX__

#include "PlanarIntersector.hxx"
#include "InterpolationUtils.hxx"
#include "TranslationRotationMatrix.hxx"

#include <iostream>
#include <limits>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  PlanarIntersector<MyMeshType,MyMatrix>::PlanarIntersector(const MyMeshType& meshT, const MyMeshType& meshS, double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf, double medianPlane, bool doRotate, int orientation, int printLevel):
    _meshT(meshT),_meshS(meshS),
    _dim_caracteristic(dimCaracteristic),_max_distance_3Dsurf_intersect(md3DSurf),_min_dot_btw_3Dsurf_intersect(minDot3DSurf),_precision(precision),_median_plane(medianPlane),
    _do_rotate(doRotate),_orientation(orientation),_print_level(printLevel)
  {
    _connectT=meshT.getConnectivityPtr();
    _connectS=meshS.getConnectivityPtr();
    _connIndexT=meshT.getConnectivityIndexPtr();
    _connIndexS=meshS.getConnectivityIndexPtr();
    _coordsT=meshT.getCoordinatesPtr();
    _coordsS=meshS.getCoordinatesPtr();
  }

  template<class MyMeshType, class MyMatrix>
  PlanarIntersector<MyMeshType,MyMatrix>::~PlanarIntersector()
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
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::createBoundingBoxes(const MyMeshType& mesh, std::vector<double>& bbox)
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
            bbox[2*SPACEDIM*ibox+2*idim]   =  std::numeric_limits<double>::max();
            bbox[2*SPACEDIM*ibox+2*idim+1] = -std::numeric_limits<double>::max();
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
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::getElemBB(double* bb, const MyMeshType& mesh, ConnType iP, ConnType nb_nodes)
  {
    const double* coords = mesh.getCoordinatesPtr();
    const ConnType* conn_index = mesh.getConnectivityIndexPtr();
    const ConnType* conn = mesh.getConnectivityPtr();
    //initializing bounding box limits
    for(int idim=0; idim<SPACEDIM; idim++)
      {
        bb[2*idim  ] =  std::numeric_limits<double>::max();
        bb[2*idim+1] = -std::numeric_limits<double>::max();
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
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::adjustBoundingBoxes(std::vector<double>& bbox, double surf3DAdjustmentEps, double surf3DAdjustmentEpsAbs)
  {
    /* We build the segment tree for locating possible matching intersections*/
  
    long size = bbox.size()/(2*SPACEDIM);
    for (int i=0; i<size; i++)
      {
        double max=- std::numeric_limits<double>::max();
        for(int idim=0; idim<SPACEDIM; idim++)
          {            
            double Dx=bbox[i*2*SPACEDIM+1+2*idim]-bbox[i*2*SPACEDIM+2*idim];
            max=(max<Dx)?Dx:max;
          }
        for(int idim=0; idim<SPACEDIM; idim++)
          {            
            bbox[i*2*SPACEDIM+2*idim  ] -= surf3DAdjustmentEps*max+surf3DAdjustmentEpsAbs;
            bbox[i*2*SPACEDIM+2*idim+1] += surf3DAdjustmentEps*max+surf3DAdjustmentEpsAbs;
          }
      }
  }

  /*!
   * @param icellT id in target mesh in format of MyMeshType.
   * @param coordsT output val that stores coordinates of the target cell automatically resized to the right length.
   */
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(ConnType icellT, std::vector<double>& coordsT)
  {
    int nbNodesT=_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)+1]-_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)];
    coordsT.resize(SPACEDIM*nbNodesT);
    for (ConnType iT=0; iT<nbNodesT; iT++)
      for(int idim=0; idim<SPACEDIM; idim++)
        coordsT[SPACEDIM*iT+idim]=_coordsT[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectT[OTT<ConnType,numPol>::conn2C(_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)]+iT)])+idim];
  }

  /*!
   * @param icellS id in source mesh in format of MyMeshType.
   * @param coordsS output val that stores coordinates of the source cell automatically resized to the right length.
   */
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(ConnType icellS, std::vector<double>& coordsS)
  {
    int nbNodesS=_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)+1]-_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)];
    coordsS.resize(SPACEDIM*nbNodesS);
    for (ConnType iS=0; iS<nbNodesS; iS++)
      for(int idim=0; idim<SPACEDIM; idim++)
        coordsS[SPACEDIM*iS+idim]=_coordsS[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectS[OTT<ConnType,numPol>::conn2C(_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)]+iS)])+idim];
  }

  /*!
   * @param icellT id in target mesh in format of MyMeshType.
   * @param offset is a value in C format that indicates the number of circular permutation.
   * @param coordsT output val that stores coordinates of the target cell automatically resized to the right length.
   */
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinatesPermute(ConnType icellT, int offset, std::vector<double>& coordsT)
  {
    int nbNodesT=_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)+1]-_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)];
    coordsT.resize(SPACEDIM*nbNodesT);
    for (ConnType iTTmp=0; iTTmp<nbNodesT; iTTmp++)
      {
        ConnType iT=(iTTmp+offset)%nbNodesT;
        for(int idim=0; idim<SPACEDIM; idim++)
          coordsT[SPACEDIM*iTTmp+idim]=_coordsT[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectT[OTT<ConnType,numPol>::conn2C(_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)]+iT)])+idim];
      }
  }

  /*!
   * @param icellS id in source mesh in format of MyMeshType.
   * @param offset is a value in C format that indicates the number of circular permutation.
   * @param coordsS output val that stores coordinates of the source cell automatically resized to the right length.
   */
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinatesPermute(ConnType icellS, int offset, std::vector<double>& coordsS)
  {
    int nbNodesS=_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)+1]-_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)];
    coordsS.resize(SPACEDIM*nbNodesS);
    for (ConnType iSTmp=0; iSTmp<nbNodesS; iSTmp++)
      {
        ConnType iS=(iSTmp+offset)%nbNodesS;
        for(int idim=0; idim<SPACEDIM; idim++)
          coordsS[SPACEDIM*iSTmp+idim]=_coordsS[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectS[OTT<ConnType,numPol>::conn2C(_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)]+iS)])+idim];
      }
  }
  
  /*!
   * @param icellT id in target mesh in format of MyMeshType.
   * @param icellS id in source mesh in format of MyMeshType.
   * @param nbNodesT nb of nodes of the target cell.
   * @param nbNodesS nb of nodes of the source cell.
   * @param coordsT output val that stores coordinates of the target cell automatically resized to the right length.
   * @param coordsS output val that stores coordinates of the source cell automatically resized to the right length.
   * @param orientation is an output value too, only set if SPACEDIM==3.
   */
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::getRealCoordinates(ConnType icellT, ConnType icellS, ConnType nbNodesT, ConnType nbNodesS, std::vector<double>& coordsT, std::vector<double>& coordsS, int& orientation)
  {
    coordsT.resize(SPACEDIM*nbNodesT);
    coordsS.resize(SPACEDIM*nbNodesS);
    for (int idim=0; idim<SPACEDIM; idim++)
      {
        for (ConnType iT=0; iT<nbNodesT; iT++)
          coordsT[SPACEDIM*iT+idim] = _coordsT[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectT[OTT<ConnType,numPol>::conn2C(_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)]+iT)])+idim];
        for (ConnType iS=0; iS<nbNodesS; iS++)
          coordsS[SPACEDIM*iS+idim] = _coordsS[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectS[OTT<ConnType,numPol>::conn2C(_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)]+iS)])+idim];
      }
    
    //project cells T and S on the median plane and rotate the median plane
    if(SPACEDIM==3)
      orientation = projectionThis(&coordsT[0], &coordsS[0], nbNodesT, nbNodesS);
    
    //DEBUG PRINTS
    if(_print_level >= 3) 
      {
        std::cout << std::endl << "Cell coordinates (possibly after projection)" << std::endl;
        std::cout << std::endl << "icellT= " << icellT << ", nb nodes T= " <<  nbNodesT << std::endl;
        for(int iT =0; iT< nbNodesT; iT++)
          {for (int idim=0; idim<SPACEDIM; idim++) std::cout << coordsT[SPACEDIM*iT+idim] << " "; std::cout << std::endl;}
        std::cout << std::endl << "icellS= " << icellS << ", nb nodes S= " <<  nbNodesS << std::endl;
        for(int iS =0; iS< nbNodesS; iS++)
          {for (int idim=0; idim<SPACEDIM; idim++) std::cout << coordsS[SPACEDIM*iS+idim]<< " "; std::cout << std::endl; }
      }
  }
  
  /*!
   * Filtering out zero surfaces and badly oriented surfaces
   * _orientation = -1,0,1,2
   * -1 : the intersection is taken into account if target and cells have different orientation
   * 0 : the intersection is always taken into account
   * 1 : the intersection is taken into account if target and cells have the same orientation
   * 2 : the absolute value of intersection is always taken into account 
   */
  template<class MyMeshType, class MyMatrix>
  double PlanarIntersector<MyMeshType,MyMatrix>::getValueRegardingOption(double val) const
  {
    if(_orientation==0)
      return val;
    if(_orientation==2)
      return fabs(val);
    if (( val > 0.0 && _orientation==1) || ( val < 0.0 && _orientation==-1 ))
      return _orientation*val;
    return 0.;
  }

  template<class MyMeshType, class MyMatrix>
  int PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(double *Coords_A, double *Coords_B, int nb_NodesA, int nb_NodesB)
  {
    return Projection(Coords_A,Coords_B,nb_NodesA,nb_NodesB,_dim_caracteristic*_precision,_max_distance_3Dsurf_intersect,_min_dot_btw_3Dsurf_intersect,_median_plane,_do_rotate);
  }

  template<class MyMeshType, class MyMatrix>
  int PlanarIntersector<MyMeshType,MyMatrix>::Projection(double *Coords_A, double *Coords_B, 
                                                         int nb_NodesA, int nb_NodesB, double epsilon, double md3DSurf, double minDot3DSurf, double median_plane, bool do_rotate)
  {
    double normal_A[3]={0,0,0};
    double normal_B[3]={0,0,0};
    double linear_comb[3];
    double proj;
    bool same_orientation;

    //Find the normal to cells A and B
    int i_A1(1);
    while(i_A1<nb_NodesA && distance2<SPACEDIM>(Coords_A,&Coords_A[SPACEDIM*i_A1])< epsilon)
      i_A1++;
    int i_A2(i_A1+1);
    crossprod<SPACEDIM>(Coords_A, &Coords_A[SPACEDIM*i_A1], &Coords_A[SPACEDIM*i_A2],normal_A);
    double normA(sqrt(dotprod<SPACEDIM>(normal_A,normal_A)));
    while(i_A2<nb_NodesA && normA < epsilon)
      {
        crossprod<SPACEDIM>(Coords_A, &Coords_A[SPACEDIM*i_A1], &Coords_A[SPACEDIM*i_A2],normal_A);
        i_A2++;
        normA = sqrt(dotprod<SPACEDIM>(normal_A,normal_A));

      }
    int i_B1(1);
    while(i_B1<nb_NodesB && distance2<SPACEDIM>(Coords_B,Coords_B+SPACEDIM*i_B1)< epsilon)
      i_B1++;
    int i_B2(i_B1+1);
    crossprod<SPACEDIM>(Coords_B, Coords_B+SPACEDIM*i_B1, Coords_B+SPACEDIM*i_B2,normal_B);
    double normB(sqrt(dotprod<SPACEDIM>(normal_B,normal_B)));
    while(i_B2<nb_NodesB && normB < epsilon)
      {
        crossprod<SPACEDIM>(Coords_B, Coords_B+SPACEDIM*i_B1, Coords_B+SPACEDIM*i_B2,normal_B);
        i_B2++;
        normB = sqrt(dotprod<SPACEDIM>(normal_B,normal_B));
      }

    //fabien option
    if(md3DSurf>0.)
      {
        double coords_GA[3];
        for(int i=0;i<3;i++)
          {
            coords_GA[i]=0.;
            for (int j=0;j<nb_NodesA;j++)
              coords_GA[i]+=Coords_A[3*j+i];
            coords_GA[i]/=nb_NodesA;
          }
        double G1[3],G2[3],G3[3];
        for(int i=0;i<3;i++)
          {
            G1[i]=Coords_B[i]-coords_GA[i];
            G2[i]=Coords_B[i+3]-coords_GA[i];
            G3[i]=Coords_B[i+6]-coords_GA[i];
          }
        double prodvect[3];
        prodvect[0]=G1[1]*G2[2]-G1[2]*G2[1];
        prodvect[1]=G1[2]*G2[0]-G1[0]*G2[2];
        prodvect[2]=G1[0]*G2[1]-G1[1]*G2[0];
        double prodscal=prodvect[0]*G3[0]+prodvect[1]*G3[1]+prodvect[2]*G3[2];
        if(fabs(prodscal)>md3DSurf)
          return 0;
      }
    if(i_A2<nb_NodesA && i_B2<nb_NodesB)
      {
        //Build the normal of the median plane

        double dotProd(dotprod<SPACEDIM>(normal_A,normal_B)/(normA*normB));

        if(fabs(dotProd)<minDot3DSurf)
          return 0;

        same_orientation=(dotProd>=0);
        
        if(!same_orientation)
          for(int idim =0; idim< SPACEDIM; idim++)
            normal_A[idim] *=-1;
        
        double normBB(sqrt(dotprod<SPACEDIM>(normal_B,normal_B)));
        
        for(int idim =0; idim< SPACEDIM; idim++)
          linear_comb[idim] = median_plane*normal_A[idim]/normA + (1-median_plane)*normal_B[idim]/normBB;
        double norm= sqrt(dotprod<SPACEDIM>(linear_comb,linear_comb));

        //Necessarily: norm>epsilon, no need to check
        for(int idim =0; idim< SPACEDIM; idim++)
          linear_comb[idim]/=norm;
        
        //Project the nodes of A and B on the median plane
        for(int i_A=0; i_A<nb_NodesA; i_A++)
          {
            proj = dotprod<SPACEDIM>(&Coords_A[SPACEDIM*i_A],linear_comb);
            for(int idim =0; idim< SPACEDIM; idim++)
              Coords_A[SPACEDIM*i_A+idim] -=  proj*linear_comb[idim];
          }
        for(int i_B=0; i_B<nb_NodesB; i_B++)
          {
            proj = dotprod<SPACEDIM>(Coords_B+SPACEDIM*i_B,linear_comb);
            for(int idim =0; idim< SPACEDIM; idim++)
              Coords_B[SPACEDIM*i_B+idim] -=  proj*linear_comb[idim];
          }
        
        //Buid the matrix sending  A into the Oxy plane and apply it to A and B  
        if(do_rotate)
          {
            TranslationRotationMatrix rotation;
            //rotate3DTriangle(Coords_A, &Coords_A[SPACEDIM*i_A1], &Coords_A[SPACEDIM*i_A2], rotation);
            Rotate3DTriangle(Coords_B, Coords_B+SPACEDIM*i_B1, Coords_B+SPACEDIM*i_B2, rotation);
            for (int i=0; i<nb_NodesA; i++)
              rotation.transform_vector(Coords_A+SPACEDIM*i);
            for (int i=0; i<nb_NodesB; i++)
              rotation.transform_vector(Coords_B+SPACEDIM*i);
          }
        return same_orientation?1:-1;
      }
    else
      {
        std::cout << " Degenerated cell " << "epsilon = " << epsilon << std::endl;
        std::cout << " i_A1= " << i_A1 << " i_A2= " << i_A2 << std::endl;
        std::cout << " distance2<SPACEDIM>(Coords_A,&Coords_A[i_A1])= " <<  distance2<SPACEDIM>(Coords_A,&Coords_A[i_A1]) << std::endl;
        std::cout << "abs(normal_A) = " << fabs(normal_A[0]) << " ; " <<fabs( normal_A[1]) << " ; " << fabs(normal_A[2]) << std::endl;
        std::cout << " i_B1= " << i_B1 << " i_B2= " << i_B2 << std::endl; 
        std::cout << " distance2<SPACEDIM>(&Coords_B[0],&Coords_B[i_B1])= " <<  distance2<SPACEDIM>(Coords_B,Coords_B+i_B1) << std::endl;
        std::cout << "normal_B = " << normal_B[0] << " ; " << normal_B[1] << " ; " << normal_B[2] << std::endl;
        return 1;
      }
  }
  
  template<class MyMeshType, class MyMatrix>
  void PlanarIntersector<MyMeshType,MyMatrix>::Rotate3DTriangle(double* PP1, double*PP2, double*PP3,
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
