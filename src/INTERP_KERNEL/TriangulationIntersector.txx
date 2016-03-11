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
#ifndef __TRIANGULATIONINTERSECTOR_TXX__
#define __TRIANGULATIONINTERSECTOR_TXX__

#include "TriangulationIntersector.hxx"
#include "PlanarIntersectorP0P0.txx"
#include "PlanarIntersectorP0P1.txx"
#include "PlanarIntersectorP1P0.txx"
#include "PlanarIntersectorP1P0Bary.txx"

#include "InterpolationUtils.hxx"
#include "PlanarIntersector.hxx"

#include <iostream>

#define TRI_INTERSECTOR TriangulationIntersector<MyMeshType,MyMatrix,InterpType>
#define TRI_INTER_TEMPLATE template<class MyMeshType, class MyMatrix, \
                    template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>

namespace INTERP_KERNEL
{
  TRI_INTER_TEMPLATE
  TRI_INTERSECTOR::TriangulationIntersector(const MyMeshType& meshT, const MyMeshType& meshS, 
                                            double DimCaracteristic, double Precision, double md3DSurf, double minDot3DSurf,
                                            double MedianPlane, int orientation, int PrintLevel)
    :InterpType<MyMeshType,MyMatrix,TRI_INTERSECTOR >(meshT,meshS,DimCaracteristic, Precision, md3DSurf, minDot3DSurf,
                                                      MedianPlane, true, orientation, PrintLevel)
  {
    if(PlanarIntersector<MyMeshType,MyMatrix>::_print_level >= 1)
      {
        std::cout << "  - intersection type = triangles " << std::endl;
        if(SPACEDIM==3) std::cout << "_do_rotate = true"<< std::endl;
      }
  }
  
  TRI_INTER_TEMPLATE
  double TRI_INTERSECTOR::intersectGeometry(ConnType icellT,   ConnType icellS,
                                            ConnType nbNodesT, ConnType nbNodesS)
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

  TRI_INTER_TEMPLATE
  double TRI_INTERSECTOR::intersectGeometryWithQuadrangle(const double             * quadrangle,
                                                          const std::vector<double>& sourceCoords,
                                                          bool                       isSourceQuad)
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

  TRI_INTER_TEMPLATE
  double TRI_INTERSECTOR::intersectGeometryGeneral(const std::vector<double>& targetCoords,
                                                   const std::vector<double>& sourceCoords)
  {
    double result = 0.;
    ConnType nbNodesS=sourceCoords.size()/SPACEDIM;
    ConnType nbNodesT=targetCoords.size()/SPACEDIM;
    //Compute the intersection area
    double area[SPACEDIM];
    for(ConnType iT = 1; iT<nbNodesT-1; iT++)
      {
        for(ConnType iS = 1; iS<nbNodesS-1; iS++)
          {
            std::vector<double> inter;
            INTERP_KERNEL::intersec_de_triangle(&targetCoords[0],&targetCoords[SPACEDIM*iT],&targetCoords[SPACEDIM*(iT+1)],
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
          }
      }
    return result;
  }

  //================================================================================
  /*!
   * \brief Intersect a triangle and a polygon for P1P0 barycentric algorithm
   *  \param targetCell - list of coordinates of target polygon in full interlace
   *  \param targetCellQuadratic - specifies if target polygon is quadratic or not
   *  \param sourceTria - list of coordinates of source triangle
   *  \param res - coefficients a,b and c associated to nodes of sourceTria
   */
  //================================================================================

  TRI_INTER_TEMPLATE
  double TRI_INTERSECTOR::intersectGeoBary(const std::vector<double>& targetCell,
                                           bool                       targetCellQuadratic,
                                           const double *             sourceTria,
                                           std::vector<double>&       res)
  {
    std::vector<const double* > sourceCell(3);
    sourceCell[0] = &sourceTria[0];
    sourceCell[1] = &sourceTria[SPACEDIM];
    sourceCell[2] = &sourceTria[SPACEDIM*2];

    //Compute the intersection area
    double inter_area[SPACEDIM], total_area = 0.;
    double total_barycenter[SPACEDIM]={0.,0.};

    const ConnType nbNodesT=targetCell.size()/SPACEDIM;
    for(ConnType iT = 1; iT<nbNodesT-1; iT++)
    {
      std::vector<double> inter;
      INTERP_KERNEL::intersec_de_triangle(&targetCell[0],&targetCell[SPACEDIM*iT],&targetCell[SPACEDIM*(iT+1)],
                                          sourceCell[0], sourceCell[1], sourceCell[2],
                                          inter, PlanarIntersector<MyMeshType,MyMatrix>::_dim_caracteristic,
                                          PlanarIntersector<MyMeshType,MyMatrix>::_precision);
      ConnType nb_inter=((ConnType)inter.size())/2;
      if(nb_inter >3) inter=reconstruct_polygon(inter);
      for(ConnType i = 1; i<nb_inter-1; i++)
      {
        INTERP_KERNEL::crossprod<2>(&inter[0],&inter[2*i],&inter[2*(i+1)],inter_area);
        inter_area[0] = 0.5 * fabs( inter_area[0] );
        total_area += inter_area[0];
        std::vector<double> inter_bary=INTERP_KERNEL::bary_poly(inter);
        total_barycenter[0] += inter_area[0] * inter_bary[0];
        total_barycenter[1] += inter_area[0] * inter_bary[1];
      }
    }
    if ( total_area > std::numeric_limits<double>::min() )
    {
      total_barycenter[0] /= total_area;
      total_barycenter[1] /= total_area;
      res.resize(3);
      barycentric_coords( sourceCell, &total_barycenter[0], &res[0]);
      res[0] *= total_area;
      res[1] *= total_area;
      res[2] *= total_area;
    }
    else
    {
      total_area = 0;
    }
    return total_area;
  }

}

#endif
