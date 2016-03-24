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
#ifndef __CONVEXINTERSECTOR_TXX__
#define __CONVEXINTERSECTOR_TXX__

#include "ConvexIntersector.hxx"
#include "PlanarIntersectorP0P0.txx"
#include "PlanarIntersectorP0P1.txx"
#include "PlanarIntersectorP1P0.txx"
#include "PlanarIntersectorP1P1.txx"
#include "PlanarIntersectorP1P0Bary.txx"
#include "PlanarIntersectorP0P1Bary.txx"

#include "PolygonAlgorithms.txx"

#include <iostream>

#define CONVINTERSECTOR_TEMPLATE template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
#define CONVEX_INTERSECTOR_ ConvexIntersector<MyMeshType,MyMatrix,InterpType>

namespace INTERP_KERNEL
{
  CONVINTERSECTOR_TEMPLATE
  CONVEX_INTERSECTOR_::ConvexIntersector(const MyMeshType& meshT, const MyMeshType& meshS, 
                                         double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf,
                                         double medianPlane, bool doRotate , int oriantation, int printLevel)
    :InterpType<MyMeshType,MyMatrix,CONVEX_INTERSECTOR_ >(meshT,meshS,dimCaracteristic, precision, md3DSurf, minDot3DSurf, medianPlane, doRotate, oriantation, printLevel),
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

  CONVINTERSECTOR_TEMPLATE
  double CONVEX_INTERSECTOR_::intersectGeometry(ConnType icellT,   ConnType icellS,
                                                ConnType nbNodesT, ConnType nbNodesS)
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

  CONVINTERSECTOR_TEMPLATE
  double CONVEX_INTERSECTOR_::intersectGeometryWithQuadrangle(const double             * quadrangle,
                                                              const std::vector<double>& sourceCoords,
                                                              bool                       isSourceQuad)
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

  CONVINTERSECTOR_TEMPLATE
  double CONVEX_INTERSECTOR_::intersectGeometryGeneral(const std::vector<double>& targetCoords,
                                                       const std::vector<double>& sourceCoords)
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

  //================================================================================
  /*!
   * \brief Intersect a triangle and a polygon for P1P0 barycentric algorithm
   *  \param targetCell - list of coordinates of target polygon in full interlace
   *  \param targetCellQuadratic - specifies if target polygon is quadratic or not
   *  \param sourceTria - list of coordinates of source triangle
   *  \param res - coefficients a,b and c associated to nodes of sourceTria
   */
  //================================================================================

  CONVINTERSECTOR_TEMPLATE
  double CONVEX_INTERSECTOR_::intersectGeoBary(const std::vector<double>& targetCell,
                                               bool                       targetCellQuadratic,
                                               const double *             sourceTria,
                                               std::vector<double>&       res)
  {
    double area = 0;
    double barycenter[SPACEDIM] = {0., 0.};
    int nbOfNodesT=targetCell.size()/SPACEDIM;

    /*** Compute the intersection area ***/
    INTERP_KERNEL::PolygonAlgorithms<SPACEDIM> P(_epsilon, PlanarIntersector<MyMeshType,MyMatrix>::_precision);
    std::deque<double> inter =  P.intersectConvexPolygons(sourceTria, &targetCell[0], 3, nbOfNodesT);
    double cross[SPACEDIM];
    int nb_inter =((int)inter.size())/SPACEDIM;
    for(int i = 1; i<nb_inter-1; i++)
    {
      INTERP_KERNEL::crossprod<SPACEDIM>(&inter[0],&inter[SPACEDIM*i],&inter[SPACEDIM*(i+1)],cross);
      area += 0.5*norm<SPACEDIM>(cross);
      barycenter[0] += inter[SPACEDIM*i];
      barycenter[1] += inter[SPACEDIM*i+1];
    }
    if ( area > std::numeric_limits<double>::min() )
    {
      barycenter[0] = ( barycenter[0] + inter[0] + inter[SPACEDIM*(nb_inter-1)]  ) / nb_inter;
      barycenter[1] = ( barycenter[1] + inter[1] + inter[SPACEDIM*(nb_inter-1)+1]) / nb_inter;
      res.resize(3);
      barycentric_coords<2>( sourceTria, &barycenter[0], &res[0]);
      res[0] *= area;
      res[1] *= area;
      res[2] *= area;
    }
    else
    {
      area = 0;
    }
    return area;
  }
}

#endif
