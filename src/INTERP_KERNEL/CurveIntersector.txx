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
#ifndef __CURVEINTERSECTOR_TXX__
#define __CURVEINTERSECTOR_TXX__

#include "CurveIntersector.hxx"
#include "InterpolationUtils.hxx"

#include <limits>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  CurveIntersector<MyMeshType,MyMatrix>
  ::CurveIntersector(const MyMeshType& meshT, const MyMeshType& meshS,
                     double precision, double tolerance, double medianLine, int printLevel):
    _meshT(meshT),
    _meshS(meshS),
    _tolerance(tolerance),
    _precision(precision),
    _median_line(medianLine),
    _print_level(printLevel)
  {
    if ( SPACEDIM != 1 && SPACEDIM != 2 )
      throw Exception("CurveIntersector(): space dimension of mesh must be 1 or 2");
    if ( MESHDIM != 1 )
      throw Exception("CurveIntersector(): mesh dimension must be 1");

    _connectT = meshT.getConnectivityPtr();
    _connectS = meshS.getConnectivityPtr();
    _connIndexT = meshT.getConnectivityIndexPtr();
    _connIndexS = meshS.getConnectivityIndexPtr();
    _coordsT = meshT.getCoordinatesPtr();
    _coordsS = meshS.getCoordinatesPtr();
  }

  template<class MyMeshType, class MyMatrix>
  CurveIntersector<MyMeshType,MyMatrix>::~CurveIntersector()
  {
  }

  //================================================================================
  /*!
    \brief creates the bounding boxes for all the cells of mesh \a mesh

    \param mesh structure pointing to the mesh
    \param bbox vector containing the bounding boxes
  */
  //================================================================================

  template<class MyMeshType, class MyMatrix>
  void CurveIntersector<MyMeshType,MyMatrix>::createBoundingBoxes (const MyMeshType&    mesh,
                                                                   std::vector<double>& bbox)
  {
    long nbelems = mesh.getNumberOfElements();
    bbox.resize(2*SPACEDIM* nbelems);
    const double* coords = mesh.getCoordinatesPtr();
    const ConnType* conn = mesh.getConnectivityPtr();
    const ConnType* conn_index = mesh.getConnectivityIndexPtr();  
    int ibox=0;
    for(long icell=0; icell<nbelems; icell++)
      {
        int nb_nodes_per_elem = conn_index[icell+1]-conn_index[icell];
        //initializing bounding box limits
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            bbox[2*SPACEDIM*ibox+2*idim]   =  std::numeric_limits<double>::max();
            bbox[2*SPACEDIM*ibox+2*idim+1] = -std::numeric_limits<double>::max();
          }
        //updating the bounding box with each node of the element
        for (int j=0; j<nb_nodes_per_elem; j++)
          {
            const double* coord_node = coords + 
              SPACEDIM*OTT<ConnType,numPol>
              ::coo2C(conn[OTT<ConnType,numPol>::conn2C(conn_index[icell]+j)]);
            for(int idim=0; idim<SPACEDIM; idim++)
              {
                double x = *(coord_node+idim);
                bbox[ibox*2*SPACEDIM + 2*idim]   =
                  ( bbox[ibox*2*SPACEDIM + 2*idim] < x ) ? bbox[ibox*2*SPACEDIM + 2*idim] : x;
                bbox[ibox*2*SPACEDIM + 2*idim+1] =
                  ( bbox[ibox*2*SPACEDIM + 2*idim+1] > x ) ? bbox[ibox*2*SPACEDIM + 2*idim+1] : x;
              }
          }
        ibox++;
      }
  }

  /*!
    Computes the bouding box of a given element. iP in numPol mode.
  */
  template<class MyMeshType, class MyMatrix>
  void CurveIntersector<MyMeshType,MyMatrix>::getElemBB (double*           bb,
                                                         const MyMeshType& mesh,
                                                         ConnType          iP,
                                                         ConnType          nb_nodes)
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
        const double* coord_node = coords +
          SPACEDIM*(OTT<ConnType,numPol>::coo2C(conn[OTT<ConnType,numPol>::conn2C(conn_index[OTT<ConnType,numPol>::ind2C(iP)]+i)]));
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            double x = *(coord_node+idim);
            bb[2*idim  ] = (x<bb[2*idim  ]) ? x : bb[2*idim  ];
            bb[2*idim+1] = (x>bb[2*idim+1]) ? x : bb[2*idim+1];
          }
      }
  }
  
  /*!
   * \param [in] startOfSeg - input coming from intersectSegments or intersectSegmentsInternal
   * \param [in] endOfSeg - input coming from intersectSegments or intersectSegmentsInternal. Assume that endOfSeg>startOfSeg.
   * \param [in] pt - position of point that the method computes the bary coords for.
   */
  template<class MyMeshType, class MyMatrix>
  bool CurveIntersector<MyMeshType,MyMatrix>::ComputeBaryCoordsOf(double startOfSeg, double endOfSeg, double pt, double& startPos, double& endPos)
  {
    double deno(endOfSeg-startOfSeg);
    startPos=(endOfSeg-pt)/deno;
    endPos=1.-startPos;
    return startPos>=0. && endPos>=0.;
  }

  /*! Readjusts a set of bounding boxes so that they are extended
    in all dimensions for avoiding missing interesting intersections

    \param bbox vector containing the bounding boxes
  */
  template<class MyMeshType, class MyMatrix>
  void CurveIntersector<MyMeshType,MyMatrix>::adjustBoundingBoxes (std::vector<double>& bbox,
                                                                   double adjustmentEpsAbs)
  {
    long size = bbox.size()/(2*SPACEDIM);
    for (int i=0; i<size; i++)
      {
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            bbox[i*2*SPACEDIM+2*idim  ] -= adjustmentEpsAbs;
            bbox[i*2*SPACEDIM+2*idim+1] += adjustmentEpsAbs;
          }
      }
  }

  /*!
   * @param icellT id in target mesh in format of MyMeshType.
   * @param coordsT output val that stores coordinates of the target cell
   * automatically resized to the right length.
   * @return true if segment is quadratic and in this case coordinates of medium node
   * are placed in the middle of coordsT
   */
  template<class MyMeshType, class MyMatrix>
  bool CurveIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(ConnType icellT, std::vector<double>& coordsT) const
  {
    int nbNodesT(_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)+1] - _connIndexT[OTT<ConnType,numPol>::ind2C(icellT)]);
    coordsT.resize(SPACEDIM*nbNodesT);
    for (ConnType iT=0; iT<nbNodesT; iT++)
      {
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            coordsT[SPACEDIM*iT+idim] =
              _coordsT[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectT[OTT<ConnType,numPol>::conn2C(_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)]+iT)])+idim];
          }
      }
    if ( nbNodesT > 2 )
      {
        for(int idim=0; idim<SPACEDIM; idim++)
          std::swap( coordsT[SPACEDIM*1+idim], coordsT[SPACEDIM*2+idim]);
        return true;
      }
    return false;
  }
  
  template<class MyMeshType, class MyMatrix>
  typename MyMeshType::MyConnType CurveIntersector<MyMeshType,MyMatrix>::getNodeIdOfTargetCellAt(ConnType icellT, ConnType nodeIdInCellT) const
  {
    int nbNodesT(_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)+1] - _connIndexT[OTT<ConnType,numPol>::ind2C(icellT)]);
    if(nodeIdInCellT>=0 && nodeIdInCellT<nbNodesT)
      return OTT<ConnType,numPol>::coo2C(_connectT[OTT<ConnType,numPol>::conn2C(_connIndexT[OTT<ConnType,numPol>::ind2C(icellT)]+nodeIdInCellT)]);
    else
      throw Exception("getNodeIdOfTargetCellAt : error in nodeId in cell");
  }

  /*!
   * @param icellS id in source mesh in format of MyMeshType.
   * @param coordsS output val that stores coordinates of the source cell automatically resized to the right length.
   * @return true if segment is quadratic and in this case coordinates of medium node
   * are placed in the middle of coordsS
   */
  template<class MyMeshType, class MyMatrix>
  bool CurveIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(ConnType icellS, std::vector<double>& coordsS) const
  {
    int nbNodesS = _connIndexS[OTT<ConnType,numPol>::ind2C(icellS)+1] - _connIndexS[OTT<ConnType,numPol>::ind2C(icellS)];
    coordsS.resize(SPACEDIM*nbNodesS);
    for(ConnType iS=0; iS<nbNodesS; iS++)
      {
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            coordsS[SPACEDIM*iS+idim] =
              _coordsS[SPACEDIM*OTT<ConnType,numPol>::coo2C(_connectS[OTT<ConnType,numPol>::conn2C(_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)]+iS)])+idim];
          }
      }
    if ( nbNodesS > 2 )
      {
        for(int idim=0; idim<SPACEDIM; idim++)
          std::swap( coordsS[SPACEDIM*1+idim], coordsS[SPACEDIM*2+idim]);
        return true;
      }
    return false;
  }

  template<class MyMeshType, class MyMatrix>
  typename MyMeshType::MyConnType CurveIntersector<MyMeshType,MyMatrix>::getNodeIdOfSourceCellAt(ConnType icellS, ConnType nodeIdInCellS) const
  {
    int nbNodesS(_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)+1] - _connIndexS[OTT<ConnType,numPol>::ind2C(icellS)]);
    if(nodeIdInCellS>=0 && nodeIdInCellS<nbNodesS)
      return OTT<ConnType,numPol>::coo2C(_connectS[OTT<ConnType,numPol>::conn2C(_connIndexS[OTT<ConnType,numPol>::ind2C(icellS)]+nodeIdInCellS)]);
    else
      throw Exception("getNodeIdOfSourceCellAt : error in nodeId in cell");
  }

  /*!
   * \brief Return dual segments of given segment
   *  \param icell - given segment in C mode
   *  \param mesh - mesh
   *  \param segments - dual segments
   */
  template<class MyMeshType, class MyMatrix>
  void CurveIntersector<MyMeshType,MyMatrix>::getDualSegments(ConnType                   icell,
                                                              const MyMeshType&          mesh,
                                                              std::vector<TDualSegment>& segments)
  {
    // get coordinates of cell nodes
    int nbNodes;
    std::vector<double> ncoords;
    std::vector<int>    nodeIds;
    {
      const ConnType *connect   = mesh.getConnectivityPtr();
      const ConnType *connIndex = mesh.getConnectivityIndexPtr();
      const double *coords      = mesh.getCoordinatesPtr();

      nbNodes = connIndex[icell+1] - connIndex[icell];

      ncoords.resize(SPACEDIM*nbNodes);
      nodeIds.resize(nbNodes);

      for(int i=0; i<nbNodes; i++)
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            nodeIds[i] = connect[OTT<ConnType,numPol>::conn2C(connIndex[OTT<ConnType,numPol>::ind2C(icell)]+i)];
            ncoords[SPACEDIM*i+idim] = coords[SPACEDIM*OTT<ConnType,numPol>::coo2C(nodeIds[i])+idim];
          }
      if ( nbNodes > 2 ) // quadratic segment, put medium node in the middle
        {
          for(int idim=0; idim<SPACEDIM; idim++)
            std::swap( ncoords[SPACEDIM*1+idim], ncoords[SPACEDIM*2+idim]);
          std::swap( nodeIds[1], nodeIds[2] );
        }
    }

    // fill segments
    segments.clear();
    segments.reserve( 2*nbNodes );
    for(int i=0; i<nbNodes-1; i++)
      {
        segments.push_back(TDualSegment());
        TDualSegment& seg1 = segments.back();
        segments.push_back(TDualSegment());
        TDualSegment& seg2 = segments.back();

        seg1._nodeId = nodeIds[i];
        seg2._nodeId = nodeIds[i+1];

        seg1._coords.resize( SPACEDIM * 2 );
        seg2._coords.resize( SPACEDIM * 2 );

        for(int idim=0; idim<SPACEDIM; idim++)
          {
            double c1 = ncoords[SPACEDIM*i+idim];
            double c2 = ncoords[SPACEDIM*(i+1)+idim];
            double m = 0.5 * ( c1 + c2 );
            seg1._coords[ idim ] = c1;
            seg1._coords[ SPACEDIM + idim ] = m;
            seg2._coords[ idim ] = m;
            seg2._coords[ SPACEDIM + idim ] = c2;
          }
      }
  }

  template<class MyMeshType, class MyMatrix>
  bool CurveIntersector<MyMeshType,MyMatrix>::projectionThis(const double *coordsT, const double *coordsS,
                                                             double& xs0, double& xs1, double& xt0, double& xt1) const
  {
    xt0 = coordsT[0]; xt1 = coordsT[1];
    xs0 = coordsS[0]; xs1 = coordsS[1];
    if ( SPACEDIM == 2 )
      {
        // Pass 2D->1D

        enum { X=0, Y };

        // check if two segments overlap in 2D within tolerance

        const double* t0 = coordsT;
        const double* t1 = coordsT + 2;
        double t01[2] = { t1[X]-t0[X], t1[Y]-t0[Y] }; // tgt segment direction
        double tSize = sqrt( t01[X]*t01[X] + t01[Y]*t01[Y] ); // tgt segment size
        if ( tSize < _precision )
          return false; // degenerated segment
        t01[X] /= tSize, t01[Y] /= tSize; // normalize t01

        const double* s0 = coordsS;
        const double* s1 = coordsS + 2;
        double t0s0[2] = { s0[X]-t0[X], s0[Y]-t0[Y] };
        double t0s1[2] = { s1[X]-t0[X], s1[Y]-t0[Y] };
        double nt01_x_t0s0 = t0s0[X] * t01[Y] - t0s0[Y] * t01[X]; // t0s0 dot norm of t01
        double nt01_x_t0s1 = t0s1[X] * t01[Y] - t0s1[Y] * t01[X]; // t0s1 dot norm of t01
        double dist_ts0 = fabs( nt01_x_t0s0 ); // dist from tgt seg to s0
        double dist_ts1 = fabs( nt01_x_t0s1 ); // dist from tgt seg to s1
        bool s0_out_of_tol = ( dist_ts0 > _tolerance );
        bool s1_out_of_tol = ( dist_ts1 > _tolerance );
        if ( nt01_x_t0s0 * nt01_x_t0s1 > 0 && ( s0_out_of_tol || s1_out_of_tol ))
          return false; // tgt segment is to far from src segment

        double S0[2] = { s0[X], s0[Y] };
        double S1[2] = { s1[X], s1[Y] };
        if ( s0_out_of_tol ) // put s0 within tolerance
          {
            double t = _tolerance * nt01_x_t0s0 / dist_ts0; // signed tolerance
            double r = ( nt01_x_t0s0 - t ) / ( nt01_x_t0s0 - nt01_x_t0s1 );
            S0[X] = s0[X] * ( 1.-r ) + s1[X] * r;
            S0[Y] = s0[Y] * ( 1.-r ) + s1[Y] * r;
          }
        if ( s1_out_of_tol ) // put s1 within tolerance
          {
            double t = _tolerance * nt01_x_t0s1 / dist_ts1; // signed tolerance
            double r = ( nt01_x_t0s1 - t ) / ( nt01_x_t0s1 - nt01_x_t0s0 );
            S1[X] = s1[X] * ( 1.-r ) + s0[X] * r;
            S1[Y] = s1[Y] * ( 1.-r ) + s0[Y] * r;
          }

        // project tgt and src segments to median line

        double s01[2] = { S1[X]-S0[X], S1[Y]-S0[Y] }; // src segment direction
        double sSize = sqrt( s01[X]*s01[X] + s01[Y]*s01[Y] ); // src segment size
        if ( sSize < _precision )
          return false; // degenerated segment
        s01[X] /= sSize, s01[Y] /= sSize; // normalize s01

        // make t01 and s01 codirected
        double t01_x_s01 = t01[X] * s01[X] + t01[Y] * s01[Y]; // t01 dot s01
        if ( t01_x_s01 < 0 )
          s01[X] = -s01[X], s01[Y] = -s01[Y];

        double medianDir[2] = {
          t01[X] * ( 1.-_median_line) + s01[X] * _median_line, 
          t01[Y] * ( 1.-_median_line) + s01[Y] * _median_line
        };
        double medianSize = sqrt( medianDir[X]*medianDir[X] + medianDir[Y]*medianDir[Y] );
        if ( medianSize < std::numeric_limits<double>::min() )
          return false; // strange...
        medianDir[X] /= medianSize, medianDir[Y] /= medianSize;

        xt0 = t0[X] * medianDir[X] + t0[Y] * medianDir[Y];
        xt1 = t1[X] * medianDir[X] + t1[Y] * medianDir[Y];
        xs0 = S0[X] * medianDir[X] + S0[Y] * medianDir[Y];
        xs1 = S1[X] * medianDir[X] + S1[Y] * medianDir[Y];

      } // if ( SPACEDIM == 2 )
    return true;
  }
  
  /*!
   * \brief Return length of intersection of two segments
   */
  template<class MyMeshType, class MyMatrix>
  double CurveIntersector<MyMeshType,MyMatrix>::intersectSegmentsInternal(const double *coordsT, const double *coordsS, double& xs0, double& xs1, double& xt0, double& xt1) const
  {
    if(!projectionThis(coordsT,coordsS,xs0,xs1,xt0,xt1))
      return 0.;
    
    if ( xt0 > xt1 ) std::swap( xt0, xt1 );
    if ( xs0 > xs1 ) std::swap( xs0, xs1 );

    double x0 = std::max( xt0, xs0 );
    double x1 = std::min( xt1, xs1 );
    return ( x0 < x1 ) ? ( x1 - x0 ) : 0.;
  }
  
  /*!
   * \brief Return length of intersection of two segments
   */
  template<class MyMeshType, class MyMatrix>
  double CurveIntersector<MyMeshType,MyMatrix>::intersectSegments(const double *coordsT, const double *coordsS) const
  {
    double xs0,xs1,xt0,xt1;
    return intersectSegmentsInternal(coordsT,coordsS,xs0,xs1,xt0,xt1);
  }

}

#endif
