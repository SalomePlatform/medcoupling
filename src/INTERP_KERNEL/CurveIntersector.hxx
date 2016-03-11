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

#ifndef __CURVEINTERSECTOR_HXX__
#define __CURVEINTERSECTOR_HXX__

#include "TargetIntersector.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  class CurveIntersector : public TargetIntersector<MyMeshType,MyMatrix>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    CurveIntersector(const MyMeshType& meshT, const MyMeshType& meshS,
                     double  precision, double adjustmentEpsAbs, double medianLine, int printLevel);
    virtual ~CurveIntersector();
    void createBoundingBoxes(const MyMeshType& mesh, std::vector<double>& bbox);
    void adjustBoundingBoxes(std::vector<double>& bbox, double adjustmentEpsAbs);
    static void getElemBB(double* bb, const MyMeshType& mesh, ConnType iP, ConnType nb_nodes);
    static bool ComputeBaryCoordsOf(double startOfSeg, double endOfSeg, double pt, double& startPos, double& endPos);
  protected :
    bool projectionThis(const double *coordsT, const double *coordsS, double& xs0, double& xs1, double& xt0, double& xt1) const;
    bool getRealTargetCoordinates(ConnType icellT, std::vector<double>& coordsT) const;
    typename MyMeshType::MyConnType getNodeIdOfTargetCellAt(ConnType icellT, ConnType nodeIdInCellT) const;
    bool getRealSourceCoordinates(ConnType icellS, std::vector<double>& coordsS) const;
    typename MyMeshType::MyConnType getNodeIdOfSourceCellAt(ConnType icellT, ConnType nodeIdInCellT) const;
    double intersectSegments(const double *coordsT, const double *coordsS) const;
    double intersectSegmentsInternal(const double *coordsT, const double *coordsS, double& xs0, double& xs1, double& xt0, double& xt1) const;
    
    struct TDualSegment
    {
      std::vector<double> _coords;
      int                 _nodeId; // in mesh mode
    };
    static void getDualSegments(ConnType                   icell,
                                const MyMeshType&          mesh,
                                std::vector<TDualSegment>& segments);

  protected:
    const ConnType *_connectT;
    const ConnType *_connectS;
    const double *_coordsT;
    const double *_coordsS;
    const ConnType *_connIndexT;
    const ConnType *_connIndexS;
    const MyMeshType& _meshT;
    const MyMeshType& _meshS;
    double _tolerance;
    double _precision;
    double _median_line;
    int _print_level;
  };
}

#endif
