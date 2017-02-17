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

#ifndef __PLANARINTERSECTOR_HXX__
#define __PLANARINTERSECTOR_HXX__

#include "TargetIntersector.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <map>
#include <set>

namespace INTERP_KERNEL
{
  class TranslationRotationMatrix;
  
  template<class MyMeshType, class MyMatrix>
  class PlanarIntersector : public TargetIntersector<MyMeshType,MyMatrix>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
    typedef typename std::map<int,std::set<int> > DuplicateFacesType;
  public:
    PlanarIntersector(const MyMeshType& meshT, const MyMeshType& meshS, double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf, double medianPlane, bool doRotate, int orientation, int printLevel);
    virtual ~PlanarIntersector();
    void createBoundingBoxes(const MyMeshType& mesh, std::vector<double>& bbox);
    void adjustBoundingBoxes(std::vector<double>& bbox, double surf3DAdjustmentEps, double surf3DAdjustmentEpsAbs);
    inline void getElemBB(double* bb, const MyMeshType& mesh, ConnType iP, ConnType nb_nodes);
    static int Projection(double *Coords_A, double *Coords_B,
                          int nb_NodesA, int nb_NodesB, double epsilon, double md3DSurf, double minDot3DSurf, double median_plane, bool do_rotate);
    virtual const DuplicateFacesType* getIntersectFaces() const
    {
      return NULL;
    }
  protected :
    int projectionThis(double *Coords_A, double *Coords_B, int nb_NodesA, int nb_NodesB);
    void getRealTargetCoordinates(ConnType icellT, std::vector<double>& coordsT);
    void getRealSourceCoordinates(ConnType icellS, std::vector<double>& coordsS);
    void getRealTargetCoordinatesPermute(ConnType icellT, int offset, std::vector<double>& coordsT);
    void getRealSourceCoordinatesPermute(ConnType icellS, int offset, std::vector<double>& coordsS);
    void getRealCoordinates(ConnType icellT, ConnType icellS, ConnType nbNodesT, ConnType nbNodesS, std::vector<double>& coordsT, std::vector<double>& coordsS, int& orientation);
    double getValueRegardingOption(double val) const;
    static void Rotate3DTriangle( double* PP1, double*PP2, double*PP3,
                                  TranslationRotationMatrix& rotation_matrix);
  protected:
    const ConnType *_connectT;
    const ConnType *_connectS;
    const double *_coordsT;
    const double *_coordsS;
    const ConnType *_connIndexT;
    const ConnType *_connIndexS;
    const MyMeshType& _meshT;
    const MyMeshType& _meshS;
    double _dim_caracteristic;
    double _max_distance_3Dsurf_intersect;
    double _min_dot_btw_3Dsurf_intersect;
    double _precision;
    double _median_plane;
    bool _do_rotate;
    int _orientation;
    int _print_level;
  };
}

#endif
