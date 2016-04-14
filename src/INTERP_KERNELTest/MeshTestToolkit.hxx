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

#ifndef __TU_MESH_TEST_TOOLKIT_HXX__
#define __TU_MESH_TEST_TOOLKIT_HXX__

#include "Interpolation3D.hxx"
#include "Interpolation3D.txx"
#include "InterpolationPlanar.hxx"

#include <vector>
#include <map>

#define ERR_TOL 1.0e-8

typedef std::vector<std::map<int,double> > IntersectionMatrix;

namespace INTERP_KERNEL
{
  class Interpolation3D;
}


namespace MEDCoupling
{
  class MEDCouplingUMesh;
}

namespace INTERP_TEST
{
  /**
   * \brief Class providing services for mesh intersection tests.
   *
   */
  template<int SPACEDIM, int MESHDIM>
  class MeshTestToolkit
  {

  public:
    double _precision;
    INTERP_KERNEL::IntersectionType _intersectionType;//Used only in the case MESHDIM==2 (planar intersections)

    MeshTestToolkit():_precision(1.e-6),_intersectionType(INTERP_KERNEL::Triangulation)  {}
  
    ~MeshTestToolkit() {}

    void intersectMeshes(const char* mesh1, const char* mesh2, const double correctVol, const double prec = 1.0e-5, bool doubleTest = true) const;

    // 1.0e-5 here is due to limited precision of "correct" volumes calculated in Salome
    void intersectMeshes(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, const double correctVol, const double prec = 1.0e-5, bool doubleTest = true) const;
  
    void dumpIntersectionMatrix(const IntersectionMatrix& m) const;

    double sumRow(const IntersectionMatrix& m, int i) const;

    double sumCol(const IntersectionMatrix& m, int i) const;

    void getVolumes(MEDCoupling::MEDCouplingUMesh& mesh, double* tab) const;

    bool testVolumes(const IntersectionMatrix& m,  MEDCoupling::MEDCouplingUMesh& sMesh,  MEDCoupling::MEDCouplingUMesh& tMesh) const;

    double sumVolume(const IntersectionMatrix& m) const;

    bool areCompatitable( const IntersectionMatrix& m1,  const IntersectionMatrix& m2) const;

    bool testTranspose(const IntersectionMatrix& m1, const IntersectionMatrix& m2) const;

    bool testDiagonal(const IntersectionMatrix& m) const;
  
    void calcIntersectionMatrix(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, IntersectionMatrix& m) const;
  
  };
}
#endif
