#ifndef __TU_MESH_TEST_TOOLKIT_HXX__
#define __TU_MESH_TEST_TOOLKIT_HXX__

#include "../Interpolation3D.hxx"
#include "../Interpolation3D.txx"
#include "../InterpolationPlanar.hxx"

#include <vector>
#include <map>

#define ERR_TOL 1.0e-8

typedef std::vector<std::map<int,double> > IntersectionMatrix;

namespace INTERP_KERNEL
{
class Interpolation3D;
}


namespace MEDMEM {
class MESH;
};

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

    void getVolumes( MEDMEM::MESH& mesh,const double*& tab) const;

    bool testVolumes(const IntersectionMatrix& m,  MEDMEM::MESH& sMesh,  MEDMEM::MESH& tMesh) const;

    double sumVolume(const IntersectionMatrix& m) const;

    bool areCompatitable( const IntersectionMatrix& m1,  const IntersectionMatrix& m2) const;

    bool testTranspose(const IntersectionMatrix& m1, const IntersectionMatrix& m2) const;

    bool testDiagonal(const IntersectionMatrix& m) const;
  
    void calcIntersectionMatrix(const char* mesh1path, const char* mesh1, const char* mesh2path, const char* mesh2, IntersectionMatrix& m) const;
	
  };
}
#endif
