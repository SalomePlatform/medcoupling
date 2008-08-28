#ifndef _PLANAR_INTERSECTOR_HXX_
#define _PLANAR_INTERSECTOR_HXX_

#include "NormalizedUnstructuredMesh.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  class TranslationRotationMatrix;
  
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  class INTERPKERNEL_EXPORT PlanarIntersector
  {
  public:
    PlanarIntersector(double dimCaracteristic, double precision, double medianPlane, bool doRotate, int printLevel);
    virtual ~ PlanarIntersector() { }
  
    //Tool for cell intersection, result is always positive
    virtual double intersectCells(int icell_A, int icell_B, int nb_NodesA, int nb_NodesB) = 0;
  
    //Tool for cell filtering
    void createBoundingBoxes(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh, 
                             std::vector<double>& bbox);
    void adjustBoundingBoxes(std::vector<double>& bbox, double Surf3DAdjustmentEps );
    inline void getElemBB(double* bb, const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh, ConnType iP, ConnType nb_nodes);
  protected :
    static void projection(std::vector< double>& Coords_A, std::vector< double>& Coords_B, 
                           int nb_NodesA, int nb_NodesB, double epsilon, double median_plane, bool do_rotate);
    static void rotate3DTriangle( double* PP1, double*PP2, double*PP3,
                                  TranslationRotationMatrix& rotation_matrix);
  protected:
    double _dimCaracteristic;
    double _precision;
    double _medianPlane;
    bool _doRotate;
    int _printLevel;
  };
}

#endif
