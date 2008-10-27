#ifndef __CELLMODEL_INTERP_KERNEL_HXX__
#define __CELLMODEL_INTERP_KERNEL_HXX__

#include "NormalizedUnstructuredMesh.hxx"

#include <map>

namespace INTERP_KERNEL
{
  /*!
   * This class descibes all static elements (different from polygons and polyhedron) 3D, 2D and 1D.
   */
  class CellModel
  {
  public:
    static const unsigned MAX_NB_OF_SONS=6;
    static const unsigned MAX_NB_OF_NODES_PER_ELEM=30;
  private:
    CellModel(NormalizedCellType type);
    static void buildUniqueInstance();
  public:
    static const CellModel& getCellModel(NormalizedCellType type);
    //! sonId is in C format.
    unsigned getDimension() const { return _dim; }
    const unsigned *getNodesConstituentTheSon(unsigned sonId) const { return _sonsCon[sonId]; }
    unsigned getNumberOfNodes() const { return _nbOfPts; }
    unsigned getNumberOfSons() const { return _nbOfSons; }
    unsigned getNumberOfNodesConstituentTheSon(unsigned sonId) const { return _nbOfSonsCon[sonId]; }
    NormalizedCellType getSonType(unsigned sonId) const { return _sonsType[sonId]; }
    void fillSonCellNodalConnectivity(int sonId, const int *nodalConn, int *sonNodalConn) const;
  private:
    unsigned _dim;
    unsigned _nbOfPts;
    unsigned _nbOfSons;
    unsigned _sonsCon[MAX_NB_OF_SONS][MAX_NB_OF_NODES_PER_ELEM];
    unsigned _nbOfSonsCon[MAX_NB_OF_SONS];
    NormalizedCellType _sonsType[MAX_NB_OF_SONS];
    static std::map<NormalizedCellType,CellModel> _mapOfUniqueInstance;
  };
}

#endif
