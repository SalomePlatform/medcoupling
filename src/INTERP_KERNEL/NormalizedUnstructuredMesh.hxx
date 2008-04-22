#ifndef __NORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __NORMALIZEDUNSTRUCTUREDMESH_HXX__

namespace INTERP_KERNEL
{
  typedef enum
    {
      ALL_C_MODE       ,
      ALL_FORTRAN_MODE
    } NumberingPolicy;


  typedef enum
    {
      NORM_TRI3    =  3,
      NORM_QUAD4   =  4,
      NORM_POLYGON =  5,
      NORM_TRI6    =  6,
      NORM_QUAD8   =  8,
      //
      NORM_TETRA4  = 14,
      NORM_PYRA5   = 15,
      NORM_PENTA6  = 16,
      NORM_HEXA8   = 18
    } NormalizedCellType;

  class GenericMesh
  {};

  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  class NormalizedUnstructuredMesh : public GenericMesh
  {
  public:
    void getBoundingBox(double *boundingBox) const { asLeaf().getBoundingBox(boundingBox); }
    NormalizedCellType getTypeOfElement(ConnType eltId) const { return asLeaf().getTypeOfElement(eltId); }
    unsigned char getNumberOfNodesOfElement(ConnType eltId) const { return asLeaf().getNumberOfNodesOfElement(eltId); }
    unsigned long getNumberOfElements() const { return asLeaf().getNumberOfElements(); }
    const ConnType *getConnectivityPtr() const { return asLeaf().getConnectivityPtr(); }
    const double *getCoordinatesPtr() const { return asLeaf().getCoordinatesPtr(); }
    const ConnType *getConnectivityIndexPtr() const { return asLeaf().getConnectivityIndexPtr(); }
    void ReleaseTempArrays() { return asLeaf().ReleaseTempArrays(); }
  protected:
    MyMeshType& asLeaf() { return (MyMeshType&)(*this); }
    const MyMeshType& asLeaf() const { return (MyMeshType&)(*this); }
  };
}

#endif
