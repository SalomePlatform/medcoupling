#ifndef __MEDNORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __MEDNORMALIZEDUNSTRUCTUREDMESH_HXX__

#include "NormalizedUnstructuredMesh.hxx"

namespace MEDMEM
{
  class MESH;
}

template<int SPACEDIM,int MESHDIM>
class MEDNormalizedUnstructuredMesh : public INTERP_KERNEL::NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,int,INTERP_KERNEL::ALL_FORTRAN_MODE,MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM> >
{
public:
  MEDNormalizedUnstructuredMesh(const MEDMEM::MESH *mesh);
  void getBoundingBox(double *boundingBox) const;
  INTERP_KERNEL::NormalizedCellType getTypeOfElement(int eltId) const;
  unsigned char getNumberOfNodesOfElement(int eltId) const;
  unsigned long getNumberOfElements() const;
  const int *getConnectivityPtr() const;
  const double *getCoordinatesPtr() const;
  const int *getConnectivityIndexPtr() const;
  void ReleaseTempArrays();
protected:
  const MEDMEM::MESH *_meshInMedMode;
};

#endif
