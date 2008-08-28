#ifndef __VTKNORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __VTKNORMALIZEDUNSTRUCTUREDMESH_HXX__

#include "NormalizedUnstructuredMesh.hxx"

#include "vtkType.h"

class vtkUnstructuredGrid;

template<int MESHDIM>
class INTERPKERNEL_EXPORT VTKNormalizedUnstructuredMesh : public INTERP_KERNEL::NormalizedUnstructuredMesh<3,MESHDIM,vtkIdType,INTERP_KERNEL::ALL_C_MODE,VTKNormalizedUnstructuredMesh<MESHDIM> >
{
public:
  VTKNormalizedUnstructuredMesh(vtkUnstructuredGrid *mesh);
  ~VTKNormalizedUnstructuredMesh();
  void getBoundingBox(double *boundingBox) const;
  NormalizedCellType getTypeOfElement(vtkIdType eltId) const;
  unsigned long getNumberOfElements() const;
  const vtkIdType *getConnectivityPtr() const;
  const double *getCoordinatesPtr() const;
  const vtkIdType *getConnectivityIndexPtr() const;
  void ReleaseTempArrays();
protected:
  void putinMEDFormat() const;
protected:
  vtkUnstructuredGrid *_meshInVtkMode;
  mutable vtkIdType *_tmpIndexArray;
};

#endif
