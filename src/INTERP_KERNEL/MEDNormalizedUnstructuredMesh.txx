#ifndef __MEDNORMALIZEDUNSTRUCTUREDMESH_TXX__
#define __MEDNORMALIZEDUNSTRUCTUREDMESH_TXX__

#include "MEDNormalizedUnstructuredMesh.hxx"

#include "MEDMEM_Mesh.hxx"

template<int SPACEDIM,int MESHDIM>
MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::MEDNormalizedUnstructuredMesh(const MEDMEM::MESH *mesh):_meshInMedMode(mesh)
{
}

template<int SPACEDIM,int MESHDIM>
void MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getBoundingBox(double *boundingBox) const
{
  vector< vector<double> > ret=_meshInMedMode->getBoundingBox();
  //low left corner
  int i;
  for(i=0;i<SPACEDIM;i++)
    boundingBox[i]=ret[0][i];
  for(i=0;i<SPACEDIM;i++)
    boundingBox[i+SPACEDIM]=ret[1][i];
}

template<int SPACEDIM,int MESHDIM>
INTERP_KERNEL::NormalizedCellType MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getTypeOfElement(int eltId) const
{
  MED_EN::medGeometryElement type=_meshInMedMode->getElementTypeWithPoly(MED_EN::MED_CELL,eltId);
  if(type==MED_EN::MED_POLYGON)
    return INTERP_KERNEL::NORM_POLYGON;
  return (INTERP_KERNEL::NormalizedCellType)(((unsigned long)type/100-2)*10+((unsigned long)type%100));
}

template<int SPACEDIM,int MESHDIM>
unsigned char MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getNumberOfNodesOfElement(int eltId) const
{
  const int *ind=_meshInMedMode->getConnectivityIndex(MED_EN::MED_NODAL, MED_EN::MED_CELL);
  return (unsigned char) (ind[eltId]-ind[eltId-1]);
}

template<int SPACEDIM,int MESHDIM>
unsigned long MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getNumberOfElements() const
{
  return _meshInMedMode->getNumberOfElements(MED_EN::MED_CELL, MED_EN::MED_ALL_ELEMENTS);
}

template<int SPACEDIM,int MESHDIM>
const int *MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getConnectivityPtr() const
{
  return _meshInMedMode->getConnectivity(MED_EN::MED_FULL_INTERLACE,MED_EN::MED_NODAL,
                                         MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);
}

template<int SPACEDIM,int MESHDIM>
const double *MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getCoordinatesPtr() const
{
  return _meshInMedMode->getCoordinates(MED_EN::MED_FULL_INTERLACE);
}

template<int SPACEDIM,int MESHDIM>
const int *MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getConnectivityIndexPtr() const
{
  return _meshInMedMode->getConnectivityIndex(MED_EN::MED_NODAL, MED_EN::MED_CELL);
}

template<int SPACEDIM,int MESHDIM>
void MEDNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::ReleaseTempArrays()
{
}

#endif
