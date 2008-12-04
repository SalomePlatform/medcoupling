//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
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
