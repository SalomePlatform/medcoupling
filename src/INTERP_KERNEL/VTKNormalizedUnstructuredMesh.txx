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
#ifndef __VTKNORMALIZEDUNSTRUCTUREDMESH_TXX__
#define __VTKNORMALIZEDUNSTRUCTUREDMESH_TXX__

#include "VTKNormalizedUnstructuredMesh.hxx"

#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"

template<int MESHDIM>
VTKNormalizedUnstructuredMesh<MESHDIM>::VTKNormalizedUnstructuredMesh(vtkUnstructuredGrid *mesh):_mesh_in_vtk_mode(mesh),
                                                                                                 _tmp_index_array(0)
{
}

template<int MESHDIM>
VTKNormalizedUnstructuredMesh<MESHDIM>::~VTKNormalizedUnstructuredMesh()
{
  _mesh_in_vtk_mode->Delete();
  releaseTempArrays();
}

template<int MESHDIM>
void VTKNormalizedUnstructuredMesh<MESHDIM>::getBoundingBox(double *boundingBox) const
{
  double tmp[6];
  _mesh_in_vtk_mode->GetBounds(tmp);
  for(unsigned i=0;i<3;i++)
    {
      boundingBox[i]=tmp[2*i];
      boundingBox[3+i]=tmp[2*i+1];
    }
}

template<int MESHDIM>
NormalizedCellType VTKNormalizedUnstructuredMesh<MESHDIM>::getTypeOfElement(vtkIdType eltId) const
{
  int cellType=_mesh_in_vtk_mode->GetCellType(eltId);
  int convTab[30]={0,0,0,0,0,(int)NORM_TRI3,0,(int)NORM_POLYGON,0,(int)NORM_QUAD4,(int)NORM_TETRA4,0,(int)NORM_HEXA8
                   0,(int)NORM_PYRA5,0,0,0,(int)NORM_TRI6,(int)NORM_QUAD8,};
}

template<int MESHDIM>
unsigned long VTKNormalizedUnstructuredMesh<MESHDIM>::getNumberOfElements() const
{
  return _mesh_in_vtk_mode->GetNumberOfCells();
}

template<int MESHDIM>
unsigned long VTKNormalizedUnstructuredMesh<MESHDIM>::getNumberOfNodes() const
{
  return _mesh_in_vtk_mode->GetNumberOfPoints();
}

template<int MESHDIM>
const vtkIdType *VTKNormalizedUnstructuredMesh<MESHDIM>::getConnectivityPtr() const
{
  vtkIdType *ret=_mesh_in_vtk_mode->GetCells()->GetPointer();
  if(_tmp_index_array)
    return ret;
  else
    {
      putinMEDFormat();
      return ret;
    }
}

template<int MESHDIM>
const double *VTKNormalizedUnstructuredMesh<MESHDIM>::getCoordinatesPtr() const
{
  return (const double *)_mesh_in_vtk_mode->GetPoints()->GetVoidPointer(0);
}

template<int MESHDIM>
const vtkIdType *VTKNormalizedUnstructuredMesh<MESHDIM>::getConnectivityIndexPtr() const
{
  if(_tmp_index_array)
    return _tmp_index_array;
  else
    {
      putinMEDFormat();
      return _tmp_index_array;
    }
}

template<int MESHDIM>
void VTKNormalizedUnstructuredMesh<MESHDIM>::putinMEDFormat() const
{
  long nbOfElem=getNumberOfElements();
  _tmp_index_array=new vtkIdType[nbOfElem+1];
  _tmp_index_array[0]=0;
  vtkIdType *coarseConn=_mesh_in_vtk_mode->GetCells()->GetPointer();
  long ptInCC=0;
  vtkIdType *finalConn=coarseConn;
  for(long i=0;i<nbOfElem;i++)
    {
      vtkIdType cellLgth=coarseConn[ptInCC];
      for(vtkIdType j=0;j<cellLgth;j++)
        *finalConn++=coarseConn[ptInCC+j+1];
      _tmp_index_array[i+1]=_tmp_index_array[i]+cellLgth;
      ptInCC+=cellLgth+1;
    }
  int gh=0;
  gh++;
}

template<int MESHDIM>
void VTKNormalizedUnstructuredMesh<MESHDIM>::releaseTempArrays()
{
  delete [] _tmp_index_array;
  _tmp_index_array=0;
}

#endif
