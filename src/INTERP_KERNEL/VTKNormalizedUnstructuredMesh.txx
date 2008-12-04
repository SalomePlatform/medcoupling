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
#ifndef __VTKNORMALIZEDUNSTRUCTUREDMESH_TXX__
#define __VTKNORMALIZEDUNSTRUCTUREDMESH_TXX__

#include "VTKNormalizedUnstructuredMesh.hxx"

#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"

template<int MESHDIM>
VTKNormalizedUnstructuredMesh<MESHDIM>::VTKNormalizedUnstructuredMesh(vtkUnstructuredGrid *mesh):_meshInVtkMode(mesh),
                                                                                                 _tmpIndexArray(0)
{
}

template<int MESHDIM>
VTKNormalizedUnstructuredMesh<MESHDIM>::~VTKNormalizedUnstructuredMesh()
{
  _meshInVtkMode->Delete();
  ReleaseTempArrays();
}

template<int MESHDIM>
void VTKNormalizedUnstructuredMesh<MESHDIM>::getBoundingBox(double *boundingBox) const
{
  double tmp[6];
  _meshInVtkMode->GetBounds(tmp);
  for(unsigned i=0;i<3;i++)
    {
      boundingBox[i]=tmp[2*i];
      boundingBox[3+i]=tmp[2*i+1];
    }
}

template<int MESHDIM>
NormalizedCellType VTKNormalizedUnstructuredMesh<MESHDIM>::getTypeOfElement(vtkIdType eltId) const
{
  int cellType=_meshInVtkMode->GetCellType(eltId);
  int convTab[30]={0,0,0,0,0,(int)NORM_TRI3,0,(int)NORM_POLYGON,0,(int)NORM_QUAD4,(int)NORM_TETRA4,0,(int)NORM_HEXA8
                   0,(int)NORM_PYRA5,0,0,0,(int)NORM_TRI6,(int)NORM_QUAD8,};
}

template<int MESHDIM>
unsigned long VTKNormalizedUnstructuredMesh<MESHDIM>::getNumberOfElements() const
{
  return _meshInVtkMode->GetNumberOfCells();
}

template<int MESHDIM>
const vtkIdType *VTKNormalizedUnstructuredMesh<MESHDIM>::getConnectivityPtr() const
{
  vtkIdType *ret=_meshInVtkMode->GetCells()->GetPointer();
  if(_tmpIndexArray)
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
  return (const double *)_meshInVtkMode->GetPoints()->GetVoidPointer(0);
}

template<int MESHDIM>
const vtkIdType *VTKNormalizedUnstructuredMesh<MESHDIM>::getConnectivityIndexPtr() const
{
  if(_tmpIndexArray)
    return _tmpIndexArray;
  else
    {
      putinMEDFormat();
      return _tmpIndexArray;
    }
}

template<int MESHDIM>
void VTKNormalizedUnstructuredMesh<MESHDIM>::putinMEDFormat() const
{
  long nbOfElem=getNumberOfElements();
  _tmpIndexArray=new vtkIdType[nbOfElem+1];
  _tmpIndexArray[0]=0;
  vtkIdType *coarseConn=_meshInVtkMode->GetCells()->GetPointer();
  long ptInCC=0;
  vtkIdType *finalConn=coarseConn;
  for(long i=0;i<nbOfElem;i++)
    {
      vtkIdType cellLgth=coarseConn[ptInCC];
      for(vtkIdType j=0;j<cellLgth;j++)
        *finalConn++=coarseConn[ptInCC+j+1];
      _tmpIndexArray[i+1]=_tmpIndexArray[i]+cellLgth;
      ptInCC+=cellLgth+1;
    }
  int gh=0;
  gh++;
}

template<int MESHDIM>
void VTKNormalizedUnstructuredMesh<MESHDIM>::ReleaseTempArrays()
{
  delete [] _tmpIndexArray;
  _tmpIndexArray=0;
}

#endif
