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
#ifndef __VTKNORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __VTKNORMALIZEDUNSTRUCTUREDMESH_HXX__

#include "NormalizedUnstructuredMesh.hxx"

#include "vtkType.h"

class vtkUnstructuredGrid;

template<int MESHDIM>
class INTERPKERNEL_EXPORT VTKNormalizedUnstructuredMesh : public INTERP_KERNEL::GenericMesh
{
public:
  static const int MY_SPACEDIM=3;
  static const int MY_MESHDIM=MESHDIM;
  typedef vtkIdType MyConnType;
  static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
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
