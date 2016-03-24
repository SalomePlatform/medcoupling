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
// File      : MEDCouplingNormalizedCartesianMesh.txx
// Created   : Mon Aug 17 12:00:38 2009
// Author    : Edward AGAPOV (eap)
//

#include "MEDCouplingNormalizedCartesianMesh.hxx"
#include "MEDCouplingCMesh.hxx"

template<int SPACEDIM>
MEDCouplingNormalizedCartesianMesh<SPACEDIM>::MEDCouplingNormalizedCartesianMesh(const MEDCoupling::MEDCouplingCMesh *mesh):_mesh(mesh)
{
  if(_mesh)
    _mesh->incrRef();
}

template<int SPACEDIM>
MEDCouplingNormalizedCartesianMesh<SPACEDIM>::~MEDCouplingNormalizedCartesianMesh()
{
  if(_mesh)
    _mesh->decrRef();
}

template<int SPACEDIM>
unsigned long MEDCouplingNormalizedCartesianMesh<SPACEDIM>::getNumberOfElements() const
{
  return _mesh->getNumberOfCells();
}

template<int SPACEDIM>
unsigned long MEDCouplingNormalizedCartesianMesh<SPACEDIM>::nbCellsAlongAxis(int axis) const
{
  return _mesh->getCoordsAt(axis)->getNumberOfTuples() - 1;
}

template<int SPACEDIM>
const double * MEDCouplingNormalizedCartesianMesh<SPACEDIM>::getCoordsAlongAxis(int axis) const
{
  return _mesh->getCoordsAt(axis)->getConstPointer();
}
