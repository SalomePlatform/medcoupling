//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#include <list>
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Exception.hxx"
#include "PointLocatorAlgos.txx"
#include "PointLocator.hxx"

namespace INTERP_KERNEL {
  PointLocator::PointLocator(const MEDMEM::MESH& mesh)
  {
    int meshdim=mesh.getMeshDimension();
    int spacedim=mesh.getSpaceDimension();
    if (meshdim != spacedim) throw MEDMEM::MEDEXCEPTION("Locator is not implemented for meshdim != spacedim");
    switch (meshdim)
      {
      case 2:
        _medmesh = new MEDNormalizedUnstructuredMesh<2,2> (&mesh);
        _point_locator=new PointLocatorAlgos<MEDNormalizedUnstructuredMesh<2,2> >(*(static_cast<MEDNormalizedUnstructuredMesh<2,2>* >(_medmesh)));
        break;
      case 3:
        _medmesh = new MEDNormalizedUnstructuredMesh<3,3> (&mesh);
        _point_locator=new PointLocatorAlgos<MEDNormalizedUnstructuredMesh<3,3> >(*(static_cast<MEDNormalizedUnstructuredMesh<3,3>* >(_medmesh)));
        break;
      }
  }
  PointLocator::~PointLocator()
  {
    delete _medmesh;
    delete _point_locator;
  }

std::list<int> PointLocator::locate(const double* x)
  {
    return _point_locator->locates(x);
  }
}
