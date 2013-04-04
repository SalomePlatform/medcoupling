// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingCMesh.hxx"

#include "MEDMeshMaker.hxx"

using namespace ParaMEDMEM;

ParaMEDMEM::MEDCouplingUMesh *MEDMeshMaker(int dim, int nbedge)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingCMesh> c=MEDCouplingCMesh::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::New();
  arr->alloc(nbedge+1,1); arr->iota(0.); arr->applyLin(1./double(nbedge),0.);
  switch(dim)
  {
  case 2:
    {
      c->setCoords(arr,arr);
      break;
    }
  case 3:
    {
      c->setCoords(arr,arr,arr);
      break;
    }
  default:
    throw INTERP_KERNEL::Exception("MEDMeshMaker : only dim 2 or 3 supported !");
  }
  return c->buildUnstructured();
}
