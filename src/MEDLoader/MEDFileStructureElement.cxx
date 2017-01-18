// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#include "MEDFileStructureElement.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileMeshLL.hxx"
#include "MEDFileSafeCaller.txx"

#include "InterpKernelAutoPtr.hxx"

using namespace MEDCoupling;

MEDFileStructureElement *MEDFileStructureElement::New(med_idt fid, int idSE)
{
  return new MEDFileStructureElement(fid,idSE);
}

MEDFileStructureElement::MEDFileStructureElement(med_idt fid, int idSE)
{
  INTERP_KERNEL::AutoPtr<char> modelName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE)),supportMeshName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  med_geometry_type sgeoType;
  med_entity_type entiyType;
  int nConsAttr(0),nVarAttr(0);
  {
    med_bool anyPfl;
    int nnode(0),ncell(0);
    MEDstructElementInfo(fid,idSE+1,modelName,&_id_type,&_dim,supportMeshName,&entiyType,&nnode,&ncell,&sgeoType,&nConsAttr,&anyPfl,&nVarAttr);
  }
  _model_name=MEDLoaderBase::buildStringFromFortran(modelName,MED_NAME_SIZE);
  _geo_type=MEDFileMesh::ConvertFromMEDFile(sgeoType);
}

MEDFileStructureElements *MEDFileStructureElements::New(med_idt fid)
{
  return new MEDFileStructureElements(fid);
}

MEDFileStructureElements *MEDFileStructureElements::New()
{
  return new MEDFileStructureElements;
}

std::vector<const BigMemoryObject *> MEDFileStructureElements::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDFileStructureElement> >::const_iterator it=_elems.begin();it!=_elems.end();it++)
    ret.push_back(*it);
  return ret;
}

std::size_t MEDFileStructureElements::getHeapMemorySizeWithoutChildren() const
{
  return _elems.capacity()*sizeof(MEDFileStructureElement);
}

void MEDFileStructureElements::writeLL(med_idt fid) const
{
}

MEDFileStructureElements::MEDFileStructureElements(med_idt fid)
{
  int nbSE(MEDnStructElement(fid));
  _elems.resize(nbSE);
  for(int i=0;i<nbSE;i++)
    _elems[i]=MEDFileStructureElement::New(fid,i);
}

MEDFileStructureElements::MEDFileStructureElements()
{
}

MEDFileStructureElements::~MEDFileStructureElements()
{
}

