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

#include "MEDFileMeshSupport.hxx"

#include "MEDLoaderBase.hxx"
#include "MEDFileMeshLL.hxx"
#include "MEDFileSafeCaller.txx"

#include "InterpKernelAutoPtr.hxx"

using namespace MEDCoupling;

MEDFileMeshSupport *MEDFileMeshSupport::New()
{
  return new MEDFileMeshSupport;
}

MEDFileMeshSupport *MEDFileMeshSupport::New(const std::string& fileName, int smid)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,smid);
}

MEDFileMeshSupport *MEDFileMeshSupport::New(med_idt fid, int smid)
{
  return new MEDFileMeshSupport(fid,smid);
}

std::vector<const BigMemoryObject *> MEDFileMeshSupport::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

std::size_t MEDFileMeshSupport::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

int MEDFileMeshSupport::getSpaceDim() const
{
  return _space_dim;
}

void MEDFileMeshSupport::setSpaceDim(int dim)
{
  _space_dim=dim;
}

int MEDFileMeshSupport::getMeshDim() const
{
  return _mesh_dim;
}

void MEDFileMeshSupport::setMeshDim(int dim)
{
  _mesh_dim=dim;
}

void MEDFileMeshSupport::writeLL(med_idt fid) const
{
}

MEDFileMeshSupport::MEDFileMeshSupport(med_idt fid, int smid)
{
  INTERP_KERNEL::AutoPtr<char> msn(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> description(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
  med_axis_type axType;
  int nAxis(MEDsupportMeshnAxis(fid,smid+1));
  INTERP_KERNEL::AutoPtr<char> axisName(new char[MED_SNAME_SIZE*nAxis+1]),axisUnit(new char[MED_SNAME_SIZE*nAxis+1]);
  MEDFILESAFECALLERRD0(MEDsupportMeshInfo,(fid,smid+1,msn,&_space_dim,&_mesh_dim,description,&axType,axisName,axisUnit));
  _name=MEDLoaderBase::buildStringFromFortran(msn,MED_NAME_SIZE);
  _description=MEDLoaderBase::buildStringFromFortran(description,MED_COMMENT_SIZE);
  _axis_type=MEDFileMeshL2::TraduceAxisType(axType);
  int nCoords(0),nmodels(0);
  {
    med_bool chgt=MED_FALSE,trsf=MED_FALSE;
    nCoords=MEDmeshnEntity(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_NODE,MED_NO_GEOTYPE,MED_COORDINATE,MED_NODAL,&chgt,&trsf);
    nmodels=MEDmeshnEntity(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_STRUCT_ELEMENT,MED_GEO_ALL,MED_CONNECTIVITY,MED_NODAL,&chgt,&trsf);
  }
  _coords=DataArrayDouble::New(); _coords->alloc(nCoords,nAxis);
  for(int i=0;i<nAxis;i++)
    {
      std::string info(DataArray::BuildInfoFromVarAndUnit(MEDLoaderBase::buildStringFromFortran(axisName+i*MED_SNAME_SIZE,MED_SNAME_SIZE),
                                                          MEDLoaderBase::buildStringFromFortran(axisUnit+i*MED_SNAME_SIZE,MED_SNAME_SIZE)));
      _coords->setInfoOnComponent(i,info);
    }
  MEDFILESAFECALLERRD0(MEDmeshNodeCoordinateRd,(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_FULL_INTERLACE,_coords->getPointer()));
  {
    med_bool chgt=MED_FALSE,trsf=MED_FALSE;
    if(MEDmeshnEntity(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_NODE,MED_NO_GEOTYPE,MED_FAMILY_NUMBER,MED_NODAL,&chgt,&trsf)>0)
      {
        _fam_coords=DataArrayInt::New();
        _fam_coords->alloc(nCoords,1);
        MEDFILESAFECALLERRD0(MEDmeshEntityFamilyNumberRd,(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_NODE,MED_NO_GEOTYPE,_fam_coords->getPointer()));
      }
    else
      _fam_coords=0;
    if(MEDmeshnEntity(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_NODE,MED_NO_GEOTYPE,MED_NUMBER,MED_NODAL,&chgt,&trsf)>0)
      {
        _num_coords=DataArrayInt::New();
        _num_coords->alloc(nCoords,1);
        MEDFILESAFECALLERRD0(MEDmeshEntityNumberRd,(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_NODE,MED_NO_GEOTYPE,_num_coords->getPointer()));
      }
    else
    _num_coords=0;
    if(MEDmeshnEntity(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_NODE,MED_NO_GEOTYPE,MED_NAME,MED_NODAL,&chgt,&trsf)>0)
      {
        _name_coords=DataArrayAsciiChar::New();
        _name_coords->alloc(nCoords+1,MED_SNAME_SIZE);//not a bug to avoid the memory corruption due to last \0 at the end
        MEDFILESAFECALLERRD0(MEDmeshEntityNameRd,(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_NODE,MED_NO_GEOTYPE,_name_coords->getPointer()));
        _name_coords->reAlloc(nCoords);//not a bug to avoid the memory corruption due to last \0 at the end
      }
    else
      _name_coords=0;
  }
    //med_bool withnodename;
  //MEDmeshNodeRd(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_FULL_INTERLACE,_coords->getPointer(),);
  std::cerr << nCoords << " ** " << nmodels << std::endl;
}

MEDFileMeshSupport::~MEDFileMeshSupport()
{
}

MEDFileMeshSupport::MEDFileMeshSupport():_space_dim(-1),_mesh_dim(-1),_axis_type(AX_CART)
{
}

MEDFileMeshSupports *MEDFileMeshSupports::New(med_idt fid)
{
  return new MEDFileMeshSupports(fid);
}

MEDFileMeshSupports *MEDFileMeshSupports::New()
{
  return new MEDFileMeshSupports;
}

MEDFileMeshSupports::MEDFileMeshSupports(med_idt fid)
{
  int nbSM(MEDnSupportMesh(fid));
  _supports.resize(nbSM);
  for(int i=0;i<nbSM;i++)
    _supports[i]=MEDFileMeshSupport::New(fid,i);
}

MEDFileMeshSupports::MEDFileMeshSupports()
{
}

MEDFileMeshSupports::~MEDFileMeshSupports()
{
}

std::vector<const BigMemoryObject *> MEDFileMeshSupports::getDirectChildrenWithNull() const
{
  std::size_t sz(_supports.size());
  std::vector<const BigMemoryObject *> ret(sz);
  for(std::size_t i=0;i<sz;i++)
    ret[i]=_supports[i];
  return ret;
}

std::size_t MEDFileMeshSupports::getHeapMemorySizeWithoutChildren() const
{
  return _supports.capacity()*sizeof(MCAuto<MEDFileMeshSupport>);
}

void MEDFileMeshSupports::writeLL(med_idt fid) const
{
  for(std::vector< MCAuto<MEDFileMeshSupport> >::const_iterator it=_supports.begin();it!=_supports.end();it++)
    if((*it).isNotNull())
      (*it)->writeLL(fid);
}
