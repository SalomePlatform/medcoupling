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

MEDFileMeshSupports *MEDFileMeshSupports::New(const std::string& fileName)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid);
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
    {
      INTERP_KERNEL::AutoPtr<char> msn(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
      INTERP_KERNEL::AutoPtr<char> description(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
      med_axis_type axType;
      int nAxis(MEDsupportMeshnAxis(fid,i+1));
      INTERP_KERNEL::AutoPtr<char> axisName(new char[MED_SNAME_SIZE*nAxis+1]),axisUnit(new char[MED_SNAME_SIZE*nAxis+1]);
      int spaceDim(0),meshDim(0);
      MEDFILESAFECALLERRD0(MEDsupportMeshInfo,(fid,i+1,msn,&spaceDim,&meshDim,description,&axType,axisName,axisUnit));
      std::string name(MEDLoaderBase::buildStringFromFortran(msn,MED_NAME_SIZE));
      _supports[i]=MEDFileUMesh::New(fid,name);
    }
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
  return _supports.capacity()*sizeof(MCAuto<MEDFileUMesh>);
}

void MEDFileMeshSupports::writeLL(med_idt fid) const
{
  for(std::vector< MCAuto<MEDFileUMesh> >::const_iterator it=_supports.begin();it!=_supports.end();it++)
    if((*it).isNotNull())
      (*it)->writeLL(fid);
}

std::vector<std::string> MEDFileMeshSupports::getSupMeshNames() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto<MEDFileUMesh> >::const_iterator it=_supports.begin();it!=_supports.end();it++)
    if((*it).isNotNull())
      ret.push_back((*it)->getName());
  return ret;
}

const MEDFileUMesh *MEDFileMeshSupports::getSupMeshWithName(const std::string& name) const
{
  std::vector<std::string> mns;
  for(std::vector< MCAuto<MEDFileUMesh> >::const_iterator it=_supports.begin();it!=_supports.end();it++)
    {
      if((*it).isNotNull())
        {
          std::string na((*it)->getName());
          if(na==name)
            return *it;
          else
            mns.push_back(na);
        }
    }
  std::ostringstream oss;
  oss << "MEDFileMeshSupports::getSupMeshWithName : no such name \"" << name << "\". Possibilitities are :";
  std::copy(mns.begin(),mns.end(),std::ostream_iterator<std::string>(oss,","));
  oss << " !";
  throw INTERP_KERNEL::Exception(oss.str());
}

int MEDFileMeshSupports::getNumberOfNodesInConnOf(TypeOfField entity, const std::string& name) const
{
  const MEDFileUMesh *sup(getSupMeshWithName(name));
  switch(entity)
    {
    case ON_NODES:
      return sup->getNumberOfNodes();
    case ON_CELLS:
      {
        std::vector<INTERP_KERNEL::NormalizedCellType> gt(sup->getAllGeoTypes());
        if(gt.size()!=1)
          throw INTERP_KERNEL::Exception("MEDFileMeshSupports::getNumberOfNodesInConnOf : on cells only one geometric type allowed !");
        const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(gt[0]));
        return sup->getNumberOfCellsAtLevel(0)*cm.getNumberOfNodes();
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileMeshSupports::getNumberOfNodesInConnOf : not recognized entity type !");
    }
}
