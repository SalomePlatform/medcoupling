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
// Author : Anthony Geay (CEA/DEN)

#include "MEDFileFieldOverView.hxx"
#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"

#include "CellModel.hxx"

using namespace ParaMEDMEM;

MEDFileMeshStruct *MEDFileMeshStruct::New(const MEDFileMesh *mesh)
{
  return new MEDFileMeshStruct(mesh);
}

std::size_t MEDFileMeshStruct::getHeapMemorySize() const
{
  return 0;
}

MEDFileMeshStruct::MEDFileMeshStruct(const MEDFileMesh *mesh):_mesh(mesh)
{
  std::vector<int> levs=mesh->getNonEmptyLevels();
  _name=mesh->getName();
  _nb_nodes=mesh->getNumberOfNodes();
  _geo_types_distrib.resize(levs.size());
  for(std::vector<int>::const_iterator lev=levs.begin();lev!=levs.end();lev++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> mLev=mesh->getGenMeshAtLevel(*lev);
      _geo_types_distrib[-(*lev)]=mLev->getDistributionOfTypes();
    }
}

int MEDFileMeshStruct::getNumberOfElemsOfGeoType(INTERP_KERNEL::NormalizedCellType t) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< std::vector<int> >::const_iterator it1=_geo_types_distrib.begin();it1!=_geo_types_distrib.end();it1++)
    {
      std::size_t sz=(*it1).size();
      if(sz%3!=0)
        throw INTERP_KERNEL::Exception("MEDFileMeshStruct::getNumberOfElemsOfGeoType : internal error in code !");
      std::size_t nbGeo=sz/3;
      for(std::size_t i=0;i<nbGeo;i++)
        if((*it1)[3*i]==(int)t)
          return (*it1)[3*i+1];
    }
  throw INTERP_KERNEL::Exception("The specified geometric type is not present in the mesh structure !");
}

//=

MEDFileField1TSStructItem2::MEDFileField1TSStructItem2(INTERP_KERNEL::NormalizedCellType a, const std::pair<int,int>& b, const std::string& c, const std::string& d):_geo_type(a),_start_end(b),_pfl(DataArrayInt::New()),_loc(d),_nb_of_entity(-1)
{
  _pfl->setName(c.c_str());
}

void MEDFileField1TSStructItem2::checkWithMeshStructForCells(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception)
{
  int nbOfEnt=mst->getNumberOfElemsOfGeoType(_geo_type);
  checkInRange(nbOfEnt,1,globs);
}

void MEDFileField1TSStructItem2::checkWithMeshStructForGaussNE(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception)
{
  int nbOfEnt=mst->getNumberOfElemsOfGeoType(_geo_type);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
  checkInRange(nbOfEnt,(int)cm.getNumberOfNodes(),globs);
}

void MEDFileField1TSStructItem2::checkWithMeshStructForGaussPT(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception)
{
  if(!globs)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkWithMeshStructForGaussPT : no globals specified !");
  if(_loc.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkWithMeshStructForGaussPT : no localization specified !");
  const MEDFileFieldLoc& loc=globs->getLocalization(_loc.c_str());
  int nbOfEnt=mst->getNumberOfElemsOfGeoType(_geo_type);
  checkInRange(nbOfEnt,loc.getNumberOfGaussPoints(),globs);
}

/*!
 * \param [in] nbOfEntity - number of entity that can be either cells or nodes. Not other possiblity.
 * \param [in] nip - number of integration points. 1 for ON_CELLS and NO_NODES
 */
void MEDFileField1TSStructItem2::checkInRange(int nbOfEntity, int nip, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception)
{
  _nb_of_entity=nbOfEntity;
  if(_pfl->getName().empty())
    {
      if(nbOfEntity!=(_start_end.second-_start_end.first)/nip)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkInRange : Mismatch between number of entities and size of node field !");
      return ;
    }
  else
    {
      if(!globs)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkInRange : Presence of a profile on field whereas no globals found in file !");
      const DataArrayInt *pfl=globs->getProfile(_pfl->getName().c_str());
      if(!pfl)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkInRange : Presence of a profile on field whereas no such profile found in file !");
      if(!pfl->checkAllIdsInRange(0,nbOfEntity))
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkInRange : The profile specified is invalid !");
    }
}

bool MEDFileField1TSStructItem2::operator==(const MEDFileField1TSStructItem2& other) const throw(INTERP_KERNEL::Exception)
{
  return _geo_type==other._geo_type && _start_end==other._start_end && _pfl->getName()==other._pfl->getName();
}

//=

MEDFileField1TSStructItem::MEDFileField1TSStructItem(TypeOfField a, const std::vector< MEDFileField1TSStructItem2 >& b):_type(a),_items(b)
{
}

void MEDFileField1TSStructItem::checkWithMeshStruct(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception)
{
  switch(_type)
    {
    case ON_NODES:
      {
        int nbOfEnt=mst->getNumberOfNodes();
        if(_items.size()!=1)
          throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::checkWithMeshStruct : for nodes field only one subdivision supported !");
        _items[0].checkInRange(nbOfEnt,1,globs);
        break ;
      }
    case ON_CELLS:
      {
        for(std::vector< MEDFileField1TSStructItem2 >::iterator it=_items.begin();it!=_items.end();it++)
          (*it).checkWithMeshStructForCells(mst,globs);
        break;
      }
    case ON_GAUSS_NE:
      {
        for(std::vector< MEDFileField1TSStructItem2 >::iterator it=_items.begin();it!=_items.end();it++)
          (*it).checkWithMeshStructForGaussNE(mst,globs);
        break;
      }
    case ON_GAUSS_PT:
      {
        for(std::vector< MEDFileField1TSStructItem2 >::iterator it=_items.begin();it!=_items.end();it++)
          (*it).checkWithMeshStructForGaussPT(mst,globs);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::checkWithMeshStruct : not managed field type !");
    }
}

bool MEDFileField1TSStructItem::operator==(const MEDFileField1TSStructItem& other) const throw(INTERP_KERNEL::Exception)
{
  if(_type!=other._type)
    return false;
  if(_items.size()!=other._items.size())
    return false;
  for(std::size_t i=0;i<_items.size();i++)
    if(!(_items[i]==other._items[i]))
      return false;
  return true;
}

bool MEDFileField1TSStructItem::isEntityCell() const
{
  if(_type==ON_NODES)
    return false;
  else
    return true;
}

class CmpGeo
{
public:
  CmpGeo(INTERP_KERNEL::NormalizedCellType geoTyp):_geo_type(geoTyp) { }
  bool operator()(const std::pair< INTERP_KERNEL::NormalizedCellType, std::vector<std::size_t> > & v) const { return _geo_type==v.first; }
private:
  INTERP_KERNEL::NormalizedCellType _geo_type;
};

MEDFileField1TSStructItem MEDFileField1TSStructItem::simplifyMeOnCellEntity(const MEDFileFieldGlobs *globs) const throw(INTERP_KERNEL::Exception)
{
  if(!isEntityCell())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::simplifyMeOnCellEntity : must be on ON_CELLS, ON_GAUSS_NE or ON_GAUSS_PT !");
  std::vector< std::pair< INTERP_KERNEL::NormalizedCellType, std::vector<std::size_t> > > m;
  std::size_t i=0;
  for(std::vector< MEDFileField1TSStructItem2 >::const_iterator it=_items.begin();it!=_items.end();it++,i++)
    {
      std::vector< std::pair< INTERP_KERNEL::NormalizedCellType, std::vector<std::size_t> > >::iterator it0(std::find_if(m.begin(),m.end(),CmpGeo((*it).getGeo())));
      if(it0==m.end())
        m.push_back(std::pair< INTERP_KERNEL::NormalizedCellType, std::vector<std::size_t> >((*it).getGeo(),std::vector<std::size_t>(1,i)));
      (*it0).second.push_back(i);
    }
  if(m.size()==_items.size())
    {
      MEDFileField1TSStructItem ret(*this);
      ret._type=ON_CELLS;
      return ret;
    }
  return *this;
}

//=

MEDFileField1TSStruct *MEDFileField1TSStruct::New(const MEDFileAnyTypeField1TS *ref) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TSStruct(ref);
}

MEDFileField1TSStruct::MEDFileField1TSStruct(const MEDFileAnyTypeField1TS *ref)
{
  _already_checked.push_back(BuildItemFrom(ref));
}

void MEDFileField1TSStruct::checkWithMeshStruct(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception)
{
  if(_already_checked.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::checkWithMeshStruct : not correctly initialized !");
  _already_checked.back().checkWithMeshStruct(mst,globs);
}

bool MEDFileField1TSStruct::isEqualConsideringThePast(const MEDFileAnyTypeField1TS *other) throw(INTERP_KERNEL::Exception)
{
  MEDFileField1TSStructItem b(BuildItemFrom(other));
  for(std::vector<MEDFileField1TSStructItem>::const_iterator it=_already_checked.begin();it!=_already_checked.end();it++)
    {
      if((*it)==b)
        return true;
    }
  return false;
}

bool MEDFileField1TSStruct::isSupportSameAs(const MEDFileAnyTypeField1TS *other) throw(INTERP_KERNEL::Exception)
{
  if(_already_checked.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::isSupportSameAs : no ref !");
  MEDFileField1TSStructItem b(BuildItemFrom(other));
  if(!_already_checked[0].isEntityCell() || !b.isEntityCell())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::isSupportSameAs : only available on cell entities !");
  return false;
}

std::size_t MEDFileField1TSStruct::getHeapMemorySize() const
{
  return 0;
}

MEDFileField1TSStructItem MEDFileField1TSStruct::BuildItemFrom(const MEDFileAnyTypeField1TS *ref)
{
  TypeOfField atype;
  std::vector< MEDFileField1TSStructItem2 > anItems;
  //
  std::vector< std::vector<std::string> > pfls,locs;
  std::vector< std::vector<TypeOfField> > typesF;
  std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes;
  std::vector< std::vector<std::pair<int,int> > > strtEnds=ref->getFieldSplitedByType(0,geoTypes,typesF,pfls,locs);
  std::size_t nbOfGeoTypes(geoTypes.size());
  if(nbOfGeoTypes==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct : not null by empty ref  !");
  bool isFirst=true;
  for(std::size_t i=0;i<nbOfGeoTypes;i++)
    {
      std::size_t sz=typesF.size();
      if(strtEnds[i].size()<1 || sz<1 || pfls[i].size()<1)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStruct : internal error #1 !");
      //
      if(isFirst)
        atype=typesF[i][0];
      isFirst=false;
      //
      for(std::size_t j=0;j<sz;j++)
        {
          if(atype!=typesF[i][j])
            anItems.push_back(MEDFileField1TSStructItem2(geoTypes[i],strtEnds[i][j],pfls[i][j],locs[i][j]));
          else
            throw INTERP_KERNEL::Exception("MEDFileField1TSStruct : can be applied only on single spatial discretization fields ! Call SplitPerDiscretization method !");
        }
    }
  return MEDFileField1TSStructItem(atype,anItems);
}

//=

MEDFileFastCellSupportComparator *MEDFileFastCellSupportComparator::New(const MEDFileMesh *m, const MEDFileAnyTypeFieldMultiTS *ref) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFastCellSupportComparator(m,ref);
}

MEDFileFastCellSupportComparator::MEDFileFastCellSupportComparator(const MEDFileMesh *m, const MEDFileAnyTypeFieldMultiTS *ref)
{
  _mesh_comp=MEDFileMeshStruct::New(m);
  int nbPts=ref->getNumberOfTS();
  _f1ts_cmps.resize(nbPts);
  for(int i=0;i<nbPts;i++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TS> elt=ref->getTimeStepAtPos(i);
      _f1ts_cmps[i]=MEDFileField1TSStruct::New(elt);
      _f1ts_cmps[i]->checkWithMeshStruct(_mesh_comp,elt->contentNotNull());
    }
}

std::size_t MEDFileFastCellSupportComparator::getHeapMemorySize() const
{
  /*std::size_t part1=sizeof(MEDFileFastCellSupportComparator)+_mesh_name.capacity()+_already_passed_code1.capacity()*sizeof(std::vector<int>)+_already_passed_code2.capacity()*sizeof(void*)+_m_geo_types_distrib.capacity()*sizeof(std::vector<int>);
  std::size_t part2=0;
  for(std::vector< std::vector<int> >::const_iterator it=_already_passed_code1.begin();it!=_already_passed_code1.end();it++)
    part2+=(*it).capacity()*(sizeof(int)+sizeof(const DataArrayInt *));
  for(std::vector< std::vector<int> >::const_iterator it2=_m_geo_types_distrib.begin();it2!=_m_geo_types_distrib.end();it2++)
    part2+=(*it2).capacity()*sizeof(int);
    return part1+part2;*/
  return 0;
}

bool MEDFileFastCellSupportComparator::isEqual(const MEDFileAnyTypeFieldMultiTS *other) throw(INTERP_KERNEL::Exception)
{
  int nbPts=other->getNumberOfTS();
  if(nbPts!=(int)_f1ts_cmps.size())
    {
      std::ostringstream oss; oss << "MEDFileFastCellSupportComparator::isEqualRegardingPast : unexpected nb of time steps in  input ! Should be " << _f1ts_cmps.size() << " it is in reality " << nbPts << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  for(int i=0;i<nbPts;i++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TS> elt=other->getTimeStepAtPos(i);
      if(!_f1ts_cmps[i]->isEqualConsideringThePast(elt))
        if(!_f1ts_cmps[i]->isSupportSameAs(elt))
          return false;
    }
  return true;
}
