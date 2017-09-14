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

#include "MEDFileJoint.hxx"
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileSafeCaller.txx"

#include "CellModel.hxx"
#include "InterpKernelAutoPtr.hxx"

extern med_geometry_type                 typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type                 typmai3[34];

using namespace MEDCoupling;

std::size_t MEDFileJointCorrespondence::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MCAuto<DataArrayInt>);
}

std::vector<const BigMemoryObject *> MEDFileJointCorrespondence::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

MEDFileJointCorrespondence::MEDFileJointCorrespondence():
  _is_nodal( true ),
  _loc_geo_type( INTERP_KERNEL::NORM_ERROR ),
  _rem_geo_type( INTERP_KERNEL::NORM_ERROR )
{
}

/*!
 * Constructor.
 *  \param [in] correspondence - correspondence.
 *  \param [in] is_nodal - is the correspondence of cells or nodes.
 *  \param [in] loc_geo_type - the local geometry type of correspondence.
 *  \param [in] rem_geo_type - the remote geometry type of correspondence.
 */
MEDFileJointCorrespondence::MEDFileJointCorrespondence(DataArrayInt* correspondence,
                                                       bool          isNodal,
                                                       INTERP_KERNEL::NormalizedCellType loc_geo_type,
                                                       INTERP_KERNEL::NormalizedCellType rem_geo_type):
  _is_nodal( isNodal ),
  _loc_geo_type( loc_geo_type ),
  _rem_geo_type( rem_geo_type )
{
  MEDFileJointCorrespondence::setCorrespondence( correspondence );
}

MEDFileJointCorrespondence* MEDFileJointCorrespondence::New(DataArrayInt* correspondence,
                                                            INTERP_KERNEL::NormalizedCellType loc_geo_type,
                                                            INTERP_KERNEL::NormalizedCellType rem_geo_type)
{
  return new MEDFileJointCorrespondence(correspondence, /*isNodal=*/false, loc_geo_type, rem_geo_type );
}

/*!
 * Returns a new MEDFileJointCorrespondence of nodes
 */
MEDFileJointCorrespondence *MEDFileJointCorrespondence::New(DataArrayInt* correspondence)
{
  return new MEDFileJointCorrespondence(correspondence);
}

/*!
 * Returns a new undefined MEDFileJointCorrespondence
 */

MEDFileJointCorrespondence *MEDFileJointCorrespondence::New()
{
  return new MEDFileJointCorrespondence();
}

/*!
 * Writes \a this joint into a MED file specified by its name.
 *  \param [in] fileName - the MED file name.
 *  \param [in] mode - the writing mode. For more on \a mode, see \ref AdvMEDLoaderBasics.
 *          - 2 - erase; an existing file is removed.
 *          - 1 - append; same data should not be present in an existing file.
 *          - 0 - overwrite; same data present in an existing file is overwritten.
 *  \param [in] order - order.
 *  \param [in] iteration - iteration.
 *  \throw If the mesh name is not set.
 *  \throw If \a mode == 1 and the same data is present in an existing file.
 */
void MEDFileJointCorrespondence::write(const std::string& fileName, int mode, const std::string& localMeshName, const std::string& jointName, int order, int iteration) const
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),medmod);

  std::ostringstream oss; oss << "MEDFileJointCorrespondence : error on attempt to write in file : \"" << fileName << "\"";
  MEDFileUtilities::CheckMEDCode(fid,fid,oss.str());

  if (( !_is_nodal ) &&
      ( _loc_geo_type == INTERP_KERNEL::NORM_ERROR ||
        _rem_geo_type == INTERP_KERNEL::NORM_ERROR ))
    {
      throw INTERP_KERNEL::Exception( "Geometric type not specified for a cell Joint" );
    }

  if ( (const DataArrayInt *)_correspondence )
    {
      writeLL(fid, localMeshName, jointName, order, iteration);
    }
  else
    {
      throw INTERP_KERNEL::Exception("MEDFileJointCorrespondence::write : correspondence array not defined");
    }
}

void MEDFileJointCorrespondence::writeLL(med_idt fid, const std::string& localMeshName, const std::string& jointName, int order, int iteration) const
{
  if ( _is_nodal )
    {
      MEDFILESAFECALLERWR0(MEDsubdomainCorrespondenceWr,(fid, localMeshName.c_str(), jointName.c_str(),
                                                         order, iteration,
                                                         MED_NODE, MED_NONE,
                                                         MED_NODE, MED_NONE,
                                                         _correspondence->getNbOfElems()/2,
                                                         _correspondence->getConstPointer()));
    }
  else
    {
      MEDFILESAFECALLERWR0(MEDsubdomainCorrespondenceWr,(fid, localMeshName.c_str(), jointName.c_str(),
                                                         order, iteration,
                                                         MED_CELL, typmai3[ _loc_geo_type ],
                                                         MED_CELL, typmai3[ _rem_geo_type ],
                                                         _correspondence->getNbOfElems()/2,
                                                         _correspondence->getConstPointer()));
    }
}

void MEDFileJointCorrespondence::setCorrespondence(DataArrayInt *corr)
{
  _correspondence=corr;
  if ( corr )
    corr->incrRef();
}

/*!
 * Checks if \a this and another mesh are equal.
 *  \param [in] other - the mesh to compare with.
 *  \return bool - \c true if the meshes are equal, \c false, else.
 */
bool MEDFileJointCorrespondence::isEqual(const MEDFileJointCorrespondence *other) const
{
  if(_is_nodal!=other->_is_nodal)
    return false;
  if(_loc_geo_type!=other->_loc_geo_type)
    return false;
  if(_rem_geo_type!=other->_rem_geo_type)
    return false;
  if(!_correspondence->isEqual(*other->_correspondence))
    return false;
  return true;
}

MEDFileJointCorrespondence *MEDFileJointCorrespondence::deepCopy() const
{
  MCAuto<MEDFileJointCorrespondence> ret=new MEDFileJointCorrespondence(*this);
  return ret.retn();
}

MEDFileJointCorrespondence *MEDFileJointCorrespondence::shallowCpy() const
{
  MCAuto<MEDFileJointCorrespondence> ret=new MEDFileJointCorrespondence(*this);
  return ret.retn();
}

/*!
 * Returns a string describing \a this mesh. This description includes the correspondence and
 * the number correspondence.
 *  \return std::string - the joint information string.
 */
std::string MEDFileJointCorrespondence::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(*************************************)\n(* JOINT_CORRESPOND INFORMATION: *)\n(*************************************)\n";
  oss << "- entity type of the correspondence : " << ( getIsNodal() ? "NODE" : "CELL" ) << "\n";
  oss << "- Local geometry type of the correspondence : " << INTERP_KERNEL::CellModel::GetCellModel( _loc_geo_type ).getRepr() << "\n";
  oss << "- Remote geometry type of the correspondence : " << INTERP_KERNEL::CellModel::GetCellModel( _rem_geo_type ).getRepr() << "\n";
  if ( (const DataArrayInt *)_correspondence )
    {
      oss << "- Number entity of the correspondence : " << getCorrespondence()->getNumberOfTuples() << "\n";

      const DataArrayInt* tmp=getCorrespondence();
      oss << "- Correspondence : <<";
      for(const int *it=tmp->begin();it!=tmp->end();it++)
        oss<< *it << " ";
    }
  else
    {
      oss << "- Number entity of the correspondence : 0\n";
    }
  oss << std::endl;
  return oss.str();
}


MEDFileJointOneStep::MEDFileJointOneStep():_order(-1),_iteration(-1)
{
}

std::size_t MEDFileJointOneStep::getHeapMemorySizeWithoutChildren() const
{
  return _correspondences.capacity()*sizeof(MCAuto<DataArrayInt>);
}

std::vector<const BigMemoryObject *> MEDFileJointOneStep::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

MEDFileJointOneStep *MEDFileJointOneStep::New(int dt, int it)
{
  MEDFileJointOneStep* j = new MEDFileJointOneStep();
  j->setOrder( dt );
  j->setIteration( it );
  return j;
}

/*!
 * Returns a new MEDFileJointOneStep.
 *  \param [in] fileName - the name of MED file to read.
 *  \param [in] mName - the name of the mesh to read.
 *  \param [in] jointName - the joint name.
 *  \param [in] num - the number of an iteration.
 *  \return MEDFileMesh * - a new instance of MEDFileJointOneStep.
 *  \throw If the file is not readable.
 *  \throw If there is no mesh with given attributes in the file.
 */
MEDFileJointOneStep *MEDFileJointOneStep::New(const std::string& fileName, const std::string& mName, const std::string& jointName, int num)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY);
  return new MEDFileJointOneStep(fid, mName, jointName, num);
}

MEDFileJointOneStep* MEDFileJointOneStep::New(med_idt fid, const std::string& mName, const std::string& jointName, int num)
{
  return new MEDFileJointOneStep( fid, mName, jointName, num);
}

MEDFileJointOneStep::MEDFileJointOneStep(med_idt fid, const std::string& mName, const std::string& jointName, int num)
{
  int order, iteration, ncorrespondence;
  MEDFILESAFECALLERRD0(MEDsubdomainComputingStepInfo,(fid, mName.c_str(), jointName.c_str(), num, &order, &iteration, &ncorrespondence));
  MEDFileJointOneStep::setOrder(order);
  MEDFileJointOneStep::setIteration(iteration);
  for ( int cur_it = 1; cur_it <= ncorrespondence; ++cur_it )
    {
      int num_entity;
      med_entity_type loc_ent_type, rem_ent_type;
      med_geometry_type loc_geo_type, rem_geo_type;
      MEDFILESAFECALLERRD0(MEDsubdomainCorrespondenceSizeInfo,(fid, mName.c_str(), jointName.c_str(), order, iteration, cur_it,
                                                               &loc_ent_type, &loc_geo_type, &rem_ent_type, &rem_geo_type, &num_entity));
      if ( num_entity > 0 )
        {
          MCAuto<DataArrayInt> correspondence=DataArrayInt::New();
          correspondence->alloc(num_entity*2, 1);
          MEDFILESAFECALLERRD0(MEDsubdomainCorrespondenceRd,(fid, mName.c_str(), jointName.c_str(), order, iteration, loc_ent_type,
                                                             loc_geo_type, rem_ent_type, rem_geo_type, correspondence->getPointer()));
          MEDFileJointCorrespondence *cor=MEDFileJointCorrespondence::New();
          cor->setIsNodal( loc_ent_type == MED_NODE );
          cor->setLocalGeometryType ( convertGeometryType( loc_geo_type ));
          cor->setRemoteGeometryType( convertGeometryType( rem_geo_type ));
          cor->setCorrespondence( correspondence );
          _correspondences.push_back(cor);
        }
    }
}

/*!
 * Writes \a this joint into a MED file specified by its name.
 *  \param [in] fileName - the MED file name.
 *  \param [in] mode - the writing mode. For more on \a mode, see \ref AdvMEDLoaderBasics.
 * - 2 - erase; an existing file is removed.
 * - 1 - append; same data should not be present in an existing file.
 * - 0 - overwrite; same data present in an existing file is overwritten.
 *  \throw If the mesh name is not set.
 *  \throw If \a mode == 1 and the same data is present in an existing file.
 */
void MEDFileJointOneStep::write(const std::string& fileName, int mode, const std::string& localMeshName, const std::string& jointName) const
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),medmod);
  std::ostringstream oss; oss << "MEDFileJointOneStep : error on attempt to write in file : \"" << fileName << "\"";
  MEDFileUtilities::CheckMEDCode(fid,fid,oss.str());
  if ( _correspondences.empty() )
    throw INTERP_KERNEL::Exception("MEDFileJointOneStep::write : no correspondences defined !");
  writeLL(fid, localMeshName, jointName);
}

void MEDFileJointOneStep::writeLL(med_idt fid, const std::string& localMeshName, const std::string& jointName) const
{
  for(std::vector< MCAuto<MEDFileJointCorrespondence> >::const_iterator it=_correspondences.begin();it!=_correspondences.end();it++)
    {
      (*it)->writeLL(fid, localMeshName, jointName, getOrder(), getIteration());
    }
}

void MEDFileJointOneStep::pushCorrespondence(MEDFileJointCorrespondence* correspondence)
{
  if(!correspondence)
    throw INTERP_KERNEL::Exception("MEDFileJointCorrespondence::pushCorrespondence : invalid input pointer ! should be different from 0 !");
  _correspondences.push_back(correspondence);
  correspondence->incrRef();
}

int MEDFileJointOneStep::getNumberOfCorrespondences() const
{
  return _correspondences.size();
}

/** Return a borrowed reference (caller is not responsible) */
MEDFileJointCorrespondence *MEDFileJointOneStep::getCorrespondenceAtPos(int i) const
{
  if(i<0 || i>=(int)_correspondences.size())
    {
      std::ostringstream oss; oss << "MEDFileJointOneStep::getCorrespondenceAtPos : invalid correspondence id given in parameter ! Should be in [0;" << _correspondences.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileJointCorrespondence* ret = _correspondences[i];
  return const_cast<MEDFileJointCorrespondence *>( ret );
}

/*!
 * Checks if \a this and another Joint are equal.
 *  \param [in] other - the Joint to compare with.
 *  \return bool - \c true if the Joints are equal, \c false, else.
 */
bool MEDFileJointOneStep::isEqual(const MEDFileJointOneStep *other) const
{
  if(_order!=other->_order)
    return false;
  if(_iteration!=other->_iteration)
    return false;
  if ( getNumberOfCorrespondences() != other->getNumberOfCorrespondences() )
    return false;

  std::vector<bool> found( getNumberOfCorrespondences(), false );
  for(int i=0; i<getNumberOfCorrespondences(); i++)
    {
      int j;
      for(j=0; j<getNumberOfCorrespondences(); j++)
        {
          if ( !found[ j ] &&
               getCorrespondenceAtPos(i)->isEqual(other->getCorrespondenceAtPos(j)))
            {
              found[ j ] = true;
              break;
            }
        }
      if ( j == getNumberOfCorrespondences() )
        return false;
    }
  return true;
}

MEDFileJointOneStep *MEDFileJointOneStep::deepCopy() const
{
  std::vector< MCAuto<MEDFileJointCorrespondence> > correspondences(_correspondences.size());
  std::size_t i=0;
  for(std::vector< MCAuto<MEDFileJointCorrespondence> >::const_iterator it=_correspondences.begin();it!=_correspondences.end();it++,i++)
    if((const MEDFileJointCorrespondence *)*it)
      correspondences[i]=(*it)->deepCopy();
  MCAuto<MEDFileJointOneStep> ret= new MEDFileJointOneStep;
  ret->_correspondences=correspondences;
  return ret.retn();
}

MEDFileJointOneStep *MEDFileJointOneStep::shallowCpy() const
{
  MCAuto<MEDFileJointOneStep> ret=new MEDFileJointOneStep(*this);
  return ret.retn();
}

/*!
 * Returns a string describing \a this Joint. This description includes the correspondence and
 * the number of correspondences.
 *  \return std::string - the joint information string.
 */
std::string MEDFileJointOneStep::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(*************************************)\n(* JOINT_ONE_STEP INFORMATION: *)\n(*************************************)\n";
  oss << "- Number of the correspondences : <<" << _correspondences.size() << ">>\n";
  for(std::vector< MCAuto<MEDFileJointCorrespondence> >::const_iterator it=_correspondences.begin();it!=_correspondences.end();it++)
    {
      oss << (*it)->simpleRepr();
    }
  return oss.str();
}

INTERP_KERNEL::NormalizedCellType MEDFileJointOneStep::convertGeometryType(med_geometry_type geotype)
{
  INTERP_KERNEL::NormalizedCellType result=INTERP_KERNEL::NORM_ERROR;
  for(int i=0; i<MED_N_CELL_FIXED_GEO; i++)
    {
      if (typmai[i]==geotype)
        {
          result=typmai2[i];
          break;
        }
    }
  return result;
}

std::size_t MEDFileJoint::getHeapMemorySizeWithoutChildren() const
{
  return _joint.capacity()*sizeof(MCAuto<MEDFileJointOneStep>);
}

std::vector<const BigMemoryObject *> MEDFileJoint::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

MEDFileJoint *MEDFileJoint::New()
{
  return new MEDFileJoint();
}

MEDFileJoint *MEDFileJoint::New(const std::string& fileName, const std::string& mName, int curJoint)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY);
  return new MEDFileJoint(fid,mName,curJoint);
}

MEDFileJoint *MEDFileJoint::New(med_idt fid, const std::string& mName, int curJoint)
{
  return new MEDFileJoint(fid,mName,curJoint);
}

MEDFileJoint *MEDFileJoint::New(const std::string& jointName, const std::string& locMeshName, const std::string& remoteMeshName, int remoteMeshNum)
{
  MEDFileJoint* j = new MEDFileJoint();
  j->setJointName( jointName );
  j->setLocalMeshName( locMeshName );
  j->setRemoteMeshName( remoteMeshName );
  j->setDomainNumber( remoteMeshNum );
  return j;
}

MEDFileJoint::MEDFileJoint()
{
}

/*!
 * Returns a new MEDFileJoint holding the mesh data that has been read from a given MED
 * file. The Joint to load is specified by mesh name and a Joint iteration.
 *  \param [in] fileName - the name of MED file to read.
 *  \param [in] mName - the name of the mesh to read.
 *  \param [in] curJoint - the iteration number of current joint.
 *  \return MEDFileJoint * - a new instance of MEDFileJoint.
 *  \throw If the file is not readable.
 *  \throw If there is no mesh with given attributes in the file.
 */
MEDFileJoint::MEDFileJoint(med_idt fid, const std::string& mName, int curJoint)
{
  INTERP_KERNEL::AutoPtr<char> joint_name=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> desc_name=MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE);
  INTERP_KERNEL::AutoPtr<char> rem_mesh_name=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  int domain_number=0, nstep=0, nocstpncorrespondence=0;
  MEDFILESAFECALLERRD0(MEDsubdomainJointInfo,(fid,mName.c_str(), curJoint, joint_name, desc_name, &domain_number,rem_mesh_name,
                                              &nstep, &nocstpncorrespondence));
  setLocalMeshName(mName);
  setRemoteMeshName(MEDLoaderBase::buildStringFromFortran(rem_mesh_name,MED_NAME_SIZE));
  setDescription(MEDLoaderBase::buildStringFromFortran(desc_name,MED_COMMENT_SIZE));
  setJointName(MEDLoaderBase::buildStringFromFortran(joint_name,MED_NAME_SIZE));
  setDomainNumber(domain_number);
  for(int cur_step=1; cur_step <= nstep; ++cur_step)
    {
      MEDFileJointOneStep *cor=MEDFileJointOneStep::New(fid, mName.c_str(), getJointName(), cur_step);
      _joint.push_back(cor);
    }
}

void MEDFileJoint::writeLL(med_idt fid) const
{
  // if ( _loc_mesh_name.empty() )
  //   throw INTERP_KERNEL::Exception("MEDFileJoint::write : name of a local mesh not defined!");
  MEDFILESAFECALLERWR0(MEDsubdomainJointCr,(fid,getLocalMeshName().c_str(),getJointName().c_str(),getDescription().c_str(),getDomainNumber(),getRemoteMeshName().c_str()));
  for(std::vector< MCAuto<MEDFileJointOneStep> >::const_iterator it=_joint.begin();it!=_joint.end();it++)
    (*it)->writeLL(fid, getLocalMeshName(),getJointName());
}

void MEDFileJoint::pushStep(MEDFileJointOneStep* step)
{
  if(!step)
    throw INTERP_KERNEL::Exception("MEDFileJoint::pushStep : invalid input pointer ! should be different from 0 !");
  _joint.push_back(step);
  step->incrRef();
}

int MEDFileJoint::getNumberOfSteps() const
{
  return _joint.size();
}

/** Return a borrowed reference (caller is not responsible) */
MEDFileJointOneStep *MEDFileJoint::getStepAtPos(int i) const
{
  if(i<0 || i>=(int)_joint.size())
    {
      std::ostringstream oss; oss << "MEDFileJoint::getStepAtPos : invalid step id given in parameter ! Should be in [0;" << _joint.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileJointOneStep* ret = _joint[i];
  return const_cast<MEDFileJointOneStep *>( ret );
}

/*!
 * Checks if \a this and another Joint are equal.
 *  \param [in] other - the Joint to compare with.
 *  \return bool - \c true if the Joints are equal, \c false, else.
 */
bool MEDFileJoint::isEqual(const MEDFileJoint *other) const
{
  if(_loc_mesh_name!=other->_loc_mesh_name)
    return false;
  if(_joint_name!=other->_joint_name)
    return false;
  if(_desc_name!=other->_desc_name)
      return false;
  if(_rem_mesh_name!=other->_rem_mesh_name)
    return false;
  if(_domain_number!=other->_domain_number)
    return false;
  int nbTS(getNumberOfSteps());
  if(nbTS!=other->getNumberOfSteps())
    return false;
  std::vector<bool> found(nbTS,false);
  for(int i=0;i<nbTS;i++)
    {
      int j;
      for(j=0;j<nbTS;j++)
        {
          if(!found[j] && getStepAtPos(i)->isEqual(other->getStepAtPos(j)))
            {
              found[j]=true;
              break;
            }
        }
      if(j==nbTS)
        return false;
    }
  return true;
}

MEDFileJoint *MEDFileJoint::deepCopy() const
{
  std::vector< MCAuto<MEDFileJointOneStep> > joint(_joint.size());
  std::size_t i=0;
  for(std::vector< MCAuto<MEDFileJointOneStep> >::const_iterator it=_joint.begin();it!=_joint.end();it++,i++)
    if((const MEDFileJointOneStep *)*it)
      joint[i]=(*it)->deepCopy();
  MCAuto<MEDFileJoint> ret=MEDFileJoint::New();
  ret->_joint=joint;
  return ret.retn();
}

MEDFileJoint *MEDFileJoint::shallowCpy() const
{
  MCAuto<MEDFileJoint> ret=new MEDFileJoint(*this);
  return ret.retn();
}

bool MEDFileJoint::changeJointNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  for(std::vector< std::pair<std::string,std::string> >::const_iterator it=modifTab.begin();it!=modifTab.end();it++)
    {
      if((*it).first==_joint_name)
        {
          _joint_name=(*it).second;
          return true;
        }
    }
  return false;
}

/*!
 * Returns a string describing \a this mesh. This description includes the correspondence and
 * the number correspondence.
 *  \return std::string - the joint information string.
 */
std::string MEDFileJoint::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(*************************************)\n(* JOINT INFORMATION: *)\n(*************************************)\n";
  oss << "- Local Mesh name : <<" << getLocalMeshName() << ">>\n";
  oss << "- Remote Mesh name : <<" << getRemoteMeshName() << ">>\n";
  oss << "- Description : <<" << getDescription() << ">>\n";
  oss << "- Joint name : <<" << getJointName() << ">>\n";
  oss << "- Domain number : " << getDomainNumber() << "\n";
  for(std::vector< MCAuto<MEDFileJointOneStep> >::const_iterator it=_joint.begin();it!=_joint.end();it++)
    {
      oss << (*it)->simpleRepr();
    }
  return oss.str();
}

MEDFileJoints *MEDFileJoints::New()
{
  return new MEDFileJoints;
}

MEDFileJoints *MEDFileJoints::New(const std::string& fileName, const std::string& meshName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY);
  return new MEDFileJoints( fid, meshName );
}

MEDFileJoints *MEDFileJoints::New(med_idt fid, const std::string& meshName)
{
  return new MEDFileJoints( fid, meshName );
}

void MEDFileJoints::writeLL(med_idt fid) const
{
  for(std::vector< MCAuto<MEDFileJoint> >::const_iterator it=_joints.begin();it!=_joints.end();it++)
    (*it)->writeLL(fid);
}

std::string MEDFileJoints::getMeshName() const
{
  for ( size_t i = 0; i <= _joints.size(); ++i )
    if ( (const MEDFileJoint*) _joints[i] )
      return _joints[i]->getLocalMeshName();

  return "";
}

int MEDFileJoints::getNumberOfJoints() const
{
  return _joints.size();
}

/** Return a borrowed reference (caller is not responsible) */
MEDFileJoint *MEDFileJoints::getJointAtPos(int i) const
{
  if(i<0 || i>=(int)_joints.size())
    {
      std::ostringstream oss; oss << "MEDFileJoints::getJointAtPos : invalid joint id given in parameter ! Should be in [0;" << _joints.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileJoint* ret = _joints[i];
  return const_cast<MEDFileJoint *>( ret );
}


/** Return a borrowed reference (caller is not responsible) */
MEDFileJoint *MEDFileJoints::getJointWithName(const std::string& jname) const
{
  std::vector<std::string> js=getJointsNames();
  std::vector<std::string>::iterator it=std::find(js.begin(),js.end(),jname);
  if(it==js.end())
    {
      std::ostringstream oss; oss << "MEDFileJoints::getJointWithName : Joint  \"" << jname << "\" does not exist in this ! Existing are : ";
      std::copy(js.begin(),js.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return getJointAtPos((int)std::distance(js.begin(),it));
}

std::vector<std::string> MEDFileJoints::getJointsNames() const
{
  std::vector<std::string> ret(_joints.size());
  int i=0;
  for(std::vector< MCAuto<MEDFileJoint> >::const_iterator it=_joints.begin();it!=_joints.end();it++,i++)
    {
      const MEDFileJoint *f=(*it);
      if(f)
        {
          ret[i]=f->getJointName();
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileJoints::getJointsNames : At rank #" << i << " joint is not defined !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret;
}

bool MEDFileJoints::changeJointNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  bool ret=false;
  for(std::vector< MCAuto<MEDFileJoint> >::iterator it=_joints.begin();it!=_joints.end();it++)
    {
      MEDFileJoint *cur(*it);
      if(cur)
        ret=cur->changeJointNames(modifTab) || ret;
    }
  return ret;
}

void MEDFileJoints::resize(int newSize)
{
  _joints.resize(newSize);
}

void MEDFileJoints::pushJoint(MEDFileJoint *joint)
{
  if(!joint)
    throw INTERP_KERNEL::Exception("MEDFileJoints::pushJoint() : invalid input pointer ! should be different from 0 !");
  if ( !_joints.empty() &&
       _joints[0]->getLocalMeshName() != joint->getLocalMeshName() )
    throw INTERP_KERNEL::Exception("MEDFileJoints::pushJoint() : different names of local meshes ! should be equal !");

  _joints.push_back(joint);
  joint->incrRef();
}

void MEDFileJoints::setJointAtPos(int i, MEDFileJoint *joint)
{
  if(!joint)
    throw INTERP_KERNEL::Exception("MEDFileJoints::setJointAtPos : invalid input pointer ! should be different from 0 !");
  if(i>=(int)_joints.size())
    _joints.resize(i+1);
  _joints[i]=joint;
  joint->incrRef();
}

void MEDFileJoints::destroyJointAtPos(int i)
{
  if(i<0 || i>=(int)_joints.size())
    {
      std::ostringstream oss; oss << "MEDFileJoints::destroyJointAtPos : Invalid given id in input (" << i << ") should be in [0," << _joints.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _joints.erase(_joints.begin()+i);
}

MEDFileJoints::MEDFileJoints()
{
}

MEDFileJoints::MEDFileJoints(med_idt fid, const std::string& meshName)
{
  int num_joint=MEDnSubdomainJoint(fid, meshName.c_str() );
  for(int i = 1; i <= num_joint; i++)
    _joints.push_back(MEDFileJoint::New(fid,meshName,i));
}

MEDFileJoints *MEDFileJoints::deepCopy() const
{
  std::vector< MCAuto<MEDFileJoint> > joints(_joints.size());
  std::size_t i=0;
  for(std::vector< MCAuto<MEDFileJoint> >::const_iterator it=_joints.begin();it!=_joints.end();it++,i++)
    if((const MEDFileJoint *)*it)
      joints[i]=(*it)->deepCopy();
  MCAuto<MEDFileJoints> ret=MEDFileJoints::New();
  ret->_joints=joints;
  return ret.retn();
}

std::size_t MEDFileJoints::getHeapMemorySizeWithoutChildren() const
{
  return _joints.capacity()*(sizeof(MCAuto<MEDFileJoint>));
}

std::vector<const BigMemoryObject *> MEDFileJoints::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDFileJoint> >::const_iterator it=_joints.begin();it!=_joints.end();it++)
    ret.push_back((const MEDFileJoint *)*it);
  return ret;
}

std::string MEDFileJoints::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(*****************)\n(* MEDFileJoints *)\n(*****************)\n\n";
  simpleReprWithoutHeader(oss);
  return oss.str();
}

void MEDFileJoints::simpleReprWithoutHeader(std::ostream& oss) const
{
  int nbOfJoints=getNumberOfJoints();
  oss << "There are " << nbOfJoints << " joints with the following names : \n";
  std::vector<std::string> jns=getJointsNames();
  for(int i=0;i<nbOfJoints;i++)
    oss << "  - #" << i << " \"" << jns[i] << "\"\n";
  for(std::vector< MCAuto<MEDFileJoint> >::const_iterator it=_joints.begin();it!=_joints.end();it++)
    {
      oss << (*it)->simpleRepr();
    }
}
