// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#include "MEDPARTITIONER_ConnectZone.hxx"

#include <map>

using namespace MEDPARTITIONER;

ConnectZone::ConnectZone():
  _name("")
  ,_description("")
  ,_distantDomainNumber(0)
  ,_localDomainNumber(0)
  ,_nodeCorresp(0)
  ,_faceCorresp(0)
{
  _entityCorresp.clear();
}

ConnectZone::~ConnectZone(){
  if (_nodeCorresp !=0) delete _nodeCorresp;
  if (_faceCorresp !=0) delete _faceCorresp;
  for (std::map < std::pair <int, int>,MEDPARTITIONER::SkyLineArray * >::iterator 
         iter = _entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      delete iter->second;
    }
}

ConnectZone::ConnectZone(const ConnectZone & myConnectZone):
  _name(myConnectZone._name)
  ,_description(myConnectZone._description)
  ,_distantDomainNumber(myConnectZone._distantDomainNumber)
  ,_localDomainNumber(myConnectZone._localDomainNumber)
  ,_nodeCorresp(myConnectZone._nodeCorresp)
  ,_faceCorresp(myConnectZone._faceCorresp)
  ,_entityCorresp(myConnectZone._entityCorresp)
{
}

std::string ConnectZone::getName() const 
{
  return _name;
}

std::string ConnectZone::getDescription() const     
{
  return _description;
}

int ConnectZone::getDistantDomainNumber() const 
{
  return _distantDomainNumber;
}

int ConnectZone::getLocalDomainNumber() const 
{
  return _localDomainNumber;
}

ParaMEDMEM::MEDCouplingUMesh* ConnectZone::getLocalMesh() const 
{
  return _localMesh;
}

ParaMEDMEM::MEDCouplingUMesh * ConnectZone::getDistantMesh() const 
{
  return _distantMesh;
}

bool ConnectZone::isEntityCorrespPresent(int localEntity,
                                         int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::SkyLineArray*>::const_iterator map_iter;

  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return true;
    }
  return false;
}

const int * ConnectZone::getNodeCorrespIndex() const
{
  return _nodeCorresp->getIndex();
}

const int * ConnectZone::getNodeCorrespValue() const
{
  return _nodeCorresp->getValue();
}

int ConnectZone::getNodeNumber() const
{
  return _nodeCorresp->getNumberOf();
}

const int * ConnectZone::getFaceCorrespIndex() const
{
  return _faceCorresp->getIndex();
}

const int * ConnectZone::getFaceCorrespValue() const
{
  return _faceCorresp->getValue();
}

int ConnectZone::getFaceNumber() const
{
  return _faceCorresp->getNumberOf();
}

const int * ConnectZone::getEntityCorrespIndex(int localEntity,
                                               int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::SkyLineArray*>::const_iterator map_iter;

  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return iter->second->getIndex();
    }
  return 0;
}

const int * ConnectZone::getEntityCorrespValue(int localEntity,
                                               int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::SkyLineArray*>::const_iterator map_iter;

  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return iter->second->getValue();
    }
  return 0;
}

int ConnectZone::getEntityCorrespNumber(int localEntity,
                                        int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::SkyLineArray*>::const_iterator map_iter;

  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return iter->second->getNumberOf();
    }
  return 0;
}

int ConnectZone::getEntityCorrespLength(int localEntity,
                                        int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::SkyLineArray*>::const_iterator map_iter;
  
  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return iter->second->getLength();
    }
  return 0;
}

void ConnectZone::setName(std::string name) 
{
  _name=name;
}

void ConnectZone::setDescription(std::string description)
{
  _description=description;
}

void ConnectZone::setDistantDomainNumber(int distantDomainNumber)
{
  _distantDomainNumber=distantDomainNumber;
}

void ConnectZone::setLocalDomainNumber(int localDomainNumber)
{
  _localDomainNumber=localDomainNumber;
}

void ConnectZone::setLocalMesh(ParaMEDMEM::MEDCouplingUMesh * localMesh)
{
  _localMesh=localMesh;
}

void ConnectZone::setDistantMesh(ParaMEDMEM::MEDCouplingUMesh * distantMesh)
{
  _distantMesh=distantMesh;
}

/*! transforms an int array containing 
 * the node-node connections
 * to a SkyLineArray
 */
void ConnectZone::setNodeCorresp(int * nodeCorresp, int nbnode)
{
  std::vector<int> index(nbnode+1),value(2*nbnode);
  for (int i=0; i<nbnode; i++)
    {
      index[i]=2*i;
      value[2*i]=nodeCorresp[2*i];
      value[2*i+1]=nodeCorresp[2*i+1];
    }
  index[nbnode]=2*nbnode;
  _nodeCorresp = new MEDPARTITIONER::SkyLineArray(index,value);
}

void ConnectZone::setNodeCorresp(MEDPARTITIONER::SkyLineArray* array)
{
  _nodeCorresp = array;
}

/*! transforms an int array containing 
 * the face-face connections
 * to a SkyLineArray
 */
void ConnectZone::setFaceCorresp(int * faceCorresp, int nbface)
{
  std::vector<int> index(nbface+1),value(2*nbface);
  for (int i=0; i<nbface; i++)
    {
      index[i]=2*i;
      value[2*i]=faceCorresp[2*i];
      value[2*i+1]=faceCorresp[2*i+1];
    }
  index[nbface]=2*nbface;
  _faceCorresp = new MEDPARTITIONER::SkyLineArray(index,value);
}

void ConnectZone::setFaceCorresp(MEDPARTITIONER::SkyLineArray* array)
{
  _faceCorresp = array;
}

/*! transforms an int array containing 
 * the entity-entity connections
 * to a SkyLineArray
 * 
 * the resulting SkyLineArray is put in the map
 */
void ConnectZone::setEntityCorresp(int localEntity,
                                   int distantEntity,
                                   int * entityCorresp, int nbentity)
{ 
  std::vector<int> index(nbentity+1),value(2*nbentity);
  for (int i=0; i<nbentity; i++)
    {
      index[i]=2*i;
      value[2*i]=entityCorresp[2*i];
      value[2*i+1]=entityCorresp[2*i+1];
    }
  index[nbentity]=2*nbentity;
  _entityCorresp[std::make_pair(localEntity,distantEntity)] = new MEDPARTITIONER::SkyLineArray(index,value);
}

void ConnectZone::setEntityCorresp(int localEntity,
                                   int distantEntity,
                                   MEDPARTITIONER::SkyLineArray* array)
{
  _entityCorresp[std::make_pair(localEntity,distantEntity)]=array;
}







