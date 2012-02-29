//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D, OPEN CASCADE
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

// few Med Memory include files
#include "MEDCouplingUMesh.hxx"
#include "MEDPARTITIONER_SkyLineArray.hxx"
#include "MEDPARTITIONER_ConnectZone.hxx"

// few STL include files
#include <map>

using namespace MEDPARTITIONER;

CONNECTZONE::CONNECTZONE():
  _name("")
  ,_description("")
  ,_distantDomainNumber(0)
  ,_localDomainNumber(0)
  ,_nodeCorresp(0)
  ,_faceCorresp(0)
{
  _entityCorresp.clear();
}

CONNECTZONE::~CONNECTZONE(){
  if (_nodeCorresp !=0) delete _nodeCorresp;
  if (_faceCorresp !=0) delete _faceCorresp;
  for (std::map < std::pair <int, int>,MEDPARTITIONER::MEDSKYLINEARRAY * >::iterator 
         iter = _entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      delete iter->second;
    }
}

CONNECTZONE::CONNECTZONE(const CONNECTZONE & myConnectZone):
  _name(myConnectZone._name)
  ,_description(myConnectZone._description)
  ,_distantDomainNumber(myConnectZone._distantDomainNumber)
  ,_localDomainNumber(myConnectZone._localDomainNumber)
  ,_nodeCorresp(myConnectZone._nodeCorresp)
  ,_faceCorresp(myConnectZone._faceCorresp)
  ,_entityCorresp(myConnectZone._entityCorresp)
{
}
std::string CONNECTZONE::getName() const 
{
  return _name;
}
std::string CONNECTZONE::getDescription() const 
{
  return _description;
}
int CONNECTZONE::getDistantDomainNumber() const 
{
  return _distantDomainNumber;
}
int CONNECTZONE::getLocalDomainNumber() const 
{
  return _localDomainNumber;
}

ParaMEDMEM::MEDCouplingUMesh* CONNECTZONE::getLocalMesh() const 
{
  return _localMesh;
}

ParaMEDMEM::MEDCouplingUMesh * CONNECTZONE::getDistantMesh() const 
{
  return _distantMesh;
}

bool CONNECTZONE::isEntityCorrespPresent(int localEntity,
                                         int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::MEDSKYLINEARRAY*>::const_iterator map_iter;
  
  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return true;
    }
  return false;
}                

const int * CONNECTZONE::getNodeCorrespIndex() const
{
  return _nodeCorresp->getIndex();
}

const int * CONNECTZONE::getNodeCorrespValue() const
{
  return _nodeCorresp->getValue();
}
int CONNECTZONE::getNodeNumber() const
{
  return _nodeCorresp->getNumberOf();
}
const int * CONNECTZONE::getFaceCorrespIndex() const
{
  return _faceCorresp->getIndex();
}

const int * CONNECTZONE::getFaceCorrespValue() const
{
  return _faceCorresp->getValue();
}
int CONNECTZONE::getFaceNumber() const
{
  return _faceCorresp->getNumberOf();
}
const int * CONNECTZONE::getEntityCorrespIndex(int localEntity,
                                               int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::MEDSKYLINEARRAY*>::const_iterator map_iter;
  
  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return iter->second->getIndex();
    }
  return 0;                       
}

const int * CONNECTZONE::getEntityCorrespValue(int localEntity,
                                               int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::MEDSKYLINEARRAY*>::const_iterator map_iter;
        
  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return iter->second->getValue();
    }
  return 0;                       
}

int CONNECTZONE::getEntityCorrespNumber(int localEntity,
                                        int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::MEDSKYLINEARRAY*>::const_iterator map_iter;
  
  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return iter->second->getNumberOf();
    }
  return 0;           
}


int CONNECTZONE::getEntityCorrespLength(int localEntity,
                                        int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDPARTITIONER::MEDSKYLINEARRAY*>::const_iterator map_iter;
  
  for (map_iter iter=_entityCorresp.begin(); iter != _entityCorresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return iter->second->getLength();
    }
  return 0;           
}

void CONNECTZONE::setName(std::string name) 
{
  _name=name;
}
void CONNECTZONE::setDescription(std::string description)
{
  _description=description;
}
void CONNECTZONE::setDistantDomainNumber(int distantDomainNumber)
{
  _distantDomainNumber=distantDomainNumber;
}
void CONNECTZONE::setLocalDomainNumber(int localDomainNumber)
{
  _localDomainNumber=localDomainNumber;
}
void CONNECTZONE::setLocalMesh(ParaMEDMEM::MEDCouplingUMesh * localMesh)
{
  _localMesh=localMesh;
}

void CONNECTZONE::setDistantMesh(ParaMEDMEM::MEDCouplingUMesh * distantMesh)
{
  _distantMesh=distantMesh;
}

/*! transforms an int array containing 
 * the node-node connections
 * to a MEDSKYLINEARRAY
 */
void CONNECTZONE::setNodeCorresp(int * nodeCorresp, int nbnode)
{
  std::vector<int> index(nbnode+1),value(2*nbnode);
  for (int i=0; i<nbnode; i++)
    {
      index[i]=2*i;
      value[2*i]=nodeCorresp[2*i];
      value[2*i+1]=nodeCorresp[2*i+1];
    }
  index[nbnode]=2*nbnode;
  _nodeCorresp = new MEDPARTITIONER::MEDSKYLINEARRAY(index,value);
}

void CONNECTZONE::setNodeCorresp(MEDPARTITIONER::MEDSKYLINEARRAY* array)
{
  _nodeCorresp = array;
}
/*! transforms an int array containing 
 * the face-face connections
 * to a MEDSKYLINEARRAY
 */
void CONNECTZONE::setFaceCorresp(int * faceCorresp, int nbface)
{
  std::vector<int> index(nbface+1),value(2*nbface);
  for (int i=0; i<nbface; i++)
    {
      index[i]=2*i;
      value[2*i]=faceCorresp[2*i];
      value[2*i+1]=faceCorresp[2*i+1];
    }
  index[nbface]=2*nbface;
  _faceCorresp = new MEDPARTITIONER::MEDSKYLINEARRAY(index,value);
}

void CONNECTZONE::setFaceCorresp(MEDPARTITIONER::MEDSKYLINEARRAY* array)
{
  _faceCorresp = array;
}

/*! transforms an int array containing 
 * the entity-entity connections
 * to a MEDSKYLINEARRAY
 * 
 * the resulting MEDSKYLINEARRAY is put in the map
 */
void CONNECTZONE::setEntityCorresp(int localEntity,
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
  _entityCorresp[std::make_pair(localEntity,distantEntity)] = new MEDPARTITIONER::MEDSKYLINEARRAY(index,value);
}


void CONNECTZONE::setEntityCorresp(int localEntity,
                                   int distantEntity,
                                   MEDPARTITIONER::MEDSKYLINEARRAY* array)
{
  _entityCorresp[std::make_pair(localEntity,distantEntity)]=array;
}
  






