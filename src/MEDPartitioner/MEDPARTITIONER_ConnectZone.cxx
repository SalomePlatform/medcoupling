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

#include "MEDPARTITIONER_ConnectZone.hxx"

#include "MEDCouplingSkyLineArray.hxx"

#include <map>

using namespace MEDCoupling;

MEDPARTITIONER::ConnectZone::ConnectZone():
  _name("")
  ,_description("")
  ,_local_domain_number(0)
  ,_distant_domain_number(0)
  ,_node_corresp(0)
  ,_face_corresp(0)
  ,_local_mesh(0)
  ,_distant_mesh(0)
{
}

MEDPARTITIONER::ConnectZone::~ConnectZone()
{
  for(std::map < std::pair <int, int>,MEDCouplingSkyLineArray * >::iterator iter=_entity_corresp.begin(); iter!=_entity_corresp.end();iter++)
    {
      iter->second->decrRef();
    }
}

MEDPARTITIONER::ConnectZone::ConnectZone(const ConnectZone & myConnectZone):
  _name(myConnectZone._name)
  ,_description(myConnectZone._description)
  ,_local_domain_number(myConnectZone._local_domain_number)
  ,_distant_domain_number(myConnectZone._distant_domain_number)
  ,_node_corresp(myConnectZone._node_corresp)
  ,_face_corresp(myConnectZone._face_corresp)
  ,_entity_corresp(myConnectZone._entity_corresp)
  ,_local_mesh(0)
  ,_distant_mesh(0)
{
}

std::string MEDPARTITIONER::ConnectZone::getName() const 
{
  return _name;
}

std::string MEDPARTITIONER::ConnectZone::getDescription() const     
{
  return _description;
}

int MEDPARTITIONER::ConnectZone::getDistantDomainNumber() const 
{
  return _distant_domain_number;
}

int MEDPARTITIONER::ConnectZone::getLocalDomainNumber() const 
{
  return _local_domain_number;
}

MEDCouplingUMesh *MEDPARTITIONER::ConnectZone::getLocalMesh() const
{
  return _local_mesh;
}

MEDCouplingUMesh *MEDPARTITIONER::ConnectZone::getDistantMesh() const
{
  return _distant_mesh;
}

bool MEDPARTITIONER::ConnectZone::isEntityCorrespPresent(int localEntity, int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDCouplingSkyLineArray*>::const_iterator map_iter;
  for(map_iter iter=_entity_corresp.begin(); iter != _entity_corresp.end(); iter++)
    {
      if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
        return true;
    }
  return false;
}

const int *MEDPARTITIONER::ConnectZone::getNodeCorrespIndex() const
{
  return _node_corresp->getIndex();
}

const int *MEDPARTITIONER::ConnectZone::getNodeCorrespValue() const
{
  return _node_corresp->getValues();
}

int MEDPARTITIONER::ConnectZone::getNodeNumber() const
{
  return _node_corresp->getNumberOf();
}

const MEDCouplingSkyLineArray * MEDPARTITIONER::ConnectZone::getNodeCorresp() const
{
  return (const MEDCouplingSkyLineArray *)_node_corresp;
}

const int *MEDPARTITIONER::ConnectZone::getFaceCorrespIndex() const
{
  return _face_corresp->getIndex();
}

const int *MEDPARTITIONER::ConnectZone::getFaceCorrespValue() const
{
  return _face_corresp->getValues();
}

int MEDPARTITIONER::ConnectZone::getFaceNumber() const
{
  return _face_corresp->getNumberOf();
}

const MEDCouplingSkyLineArray * MEDPARTITIONER::ConnectZone::getFaceCorresp() const
{
  return _face_corresp;
}

const int *MEDPARTITIONER::ConnectZone::getEntityCorrespIndex(int localEntity,
                                                              int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDCouplingSkyLineArray*>::const_iterator map_iter;

  for(map_iter iter=_entity_corresp.begin();iter!=_entity_corresp.end();iter++)
  {
    if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
      return iter->second->getIndex();
  }
  return 0;
}

const int *MEDPARTITIONER::ConnectZone::getEntityCorrespValue(int localEntity,
                                                              int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDCouplingSkyLineArray*>::const_iterator map_iter;

  for (map_iter iter=_entity_corresp.begin();iter!=_entity_corresp.end();iter++)
  {
    if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
      return iter->second->getValues();
  }
  return 0;
}

int MEDPARTITIONER::ConnectZone::getEntityCorrespNumber(int localEntity,
                                                        int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDCouplingSkyLineArray*>::const_iterator map_iter;

  for(map_iter iter=_entity_corresp.begin();iter!=_entity_corresp.end();iter++)
  {
    if((iter->first).first==localEntity && (iter->first).second==distantEntity)
      return iter->second->getNumberOf();
  }
  return 0;
}

int MEDPARTITIONER::ConnectZone::getEntityCorrespLength(int localEntity,
                                                        int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDCouplingSkyLineArray*>::const_iterator map_iter;

  for (map_iter iter=_entity_corresp.begin(); iter != _entity_corresp.end(); iter++)
  {
    if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
      return iter->second->getLength();
  }
  return 0;
}

const MEDCouplingSkyLineArray *
MEDPARTITIONER::ConnectZone::getEntityCorresp(int localEntity, int distantEntity) const
{
  typedef std::map<std::pair<int,int>, MEDCouplingSkyLineArray*>::const_iterator map_iter;

  for (map_iter iter=_entity_corresp.begin(); iter != _entity_corresp.end(); iter++)
  {
    if ((iter->first).first==localEntity && (iter->first).second==distantEntity)
      return iter->second;
  }
  return 0;
}

std::vector< std::pair< int,int > > MEDPARTITIONER::ConnectZone::getEntities() const
{
  std::vector< std::pair< int,int > > types;

  std::map<std::pair<int,int>, MEDCouplingSkyLineArray*>::const_iterator
    iter = _entity_corresp.begin();
  for ( ; iter != _entity_corresp.end(); iter++)
    {
      types.push_back( iter->first );
    }

  return types;
}

void MEDPARTITIONER::ConnectZone::setName(const std::string& name) 
{
  _name=name;
}

void MEDPARTITIONER::ConnectZone::setDescription(const std::string& description)
{
  _description=description;
}

void MEDPARTITIONER::ConnectZone::setDistantDomainNumber(int distantDomainNumber)
{
  _distant_domain_number=distantDomainNumber;
}

void MEDPARTITIONER::ConnectZone::setLocalDomainNumber(int localDomainNumber)
{
  _local_domain_number=localDomainNumber;
}

void MEDPARTITIONER::ConnectZone::setLocalMesh(MEDCouplingUMesh * localMesh)
{
  _local_mesh=localMesh;
}

void MEDPARTITIONER::ConnectZone::setDistantMesh(MEDCouplingUMesh * distantMesh)
{
  _distant_mesh=distantMesh;
}

/*! transforms an int array containing 
 * the node-node connections
 * to a MEDCouplingSkyLineArray
 */
void MEDPARTITIONER::ConnectZone::setNodeCorresp(const int * nodeCorresp, int nbnode)
{
  MCAuto<DataArrayInt> indexArr( DataArrayInt::New() );
  MCAuto<DataArrayInt> valueArr( DataArrayInt::New() );
  indexArr->alloc( nbnode+1 );
  valueArr->alloc( 2*nbnode );
  int * index = indexArr->getPointer();
  int * value = valueArr->getPointer();
  for (int i=0; i<nbnode; i++)
    {
      index[i]=2*i;
      value[2*i  ]=nodeCorresp[2*i];
      value[2*i+1]=nodeCorresp[2*i+1];
    }
  index[nbnode]=2*nbnode;
  setNodeCorresp( MEDCouplingSkyLineArray::New( indexArr, valueArr ));
}

void MEDPARTITIONER::ConnectZone::setNodeCorresp(MEDCouplingSkyLineArray* array)
{
  MCAuto<MEDCouplingSkyLineArray> arr(array);
  _node_corresp = arr;
}

/*! transforms an int array containing 
 * the face-face connections
 * to a MEDCouplingSkyLineArray
 */
void MEDPARTITIONER::ConnectZone::setFaceCorresp(const int * faceCorresp, int nbface)
{
  MCAuto<DataArrayInt> indexArr( DataArrayInt::New() );
  MCAuto<DataArrayInt> valueArr( DataArrayInt::New() );
  indexArr->alloc( nbface+1 );
  valueArr->alloc( 2*nbface );
  int * index = indexArr->getPointer();
  int * value = valueArr->getPointer();
  for (int i=0; i<nbface; i++)
    {
      index[i]=2*i;
      value[2*i]=faceCorresp[2*i];
      value[2*i+1]=faceCorresp[2*i+1];
    }
  index[nbface]=2*nbface;
  setFaceCorresp( MEDCouplingSkyLineArray::New( indexArr, valueArr ));
}

void MEDPARTITIONER::ConnectZone::setFaceCorresp(MEDCouplingSkyLineArray* array)
{
  MCAuto<MEDCouplingSkyLineArray> arr (array);
  _face_corresp = arr;
}

/*! transforms an int array containing 
 * the entity-entity connections
 * to a MEDCouplingSkyLineArray
 * 
 * the resulting MEDCouplingSkyLineArray is put in the map
 */
void MEDPARTITIONER::ConnectZone::setEntityCorresp(int localEntity, int distantEntity,
                                                   const int *entityCorresp, int nbentity)
{ 
  MCAuto<DataArrayInt> indexArr( DataArrayInt::New() );
  MCAuto<DataArrayInt> valueArr( DataArrayInt::New() );
  indexArr->alloc( nbentity+1 );
  valueArr->alloc( 2*nbentity );
  int * index = indexArr->getPointer();
  int * value = valueArr->getPointer();
  for (int i=0; i<nbentity; i++)
    {
      index[i]=2*i;
      value[2*i  ]=entityCorresp[2*i];
      value[2*i+1]=entityCorresp[2*i+1];
    }
  index[nbentity]=2*nbentity;
  setEntityCorresp( localEntity, distantEntity, MEDCouplingSkyLineArray::New(indexArr,valueArr));
}

void MEDPARTITIONER::ConnectZone::setEntityCorresp(int localEntity, int distantEntity,
                                                   MEDCouplingSkyLineArray *array)
{
  MEDCouplingSkyLineArray * nullArray = 0;
  std::map < std::pair <int,int>, MEDCouplingSkyLineArray * >::iterator it;
  it = _entity_corresp.insert
    ( std::make_pair( std::make_pair(localEntity,distantEntity), nullArray )).first;
  if ( it->second != nullArray ) it->second->decrRef();
  it->second = array;
}
