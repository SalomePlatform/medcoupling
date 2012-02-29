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

// File: MEDMEM_ConnectZone.hxx  Author: Nadir Bouhamou (CEA Saclay)
//
# ifndef __MEDPARTITIONER_CONNECTZONE_HXX__
# define __MEDPARTITIONER_CONNECTZONE_HXX__

#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_SkyLineArray.hxx"

// few STL include files
#include <map>
#include <string>

namespace MEDPARTITIONER {
class MEDPARTITIONER_EXPORT CONNECTZONE
{
private :
  std::string _name;
  std::string _description;
  int _localDomainNumber;
  int _distantDomainNumber;

  ParaMEDMEM::MEDCouplingUMesh * _localMesh;
  ParaMEDMEM::MEDCouplingUMesh * _distantMesh;

  MEDPARTITIONER::MEDSKYLINEARRAY * _nodeCorresp;
  MEDPARTITIONER::MEDSKYLINEARRAY * _faceCorresp;
  
  std::map < std::pair <int,int>, MEDPARTITIONER::MEDSKYLINEARRAY * > _entityCorresp;

public :
  CONNECTZONE();
  ~CONNECTZONE();
  CONNECTZONE(const CONNECTZONE & myConnectZone);

  std::string getName() const ;
  std::string getDescription() const ;
  int getDistantDomainNumber() const ;
  int getLocalDomainNumber() const ;
  ParaMEDMEM::MEDCouplingUMesh * getLocalMesh() const ;
  ParaMEDMEM::MEDCouplingUMesh * getDistantMesh() const ;

  bool isEntityCorrespPresent(int localEntity,int distantEntity) const;
  const int * getNodeCorrespIndex() const;
  const int * getNodeCorrespValue() const;
  int getNodeNumber() const;
  const int * getFaceCorrespIndex() const;
  const int * getFaceCorrespValue() const;
  int getFaceNumber() const;
  const int * getEntityCorrespIndex(int localEntity,
                                    int distantEntity) const;
  const int * getEntityCorrespValue(int localEntity,
                                    int distantEntity) const;
  int getEntityCorrespNumber(int localEntity,
                  int distantEntity) const;
  int getEntityCorrespLength(int localEntity,
                  int distantEntity) const;
  void setName(std::string name) ;
  void setDescription(std::string description) ;
  void setDistantDomainNumber(int distantDomainNumber) ;
  void setLocalDomainNumber(int distantDomainNumber) ;
  void setLocalMesh(ParaMEDMEM::MEDCouplingUMesh * localMesh) ;
  void setDistantMesh(ParaMEDMEM::MEDCouplingUMesh * distantMesh) ;

  void setNodeCorresp(int * nodeCorresp, int nbnode);
  void setNodeCorresp(MEDPARTITIONER::MEDSKYLINEARRAY* array);
  void setFaceCorresp(int * faceCorresp, int nbface);
  void setFaceCorresp(MEDPARTITIONER::MEDSKYLINEARRAY* array);
  void setEntityCorresp(int localEntity,
                        int distantEntity,
                        int * entityCorresp, int nbentity);
  void setEntityCorresp(int localEntity,
                        int distantEntity,
                        MEDPARTITIONER::MEDSKYLINEARRAY* array);
};
}
# endif
