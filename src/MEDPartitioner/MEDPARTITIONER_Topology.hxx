// Copyright (C) 2007-2025  CEA, EDF
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

#ifndef __MEDPARTITIONER_TOPOLOGY_HXX__
#define __MEDPARTITIONER_TOPOLOGY_HXX__

#include "MEDPARTITIONER.hxx"
#include "MCType.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class MEDCouplingUMesh;
}

namespace MEDPARTITIONER
{
  class Graph;
  class ConnectZone;
  class MeshCollection;
  class MEDPARTITIONER_FaceModel;

  class MEDPARTITIONER_EXPORT Topology
  {
  public:
    Topology() { }
    Topology(std::vector<MEDCoupling::MEDCouplingUMesh*>, std::vector<MEDPARTITIONER::ConnectZone*>) { }
    virtual ~Topology() { }

    /*! converts a list of global cell numbers
     *  to a distributed array with local cell numbers
     */
    virtual void convertGlobalNodeList(const mcIdType *list, mcIdType nb, mcIdType *local, int*ip) = 0;
    virtual void convertGlobalNodeList(const mcIdType *list, mcIdType nb, mcIdType *local, int ip) = 0;
    //converts a list of global node numbers
    /*! to a distributed array with local cell numbers */
    virtual void convertGlobalCellList(const mcIdType*list , mcIdType nb, mcIdType *local, int*ip) = 0;

    /*! converts a list of global face numbers
     *  to a distributed array with local face numbers
     */
     virtual void convertGlobalFaceList(const mcIdType*list , mcIdType nb, mcIdType* local, int*ip) = 0;
    virtual void convertGlobalFaceList(const mcIdType*list , mcIdType nb, mcIdType* local, int ip) = 0;
    virtual void convertGlobalFaceListWithTwins(const mcIdType *face_list, mcIdType nbface, mcIdType*& local, int*& ip, mcIdType*& full_array, mcIdType& size) = 0;
    virtual void convertGlobalNodeListWithTwins(const mcIdType *face_list, mcIdType nbnode, mcIdType*& local, int*& ip, mcIdType*& full_array, mcIdType& size) = 0;
    /*! number of doamins */
    virtual int nbDomain() const = 0;
    /*! number of cells */
    virtual mcIdType nbCells() const = 0;
    /*! number of nodes */
    virtual mcIdType nbNodes() const = 0;
    /*! number of cells on a specific domain */
    virtual mcIdType nbCells(int idomain) const = 0;
    /*! converting node global numberings to local numberings */
    virtual void convertToLocal2ndVersion(mcIdType*,mcIdType,int) = 0;
    virtual mcIdType convertNodeToGlobal(int ip,mcIdType icell) const = 0;
    virtual mcIdType convertFaceToGlobal(int ip,mcIdType icell) const = 0;
    virtual mcIdType convertCellToGlobal(int ip,mcIdType icell) const = 0;
    virtual void convertNodeToGlobal(int ip,const mcIdType *local, mcIdType n, mcIdType *global) const = 0;
    virtual void convertCellToGlobal(int ip,const mcIdType *local, mcIdType n, mcIdType *global) const = 0;
    virtual void convertFaceToGlobal(int ip,const mcIdType *local, mcIdType n, mcIdType *global) const = 0;
    /*! retrieving number of nodes */
    virtual mcIdType getNodeNumber(int idomain) const = 0;
    virtual mcIdType getNodeNumber() const = 0;
    /*! retrieving list of nodes */
    virtual void getNodeList(int idomain, mcIdType *list) const = 0;
    virtual std::vector<mcIdType> & getFusedCellNumbers(int idomain) = 0;
    virtual const std::vector<mcIdType> & getFusedCellNumbers(int idomain) const = 0;
    virtual std::vector<mcIdType> & getFusedFaceNumbers(int idomain) = 0;
    virtual const std::vector<mcIdType> & getFusedFaceNumbers(int idomain) const = 0;
    /*! retrieving number of nodes */
    virtual mcIdType getCellNumber(int idomain) const = 0;
    /*! retrieving list of nodes */
    virtual void getCellList(int idomain, mcIdType *list) const = 0;
    /*! retrieving number of faces */
    virtual mcIdType getFaceNumber(int idomain) const = 0;
    virtual mcIdType getFaceNumber() const = 0;
    /*! retrieving list of nodes */
    virtual void getFaceList(int idomain, mcIdType *list) const = 0;
    /*! adding a face to the mapping */
    virtual void appendFace(int idomain, mcIdType ilocal, mcIdType iglobal) = 0;
    /*! returns max global face number */
    virtual mcIdType getMaxGlobalFace() const = 0;
    /*! converting a global cell number to a local representation */
    virtual std::pair<int,mcIdType> convertGlobalCell(mcIdType iglobal) const = 0;
    /*! converting a global face number to a local representation */
    virtual mcIdType convertGlobalFace(mcIdType iglobal, int idomain) = 0;
    /*! converting a global node number to a local representation */
    virtual mcIdType convertGlobalNode(mcIdType iglobal, int idomain) = 0;
    /*! getting a reference to connect zones vector */
    virtual std::vector<MEDPARTITIONER::ConnectZone*>& getCZ() = 0;
  };
}

#endif
