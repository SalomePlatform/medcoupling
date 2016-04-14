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

#ifndef __MEDPARTITIONER_TOPOLOGY_HXX__
#define __MEDPARTITIONER_TOPOLOGY_HXX__

#include "MEDPARTITIONER.hxx"

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
    virtual void convertGlobalNodeList(const int *list, int nb, int *local, int*ip) = 0;
    virtual void convertGlobalNodeList(const int *list, int nb, int *local, int ip) = 0;
    //converts a list of global node numbers
    /*! to a distributed array with local cell numbers */
    virtual void convertGlobalCellList(const int*list , int nb, int *local, int*ip) = 0;
    
    /*! converts a list of global face numbers
     *  to a distributed array with local face numbers
     */
     virtual void convertGlobalFaceList(const int*list , int nb, int* local, int*ip) = 0;
    virtual void convertGlobalFaceList(const int*list , int nb, int* local, int ip) = 0;
    virtual void convertGlobalFaceListWithTwins(const int *face_list, int nbface, int*& local, int*& ip, int*& full_array, int& size) = 0;
    virtual void convertGlobalNodeListWithTwins(const int *face_list, int nbnode, int*& local, int*& ip, int*& full_array, int& size) = 0;
    /*! number of doamins */
    virtual int nbDomain() const = 0;
    /*! number of cells */
    virtual int nbCells() const = 0;
    /*! number of nodes */
    virtual int nbNodes() const = 0;
    /*! number of cells on a specific domain */
    virtual int nbCells(int idomain) const = 0;
    /*! converting node global numberings to local numberings */
    virtual void convertToLocal2ndVersion(int*,int,int) = 0;
    virtual int convertNodeToGlobal(int ip,int icell) const = 0;
    virtual int convertFaceToGlobal(int ip,int icell) const = 0;
    virtual int convertCellToGlobal(int ip,int icell) const = 0;
    virtual void convertNodeToGlobal(int ip,const int *local, int n, int *global) const = 0;
    virtual void convertCellToGlobal(int ip,const int *local, int n, int *global) const = 0;
    virtual void convertFaceToGlobal(int ip,const int *local, int n, int *global) const = 0;
    /*! retrieving number of nodes */
    virtual int getNodeNumber(int idomain) const = 0;
    virtual int getNodeNumber() const = 0;
    /*! retrieving list of nodes */
    virtual void getNodeList(int idomain, int *list) const = 0;
    virtual std::vector<int> & getFusedCellNumbers(int idomain) = 0;
    virtual const std::vector<int> & getFusedCellNumbers(int idomain) const = 0;
    virtual std::vector<int> & getFusedFaceNumbers(int idomain) = 0;
    virtual const std::vector<int> & getFusedFaceNumbers(int idomain) const = 0;
    /*! retrieving number of nodes */
    virtual int getCellNumber(int idomain) const = 0;
    /*! retrieving list of nodes */
    virtual void getCellList(int idomain, int *list) const = 0;
    /*! retrieving number of faces */
    virtual int getFaceNumber(int idomain) const = 0;
    virtual int getFaceNumber() const = 0;
    /*! retrieving list of nodes */
    virtual void getFaceList(int idomain, int *list) const = 0;
    /*! adding a face to the mapping */
    virtual void appendFace(int idomain, int ilocal, int iglobal) = 0;
    /*! returns max global face number */
    virtual int getMaxGlobalFace() const = 0;
    /*! converting a global cell number to a local representation */
    virtual std::pair<int,int> convertGlobalCell(int iglobal) const = 0;
    /*! converting a global face number to a local representation */
    virtual int convertGlobalFace(int iglobal, int idomain) = 0;
    /*! converting a global node number to a local representation */
    virtual int convertGlobalNode(int iglobal, int idomain) = 0;
    /*! getting a reference to connect zones vector */
    virtual std::vector<MEDPARTITIONER::ConnectZone*>& getCZ() = 0;
  };
}

#endif
