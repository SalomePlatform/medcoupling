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

#ifndef __MEDPARTITIONER_PARALLELTOPOLOGY_HXX__
#define __MEDPARTITIONER_PARALLELTOPOLOGY_HXX__

#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_Topology.hxx"

#include "InterpKernelHashMap.hxx"

#include <set>
#include <vector>

namespace MEDPARTITIONER
{
  class Graph;
  class MeshCollection;
  class ParaDomainSelector;

  class MEDPARTITIONER_EXPORT ParallelTopology : public Topology
  {

  public:

    ParallelTopology();
    ParallelTopology(const std::vector<MEDCoupling::MEDCouplingUMesh*>&);
    ParallelTopology(const std::vector<MEDCoupling::MEDCouplingUMesh*>&,
                     const std::vector<MEDPARTITIONER::ConnectZone*>&,
                     std::vector<int*>&,
                     std::vector<int*>&,
                     std::vector<int*>&);
    ParallelTopology(Graph* graph, Topology* oldTopology, int nbdomain, int mesh_dimension);
    ~ParallelTopology();
    
    void setGlobalNumerotationDefault(ParaDomainSelector* domainSelector);

    /*! converts a list of global cell numbers
     * to a distributed array with local cell numbers
     */
    void convertGlobalNodeList(const int*, int,int*,int*);
    void convertGlobalNodeList(const int*, int,int*,int);
    void convertGlobalNodeListWithTwins(const int* face_list, int nbnode, int*& local, int*& ip, int*& full_array, int& size);

    /*! converts a list of global node numbers
     * to a distributed array with local cell numbers
     */
    void convertGlobalCellList(const int*, int , int*, int *);

    /*! converts a list of global face numbers
     *  to a distributed array with local face numbers
     */
    void convertGlobalFaceList(const int*, int , int*, int *);  
    void convertGlobalFaceList(const int*, int , int*, int);  
    void convertGlobalFaceListWithTwins(const int* face_list, int nbface, int*& local, int*& ip, int*& full_array,int& size);

    /*! converting node global numberings to local numberings */
    void convertToLocal2ndVersion(int* nodes, int nbnodes, int idomain);

    /*! converting node local numbering to global */
    int convertNodeToGlobal(int ip, int icell) const { return _node_loc_to_glob[ip][icell]; }

    /*! converting face local numbering to global */
    int convertFaceToGlobal(int ip, int iface) const { return _face_loc_to_glob[ip][iface]; }

    /*! converting cell global numbering to local */
    int convertCellToGlobal(int ip, int icell) const { return _loc_to_glob[ip][icell]; }

    void convertNodeToGlobal(int ip, const int* local, int n, int *global) const
    {
      for (int i=0; i<n; i++)
        global[i]=_node_loc_to_glob[ip][local[i]];
    }

    void convertCellToGlobal(int ip, const int* local, int n, int *global) const
    {
      for (int i=0; i<n; i++)
        global[i]=_loc_to_glob[ip][local[i]];  
    }

    void convertFaceToGlobal(int ip, const int* local, int n, int *global) const
    {
      for (int i=0; i<n; i++) 
        global[i]=_face_loc_to_glob[ip][local[i]];
    }

    int nbDomain() const { return _nb_domain; }

    int nbCells() const { return _nb_total_cells; }
    
    int nbNodes() const { return _nb_total_nodes; }

    int nbCells( int idomain) const { return _nb_cells[idomain]; }

    /*! retrieving number of nodes */
    int getNodeNumber(int idomain) const { return _nb_nodes[idomain]; }

    int getNodeNumber() const;

    void getNodeList(int idomain, int* list) const;

    /*! retrieving cell numbers after merging in parallel mode */
    std::vector<int> & getFusedCellNumbers(int idomain) { return _cell_loc_to_glob_fuse[idomain]; }
    
    const std::vector<int>& getFusedCellNumbers(int idomain) const { return _cell_loc_to_glob_fuse[idomain]; }

    /*! retrieving face numbers after merging in parallel mode */
    std::vector<int> & getFusedFaceNumbers(int idomain) { return _face_loc_to_glob_fuse[idomain]; }

    const std::vector<int>& getFusedFaceNumbers(int idomain) const { return _face_loc_to_glob_fuse[idomain]; }

    /*! retrieving number of nodes */
    int getCellNumber(int idomain) const { return _nb_cells[idomain]; }

    int getCellDomainNumber(int global) const { return (_glob_to_loc.find(global)->second).first; }

    void getCellList(int idomain, int* list) const;

    int getFaceNumber(int idomain) const { return _nb_faces[idomain]; }

    int getFaceNumber() const;

    void getFaceList(int idomain, int* list) const;

    /*! converting a global cell number to a local representation (domain + local number) */
    std::pair<int,int> convertGlobalCell(int iglobal) const { return _glob_to_loc.find(iglobal)->second; }

    int convertGlobalFace(int iglobal, int idomain);

    int convertGlobalNode(int iglobal, int idomain);
    
    std::vector<MEDPARTITIONER::ConnectZone*>& getCZ();

    //adding a face to the topology
    void appendFace(int idomain, int ilocal, int iglobal);

    //return max global face number
    int getMaxGlobalFace() const;

  private:
    bool hasCellWithNodes( const MeshCollection&, int dom, const std::set<int>& nodes );

  private:
    //mapping global -> local
    typedef INTERP_KERNEL::HashMultiMap<int,std::pair<int,int> > TGlob2DomainLoc;

    TGlob2DomainLoc _glob_to_loc;
    TGlob2DomainLoc _node_glob_to_loc;

    //mapping local -> global
    std::vector<std::vector<int> >  _loc_to_glob;
    std::vector<std::vector <int> > _node_loc_to_glob;

    // global numbers in parallel mode
    std::vector<std::vector <int> > _cell_loc_to_glob_fuse; // glob nums after merging
    std::vector<std::vector <int> > _face_loc_to_glob_fuse; // glob nums after merging

    //mapping global -> local
    typedef INTERP_KERNEL::HashMultiMap<int,std::pair<int,int> > TGlob2LocsMap;
    TGlob2LocsMap _face_glob_to_loc;

    //mapping local -> global
    std::vector<std::vector <int> > _face_loc_to_glob;
    std::vector<int> _nb_cells;
    std::vector<int> _nb_nodes;
    std::vector<int> _nb_faces;
    int _nb_total_cells;
    int _nb_total_nodes;
    int _nb_total_faces;
    int _nb_domain;
    int _mesh_dimension;

    //links to connectzones
    std::vector<MEDPARTITIONER::ConnectZone*> _connect_zones;

  };
}
#endif
