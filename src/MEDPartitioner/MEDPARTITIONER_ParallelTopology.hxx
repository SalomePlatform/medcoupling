// Copyright (C) 2007-2021  CEA/DEN, EDF R&D
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
                     std::vector<mcIdType*>&,
                     std::vector<mcIdType*>&,
                     std::vector<mcIdType*>&);
    ParallelTopology(Graph* graph, Topology* oldTopology, int nbdomain, int mesh_dimension);
    ~ParallelTopology();
    
    void setGlobalNumerotationDefault(ParaDomainSelector* domainSelector);

    /*! converts a list of global cell numbers
     * to a distributed array with local cell numbers
     */
    void convertGlobalNodeList(const mcIdType*, mcIdType,mcIdType*,int*);
    void convertGlobalNodeList(const mcIdType*, mcIdType,mcIdType*,int);
    void convertGlobalNodeListWithTwins(const mcIdType* face_list, mcIdType nbnode, mcIdType*& local, int*& ip, mcIdType*& full_array, mcIdType& size);

    /*! converts a list of global node numbers
     * to a distributed array with local cell numbers
     */
    void convertGlobalCellList(const mcIdType*, mcIdType , mcIdType*, int *);

    /*! converts a list of global face numbers
     *  to a distributed array with local face numbers
     */
    void convertGlobalFaceList(const mcIdType*, mcIdType , mcIdType*, int *);  
    void convertGlobalFaceList(const mcIdType*, mcIdType , mcIdType*, int);  
    void convertGlobalFaceListWithTwins(const mcIdType* face_list, mcIdType nbface, mcIdType*& local, int*& ip, mcIdType*& full_array,mcIdType& size);

    /*! converting node global numberings to local numberings */
    void convertToLocal2ndVersion(mcIdType* nodes, mcIdType nbnodes, int idomain);

    /*! converting node local numbering to global */
    mcIdType convertNodeToGlobal(int ip, mcIdType icell) const { return _node_loc_to_glob[ip][icell]; }

    /*! converting face local numbering to global */
    mcIdType convertFaceToGlobal(int ip, mcIdType iface) const { return _face_loc_to_glob[ip][iface]; }

    /*! converting cell global numbering to local */
    mcIdType convertCellToGlobal(int ip, mcIdType icell) const { return _loc_to_glob[ip][icell]; }

    void convertNodeToGlobal(int ip, const mcIdType* local, mcIdType n, mcIdType *global) const
    {
      for (mcIdType i=0; i<n; i++)
        global[i]=_node_loc_to_glob[ip][local[i]];
    }

    void convertCellToGlobal(int ip, const mcIdType* local, mcIdType n, mcIdType *global) const
    {
      for (mcIdType i=0; i<n; i++)
        global[i]=_loc_to_glob[ip][local[i]];  
    }

    void convertFaceToGlobal(int ip, const mcIdType* local, mcIdType n, mcIdType *global) const
    {
      for (mcIdType i=0; i<n; i++) 
        global[i]=_face_loc_to_glob[ip][local[i]];
    }

    int nbDomain() const { return _nb_domain; }

    mcIdType nbCells() const { return _nb_total_cells; }
    
    mcIdType nbNodes() const { return _nb_total_nodes; }

    mcIdType nbCells( int idomain) const { return _nb_cells[idomain]; }

    /*! retrieving number of nodes */
    mcIdType getNodeNumber(int idomain) const { return _nb_nodes[idomain]; }

    mcIdType getNodeNumber() const;

    void getNodeList(int idomain, mcIdType* list) const;

    /*! retrieving cell numbers after merging in parallel mode */
    std::vector<mcIdType> & getFusedCellNumbers(int idomain) { return _cell_loc_to_glob_fuse[idomain]; }
    
    const std::vector<mcIdType>& getFusedCellNumbers(int idomain) const { return _cell_loc_to_glob_fuse[idomain]; }

    /*! retrieving face numbers after merging in parallel mode */
    std::vector<mcIdType> & getFusedFaceNumbers(int idomain) { return _face_loc_to_glob_fuse[idomain]; }

    const std::vector<mcIdType>& getFusedFaceNumbers(int idomain) const { return _face_loc_to_glob_fuse[idomain]; }

    /*! retrieving number of nodes */
    mcIdType getCellNumber(int idomain) const { return _nb_cells[idomain]; }

    mcIdType getCellDomainNumber(int global) const { return (_glob_to_loc.find(global)->second).first; }

    void getCellList(int idomain, mcIdType* list) const;

    mcIdType getFaceNumber(int idomain) const { return _nb_faces[idomain]; }

    mcIdType getFaceNumber() const;

    void getFaceList(int idomain, mcIdType* list) const;

    /*! converting a global cell number to a local representation (domain + local number) */
    std::pair<int,mcIdType> convertGlobalCell(mcIdType iglobal) const { return _glob_to_loc.find(iglobal)->second; }

    mcIdType convertGlobalFace(mcIdType iglobal, int idomain);

    mcIdType convertGlobalNode(mcIdType iglobal, int idomain);
    
    std::vector<MEDPARTITIONER::ConnectZone*>& getCZ();

    //adding a face to the topology
    void appendFace(int idomain, mcIdType ilocal, mcIdType iglobal);

    //return max global face number
    mcIdType getMaxGlobalFace() const;

  private:
    bool hasCellWithNodes( const MeshCollection&, int dom, const std::set<mcIdType>& nodes );

  private:
    //mapping global -> local
    typedef INTERP_KERNEL::HashMultiMap<mcIdType,std::pair<int,mcIdType> > TGlob2DomainLoc;

    TGlob2DomainLoc _glob_to_loc;
    TGlob2DomainLoc _node_glob_to_loc;

    //mapping local -> global
    std::vector<std::vector<mcIdType> >  _loc_to_glob;
    std::vector<std::vector <mcIdType> > _node_loc_to_glob;

    // global numbers in parallel mode
    std::vector<std::vector <mcIdType> > _cell_loc_to_glob_fuse; // glob nums after merging
    std::vector<std::vector <mcIdType> > _face_loc_to_glob_fuse; // glob nums after merging

    //mapping global -> local
    typedef INTERP_KERNEL::HashMultiMap<mcIdType,std::pair<int,mcIdType> > TGlob2LocsMap;
    TGlob2LocsMap _face_glob_to_loc;

    //mapping local -> global
    std::vector<std::vector <mcIdType> > _face_loc_to_glob;
    std::vector<mcIdType> _nb_cells;
    std::vector<mcIdType> _nb_nodes;
    std::vector<mcIdType> _nb_faces;
    mcIdType _nb_total_cells;
    mcIdType _nb_total_nodes;
    mcIdType _nb_total_faces;
    int _nb_domain;
    int _mesh_dimension;

    //links to connectzones
    std::vector<MEDPARTITIONER::ConnectZone*> _connect_zones;

  };
}
#endif
