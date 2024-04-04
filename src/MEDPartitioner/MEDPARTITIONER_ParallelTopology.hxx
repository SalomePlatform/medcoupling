// Copyright (C) 2007-2024  CEA, EDF
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

#include "MCIdType.hxx"
#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_Topology.hxx"

#include "InterpKernelHashMap.hxx"

#include <set>
#include <utility>
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
    ~ParallelTopology() override;
    
    void setGlobalNumerotationDefault(ParaDomainSelector* domainSelector);

    /*! converts a list of global cell numbers
     * to a distributed array with local cell numbers
     */
    void convertGlobalNodeList(const mcIdType*, mcIdType,mcIdType*,int*) override;
    void convertGlobalNodeList(const mcIdType*, mcIdType,mcIdType*,int) override;
    void convertGlobalNodeListWithTwins(const mcIdType* face_list, mcIdType nbnode, mcIdType*& local, int*& ip, mcIdType*& full_array, mcIdType& size) override;

    /*! converts a list of global node numbers
     * to a distributed array with local cell numbers
     */
    void convertGlobalCellList(const mcIdType*, mcIdType , mcIdType*, int *) override;

    /*! converts a list of global face numbers
     *  to a distributed array with local face numbers
     */
    void convertGlobalFaceList(const mcIdType*, mcIdType , mcIdType*, int *) override;  
    void convertGlobalFaceList(const mcIdType*, mcIdType , mcIdType*, int) override;  
    void convertGlobalFaceListWithTwins(const mcIdType* face_list, mcIdType nbface, mcIdType*& local, int*& ip, mcIdType*& full_array,mcIdType& size) override;

    /*! converting node global numberings to local numberings */
    void convertToLocal2ndVersion(mcIdType* nodes, mcIdType nbnodes, int idomain) override;

    /*! converting node local numbering to global */
    mcIdType convertNodeToGlobal(int ip, mcIdType icell) const override { return _node_loc_to_glob[ip][icell]; }

    /*! converting face local numbering to global */
    mcIdType convertFaceToGlobal(int ip, mcIdType iface) const override { return _face_loc_to_glob[ip][iface]; }

    /*! converting cell global numbering to local */
    mcIdType convertCellToGlobal(int ip, mcIdType icell) const override { return _loc_to_glob[ip][icell]; }

    void convertNodeToGlobal(int ip, const mcIdType* local, mcIdType n, mcIdType *global) const override
    {
      for (mcIdType i=0; i<n; i++)
        global[i]=_node_loc_to_glob[ip][local[i]];
    }

    void convertCellToGlobal(int ip, const mcIdType* local, mcIdType n, mcIdType *global) const override
    {
      for (mcIdType i=0; i<n; i++)
        global[i]=_loc_to_glob[ip][local[i]];  
    }

    void convertFaceToGlobal(int ip, const mcIdType* local, mcIdType n, mcIdType *global) const override
    {
      for (mcIdType i=0; i<n; i++) 
        global[i]=_face_loc_to_glob[ip][local[i]];
    }

    int nbDomain() const override { return _nb_domain; }

    mcIdType nbCells() const override { return _nb_total_cells; }
    
    mcIdType nbNodes() const override { return _nb_total_nodes; }

    mcIdType nbCells( int idomain) const override { return _nb_cells[idomain]; }

    /*! retrieving number of nodes */
    mcIdType getNodeNumber(int idomain) const override { return _nb_nodes[idomain]; }

    mcIdType getNodeNumber() const override;

    void getNodeList(int idomain, mcIdType* list) const override;

    /*! retrieving cell numbers after merging in parallel mode */
    std::vector<mcIdType> & getFusedCellNumbers(int idomain) override { return _cell_loc_to_glob_fuse[idomain]; }
    
    const std::vector<mcIdType>& getFusedCellNumbers(int idomain) const override { return _cell_loc_to_glob_fuse[idomain]; }

    /*! retrieving face numbers after merging in parallel mode */
    std::vector<mcIdType> & getFusedFaceNumbers(int idomain) override { return _face_loc_to_glob_fuse[idomain]; }

    const std::vector<mcIdType>& getFusedFaceNumbers(int idomain) const override { return _face_loc_to_glob_fuse[idomain]; }

    /*! retrieving number of nodes */
    mcIdType getCellNumber(int idomain) const override { return _nb_cells[idomain]; }

    mcIdType getCellDomainNumber(int global) const { return (_glob_to_loc.find(global)->second).first; }

    void getCellList(int idomain, mcIdType* list) const override;

    mcIdType getFaceNumber(int idomain) const override { return _nb_faces[idomain]; }

    mcIdType getFaceNumber() const override;

    void getFaceList(int idomain, mcIdType* list) const override;

    /*! converting a global cell number to a local representation (domain + local number) */
    std::pair<int,mcIdType> convertGlobalCell(mcIdType iglobal) const override { return _glob_to_loc.find(iglobal)->second; }

    mcIdType convertGlobalFace(mcIdType iglobal, int idomain) override;

    mcIdType convertGlobalNode(mcIdType iglobal, int idomain) override;
    
    std::vector<MEDPARTITIONER::ConnectZone*>& getCZ() override;

    //adding a face to the topology
    void appendFace(int idomain, mcIdType ilocal, mcIdType iglobal) override;

    //return max global face number
    mcIdType getMaxGlobalFace() const override;

  private:
    bool hasCellWithNodes( const MeshCollection&, int dom, const std::set<mcIdType>& nodes );

  private:
    //mapping global -> local
    using TGlob2DomainLoc = INTERP_KERNEL::HashMultiMap<mcIdType, std::pair<int, mcIdType>>;

    TGlob2DomainLoc _glob_to_loc;
    TGlob2DomainLoc _node_glob_to_loc;

    //mapping local -> global
    std::vector<std::vector<mcIdType> >  _loc_to_glob;
    std::vector<std::vector <mcIdType> > _node_loc_to_glob;

    // global numbers in parallel mode
    std::vector<std::vector <mcIdType> > _cell_loc_to_glob_fuse; // glob nums after merging
    std::vector<std::vector <mcIdType> > _face_loc_to_glob_fuse; // glob nums after merging

    //mapping global -> local
    using TGlob2LocsMap = INTERP_KERNEL::HashMultiMap<mcIdType, std::pair<int, mcIdType>>;
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
