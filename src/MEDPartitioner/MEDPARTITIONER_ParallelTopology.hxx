//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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
#ifndef __PARALLELTOPOLOGY_HXX__
#define __PARALLELTOPOLOGY_HXX__

#include <set>
#include <vector>
#include "InterpKernelHashMap.hxx"

#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"

/*
namespace INTERP_KERNEL
{
  template<> struct hash< std::pair<int,int> >
  {
    size_t operator()( const std::pair<int,int>& x ) const
    {
      return hash< int >()( x.first*1000000+x.second );
    }
  };
}*/

namespace MEDPARTITIONER {

  class Graph;
  class MESHCollection;
  class MEDPARTITIONER_FaceModel;

  class ParallelTopology:public Topology
  {

  public:

    ParallelTopology();

    ParallelTopology(const std::vector<ParaMEDMEM::MEDCouplingUMesh*>&);
    ParallelTopology(const std::vector<ParaMEDMEM::MEDCouplingUMesh*>&,
                     const std::vector<MEDPARTITIONER::CONNECTZONE*>&,
                     std::vector<int*>&,
                     std::vector<int*>&,
                     std::vector<int*>&);

    ParallelTopology(Graph* graph, Topology* oldTopology, int nbdomain, int mesh_dimension);

    ~ParallelTopology();
    
    void setGlobalNumerotationDefault(ParaDomainSelector* domainSelector);

    //!converts a list of global cell numbers
    //!to a distributed array with local cell numbers
    void convertGlobalNodeList(const int*, int,int*,int*);
    void convertGlobalNodeList(const int*, int,int*,int);
    void convertGlobalNodeListWithTwins(const int* face_list, int nbnode, int*& local, int*& ip, int*& full_array, int& size);

    //!converts a list of global node numbers
    //!to a distributed array with local cell numbers
    void convertGlobalCellList(const int*, int , int*, int *);

    //!converts a list of global face numbers
    //!to a distributed array with local face numbers
    void convertGlobalFaceList(const int*, int , int*, int *);  
    void convertGlobalFaceList(const int*, int , int*, int);  
    void convertGlobalFaceListWithTwins(const int* face_list, int nbface, int*& local, int*& ip, int*& full_array,int& size);

//     //!converting node global numberings to local numberings
    void convertToLocal2ndVersion(int* nodes, int nbnodes, int idomain);

    //!converting node local numbering to global
    inline  int convertNodeToGlobal(int ip,int icell) const
    {
      //return _node_loc_to_glob.find(make_pair(ip,icell))->second;
      return _node_loc_to_glob[ip][icell];
    }

    //!converting face local numbering to global
    inline  int convertFaceToGlobal(int ip,int iface) const
    {
      //     if (_face_loc_to_glob.find(make_pair(ip,icell))==_face_loc_to_glob.end())
      //       return -1;
      //     else
      //return _face_loc_to_glob.find(make_pair(ip,icell))->second;
      return _face_loc_to_glob[ip][iface];
    }

    //converting cell global numbering to local
    inline  int convertCellToGlobal(int ip,int icell) const
    {
      //     if (_loc_to_glob.find(make_pair(ip,icell))==_loc_to_glob.end())
      //       return -1;
      //     else
      //return _loc_to_glob.find(make_pair(ip,icell))->second;
      return _loc_to_glob[ip][icell];
    }

    inline  void convertNodeToGlobal(int ip, const int* local, int n, int* global)const
    {
      for (int i=0; i<n; i++)
        global[i]=_node_loc_to_glob[ip][local[i]];
    }

    inline  void convertCellToGlobal(int ip, const int* local, int n, int* global)const
    {
      for (int i=0; i<n; i++)
        global[i]=_loc_to_glob[ip][local[i]];  
    }

    inline  void convertFaceToGlobal(int ip, const int* local, int n, int* global)const
    {
      for (int i=0; i<n; i++) 
        global[i]=_face_loc_to_glob[ip][local[i]];
    }

    inline  int nbDomain() const
    {
      return _nb_domain;
    }

    int nbCells() const
    {
      return _nb_total_cells;
    }
    
    int nbNodes() const
    {
      return _nb_total_nodes;
    }

    inline  int nbCells( int idomain) const
    {
      return _nb_cells[idomain];
    }

    //!retrieving number of nodes
    inline  int getNodeNumber(int idomain) const
    {
      return _nb_nodes[idomain];
    }

    inline  int getNodeNumber() const
    {
      if (_node_glob_to_loc.empty()) return 0;
      std::set <int> keys;
      for (INTERP_KERNEL::HashMultiMap<int, std::pair<int,int> >::const_iterator iter= _node_glob_to_loc.begin();
           iter!=_node_glob_to_loc.end();
           iter++) {
        keys.insert(iter->first);
      }
      return keys.size();
    }

    //!retrieving list of nodes in global numbers
    inline  void getNodeList(int idomain, int* list) const
    {
      for (int i=0; i<_nb_nodes[idomain]; i++) 
        list[i]=_node_loc_to_glob[idomain][i];
    }

    //!< retrieving cell numbers after fusing in parallel mode
    std::vector<int> & getFusedCellNumbers(int idomain)
    {
      return _cell_loc_to_glob_fuse[idomain];
    }
    
    const std::vector<int> & getFusedCellNumbers(int idomain) const
    {
      return _cell_loc_to_glob_fuse[idomain];
    }

    //!< retrieving face numbers after fusing in parallel mode
    std::vector<int> & getFusedFaceNumbers(int idomain)
    {
      return _face_loc_to_glob_fuse[idomain];
    }
    const std::vector<int> & getFusedFaceNumbers(int idomain) const
    {
      return _face_loc_to_glob_fuse[idomain];
    }


    //!retrieving number of nodes
    inline  int getCellNumber(int idomain) const
    {
      return _nb_cells[idomain];
    }

    inline  int getCellDomainNumber(int global) const
    {
      return (_glob_to_loc.find(global)->second).first;
    }

    //!retrieving list of nodes in global numbers
    inline  void getCellList(int idomain, int* list) const
    {
      for (int i=0; i<_nb_cells[idomain];i++)
        list[i]=_loc_to_glob[idomain][i];
    }

    inline int getFaceNumber(int idomain) const
    {
      return _nb_faces[idomain];
    }

    inline  int getFaceNumber() const
    {
      if (_face_glob_to_loc.empty()) return 0;
      std::set <int> keys;
      for (INTERP_KERNEL::HashMultiMap<int, std::pair<int,int> >::const_iterator iter= _face_glob_to_loc.begin();
           iter!=_face_glob_to_loc.end();
           iter++) {
        keys.insert(iter->first);
      }
      return keys.size();
    }


    //!retrieving list of faces in global numbers
    inline  void getFaceList(int idomain, int* list) const
    {
      for (int i=0; i<_nb_faces[idomain];i++)   
        list[i]=_face_loc_to_glob[idomain][i];
    }

    //! converting a global cell number to a local representation (domain + local number)
    inline std::pair<int,int> convertGlobalCell(int iglobal) const
    {
      return _glob_to_loc.find(iglobal)->second;
    }

    inline int convertGlobalFace(int iglobal, int idomain)
    {
      typedef INTERP_KERNEL::HashMultiMap<int, std::pair<int,int> >::const_iterator MMiter;
      std::pair<MMiter,MMiter> eq = _face_glob_to_loc.equal_range(iglobal);
      for (MMiter it=eq.first; it != eq.second; it++) 
        if (it->second.first == idomain) return it->second.second;   
      return -1;
    }

    inline int convertGlobalNode(int iglobal, int idomain)
    {
      typedef INTERP_KERNEL::HashMultiMap<int, std::pair<int,int> >::const_iterator MMiter;
      std::pair<MMiter,MMiter> eq = _node_glob_to_loc.equal_range(iglobal);
      for (MMiter it=eq.first; it != eq.second; it++)
      {
        if (it->second.first == idomain) return it->second.second;
      }
      return -1;
    }
    
    //!adding a face to the topology
    inline void appendFace(int idomain, int ilocal, int iglobal)
    {
      _face_loc_to_glob[idomain].push_back(iglobal);
      _face_glob_to_loc.insert(std::make_pair(iglobal,std::make_pair(idomain,ilocal)));
    }

    //return max global face number
    int getMaxGlobalFace() const;

  private:
    bool hasCellWithNodes( const MESHCollection&, int dom, const std::set<int>& nodes );

  private:
    //!mapping global -> local
    typedef INTERP_KERNEL::HashMultiMap<int,std::pair<int,int> > TGlob2DomainLoc;

    TGlob2DomainLoc _glob_to_loc;

    std::vector<std::vector<int> >  _loc_to_glob;

    INTERP_KERNEL::HashMultiMap<int,std::pair<int,int> > _node_glob_to_loc;

    //!mapping local -> global
    std::vector<std::vector <int> > _node_loc_to_glob;

    // global numbers in parallel mode
    std::vector<std::vector <int> > _cell_loc_to_glob_fuse; // glob nums after fusing
    std::vector<std::vector <int> > _face_loc_to_glob_fuse; // glob nums after fusing

    //!mapping global -> local
    typedef INTERP_KERNEL::HashMultiMap<int,std::pair<int,int> > TGlob2LocsMap;
    TGlob2LocsMap _face_glob_to_loc;

    //!mapping local -> global
    std::vector<std::vector <int> > _face_loc_to_glob;
    std::vector<int> _nb_cells;
    std::vector<int> _nb_nodes;
    std::vector<int> _nb_faces;
    int _nb_total_cells;
    int _nb_total_nodes;
    int _nb_total_faces;
    int _nb_domain;
    int _mesh_dimension;
  };


}
#endif /*PARALLELTOPOLOGY_HXX_*/
