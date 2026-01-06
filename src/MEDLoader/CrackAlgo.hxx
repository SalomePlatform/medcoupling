// Copyright (C) 2007-2026  CEA, EDF
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
// See http://www.salome-platform.org/ or email :
// webmaster.salome@opencascade.com
//
// Author : Aymeric SONOLET (CEA/DES)

#ifndef SRC_MEDLOADER_CRACKALGO_HXX_
#define SRC_MEDLOADER_CRACKALGO_HXX_

#include <map>
#include <vector>
#include <memory>
#include <string>
#include <set>
#include <utility>

#include <MCType.hxx>

namespace MEDCoupling
{

class MEDFileUMesh;
class DataArrayIdType;
class MEDCouplingUMesh;

class CrackAlgo
{
   public:
    using Set = std::set<mcIdType>;
    using Map = std::map<mcIdType, mcIdType>;
    using Graph = std::map<mcIdType, Set>;
    using Map2Set = std::map<mcIdType, Set>;
    using Map2Map = std::map<mcIdType, Map>;

    CrackAlgo() {}
    ~CrackAlgo() {}
    static Map2Map Compute(MEDFileUMesh *mm, const std::string &grp_name, bool grpMustBeFullyDup = true);

    static void OpenCrack(MEDFileUMesh *mf, const Map2Map &cellOld2NewNode, const double &factor = 0.9);

   private:
    /* Find connected components using a node to node graph.
     */
    static std::vector<std::shared_ptr<std::vector<mcIdType>>> FindConnectedComponents(const Graph &graph);

    /* Depth First Search function to find connected components.
     */
    static void Dfs(
        const Graph &graph,
        const mcIdType &node,
        const std::size_t &componentId,
        std::map<mcIdType, bool> &visited,
        std::map<mcIdType, std::size_t> &componentMap
    );

    /* Converts DataArrayIdType to Set.
     *
     * Usefull for doing multiple search in the Set efficiently.
     */
    static Set DataArrayToSet(const DataArrayIdType &da);

    static DataArrayIdType *SetToDataArray(const Set &s);

    static Set GetCellsTouchingNodesToDup(
        const MEDCouplingUMesh *mf,
        const DataArrayIdType *n2cIdx,
        const DataArrayIdType *n2c,
        const DataArrayIdType *f2dup
    );

    static Map2Set GetNode2CellMap(
        const MEDCouplingUMesh *mf,
        const DataArrayIdType *n2cIdx,
        const DataArrayIdType *n2c,
        const DataArrayIdType *f2dup
    );

    /* Building the cell to cell graph.
     *
     * This graph concerns only cells touching nodes touching faces to duplicate.
     * The connection between cells sharing a face to dup is cut in this graph.
     */
    static Graph BuildCutC2CGraph(
        const DataArrayIdType *c2fIdx,
        const DataArrayIdType *c2f,
        const DataArrayIdType *f2cIdx,
        const DataArrayIdType *f2c,
        const Set &cTouchingN_dup,
        const DataArrayIdType *f2dup
    );

    static Map2Map CreateNewNodesInTopLevelMesh(const Map2Set &n2c_dup, const Graph &c2c, MEDCouplingUMesh *m0);

    static void AddMissingElementsOnLevelM1AndChangeConnectivity(
        const DataArrayIdType *f2cIdx,
        const DataArrayIdType *f2c,
        const DataArrayIdType *f2dupIdInM1,
        const DataArrayIdType *f2dupIdInMf,
        const Map2Map &cellOld2NewNode,
        MEDCouplingUMesh *m1,
        bool grpMustBeFullyDup = true
    );

    /* Create new family array for level -1, which is the extended copy of the
     * original as new elements are appended. Caller owns newFamily DAI.
     * Supposes that the elements to duplicate are all duplicated and appended
     * in the order at the end of mf.
     */
    static DataArrayIdType *CopyAndCompleteFamilyArrAtLevelM1(
        const MEDFileUMesh *mm, const MEDCouplingUMesh *mf, const DataArrayIdType *f2dup
    );

    static DataArrayIdType *CopyFamilyArrAtLev0(const MEDFileUMesh *mm);

    /* Manage node families.
     *
     * Extend the family size inplace and set new nodes family to their
     * predecessor.
     */
    static void CompleteFamilyArrAtNodeLevel(const Map2Set &addedNodes, MEDFileUMesh *mm);

    /* Aggregate CellOld2NewNodes into oldNode to newNodes.
     */
    static Map2Set BuildMap2Set(const Map2Map &cellOld2NewNode);

    /* Remomves cells from crackingMesh which are part of m skin.
     *
     * Returns a new MEDCouplingUMesh. User is responsible to delete it.
     */
    static MEDCouplingUMesh *CleanM1Mesh(const MEDCouplingUMesh &m, const MEDCouplingUMesh &crackingMesh);

    static std::pair<DataArrayIdType *, DataArrayIdType *> GetFacesInM1TouchingDuplicatedNodes(
        const Map2Set &n2c_dup,
        const DataArrayIdType *f2dupIdInM1,
        const MEDCouplingUMesh &mf,
        const MEDCouplingUMesh &m1
    );

    static DataArrayIdType *GetFacesToDupInM1(const MEDCouplingUMesh &crackMesh, const MEDCouplingUMesh &m1);

    static void ChangeConnectivityOfM1Elements(
        const DataArrayIdType *f2changeIdInM1,
        const DataArrayIdType *f2changeIdInMf,
        const Map2Map &cellOld2New,
        const DataArrayIdType *f2cIdx,
        const DataArrayIdType *f2c,
        MEDCouplingUMesh *m1
    );
};
}  // namespace MEDCoupling
#endif  // SRC_MEDLOADER_CRACKALGO_HXX_
