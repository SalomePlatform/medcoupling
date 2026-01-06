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

#include "CrackAlgo.hxx"

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <memory>
#include <utility>

#include "./MEDFileMesh.hxx"

#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingMemArray.hxx>
#include <InterpKernelException.hxx>
#include <MCType.hxx>
#include <MCAuto.hxx>

using namespace MEDCoupling;
using namespace std;
using MCU = MCAuto<MEDCouplingUMesh>;
using DAI = MCAuto<DataArrayIdType>;
using DAD = MCAuto<DataArrayDouble>;

CrackAlgo::Map2Map
CrackAlgo::Compute(MEDFileUMesh *mm, const string &grp_name, bool grpMustBeFullyDup)
{
    vector<int> levs = mm->getNonEmptyLevels();
    if (find(levs.begin(), levs.end(), 0) == levs.end() || find(levs.begin(), levs.end(), -1) == levs.end())
        throw INTERP_KERNEL::Exception(
            "MEDFileUMesh::duplicateFaces : This method "
            "works only for mesh defined on level 0 and -1 !"
        );

    MCU m0 = mm->getMeshAtLevel(0);

    MCU crackingMesh = mm->getGroup(-1, grp_name);
    MCU cleanCrackingMesh = CleanM1Mesh(*m0, *crackingMesh);

    DAI c2f(DataArrayIdType::New()), c2fIdx(DataArrayIdType::New()), f2c(DataArrayIdType::New()),
        f2cIdx(DataArrayIdType::New());
    const MCU mf = m0->buildDescendingConnectivity(c2f, c2fIdx, f2c, f2cIdx);

    DataArrayIdType *crackToMfP;
    const bool crackCellsInMf = mf->areCellsIncludedIn(cleanCrackingMesh, 2, crackToMfP);
    DAI f2dupIdInMf(crackToMfP);
    if (!crackCellsInMf)
        throw INTERP_KERNEL::Exception("crackAlong: All cells in crack are not part of Mf.");

    DAI n2c(DataArrayIdType::New()), n2cIdx(DataArrayIdType::New());
    m0->getReverseNodalConnectivity(n2c, n2cIdx);

    const Map2Set n2c_dup = GetNode2CellMap(mf, n2cIdx, n2c, f2dupIdInMf);

    const Set cTouchingN_dup = GetCellsTouchingNodesToDup(mf, n2cIdx, n2c, f2dupIdInMf);

    const Graph c2c = BuildCutC2CGraph(c2fIdx, c2f, f2cIdx, f2c, cTouchingN_dup, f2dupIdInMf);

    Map2Map cellOld2NewNode = CreateNewNodesInTopLevelMesh(n2c_dup, c2c, m0);
    // End of node creation, separation of M0

    // Faces in M1 must be added/updated

    MCU m1 = mm->getMeshAtLevel(-1);

    // Before changing M1
    const DAI f2dupIdInM1 = GetFacesToDupInM1(*cleanCrackingMesh, *m1);

    const auto f2changeM1Mf = GetFacesInM1TouchingDuplicatedNodes(n2c_dup, f2dupIdInM1, *mf, *m1);
    const DAI f2changeIdInM1 = f2changeM1Mf.first;
    const DAI f2changeIdInMf = f2changeM1Mf.second;

    AddMissingElementsOnLevelM1AndChangeConnectivity(
        f2cIdx, f2c, f2dupIdInM1, f2dupIdInMf, cellOld2NewNode, m1, grpMustBeFullyDup
    );

    ChangeConnectivityOfM1Elements(f2changeIdInM1, f2changeIdInMf, cellOld2NewNode, f2cIdx, f2c, m1);

    // TODO(aymeric): If one wants to implement
    // AddMissingElementsOnLevelM2AndChangeConnectivity and
    // ChangeConnectivityOfM2Elements it should be there.

    DAI famLev0 = CopyFamilyArrAtLev0(mm);
    DAI newFam = CopyAndCompleteFamilyArrAtLevelM1(mm, m1, f2dupIdInM1);

    mm->setMeshAtLevel(0, m0);
    mm->setMeshAtLevel(-1, m1);
    // mm->setMeshAtLevel(-2, m2);

    mm->setFamilyFieldArr(0, famLev0);
    mm->setFamilyFieldArr(-1, newFam);
    // mm->setFamilyFieldArr(-2, newFam2);

    const Map2Set addedNodes = BuildMap2Set(cellOld2NewNode);
    CompleteFamilyArrAtNodeLevel(addedNodes, mm);

    return cellOld2NewNode;
}

CrackAlgo::Map2Set
CrackAlgo::BuildMap2Set(const Map2Map &cellOld2NewNode)
{
    Map2Set res;
    for (const auto &keyElem : cellOld2NewNode)
    {
        const auto &old2NewMap = keyElem.second;
        for (const auto &old2New : old2NewMap)
        {
            res[old2New.first].insert(old2New.second);
        }
    }
    return res;
}

DataArrayIdType *
CrackAlgo::CopyAndCompleteFamilyArrAtLevelM1(
    const MEDFileUMesh *mm, const MEDCouplingUMesh *mf, const DataArrayIdType *f2dup
)
{
    DataArrayIdType *newFam = nullptr;
    const DataArrayIdType *fam = mm->getFamilyFieldAtLevel(-1);
    if (fam != nullptr)
    {
        newFam = DataArrayIdType::New();
        const auto fam_size = static_cast<size_t>(fam->getNumberOfTuples());
        newFam->alloc(static_cast<size_t>(mf->getNumberOfCells()));
        copy(fam->begin(), fam->end(), newFam->getPointer());
        mcIdType *fam_ptr = newFam->getPointer();
        size_t i_elem_added = 0;
        for (const auto &face : *f2dup)
        {
            // NOTE(aymeric): assumes all faces are duplicated
            fam_ptr[fam_size + i_elem_added] = fam_ptr[face];
            i_elem_added++;
        }
    }
    return newFam;
}

DataArrayIdType *
CrackAlgo::CopyFamilyArrAtLev0(const MEDFileUMesh *mm)
{
    DataArrayIdType *famLev0 = nullptr;
    const DataArrayIdType *fam = mm->getFamilyFieldAtLevel(0);
    if (fam != nullptr)
    {
        famLev0 = DataArrayIdType::New();
        const auto fam_size = static_cast<size_t>(fam->getNumberOfTuples());
        famLev0->alloc(fam_size);
        copy(fam->begin(), fam->end(), famLev0->getPointer());
    }
    return famLev0;
}

void
CrackAlgo::CompleteFamilyArrAtNodeLevel(const Map2Set &addedNodes, MEDFileUMesh *mm)
{
    DataArrayIdType *node_fam = mm->getFamilyFieldAtLevel(1);
    if (node_fam != nullptr)
    {
        node_fam->reAlloc(static_cast<size_t>(mm->getCoords()->getNumberOfTuples()));
        mcIdType *node_fam_ptr = node_fam->getPointer();
        for (const auto &old2NewPair : addedNodes)
        {
            const auto &oldNode = old2NewPair.first;
            const auto &newNodes = old2NewPair.second;
            for (const auto &newNode : newNodes)
            {
                node_fam_ptr[newNode] = node_fam_ptr[oldNode];
            }
        }
    }
}

MEDCouplingUMesh *
CrackAlgo::CleanM1Mesh(const MEDCouplingUMesh &m, const MEDCouplingUMesh &crackingMesh)
{
    MCU m0skin = m.computeSkin();
    DataArrayIdType *idsToKeepP;
    m0skin->areCellsIncludedIn(&crackingMesh, 2, idsToKeepP);
    DAI idsToKeep(idsToKeepP);
    // discard cells on the skin of M0
    DAI ids2 = idsToKeep->findIdsNotInRange(0, m0skin->getNumberOfCells());
    MCU otherDimM1OnSameCoords = crackingMesh.buildPartOfMySelf(ids2->begin(), ids2->end(), true);
    return otherDimM1OnSameCoords.retn();
}

void
CrackAlgo::Dfs(
    const Graph &graph,
    const mcIdType &node,
    const size_t &componentId,
    map<mcIdType, bool> &visited,
    map<mcIdType, size_t> &componentMap
)
{
    visited[node] = true;
    // Assign a unique component ID to the node
    componentMap[node] = componentId;

    for (const auto &neighbor : graph.at(node))
    {
        if (!visited[neighbor])
        {
            Dfs(graph, neighbor, componentId, visited, componentMap);
        }
    }
}

// Function to find connected components
vector<shared_ptr<vector<mcIdType>>>
CrackAlgo::FindConnectedComponents(const Graph &graph)
{
    map<mcIdType, bool> visited;
    map<mcIdType, size_t> componentMap;
    size_t componentId = 0;

    for (const auto &nodePair : graph) visited[nodePair.first] = false;

    for (const auto &nodePair : graph)
    {
        if (!visited[nodePair.first])
        {
            Dfs(graph, nodePair.first, componentId, visited, componentMap);
            componentId++;
        }
    }

    using spv = shared_ptr<vector<mcIdType>>;
    vector<spv> components;
    for (size_t i = 0; i < componentId; i++) components.push_back(make_shared<vector<mcIdType>>());
    for (const auto &nodePair : componentMap)
    {
        components[nodePair.second]->push_back(nodePair.first);
    }

    return components;
}

CrackAlgo::Set
CrackAlgo::DataArrayToSet(const DataArrayIdType &da)
{
    Set res;
    for (const auto &elem : da)
    {
        res.insert(elem);
    }
    return res;
}

DataArrayIdType *
CrackAlgo::SetToDataArray(const Set &s)
{
    DataArrayIdType *da(DataArrayIdType::New());
    da->alloc(s.size());
    mcIdType *da_p = da->rwBegin();
    size_t i = 0;
    for (const auto &elem : s)
    {
        da_p[i] = elem;
        i++;
    }
    return da;
}

CrackAlgo::Graph
CrackAlgo::BuildCutC2CGraph(
    const DataArrayIdType *c2fIdx,
    const DataArrayIdType *c2f,
    const DataArrayIdType *f2cIdx,
    const DataArrayIdType *f2c,
    const Set &cTouchingN_dup,
    const DataArrayIdType *f2dup
)
{
    Set f2dup_set = DataArrayToSet(*f2dup);
    Graph c2c;

    const mcIdType *c2fIdx_p{c2fIdx->begin()};
    const mcIdType *c2f_p{c2f->begin()};
    const mcIdType *f2cIdx_p{f2cIdx->begin()};
    const mcIdType *f2c_p{f2c->begin()};

    for (const auto &cell : cTouchingN_dup)
    {
        c2c[cell] = set<mcIdType>();
        for (auto faceIdx = c2fIdx_p[cell]; faceIdx < c2fIdx_p[cell + 1]; faceIdx++)
        {
            const mcIdType &face = c2f_p[faceIdx];
            // face is not in face to duplicate
            if (f2dup_set.find(face) == f2dup_set.end())
            {
                for (mcIdType cellIdx{f2cIdx_p[face]}; cellIdx < f2cIdx_p[face + 1]; cellIdx++)
                {
                    // getting the id of the other cell (on face is
                    // connected to two cells)
                    const auto &cell2{f2c_p[cellIdx]};
                    if ((cell != cell2) & (cTouchingN_dup.find(cell2) != cTouchingN_dup.end()))
                    {
                        c2c.at(cell).insert(cell2);
                    }
                }
            }
        }
    }
    return c2c;
}

CrackAlgo::Map2Map
CrackAlgo::CreateNewNodesInTopLevelMesh(const Map2Set &n2c_dup, const Graph &c2c, MEDCouplingUMesh *m0)
{
    DataArrayIdType *c2n = m0->getNodalConnectivity();
    DataArrayIdType *c2nIdx = m0->getNodalConnectivityIndex();
    mcIdType *c2n_ptr = c2n->rwBegin();
    const mcIdType *c2nIdx_ptr = c2nIdx->begin();

    DataArrayDouble *coords = m0->getCoords();
    mcIdType i = coords->getNumberOfTuples();
    const size_t coordsDim = coords->getNumberOfComponents();

    Map2Map cellOld2NewNode;

    for (const auto &nc_pair : n2c_dup)
    {
        const auto &node = nc_pair.first;
        const auto &cells = nc_pair.second;

        // Building local connection graph
        Graph c2c_local;
        for (const auto &cell : cells)
        {
            auto &val = c2c_local[cell];
            if (c2c.find(cell) == c2c.end())
                throw INTERP_KERNEL::Exception("A cell touching a node is not part of c2c.");
            for (const auto &neighbor : c2c.at(cell))
            {
                if (cells.find(neighbor) != cells.end())
                    val.insert(neighbor);
            }
        }

        const vector<shared_ptr<vector<mcIdType>>> compo_connex = FindConnectedComponents(c2c_local);

        // Duplicate node for compo[1:]
        const size_t i_nodes_to_add = compo_connex.size() - 1;
        coords->reAlloc(static_cast<size_t>(i) + i_nodes_to_add);

        int i_compo = 0;
        for (const auto &compo : compo_connex)
        {
            if (i_compo > 0)
            {
                // Coords are copied at the end of the vector of coords, ie the
                // node is duplicated
                coords->getTuple(node, coords->getPointer() + static_cast<size_t>(i) * coordsDim);
                for (const auto &cell : *compo)
                {
                    // The node number is replaced in all corresponding cells
                    replace(&c2n_ptr[c2nIdx_ptr[cell] + 1], &c2n_ptr[c2nIdx_ptr[cell + 1]], node, i);
                    // This map is build in order to assign later new nodes to
                    // the corresponding faces
                    cellOld2NewNode[cell][node] = i;
                }
                i++;
            }
            i_compo++;
        }
    }
    return cellOld2NewNode;
}

void
CrackAlgo::AddMissingElementsOnLevelM1AndChangeConnectivity(
    const DataArrayIdType *f2cIdx,
    const DataArrayIdType *f2c,
    const DataArrayIdType *f2dupIdInM1,
    const DataArrayIdType *f2dupIdInMf,
    const Map2Map &cellOld2NewNode,
    MEDCouplingUMesh *m1,
    bool grpMustBeFullyDup
)
{
    DataArrayIdType *f2nIdx = m1->getNodalConnectivityIndex();
    DataArrayIdType *f2n = m1->getNodalConnectivity();

    // Change connectivity index size
    const auto f2nIdx_size = static_cast<size_t>(f2nIdx->getNumberOfTuples());  // = nb_face + 1
    const auto f2dup_size = static_cast<size_t>(f2dupIdInM1->getNumberOfTuples());
    if (f2dup_size != static_cast<size_t>(f2dupIdInMf->getNumberOfTuples()))
        throw INTERP_KERNEL::Exception("There is not as many faces to dup looking into Mf and into M1.");

    f2nIdx->reAlloc(f2nIdx_size + f2dup_size);

    // Change connectivity size
    const auto f2n_size = static_cast<size_t>(f2n->getNumberOfTuples());
    const mcIdType *const f2nIdx_p = f2nIdx->begin();
    {
        size_t new_connectivity_size = f2n_size;
        for (const auto &face : *f2dupIdInM1)
            new_connectivity_size += static_cast<size_t>(f2nIdx_p[face + 1] - f2nIdx_p[face]);
        f2n->reAlloc(new_connectivity_size);
    }

    mcIdType *const f2nIdx_pw = f2nIdx->rwBegin();
    const mcIdType *const f2n_p = f2n->begin();
    mcIdType *const f2n_pw = f2n->rwBegin();
    size_t lastFace = f2nIdx_size - 1;

    const mcIdType *const f2cIdx_p = f2cIdx->begin();
    const mcIdType *const f2c_p = f2c->begin();

    for (size_t iFace = 0; iFace < f2dup_size; iFace++)
    {
        const auto &faceM1 = f2dupIdInM1->begin()[iFace];
        const auto &faceMf = f2dupIdInMf->begin()[iFace];

        // Copy the face appending it
        //
        //   1. Append last connectivity index number
        const mcIdType face_start = f2nIdx_p[faceM1];
        const mcIdType face_end = f2nIdx_p[faceM1 + 1];
        const mcIdType lenFace = face_end - face_start;
        f2nIdx_pw[lastFace + 1] = f2nIdx_p[lastFace] + lenFace;

        //   2. Append last face connectivity
        mcIdType *const new_f2n_p_start = f2n_pw + f2nIdx_p[lastFace];
        copy(f2n_p + face_start, f2n_p + face_end, new_f2n_p_start);

        // Change inplace the connectivity
        // Check that this face is an inner face, ie it is connected to two
        // cells
        const mcIdType d = f2cIdx_p[faceMf + 1] - f2cIdx_p[faceMf];
        if (d != 2)
            throw INTERP_KERNEL::Exception(
                "MEDFileMesh::duplicateFaces, the face to cell (or node to"
                "segment) DataArray does not always adress two cells."
            );

        bool noNodesAreChanged = true;

        //   1. In original face, attaching it to cell0
        const auto &cell0 = f2c_p[f2cIdx_p[faceMf]];
        // If nodes where changed in cell0
        if (cellOld2NewNode.find(cell0) != cellOld2NewNode.end())
        {
            const auto &mapO2N0 = cellOld2NewNode.at(cell0);
            for (mcIdType i_node = f2nIdx_p[faceM1] + 1; i_node < f2nIdx_p[faceM1 + 1]; i_node++)
            {
                const mcIdType node = f2n_p[i_node];
                // If this node was modified
                if (mapO2N0.find(node) != mapO2N0.end())
                {
                    const auto &new_node = mapO2N0.at(node);
                    f2n_pw[i_node] = new_node;
                    noNodesAreChanged = false;
                }
            }
        }

        //   2. In newly created face, attaching it to cell1
        const auto &cell1 = f2c_p[f2cIdx_p[faceMf] + 1];
        // If nodes where changed in cell1
        if (cellOld2NewNode.find(cell1) != cellOld2NewNode.end())
        {
            const auto &mapO2N1 = cellOld2NewNode.at(cell1);
            for (mcIdType i_node = f2nIdx_p[lastFace] + 1; i_node < f2nIdx_p[lastFace + 1]; i_node++)
            {
                const mcIdType node = f2n_p[i_node];
                // If this node was modified
                if (mapO2N1.find(node) != mapO2N1.end())
                {
                    const auto &new_node = mapO2N1.at(node);
                    f2n_pw[i_node] = new_node;
                    noNodesAreChanged = false;
                }
            }
        }

        if (noNodesAreChanged && grpMustBeFullyDup)
            throw INTERP_KERNEL::Exception(
                "duplicateFaces: A face was supposed to be duplicated but could"
                "not be. Please be aware that all M1 groups cannot be"
                "duplicated, at least two adjecent inner faces are needed to"
                "duplicate a node."
            );

        lastFace++;
    }
}

CrackAlgo::Set
CrackAlgo::GetCellsTouchingNodesToDup(
    const MEDCouplingUMesh *mf, const DataArrayIdType *n2cIdx, const DataArrayIdType *n2c, const DataArrayIdType *f2dup
)
{
    Set res{};

    const DataArrayIdType *f2nIdx = mf->getNodalConnectivityIndex();
    const DataArrayIdType *f2n = mf->getNodalConnectivity();

    const mcIdType *f2nIdx_p{f2nIdx->begin()};
    const mcIdType *f2n_p{f2n->begin()};
    const mcIdType *n2cIdx_p{n2cIdx->begin()};
    const mcIdType *n2c_p{n2c->begin()};

    for (const auto face : *f2dup)
        for (mcIdType nodeIdx{f2nIdx_p[face] + 1}; nodeIdx < f2nIdx_p[face + 1]; nodeIdx++)
        {
            const auto &node = f2n_p[nodeIdx];
            for (mcIdType cellIdx = n2cIdx_p[node]; cellIdx < n2cIdx_p[node + 1]; cellIdx++)
            {
                const auto &cell = n2c_p[cellIdx];
                res.insert(cell);
            }
        }
    return res;
}

CrackAlgo::Map2Set
CrackAlgo::GetNode2CellMap(
    const MEDCouplingUMesh *mf, const DataArrayIdType *n2cIdx, const DataArrayIdType *n2c, const DataArrayIdType *f2dup
)
{
    const DataArrayIdType *f2nIdx = mf->getNodalConnectivityIndex();
    const DataArrayIdType *f2n = mf->getNodalConnectivity();

    Map2Set n2c_dup{};
    const mcIdType *f2nIdx_p{f2nIdx->begin()};
    const mcIdType *f2n_p{f2n->begin()};
    const mcIdType *n2cIdx_p{n2cIdx->begin()};
    const mcIdType *n2c_p{n2c->begin()};

    for (const auto face : *f2dup)
        for (mcIdType nodeIdx{f2nIdx_p[face] + 1}; nodeIdx < f2nIdx_p[face + 1]; nodeIdx++)
        {
            const auto &node = f2n_p[nodeIdx];
            for (mcIdType cellIdx{n2cIdx_p[node]}; cellIdx < n2cIdx_p[node + 1]; cellIdx++)
            {
                const auto &cell = n2c_p[cellIdx];
                n2c_dup[node].insert(cell);
            }
        }
    return n2c_dup;
}

void
CrackAlgo::OpenCrack(MEDFileUMesh *mf, const Map2Map &cellOld2NewNode, const double &factor)
{
    if ((factor <= 0.0) || (factor >= 1.0))
        throw INTERP_KERNEL::Exception("factor should be between 0.0 and 1.0");
    DataArrayDouble *coords = mf->getCoords();
    const auto dim = static_cast<mcIdType>(coords->getNumberOfComponents());
    const double *coords_p = coords->begin();
    MCU m0 = mf->getMeshAtLevel(0);
    DAD barys = m0->computeCellCenterOfMass();
    const double *barys_p = barys->begin();
    for (const auto &cellPair : cellOld2NewNode)
    {
        const auto &cell = cellPair.first;
        const auto &mapO2N = cellPair.second;
        for (const auto &nodeO2N : mapO2N)
        {
            const auto &oldN = nodeO2N.first;
            const auto &newN = nodeO2N.second;
            double *pos_new_p = coords->rwBegin() + dim * newN;
            for (int i = 0; i < dim; i++)
            {
                pos_new_p[i] += (1.0 - factor) * (barys_p[cell * dim + i] - coords_p[dim * oldN + i]);
            }
        }
    }
}

DataArrayIdType *
CrackAlgo::GetFacesToDupInM1(const MEDCouplingUMesh &crackMesh, const MEDCouplingUMesh &m1)
{
    DataArrayIdType *f2dupIdInM1P{};
    const bool allIn = m1.areCellsIncludedIn(&crackMesh, 2, f2dupIdInM1P);
    if (!allIn)
        throw INTERP_KERNEL::Exception("crackAlong: cleanCrackMesh is not part of m1 mesh.");
    return f2dupIdInM1P;
}

pair<DataArrayIdType *, DataArrayIdType *>
CrackAlgo::GetFacesInM1TouchingDuplicatedNodes(
    const Map2Set &n2c_dup, const DataArrayIdType *f2dupIdInM1, const MEDCouplingUMesh &mf, const MEDCouplingUMesh &m1
)
{
    Set f2dup_set = DataArrayToSet(*f2dupIdInM1);

    DAI n2f(DataArrayIdType::New());
    DAI n2fIdx(DataArrayIdType::New());
    m1.getReverseNodalConnectivity(n2f, n2fIdx);

    // 1. I retrieve all the faces of m1 which have a potentially duplicated
    // node, I just use the reverse nodal connectivity of m1 for this.
    Set f2changeIdInM1_set{};

    for (const auto &nodeP : n2c_dup)
    {
        const auto &node = nodeP.first;
        for (mcIdType fId = n2fIdx->begin()[node]; fId < n2fIdx->begin()[node + 1]; fId++)
        {
            const auto &f = n2f->begin()[fId];
            // f is not a face to duplicate
            if (f2dup_set.find(f) == f2dup_set.end())
                f2changeIdInM1_set.insert(f);
        }
    }
    // 2. I retrieve M1 sub mesh
    DataArrayIdType *f2changeIdInM1 = SetToDataArray(f2changeIdInM1_set);
    MCU part = m1.buildPartOfMySelf(f2changeIdInM1->begin(), f2changeIdInM1->end());
    // 3. I search for the ids in mf
    DataArrayIdType *f2changeIdInMf{};
    const bool ok = mf.areCellsIncludedIn(part, 2, f2changeIdInMf);
    if (!ok)
        throw INTERP_KERNEL::Exception(
            "Some faces in -1 mesh next to the crack are not included in "
            "the face mesh issued from m0."
        );

    // 4. I return the 2 DAI
    return pair<DataArrayIdType *, DataArrayIdType *>{f2changeIdInM1, f2changeIdInMf};
}

void
CrackAlgo::ChangeConnectivityOfM1Elements(
    const DataArrayIdType *f2changeIdInM1,
    const DataArrayIdType *f2changeIdInMf,
    const Map2Map &cellOld2NewNode,
    const DataArrayIdType *f2cIdx,
    const DataArrayIdType *f2c,
    MEDCouplingUMesh *m1
)
{
    // 1. I go through the faces in question. I go through all the nodes on the
    // face. I will look in neighboring cells to see if there have been changes
    // in connectivity on the nodes of the face.
    // 2. if the changes are consistent, I change in place.
    // 3. otherwise I return a consistency error.

    DataArrayIdType *f2nIdx = m1->getNodalConnectivityIndex();
    DataArrayIdType *f2n = m1->getNodalConnectivity();

    // Change connectivity size
    const mcIdType *const f2nIdx_p = f2nIdx->begin();
    const mcIdType *const f2n_p = f2n->begin();
    mcIdType *const f2n_pw = f2n->rwBegin();

    const mcIdType *const f2cIdx_p{f2cIdx->begin()};
    const mcIdType *const f2c_p{f2c->begin()};

    const auto f2change_size = static_cast<size_t>(f2changeIdInM1->getNumberOfTuples());

    for (size_t iFace = 0; iFace < f2change_size; iFace++)
    {
        const auto &faceM1 = f2changeIdInM1->begin()[iFace];
        const auto &faceMf = f2changeIdInMf->begin()[iFace];

        // Change inplace the connectivity
        // Check that this face is an inner face, ie it is connected to two
        // cells
        const mcIdType d = f2cIdx_p[faceMf + 1] - f2cIdx_p[faceMf];

        //   1. In original face, attaching it to cell0
        const auto &cell0 = f2c_p[f2cIdx_p[faceMf]];
        // const auto& cell1 = f2c_p[f2cIdx_p[faceMf] + 1];

        // If nodes where changed in cell0
        if (cellOld2NewNode.find(cell0) != cellOld2NewNode.end())
        {
            const auto &mapO2N0 = cellOld2NewNode.at(cell0);
            for (mcIdType i_node = f2nIdx_p[faceM1] + 1; i_node < f2nIdx_p[faceM1 + 1]; i_node++)
            {
                const mcIdType node = f2n_p[i_node];
                // If this node was modified
                if (mapO2N0.find(node) != mapO2N0.end())
                {
                    const auto &new_node = mapO2N0.at(node);
                    if ((d == 2) && ((cellOld2NewNode.find(f2c_p[f2cIdx_p[faceMf] + 1]) == cellOld2NewNode.end()) ||
                                     (cellOld2NewNode.at(f2c_p[f2cIdx_p[faceMf] + 1]).find(node) ==
                                      cellOld2NewNode.at(f2c_p[f2cIdx_p[faceMf] + 1]).end()) ||
                                     (cellOld2NewNode.at(f2c_p[f2cIdx_p[faceMf] + 1]).at(node) != new_node)))
                        throw INTERP_KERNEL::Exception(
                            "There is an"
                            "incoherent change of connectivity between both"
                            "cell surrounding a face touching the crack but"
                            "not duplicated."
                        );
                    f2n_pw[i_node] = new_node;
                }
            }
        }
    }
}
