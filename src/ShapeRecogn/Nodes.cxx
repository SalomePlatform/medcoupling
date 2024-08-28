// Copyright (C) 2024  CEA, EDF
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

#include "Nodes.hxx"

using namespace MEDCoupling;

Nodes::Nodes(
    const MEDCouplingUMesh *mesh,
    const DataArrayInt64 *neighbors,
    const DataArrayInt64 *neighborsIdx)
    : mesh(mesh), coords(mesh->getCoords()), neighbors(neighbors), neighborsIdx(neighborsIdx)
{
}

mcIdType Nodes::getNbNodes() const
{
    return coords->getNumberOfTuples();
}

const std::vector<double> &Nodes::getK1() const
{
    return k1;
}

const std::vector<double> &Nodes::getK2() const
{
    return k2;
}

const std::vector<PrimitiveType> &Nodes::getPrimitiveType() const
{
    return primitives;
}

const std::vector<double> &Nodes::getNormals() const
{
    return normals;
}

const std::vector<double> &Nodes::getWeakDirections() const
{
    return weakDirections;
}

const std::vector<double> &Nodes::getMainDirections() const
{
    return mainDirections;
}

std::array<double, 3> Nodes::getNormal(mcIdType nodeId) const
{
    return {normals[3 * nodeId], normals[3 * nodeId + 1], normals[3 * nodeId + 2]};
}

double Nodes::getK1(mcIdType nodeId) const
{
    return k1[nodeId];
}

double Nodes::getK2(mcIdType nodeId) const
{
    return k2[nodeId];
}

double Nodes::getKdiff0(mcIdType nodeId) const
{
    return fabs(k1[nodeId]) > fabs(k2[nodeId]) ? k1[nodeId] : k2[nodeId];
}

double Nodes::getAdimK1(mcIdType nodeId) const
{
    return adimK1[nodeId];
}

double Nodes::getAdimK2(mcIdType nodeId) const
{
    return adimK2[nodeId];
}

double Nodes::getAdimKdiff0(mcIdType nodeId) const
{
    return fabs(adimK1[nodeId]) > fabs(adimK2[nodeId]) ? adimK1[nodeId] : adimK2[nodeId];
}

const std::vector<mcIdType> Nodes::getNeighbors(mcIdType nodeId) const
{
    mcIdType start = neighborsIdx->getIJ(nodeId, 0);
    mcIdType nbNeighbors = neighborsIdx->getIJ(nodeId + 1, 0) - start;

    std::vector<mcIdType> neighborNodes(nbNeighbors, 0);
    for (mcIdType i = 0; i < nbNeighbors; i++)
        neighborNodes[i] = neighbors->getIJ(start + i, 0);
    return neighborNodes;
}

PrimitiveType Nodes::getPrimitiveType(mcIdType nodeId) const
{
    return primitives[nodeId];
}

std::array<double, 3> Nodes::getWeakDirection(mcIdType nodeId) const
{
    return {weakDirections[3 * nodeId],
            weakDirections[3 * nodeId + 1],
            weakDirections[3 * nodeId + 2]};
}

std::array<double, 3> Nodes::getMainDirection(mcIdType nodeId) const
{
    return {mainDirections[3 * nodeId],
            mainDirections[3 * nodeId + 1],
            mainDirections[3 * nodeId + 2]};
}

std::array<double, 3> Nodes::getCoordinates(mcIdType nodeId) const
{
    std::array<double, 3> nodeCoords;
    coords->getTuple(nodeId, nodeCoords.data());
    return nodeCoords;
}
