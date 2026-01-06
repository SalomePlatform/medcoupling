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
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// Author : Anthony Geay (CEA/DEN)

#include "InterpKernelCellSimplify.hxx"
#include "CellModel.hxx"

#include <functional>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <numeric>
#include <cstring>
#include <limits>
#include <vector>
#include <list>
#include <set>

using namespace INTERP_KERNEL;

/*!
 * This method takes as input a cell with type 'type' and whose connectivity is defined by (conn,lgth)
 * It retrieves the same cell with a potentially different type (in return) whose connectivity is defined by
 * (retConn,retLgth)
 * \b WARNING for optimization reason the arrays 'retConn' and 'conn' can overlapped !
 */
INTERP_KERNEL::NormalizedCellType
CellSimplify::simplifyDegeneratedCell(
    INTERP_KERNEL::NormalizedCellType type, const mcIdType *conn, mcIdType lgth, mcIdType *retConn, mcIdType &retLgth
)
{
    const INTERP_KERNEL::CellModel &cm = INTERP_KERNEL::CellModel::GetCellModel(type);
    std::set<mcIdType> c(conn, conn + lgth);
    c.erase(-1);
    bool isObviousNonDegeneratedCell = (ToIdType(c.size()) == lgth);
    if ((cm.getDimension() == 3 && cm.isQuadratic()) || isObviousNonDegeneratedCell)
    {  // quadratic 3D, do nothing for the moment.
        retLgth = lgth;
        mcIdType *tmp = new mcIdType[lgth];  // no direct std::copy ! overlapping of conn and retConn !
        std::copy(conn, conn + lgth, tmp);
        std::copy(tmp, tmp + lgth, retConn);
        delete[] tmp;
        return type;
    }
    if (cm.getDimension() == 2)
    {
        mcIdType *tmp = new mcIdType[lgth];
        int newPos = 0;
        if (!cm.isQuadratic())
        {
            for (int i = 0; i < lgth; i++)
                if (conn[i] != conn[(i + 1) % lgth])  // zip nul segments/arcs
                    tmp[newPos++] = conn[i];
        }
        else
        {
            mcIdType quadOff = lgth / 2;
            mcIdType *tmpQuad = new mcIdType[quadOff];
            for (int i = 0; i < quadOff; i++)
                if (conn[i] != conn[(i + 1) % quadOff] ||
                    conn[i] != conn[i + quadOff])  // zip nul segments/arcs (quad point must match too)
                {
                    tmp[newPos] = conn[i];
                    tmpQuad[newPos++] = conn[(i + quadOff) % lgth];
                }
            // Merge linear and quad points into tmp
            std::copy(tmpQuad, tmpQuad + newPos, tmp + newPos);
            delete[] tmpQuad;
            newPos *= 2;  // take in quad points in the final length
        }
        INTERP_KERNEL::NormalizedCellType ret = tryToUnPoly2D(cm.isQuadratic(), tmp, newPos, retConn, retLgth);
        delete[] tmp;
        return ret;
    }
    if (cm.getDimension() == 3)
    {
        mcIdType nbOfFaces, lgthOfPolyhConn;
        mcIdType *zipFullReprOfPolyh = getFullPolyh3DCell(type, conn, lgth, nbOfFaces, lgthOfPolyhConn);
        INTERP_KERNEL::NormalizedCellType ret =
            tryToUnPoly3D(zipFullReprOfPolyh, nbOfFaces, lgthOfPolyhConn, retConn, retLgth);
        delete[] zipFullReprOfPolyh;
        return ret;
    }
    throw INTERP_KERNEL::Exception("CellSimplify::simplifyDegeneratedCell : works only with 2D and 3D cell !");
}

/*!
 * This static method tries to unpolygonize a cell whose connectivity is given by 'conn' and 'lgth'.
 * Contrary to INTERP_KERNEL::CellSimplify::simplifyDegeneratedCell method 'conn' and 'retConn' do not overlap.
 */
INTERP_KERNEL::NormalizedCellType
CellSimplify::tryToUnPoly2D(bool isQuad, const mcIdType *conn, mcIdType lgth, mcIdType *retConn, mcIdType &retLgth)
{
    retLgth = lgth;
    std::copy(conn, conn + lgth, retConn);
    if (!isQuad)
    {
        switch (lgth)
        {
            case 3:
                return INTERP_KERNEL::NORM_TRI3;
            case 4:
                return INTERP_KERNEL::NORM_QUAD4;
            default:
                return INTERP_KERNEL::NORM_POLYGON;
        }
    }
    else
    {
        switch (lgth)
        {
            case 6:
                return INTERP_KERNEL::NORM_TRI6;
            case 8:
                return INTERP_KERNEL::NORM_QUAD8;
            default:
                return INTERP_KERNEL::NORM_QPOLYG;
        }
    }
}

/*!
 * This method takes as input a 3D linear cell and put its representation in returned array. Warning the returned array
 * has to be deallocated. The length of the returned array is specified by out parameter The format of output array is
 * the following : 1,2,3,-1,3,4,2,-1,3,4,1,-1,1,2,4,NORM_TRI3,NORM_TRI3,NORM_TRI3 (faces type at the end of classical
 * polyhedron nodal description)
 */
mcIdType *
CellSimplify::getFullPolyh3DCell(
    INTERP_KERNEL::NormalizedCellType type,
    const mcIdType *conn,
    mcIdType lgth,
    mcIdType &retNbOfFaces,
    mcIdType &retLgth
)
{
    const INTERP_KERNEL::CellModel &cm = INTERP_KERNEL::CellModel::GetCellModel(type);
    unsigned nbOfFaces = cm.getNumberOfSons2(conn, lgth);
    mcIdType *tmp = new mcIdType[nbOfFaces * (lgth + 1)];
    mcIdType *work = tmp;
    std::vector<mcIdType> faces;
    for (unsigned j = 0; j < nbOfFaces; j++)
    {
        INTERP_KERNEL::NormalizedCellType type2;
        unsigned offset = cm.fillSonCellNodalConnectivity2(j, conn, lgth, work, type2);
        //
        mcIdType *tmp2 = new mcIdType[offset];
        tmp2[0] = work[0];
        mcIdType newPos = 1;
        for (unsigned k = 1; k < offset; k++)
            if (std::find(tmp2, tmp2 + newPos, work[k]) == tmp2 + newPos)
                tmp2[newPos++] = work[k];
        if (newPos < 3)
        {
            delete[] tmp2;
            continue;
        }
        mcIdType tmp3;
        faces.push_back(tryToUnPoly2D(CellModel::GetCellModel(type2).isQuadratic(), tmp2, newPos, work, tmp3));
        delete[] tmp2;
        //
        work += newPos;
        *work++ = -1;
    }
    std::copy(faces.begin(), faces.end(), --work);
    retNbOfFaces = (int)faces.size();
    retLgth = (int)std::distance(tmp, work);
    return tmp;
}

/*!
 * This static method tries to unpolygonize a cell whose connectivity is given by 'conn' (format is the same as
 * specified in method INTERP_KERNEL::CellSimplify::getFullPolyh3DCell ) and 'lgth'+'nbOfFaces'. Contrary to
 * INTERP_KERNEL::CellSimplify::simplifyDegeneratedCell method 'conn' and 'retConn' do not overlap.
 */
INTERP_KERNEL::NormalizedCellType
CellSimplify::tryToUnPoly3D(
    const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType &retLgth
)
{
    std::set<mcIdType> nodes(conn, conn + lgth);
    nodes.erase(-1);
    std::size_t nbOfNodes = nodes.size();
    std::size_t magicNumber = 100 * nbOfNodes + nbOfFaces;
    switch (magicNumber)
    {
        case 806:
            return tryToUnPolyHex8(conn, nbOfFaces, lgth, retConn, retLgth);
        case 1208:
            return tryToUnPolyHexp12(conn, nbOfFaces, lgth, retConn, retLgth);
        case 605:
            return tryToUnPolyPenta6(conn, nbOfFaces, lgth, retConn, retLgth);
        case 505:
            return tryToUnPolyPyra5(conn, nbOfFaces, lgth, retConn, retLgth);
        case 404:
            return tryToUnPolyTetra4(conn, nbOfFaces, lgth, retConn, retLgth);
        default:
            retLgth = lgth;
            std::copy(conn, conn + lgth, retConn);
            return INTERP_KERNEL::NORM_POLYHED;
    }
}

bool
CellSimplify::orientOppositeFace(
    const mcIdType *baseFace, mcIdType *retConn, const mcIdType *sideFace, mcIdType lgthBaseFace
)
{
    std::vector<mcIdType> tmp2;
    std::set<mcIdType> bases(baseFace, baseFace + lgthBaseFace);
    std::set<mcIdType> sides(sideFace, sideFace + 4);
    std::set_intersection(
        bases.begin(), bases.end(), sides.begin(), sides.end(), std::back_insert_iterator<std::vector<mcIdType> >(tmp2)
    );
    if (tmp2.size() != 2)
        return false;
    std::vector<std::pair<mcIdType, mcIdType> > baseEdges(lgthBaseFace);
    std::vector<std::pair<mcIdType, mcIdType> > oppEdges(lgthBaseFace);
    std::vector<std::pair<mcIdType, mcIdType> > sideEdges(4);
    for (mcIdType i = 0; i < lgthBaseFace; i++)
    {
        baseEdges[i] = std::pair<mcIdType, mcIdType>(baseFace[i], baseFace[(i + 1) % lgthBaseFace]);
        oppEdges[i] = std::pair<mcIdType, mcIdType>(retConn[i], retConn[(i + 1) % lgthBaseFace]);
    }
    for (int i = 0; i < 4; i++) sideEdges[i] = std::pair<mcIdType, mcIdType>(sideFace[i], sideFace[(i + 1) % 4]);
    std::vector<std::pair<mcIdType, mcIdType> > tmp;
    std::set<std::pair<mcIdType, mcIdType> > baseEdgesS(baseEdges.begin(), baseEdges.end());
    std::set<std::pair<mcIdType, mcIdType> > sideEdgesS(sideEdges.begin(), sideEdges.end());
    std::set_intersection(
        baseEdgesS.begin(),
        baseEdgesS.end(),
        sideEdgesS.begin(),
        sideEdgesS.end(),
        std::back_insert_iterator<std::vector<std::pair<mcIdType, mcIdType> > >(tmp)
    );
    if (tmp.empty())
    {
        // reverse sideFace
        for (int i = 0; i < 4; i++)
        {
            std::pair<mcIdType, mcIdType> p = sideEdges[i];
            std::pair<mcIdType, mcIdType> r(p.second, p.first);
            sideEdges[i] = r;
        }
        // end reverse sideFace
        std::set<std::pair<mcIdType, mcIdType> > baseEdgesS2(baseEdges.begin(), baseEdges.end());
        std::set<std::pair<mcIdType, mcIdType> > sideEdgesS2(sideEdges.begin(), sideEdges.end());
        std::set_intersection(
            baseEdgesS2.begin(),
            baseEdgesS2.end(),
            sideEdgesS2.begin(),
            sideEdgesS2.end(),
            std::back_insert_iterator<std::vector<std::pair<mcIdType, mcIdType> > >(tmp)
        );
        if (tmp.empty())
            return false;
    }
    if (tmp.size() != 1)
        return false;
    bool found = false;
    std::pair<mcIdType, mcIdType> pInOpp;
    for (int i = 0; i < 4 && !found; i++)
    {  // finding the pair(edge) in sideFace that do not include any node of tmp[0] edge
        found =
            (tmp[0].first != sideEdges[i].first && tmp[0].first != sideEdges[i].second &&
             tmp[0].second != sideEdges[i].first && tmp[0].second != sideEdges[i].second);
        if (found)
        {  // found ! reverse it
            pInOpp.first = sideEdges[i].second;
            pInOpp.second = sideEdges[i].first;
        }
    }
    if (!found)
        return false;
    int pos = (int)std::distance(baseEdges.begin(), std::find(baseEdges.begin(), baseEdges.end(), tmp[0]));
    std::vector<std::pair<mcIdType, mcIdType> >::iterator it = std::find(oppEdges.begin(), oppEdges.end(), pInOpp);
    if (it == oppEdges.end())  // the opposite edge of side face is not found opposite face ... maybe problem of
                               // orientation of polyhedron
        return false;
    mcIdType pos2 = ToIdType(std::distance(oppEdges.begin(), it));
    mcIdType offset = pos - pos2;
    if (offset < 0)
        offset += lgthBaseFace;
    // this is the end copy the result
    mcIdType *tmp3 = new mcIdType[lgthBaseFace];
    for (int i = 0; i < lgthBaseFace; i++) tmp3[(offset + i) % lgthBaseFace] = oppEdges[i].first;
    std::copy(tmp3, tmp3 + lgthBaseFace, retConn);
    delete[] tmp3;
    return true;
}

bool
CellSimplify::isWellOriented(
    const mcIdType *baseFace, mcIdType *retConn, const mcIdType *sideFace, mcIdType lgthBaseFace
)
{
    return true;
}

/*!
 * This method is trying to permute the connectivity of 'oppFace' face so that the k_th node of 'baseFace' is associated
 * to the k_th node in retConnOfOppFace. Excluded faces 'baseFace' and 'oppFace' all the other faces in 'conn' must be
 * QUAD4 faces. If the arrangement process succeeds true is returned and retConnOfOppFace is filled.
 */
bool
CellSimplify::tryToArrangeOppositeFace(
    const mcIdType *conn,
    mcIdType lgth,
    mcIdType lgthBaseFace,
    const mcIdType *baseFace,
    const mcIdType *oppFace,
    mcIdType nbOfFaces,
    mcIdType *retConnOfOppFace
)
{
    retConnOfOppFace[0] = oppFace[0];
    for (mcIdType j = 1; j < lgthBaseFace; j++) retConnOfOppFace[j] = oppFace[lgthBaseFace - j];
    const mcIdType *curFace = conn;
    int sideFace = 0;
    bool ret = true;
    for (int i = 0; i < nbOfFaces && ret; i++)
    {
        if (curFace != baseFace && curFace != oppFace)
        {
            if (sideFace == 0)
                ret = orientOppositeFace(baseFace, retConnOfOppFace, curFace, lgthBaseFace);
            else
                ret = isWellOriented(baseFace, retConnOfOppFace, curFace, lgthBaseFace);
            sideFace++;
        }
        curFace = std::find(curFace, conn + lgth, -1);
        curFace++;
    }
    return ret;
}

/*!
 * Cell with 'conn' connectivity has been detected as a good candidate. Full check of this. If yes NORM_HEXA8 is
 * returned. This method is only callable if in 'conn' there is 8 nodes and 6 faces. If fails a POLYHED is returned.
 */
INTERP_KERNEL::NormalizedCellType
CellSimplify::tryToUnPolyHex8(
    const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType &retLgth
)
{
    if (std::find_if(
            conn + lgth,
            conn + lgth + nbOfFaces,
            std::bind(std::not_equal_to<mcIdType>(), std::placeholders::_1, ToIdType(INTERP_KERNEL::NORM_QUAD4))
        ) == conn + lgth + nbOfFaces)
    {  // 6 faces are QUAD4.
        int oppositeFace = -1;
        std::set<mcIdType> conn1(conn, conn + 4);
        for (int i = 1; i < 6 && oppositeFace < 0; i++)
        {
            std::vector<mcIdType> tmp;
            std::set<mcIdType> conn2(conn + 5 * i, conn + 5 * i + 4);
            std::set_intersection(
                conn1.begin(),
                conn1.end(),
                conn2.begin(),
                conn2.end(),
                std::back_insert_iterator<std::vector<mcIdType> >(tmp)
            );
            if (tmp.empty())
                oppositeFace = i;
        }
        if (oppositeFace >= 1)
        {  // oppositeFace of face#0 found.
            mcIdType tmp2[4];
            if (tryToArrangeOppositeFace(conn, lgth, 4, conn, conn + 5 * oppositeFace, 6, tmp2))
            {
                std::copy(conn, conn + 4, retConn);
                std::copy(tmp2, tmp2 + 4, retConn + 4);
                retLgth = 8;
                return INTERP_KERNEL::NORM_HEXA8;
            }
        }
    }
    retLgth = lgth;
    std::copy(conn, conn + lgth, retConn);
    return INTERP_KERNEL::NORM_POLYHED;
}

INTERP_KERNEL::NormalizedCellType
CellSimplify::tryToUnPolyHexp12(
    const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType &retLgth
)
{
    std::size_t nbOfHexagon = std::count(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_POLYGON));
    std::size_t nbOfQuad = std::count(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_QUAD4));
    if (nbOfQuad == 6 && nbOfHexagon == 2)
    {
        const mcIdType *hexag0 = std::find(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_POLYGON));
        std::size_t hexg0Id = std::distance(conn + lgth, hexag0);
        const mcIdType *hexag1 = std::find(hexag0 + 1, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_POLYGON));
        std::size_t hexg1Id = std::distance(conn + lgth, hexag1);
        const mcIdType *connHexag0 = conn + 5 * hexg0Id;
        std::size_t lgthH0 = std::distance(connHexag0, std::find(connHexag0, conn + lgth, -1));
        if (lgthH0 == 6)
        {
            const mcIdType *connHexag1 = conn + 5 * hexg0Id + 7 + (hexg1Id - hexg0Id - 1) * 5;
            std::size_t lgthH1 = std::distance(connHexag1, std::find(connHexag1, conn + lgth, -1));
            if (lgthH1 == 6)
            {
                std::vector<mcIdType> tmp;
                std::set<mcIdType> conn1(connHexag0, connHexag0 + 6);
                std::set<mcIdType> conn2(connHexag1, connHexag1 + 6);
                std::set_intersection(
                    conn1.begin(),
                    conn1.end(),
                    conn2.begin(),
                    conn2.end(),
                    std::back_insert_iterator<std::vector<mcIdType> >(tmp)
                );
                if (tmp.empty())
                {
                    mcIdType tmp2[6];
                    if (tryToArrangeOppositeFace(conn, lgth, 6, connHexag0, connHexag1, 8, tmp2))
                    {
                        std::copy(connHexag0, connHexag0 + 6, retConn);
                        std::copy(tmp2, tmp2 + 6, retConn + 6);
                        retLgth = 12;
                        return INTERP_KERNEL::NORM_HEXGP12;
                    }
                }
            }
        }
    }
    retLgth = lgth;
    std::copy(conn, conn + lgth, retConn);
    return INTERP_KERNEL::NORM_POLYHED;
}

/*!
 * Cell with 'conn' connectivity has been detected as a good candidate. Full check of this. If yes NORM_PENTA6 is
 * returned. If fails a POLYHED is returned.
 */
INTERP_KERNEL::NormalizedCellType
CellSimplify::tryToUnPolyPenta6(
    const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType &retLgth
)
{
    std::size_t nbOfTriFace = std::count(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_TRI3));
    std::size_t nbOfQuadFace = std::count(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_QUAD4));
    if (nbOfTriFace == 2 && nbOfQuadFace == 3)
    {
        std::size_t tri3_0 = std::distance(
            conn + lgth, std::find(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_TRI3))
        );
        std::size_t tri3_1 = std::distance(
            conn + lgth,
            std::find(conn + lgth + tri3_0 + 1, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_TRI3))
        );
        const mcIdType *tri_0 = 0, *tri_1 = 0;
        const mcIdType *w = conn;
        for (std::size_t i = 0; i < 5; i++)
        {
            if (i == tri3_0)
                tri_0 = w;
            if (i == tri3_1)
                tri_1 = w;
            w = std::find(w, conn + lgth, -1);
            w++;
        }
        std::vector<mcIdType> tmp;
        std::set<mcIdType> conn1(tri_0, tri_0 + 3);
        std::set<mcIdType> conn2(tri_1, tri_1 + 3);
        std::set_intersection(
            conn1.begin(),
            conn1.end(),
            conn2.begin(),
            conn2.end(),
            std::back_insert_iterator<std::vector<mcIdType> >(tmp)
        );
        if (tmp.empty())
        {
            mcIdType tmp2[3];
            if (tryToArrangeOppositeFace(conn, lgth, 3, tri_0, tri_1, 5, tmp2))
            {
                std::copy(tri_0, tri_0 + 3, retConn);
                std::copy(tmp2, tmp2 + 3, retConn + 3);
                retLgth = 6;
                return INTERP_KERNEL::NORM_PENTA6;
            }
        }
    }
    retLgth = lgth;
    std::copy(conn, conn + lgth, retConn);
    return INTERP_KERNEL::NORM_POLYHED;
}

/*!
 * Cell with 'conn' connectivity has been detected as a good candidate. Full check of this. If yes NORM_PYRA5 is
 * returned. If fails a POLYHED is returned.
 */
INTERP_KERNEL::NormalizedCellType
CellSimplify::tryToUnPolyPyra5(
    const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType &retLgth
)
{
    std::size_t nbOfTriFace = std::count(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_TRI3));
    std::size_t nbOfQuadFace = std::count(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_QUAD4));
    if (nbOfTriFace == 4 && nbOfQuadFace == 1)
    {
        std::size_t quad4_pos = std::distance(
            conn + lgth, std::find(conn + lgth, conn + lgth + nbOfFaces, ToIdType(INTERP_KERNEL::NORM_QUAD4))
        );
        const mcIdType *quad4 = 0;
        const mcIdType *w = conn;
        for (std::size_t i = 0; i < 5 && quad4 == 0; i++)
        {
            if (i == quad4_pos)
                quad4 = w;
            w = std::find(w, conn + lgth, -1);
            w++;
        }
        std::set<mcIdType> quad4S(quad4, quad4 + 4);
        w = conn;
        bool ok = true;
        mcIdType point = -1;
        for (std::size_t i = 0; i < 5 && ok; i++)
        {
            if (i != quad4_pos)
            {
                std::vector<mcIdType> tmp;
                std::set<mcIdType> conn2(w, w + 3);
                std::set_intersection(
                    conn2.begin(),
                    conn2.end(),
                    quad4S.begin(),
                    quad4S.end(),
                    std::back_insert_iterator<std::vector<mcIdType> >(tmp)
                );
                ok = tmp.size() == 2;
                tmp.clear();
                std::set_difference(
                    conn2.begin(),
                    conn2.end(),
                    quad4S.begin(),
                    quad4S.end(),
                    std::back_insert_iterator<std::vector<mcIdType> >(tmp)
                );
                ok = ok && tmp.size() == 1;
                if (ok)
                {
                    if (point >= 0)
                        ok = point == tmp[0];
                    else
                        point = tmp[0];
                }
            }
            w = std::find(w, conn + lgth, -1);
            w++;
        }
        if (ok && point >= 0)
        {
            std::copy(quad4, quad4 + 4, retConn);
            retConn[4] = point;
            retLgth = 5;
            return INTERP_KERNEL::NORM_PYRA5;
        }
    }
    retLgth = lgth;
    std::copy(conn, conn + lgth, retConn);
    return INTERP_KERNEL::NORM_POLYHED;
}

/*!
 * Cell with 'conn' connectivity has been detected as a good candidate. Full check of this. If yes NORM_TETRA4 is
 * returned. If fails a POLYHED is returned.
 */
INTERP_KERNEL::NormalizedCellType
CellSimplify::tryToUnPolyTetra4(
    const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType &retLgth
)
{
    if (std::find_if(
            conn + lgth,
            conn + lgth + nbOfFaces,
            std::bind(std::not_equal_to<mcIdType>(), std::placeholders::_1, ToIdType(INTERP_KERNEL::NORM_TRI3))
        ) == conn + lgth + nbOfFaces)
    {
        std::set<mcIdType> tribase(conn, conn + 3);
        mcIdType point = -1;
        bool ok = true;
        for (int i = 1; i < 4 && ok; i++)
        {
            std::vector<mcIdType> tmp;
            std::set<mcIdType> conn2(conn + i * 4, conn + 4 * i + 3);
            std::set_intersection(
                conn2.begin(),
                conn2.end(),
                tribase.begin(),
                tribase.end(),
                std::back_insert_iterator<std::vector<mcIdType> >(tmp)
            );
            ok = tmp.size() == 2;
            tmp.clear();
            std::set_difference(
                conn2.begin(),
                conn2.end(),
                tribase.begin(),
                tribase.end(),
                std::back_insert_iterator<std::vector<mcIdType> >(tmp)
            );
            ok = ok && tmp.size() == 1;
            if (ok)
            {
                if (point >= 0)
                    ok = point == tmp[0];
                else
                    point = tmp[0];
            }
        }
        if (ok && point >= 0)
        {
            std::copy(conn, conn + 3, retConn);
            retConn[3] = point;
            retLgth = 4;
            return INTERP_KERNEL::NORM_TETRA4;
        }
    }
    retLgth = lgth;
    std::copy(conn, conn + lgth, retConn);
    return INTERP_KERNEL::NORM_POLYHED;
}

/*!
 * Tell whether a cell is exactly flat.
 * For the moment only handle:
 *  - fully degenerated polygons (polygon with 1 point, or 2 if quadratic)
 *  - quad polygon with 2 points and two identical quad points
 */
bool
CellSimplify::isFlatCell(const mcIdType *conn, mcIdType pos, mcIdType lgth, NormalizedCellType type)
{
    const INTERP_KERNEL::CellModel &cm = INTERP_KERNEL::CellModel::GetCellModel(type);
    if (lgth <=
        2)  // a polygon with a single, or two points has been returned. This check also captures degenerated quadratics
        return true;
    if (cm.isQuadratic() && lgth == 4)  // test for flat quadratic polygon with 2 edges ...
        if (conn[pos + 1 + lgth / 2] == conn[pos + 1 + lgth / 2 + 1])  // the only 2 quad points are equal
            return true;
    return false;
}
