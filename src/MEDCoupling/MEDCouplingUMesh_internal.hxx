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
// Author : Anthony Geay (EdF)

using namespace MEDCoupling;

class MinusOneSonsGenerator
{
   public:
    MinusOneSonsGenerator(const INTERP_KERNEL::CellModel &cm) : _cm(cm) {}
    unsigned getNumberOfSons2(const mcIdType *conn, mcIdType lgth) const { return _cm.getNumberOfSons2(conn, lgth); }
    unsigned fillSonCellNodalConnectivity2(
        int sonId,
        const mcIdType *nodalConn,
        mcIdType lgth,
        mcIdType *sonNodalConn,
        INTERP_KERNEL::NormalizedCellType &typeOfSon
    ) const
    {
        return _cm.fillSonCellNodalConnectivity2(sonId, nodalConn, lgth, sonNodalConn, typeOfSon);
    }
    static const int DELTA = 1;

   private:
    const INTERP_KERNEL::CellModel &_cm;
};

class MinusOneSonsGeneratorBiQuadratic
{
   public:
    MinusOneSonsGeneratorBiQuadratic(const INTERP_KERNEL::CellModel &cm) : _cm(cm) {}
    unsigned getNumberOfSons2(const mcIdType *conn, mcIdType lgth) const { return _cm.getNumberOfSons2(conn, lgth); }
    unsigned fillSonCellNodalConnectivity2(
        int sonId,
        const mcIdType *nodalConn,
        mcIdType lgth,
        mcIdType *sonNodalConn,
        INTERP_KERNEL::NormalizedCellType &typeOfSon
    ) const
    {
        return _cm.fillSonCellNodalConnectivity4(sonId, nodalConn, lgth, sonNodalConn, typeOfSon);
    }
    static const int DELTA = 1;

   private:
    const INTERP_KERNEL::CellModel &_cm;
};

class MinusTwoSonsGenerator
{
   public:
    MinusTwoSonsGenerator(const INTERP_KERNEL::CellModel &cm) : _cm(cm) {}
    unsigned getNumberOfSons2(const mcIdType *conn, mcIdType lgth) const
    {
        return _cm.getNumberOfEdgesIn3D(conn, lgth);
    }
    unsigned fillSonCellNodalConnectivity2(
        int sonId,
        const mcIdType *nodalConn,
        mcIdType lgth,
        mcIdType *sonNodalConn,
        INTERP_KERNEL::NormalizedCellType &typeOfSon
    ) const
    {
        return _cm.fillSonEdgesNodalConnectivity3D(sonId, nodalConn, lgth, sonNodalConn, typeOfSon);
    }
    static const int DELTA = 2;

   private:
    const INTERP_KERNEL::CellModel &_cm;
};

class MicroEdgesGenerator2D
{
   public:
    MicroEdgesGenerator2D(const INTERP_KERNEL::CellModel &cm) : _cm(cm) {}
    unsigned getNumberOfSons2(const mcIdType *conn, mcIdType lgth) const { return _cm.getNumberOfMicroEdges(); }
    unsigned fillSonCellNodalConnectivity2(
        int sonId,
        const mcIdType *nodalConn,
        mcIdType lgth,
        mcIdType *sonNodalConn,
        INTERP_KERNEL::NormalizedCellType &typeOfSon
    ) const
    {
        return _cm.fillMicroEdgeNodalConnectivity(sonId, nodalConn, sonNodalConn, typeOfSon);
    }
    static const int DELTA = 1;

   private:
    const INTERP_KERNEL::CellModel &_cm;
};

class MicroEdgesGenerator3D
{
   public:
    MicroEdgesGenerator3D(const INTERP_KERNEL::CellModel &cm) : _cm(cm) {}
    unsigned getNumberOfSons2(const mcIdType *conn, mcIdType lgth) const { return _cm.getNumberOfMicroEdges(); }
    unsigned fillSonCellNodalConnectivity2(
        int sonId,
        const mcIdType *nodalConn,
        mcIdType lgth,
        mcIdType *sonNodalConn,
        INTERP_KERNEL::NormalizedCellType &typeOfSon
    ) const
    {
        return _cm.fillMicroEdgeNodalConnectivity(sonId, nodalConn, sonNodalConn, typeOfSon);
    }
    static const int DELTA = 2;

   private:
    const INTERP_KERNEL::CellModel &_cm;
};

mcIdType
MEDCouplingFastNbrer(
    mcIdType id,
    mcIdType nb,
    const INTERP_KERNEL::CellModel &cm,
    bool compute,
    const mcIdType *conn1,
    const mcIdType *conn2
);
mcIdType
MEDCouplingOrientationSensitiveNbrer(
    mcIdType id,
    mcIdType nb,
    const INTERP_KERNEL::CellModel &cm,
    bool compute,
    const mcIdType *conn1,
    const mcIdType *conn2
);

namespace MEDCoupling
{
template <const int SPACEDIMM>
class DummyClsMCUG
{
   public:
    static const int MY_SPACEDIM = SPACEDIMM;
    static const int MY_MESHDIM = 8;
    typedef mcIdType MyConnType;
    static const INTERP_KERNEL::NumberingPolicy My_numPol = INTERP_KERNEL::ALL_C_MODE;
    // begin
    // useless, but for windows compilation ...
    const double *getCoordinatesPtr() const { return 0; }
    const MyConnType *getConnectivityPtr() const { return 0; }
    const MyConnType *getConnectivityIndexPtr() const { return 0; }
    INTERP_KERNEL::NormalizedCellType getTypeOfElement(MyConnType) const
    {
        return (INTERP_KERNEL::NormalizedCellType)0;
    }
    // end
};
}  // namespace MEDCoupling

template <int SPACEDIM>
void
MEDCouplingUMesh::getCellsContainingPointsAlg(
    const double *coords,
    const double *pos,
    mcIdType nbOfPoints,
    double eps,
    MCAuto<DataArrayIdType> &elts,
    MCAuto<DataArrayIdType> &eltsIndex,
    std::function<bool(INTERP_KERNEL::NormalizedCellType, int)> sensibilityTo2DQuadraticLinearCellsFunc
) const
{
    // Override precision for this method only:
    INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);

    elts = DataArrayIdType::New();
    eltsIndex = DataArrayIdType::New();
    eltsIndex->alloc(nbOfPoints + 1, 1);
    eltsIndex->setIJ(0, 0, 0);
    elts->alloc(0, 1);
    mcIdType *eltsIndexPtr(eltsIndex->getPointer());
    MCAuto<DataArrayDouble> bboxArr(getBoundingBoxForBBTree(eps));
    const double *bbox(bboxArr->begin());
    mcIdType nbOfCells = getNumberOfCells();
    const mcIdType *conn = _nodal_connec->getConstPointer();
    const mcIdType *connI = _nodal_connec_index->getConstPointer();
    double bb[2 * SPACEDIM];
    BBTree<SPACEDIM, mcIdType> myTree(&bbox[0], 0, 0, nbOfCells, -eps);
    for (mcIdType i = 0; i < nbOfPoints; i++)
    {
        eltsIndexPtr[i + 1] = eltsIndexPtr[i];
        for (int j = 0; j < SPACEDIM; j++)
        {
            bb[2 * j] = pos[SPACEDIM * i + j];
            bb[2 * j + 1] = pos[SPACEDIM * i + j];
        }
        std::vector<mcIdType> candidates;
        myTree.getIntersectingElems(bb, candidates);
        for (std::vector<mcIdType>::const_iterator iter = candidates.begin(); iter != candidates.end(); iter++)
        {
            mcIdType sz(connI[(*iter) + 1] - connI[*iter] - 1);
            INTERP_KERNEL::NormalizedCellType ct((INTERP_KERNEL::NormalizedCellType)conn[connI[*iter]]);
            bool status(false);
            // [ABN] : point locator algorithms are not impl. for POLY or QPOLY in spaceDim3
            if (SPACEDIM != 2 &&
                (ct == INTERP_KERNEL::NORM_POLYGON || sensibilityTo2DQuadraticLinearCellsFunc(ct, _mesh_dim)))
                throw INTERP_KERNEL::Exception(
                    "MEDCouplingUMesh::getCellsContainingPointsAlg : not implemented yet for POLYGON and QPOLYGON in "
                    "spaceDim 3 !"
                );
            // Keep calling simple algorithm when this is desired and simple for speed reasons:
            if (SPACEDIM == 2 && ct != INTERP_KERNEL::NORM_POLYGON &&
                !sensibilityTo2DQuadraticLinearCellsFunc(ct, _mesh_dim))
                status = INTERP_KERNEL::PointLocatorAlgos<DummyClsMCUG<2> >::isElementContainsPointAlgo2DSimple2(
                    pos + i * SPACEDIM, ct, coords, conn + connI[*iter] + 1, sz, eps
                );
            else
                status = INTERP_KERNEL::PointLocatorAlgos<DummyClsMCUG<SPACEDIM> >::isElementContainsPoint(
                    pos + i * SPACEDIM, ct, coords, conn + connI[*iter] + 1, sz, eps
                );
            if (status)
            {
                eltsIndexPtr[i + 1]++;
                elts->pushBackSilent(*iter);
            }
        }
    }
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 */
template <class SonsGenerator>
MEDCouplingUMesh *
MEDCouplingUMesh::buildDescendingConnectivityGen(
    DataArrayIdType *desc,
    DataArrayIdType *descIndx,
    DataArrayIdType *revDesc,
    DataArrayIdType *revDescIndx,
    DimM1DescNbrer nbrer
) const
{
    if (!desc || !descIndx || !revDesc || !revDescIndx)
        throw INTERP_KERNEL::Exception(
            "MEDCouplingUMesh::buildDescendingConnectivityGen : present of a null pointer in input !"
        );
    checkConnectivityFullyDefined();
    mcIdType nbOfCells = getNumberOfCells();
    mcIdType nbOfNodes = getNumberOfNodes();
    MCAuto<DataArrayIdType> revNodalIndx = DataArrayIdType::New();
    revNodalIndx->alloc(nbOfNodes + 1, 1);
    revNodalIndx->fillWithZero();
    mcIdType *revNodalIndxPtr = revNodalIndx->getPointer();
    const mcIdType *conn = _nodal_connec->getConstPointer();
    const mcIdType *connIndex = _nodal_connec_index->getConstPointer();
    std::string name = "Mesh constituent of ";
    name += getName();
    MCAuto<MEDCouplingUMesh> ret = MEDCouplingUMesh::New(name, getMeshDimension() - SonsGenerator::DELTA);
    ret->setCoords(getCoords());
    ret->allocateCells(2 * nbOfCells);
    descIndx->alloc(nbOfCells + 1, 1);
    MCAuto<DataArrayIdType> revDesc2(DataArrayIdType::New());
    revDesc2->reserve(2 * nbOfCells);
    mcIdType *descIndxPtr = descIndx->getPointer();
    *descIndxPtr++ = 0;
    for (mcIdType eltId = 0; eltId < nbOfCells; eltId++, descIndxPtr++)
    {
        mcIdType pos = connIndex[eltId];
        mcIdType posP1 = connIndex[eltId + 1];
        const INTERP_KERNEL::CellModel &cm =
            INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[pos]);
        SonsGenerator sg(cm);
        unsigned nbOfSons = sg.getNumberOfSons2(conn + pos + 1, posP1 - pos - 1);
        INTERP_KERNEL::AutoPtr<mcIdType> tmp = new mcIdType[posP1 - pos];
        for (unsigned i = 0; i < nbOfSons; i++)
        {
            INTERP_KERNEL::NormalizedCellType cmsId;
            unsigned nbOfNodesSon = sg.fillSonCellNodalConnectivity2(i, conn + pos + 1, posP1 - pos - 1, tmp, cmsId);
            for (unsigned k = 0; k < nbOfNodesSon; k++)
                if (tmp[k] >= 0)
                    revNodalIndxPtr[tmp[k] + 1]++;
            ret->insertNextCell(cmsId, nbOfNodesSon, tmp);
            revDesc2->pushBackSilent(eltId);
        }
        descIndxPtr[0] = descIndxPtr[-1] + ToIdType(nbOfSons);
    }
    mcIdType nbOfCellsM1 = ret->getNumberOfCells();
    std::transform(
        revNodalIndxPtr + 1,
        revNodalIndxPtr + nbOfNodes + 1,
        revNodalIndxPtr,
        revNodalIndxPtr + 1,
        std::plus<mcIdType>()
    );
    MCAuto<DataArrayIdType> revNodal = DataArrayIdType::New();
    revNodal->alloc(revNodalIndx->back(), 1);
    std::fill(revNodal->getPointer(), revNodal->getPointer() + revNodalIndx->back(), -1);
    mcIdType *revNodalPtr = revNodal->getPointer();
    const mcIdType *connM1 = ret->getNodalConnectivity()->getConstPointer();
    const mcIdType *connIndexM1 = ret->getNodalConnectivityIndex()->getConstPointer();
    for (mcIdType eltId = 0; eltId < nbOfCellsM1; eltId++)
    {
        const mcIdType *strtNdlConnOfCurCell = connM1 + connIndexM1[eltId] + 1;
        const mcIdType *endNdlConnOfCurCell = connM1 + connIndexM1[eltId + 1];
        for (const mcIdType *iter = strtNdlConnOfCurCell; iter != endNdlConnOfCurCell; iter++)
            if (*iter >= 0)  // for polyhedrons
                *std::find_if(
                    revNodalPtr + revNodalIndxPtr[*iter],
                    revNodalPtr + revNodalIndxPtr[*iter + 1],
                    std::bind(std::equal_to<mcIdType>(), std::placeholders::_1, -1)
                ) = eltId;
    }
    //
    DataArrayIdType *commonCells = 0, *commonCellsI = 0;
    FindCommonCellsAlg(
        3,
        0,
        ret->getNodalConnectivity(),
        ret->getNodalConnectivityIndex(),
        revNodal,
        revNodalIndx,
        commonCells,
        commonCellsI
    );
    MCAuto<DataArrayIdType> commonCellsTmp(commonCells), commonCellsITmp(commonCellsI);
    const mcIdType *commonCellsPtr(commonCells->getConstPointer()), *commonCellsIPtr(commonCellsI->getConstPointer());
    mcIdType newNbOfCellsM1 = -1;
    MCAuto<DataArrayIdType> o2nM1 = DataArrayIdType::ConvertIndexArrayToO2N(
        nbOfCellsM1, commonCells->begin(), commonCellsI->begin(), commonCellsI->end(), newNbOfCellsM1
    );
    std::vector<bool> isImpacted(nbOfCellsM1, false);
    for (const mcIdType *work = commonCellsI->begin(); work != commonCellsI->end() - 1; work++)
        for (mcIdType work2 = work[0]; work2 != work[1]; work2++) isImpacted[commonCellsPtr[work2]] = true;
    const mcIdType *o2nM1Ptr = o2nM1->getConstPointer();
    MCAuto<DataArrayIdType> n2oM1 = o2nM1->invertArrayO2N2N2OBis(newNbOfCellsM1);
    const mcIdType *n2oM1Ptr = n2oM1->getConstPointer();
    MCAuto<MEDCouplingUMesh> ret2 =
        static_cast<MEDCouplingUMesh *>(ret->buildPartOfMySelf(n2oM1->begin(), n2oM1->end(), true));
    ret2->copyTinyInfoFrom(this);
    desc->alloc(descIndx->back(), 1);
    mcIdType *descPtr = desc->getPointer();
    const INTERP_KERNEL::CellModel &cmsDft = INTERP_KERNEL::CellModel::GetCellModel(INTERP_KERNEL::NORM_POINT1);
    for (mcIdType i = 0; i < nbOfCellsM1; i++, descPtr++)
    {
        if (!isImpacted[i])
            *descPtr = nbrer(o2nM1Ptr[i], 0, cmsDft, false, 0, 0);
        else
        {
            if (i != n2oM1Ptr[o2nM1Ptr[i]])
            {
                const INTERP_KERNEL::CellModel &cms =
                    INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)connM1[connIndexM1[i]]);
                *descPtr = nbrer(
                    o2nM1Ptr[i],
                    connIndexM1[i + 1] - connIndexM1[i] - 1,
                    cms,
                    true,
                    connM1 + connIndexM1[n2oM1Ptr[o2nM1Ptr[i]]] + 1,
                    connM1 + connIndexM1[i] + 1
                );
            }
            else
                *descPtr = nbrer(o2nM1Ptr[i], 0, cmsDft, false, 0, 0);
        }
    }
    revDesc->reserve(newNbOfCellsM1);
    revDescIndx->alloc(newNbOfCellsM1 + 1, 1);
    mcIdType *revDescIndxPtr = revDescIndx->getPointer();
    *revDescIndxPtr++ = 0;
    const mcIdType *revDesc2Ptr = revDesc2->getConstPointer();
    for (mcIdType i = 0; i < newNbOfCellsM1; i++, revDescIndxPtr++)
    {
        mcIdType oldCellIdM1 = n2oM1Ptr[i];
        if (!isImpacted[oldCellIdM1])
        {
            revDesc->pushBackSilent(revDesc2Ptr[oldCellIdM1]);
            revDescIndxPtr[0] = revDescIndxPtr[-1] + 1;
        }
        else
        {
            for (mcIdType j = commonCellsIPtr[0]; j < commonCellsIPtr[1]; j++)
                revDesc->pushBackSilent(revDesc2Ptr[commonCellsPtr[j]]);
            revDescIndxPtr[0] = revDescIndxPtr[-1] + commonCellsIPtr[1] - commonCellsIPtr[0];
            commonCellsIPtr++;
        }
    }
    //
    return ret2.retn();
}
