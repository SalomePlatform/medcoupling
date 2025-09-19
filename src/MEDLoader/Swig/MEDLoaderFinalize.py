#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2023-2025  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

import logging


def MEDFileUMeshFuseNodesAndCells(
    self, compType=2, eps=1e-6, logLev=logging.INFO, infoWrapNodes=None
):
    """
    Method fusing nodes in this, then fusing cells in this. Fusion is done following eps and compType.

    :param compType : see MEDCouplingPointSet.zipConnectivityTraducer method for explanations
    :param eps: see DataArrayDouble.findCommonTuples for explanations.
    :param logLev: Integer specifying log level
    :param n2oHolder: Optional output param storing for each level ext n2o conversions applied during transformation. The storage should follow pydict concept. ket is levelext. Value is n2o array associated.
    :param infoWrapNodes: Optional input. If precised contains tuple of size 4 ( cNode, ciNodes, o2nNodes, n2oNodes ) fully defining the merge of nodes.
    :return: MEDFileUMesh instance containing the result of nodes and cells fusion
    """
    mmOut, _ = MEDFileUMeshFuseNodesAndCellsAdv(
        self, compType, eps, logLev, infoWrapNodes
    )
    return mmOut


def getLogger(level=logging.INFO):
    FORMAT = "%(levelname)s : %(asctime)s : [%(filename)s:%(funcName)s:%(lineno)s] : %(message)s"
    logging.basicConfig(format=FORMAT, level=level)
    return logging.getLogger()


def MEDFileUMeshFuseNodesAndCellsAdv(
    self, compType=2, eps=1e-6, logLev=logging.INFO, infoWrapNodes=None
):
    """
    Same than MEDFileUMeshfuseNodesAndCells except that

    WARNING !!!! self is expected to follow MED file convention. Positive fam values for nodes, negative for cells.

    :return: a tuple a size 2. First element is MEDFileUMesh instance containing the result of nodes and cells fusion, 2nd element is dict storing for each level ext n2o conversions applied during transformation. key is levelext. Value is n2o array associated.
    """
    import MEDLoader as ml

    def updateMap(
        mm: ml.MEDFileUMesh,
        lev: int,
        famMap: ml.DataArrayInt,
        famMapI: ml.DataArrayInt,
        famIdOffset: int,
    ):
        """
        mm instance to be updated
        """

        def famIdManager(lev, famId):
            if lev <= 0:
                return -famId
            else:
                return famId

        nbOfPartSetToBeUpdated = len(famMapI) - 1
        for partSetId in range(nbOfPartSetToBeUpdated):
            newFamId = famIdManager(lev, famMap[famMapI[partSetId]] + famIdOffset)
            newFamName = f"Family_{newFamId}"
            logger.debug(f"For level {lev} new family : {newFamId}")
            mm.addFamily(newFamName, newFamId)
            for famId in famMap[famMapI[partSetId] + 1 : famMapI[partSetId + 1]]:
                zeFamId = famIdManager(lev, int(famId))
                if not mm.existsFamily(zeFamId):
                    continue
                grpsToBeUpdated = mm.getGroupsOnFamily(mm.getFamilyNameGivenId(zeFamId))
                for grpToBeUpdated in grpsToBeUpdated:
                    mm.addFamilyOnGrp(grpToBeUpdated, newFamName)
        pass

    famIdZeroForNewFamilies = self.getMaxAbsFamilyIdInArrays() + 1
    n2oHolder = {}
    logger = getLogger(level=logLev)
    initNbNodes = len(self.getCoords())
    if infoWrapNodes is None:
        logger.info(f"No n2onodes given. Trying to compute it using given eps = {eps}")
        logger.info(f"Begin merging nodes with eps = {eps}")
        cNode, ciNodes = self.getCoords().findCommonTuples(eps)
        logger.info(
            f"End of merging nodes with eps = {eps} : Nb of nodes groups to be merged : {len(ciNodes) - 1} / {self.getNumberOfNodes()}"
        )
        o2nNodes, newNbNodes = ml.DataArrayInt.ConvertIndexArrayToO2N(
            initNbNodes, cNode, ciNodes
        )
        n2oNodes = o2nNodes.invertArrayO2N2N2O(newNbNodes)
    else:
        cNode, ciNodes, o2nNodes, n2oNodes = infoWrapNodes
    newCoords = self.getCoords()[n2oNodes]
    n2oHolder[1] = n2oNodes
    # creation of
    mmOut = ml.MEDFileUMesh()
    mmOut.copyFamGrpMapsFrom(self)

    for lev in self.getNonEmptyLevels():
        logger.debug(f"Begin level {lev}")
        m1 = self[lev].deepCopy()
        logger.debug(f"Begin renumbering connectivity of level {lev}")
        m1.renumberNodesInConn(o2nNodes)
        logger.debug(f"End renumbering connectivity of level {lev}")
        m1.setCoords(newCoords)
        logger.info(f"Begin of finding of same cells of level {lev}")
        cce, ccei = m1.findCommonCells(compType, 0)
        logger.info(
            f"End of finding of same cells of level {lev} : Nb of cells groups to be merged : {len(ccei) - 1} / {m1.getNumberOfCells()}"
        )
        famsCell = self.getFamilyFieldAtLevel(lev)
        if famsCell:
            famsCell = -famsCell
            localFamRef = famsCell.getMaxAbsValueInArray() + 1
            famsMergedCell, famMap, famMapI = famsCell.forThisAsPartitionBuildReduction(
                cce, ccei
            )  # <- method updating family field array
            nbOfNewFamsToCreate = famMapI.getNumberOfTuples() - 1
            famsMergedCell[famsMergedCell.findIdsGreaterOrEqualTo(localFamRef)] += (
                famIdZeroForNewFamilies - localFamRef
            )
            updateMap(
                mmOut, lev, famMap, famMapI, famIdZeroForNewFamilies - localFamRef
            )
            famsMergedCell = -famsMergedCell
            famIdZeroForNewFamilies += nbOfNewFamsToCreate
        o2nCells, newNbCells = ml.DataArrayInt.ConvertIndexArrayToO2N(
            m1.getNumberOfCells(), cce, ccei
        )
        n2oCells = o2nCells.invertArrayO2N2N2O(newNbCells)
        n2oHolder[lev] = n2oCells
        m1 = m1[n2oCells]
        m1.setCoords(newCoords)
        m1.setName(self.getName())
        mmOut[lev] = m1
        if famsCell:
            mmOut.setFamilyFieldArr(lev, famsMergedCell)

    famsNode = self.getFamilyFieldAtLevel(1)
    if famsNode:
        localFamRef = famsNode.getMaxAbsValueInArray() + 1
        famsMergedNode, famMap, famMapI = famsNode.forThisAsPartitionBuildReduction(
            cNode, ciNodes
        )
        famsMergedNode[famsMergedNode.findIdsGreaterOrEqualTo(localFamRef)] += (
            famIdZeroForNewFamilies - localFamRef
        )
        updateMap(mmOut, 1, famMap, famMapI, famIdZeroForNewFamilies - localFamRef)
        mmOut.setFamilyFieldArr(1, famsMergedNode)
    return mmOut, n2oHolder


def FindIdFromPathAndPattern(fname, pat):
    import re
    from pathlib import Path

    patRe = Path(pat).name.replace("*", "([\d]+)")
    patReRe = re.compile(patRe)
    m = patReRe.match(Path(fname).name)
    if not m:
        raise RuntimeError("Unrecognized pattern {} in file {}".format(pat, fname))
    return int(m.group(1))


def GetNodesFusionInfoFromJointsOf(pat: str):
    """
    [EDF32671]. This method expects that each MED file fitting pat pattern contains joints with correspondance on NODES. If yes a n2o conversion array will be computed and returned
    This output may be used by MEDFileUMesh.fuseNodesAndCells.


    :param pat : Pattern pointing to MED files candidates of fusion.
    :return: tuple of size 4 ( cNode, ciNodes, o2nNodes, n2oNodes ) fully defining the merge of nodes ( may be useful for MEDFileUMesh.fuseNodesAndCells )
    """

    import re
    from glob import glob
    import MEDLoader as ml

    def GetAllCommonNodesRegardingJoints(iPart, fileToMerge):
        """
        [EDF32671] : Return list of 2 components arrays giving correspondances. Name of components returns the rank of proc attached.
        """

        def RetriveCorrespondanceForOneJoint(iPart, joint):
            """
            Returns for one joint the
            """
            if joint.getNumberOfSteps() != 1:
                raise NotImplementedError("Juste single timestep joint supported")
            # p0 is for receiving proc. p1 is for sending proc
            p0, p1 = [int(p) for p in re.split("[\s]+", joint.getJointName())]
            if (p0 != iPart) and (p1 != iPart):
                raise RuntimeError(
                    "Unexpected joint name {!r} in proc {}".format(
                        joint.getJointName(), iPart
                    )
                )
            # Find correspondance on NODES. Expected to have exactly one
            cors = [
                cor
                for cor in [
                    joint[0].getCorrespondenceAtPos(i)
                    for i in range(joint[0].getNumberOfCorrespondences())
                ]
                if cor.getLocalGeometryType() == ml.NORM_ERROR
            ]
            if len(cors) != 1:
                raise RuntimeError(
                    "No correspondances lying on NODES in {}".format(fileToMerge)
                )
            cor = cors[0].getCorrespondence()
            # put correspondence cor in right shape. 2 components and in C format
            cor = cor - 1
            cor.rearrange(2)
            #
            if p0 == iPart:
                # receiving
                pOther = p1
            else:
                # sending
                pOther = p0
            if pOther < iPart:
                return None
            cor.setInfoOnComponent(0, str(iPart))
            cor.setInfoOnComponent(1, str(pOther))
            return cor

        allMeshNames = ml.GetMeshNames(fileToMerge)
        if len(allMeshNames) != 1:
            raise NotImplementedError(
                "{} contains not exactly one mesh".format(fileToMerge)
            )
        joints = ml.MEDFileJoints(fileToMerge, allMeshNames[0])

        ret = [
            RetriveCorrespondanceForOneJoint(iPart, joints[iJoint])
            for iJoint in range(joints.getNumberOfJoints())
        ]
        return [elt for elt in ret if elt is not None]

    filesToMerge = sorted(glob(pat), key=lambda x: FindIdFromPathAndPattern(x, pat))

    nbNodesPerProc = []
    allNodesCorr = []

    for fileToMerge in filesToMerge:
        iPart = FindIdFromPathAndPattern(fileToMerge, pat)
        allMeshNames = ml.GetMeshNames(fileToMerge)
        if len(allMeshNames) != 1:
            raise NotImplementedError(
                "{} contains not exactly one mesh".format(fileToMerge)
            )
        _, _, _, curNbNodes = ml.GetUMeshGlobalInfo(fileToMerge, allMeshNames[0])
        curNodeCorr = GetAllCommonNodesRegardingJoints(iPart, fileToMerge)
        allNodesCorr += curNodeCorr
        nbNodesPerProc.append(curNbNodes)

    # apply node offsets

    nodeOffsets = ml.DataArrayInt(nbNodesPerProc)
    nodeOffsets.computeOffsetsFull()
    for elt in allNodesCorr:
        elt[:, 0] += nodeOffsets[int(elt.getInfoOnComponent(0))]
        elt[:, 1] += nodeOffsets[int(elt.getInfoOnComponent(1))]
    c, ci = ml.DataArrayInt.Aggregate(allNodesCorr).fromListOfPairsToIndexArray()
    totalNbOfNodesWithDup = sum(nbNodesPerProc)
    o2nNodes, newNbNodes = ml.DataArrayInt.ConvertIndexArrayToO2N(
        totalNbOfNodesWithDup, c, ci
    )
    n2oNodes = o2nNodes.invertArrayO2N2N2O(newNbNodes)
    return c, ci, o2nNodes, n2oNodes


def AggregateMEDFilesNoProfilesNoFusion(pat: str, fnameOut: str, logLev=logging.INFO):
    """
    This method is useful to aggregate split MED files into a single one.

    This method fuse content of MED files in pat and put the result into fnameOut. For the moment profiles are not managed. And only one mesh is supported. All MED files
    for pat are expected to have the same structure.

    Pay attention nodes and cells may be duplicated by this method. To remove cells/nodes duplication call fuseCellsAndNodes

    :param pat: pattern of MED files to be aggregated
    :param fnameOut: output file storing the result
    """
    import MEDLoader as ml
    import contextlib
    from glob import glob
    from distutils.version import StrictVersion

    logger = getLogger(logLev)
    filesToMerge = sorted(glob(pat), key=lambda x: FindIdFromPathAndPattern(x, pat))
    inpVersion = StrictVersion(ml.MEDFileVersionOfFileStr(filesToMerge[0])).version
    meshes = [ml.MEDFileMesh.New(elt) for elt in filesToMerge]
    mm = ml.MEDFileUMesh.Aggregate(meshes)
    mm.writeXX(fnameOut, 2, *inpVersion)
    allFields = ml.GetAllFieldNames(filesToMerge[0])
    ## Trés important on vérifie l'absence de profile
    for elt in allFields:
        f1ts = ml.MEDFileField1TS(filesToMerge[0], elt)
        assert len(f1ts.getPflsReallyUsed()) == 0
    ##

    fmts = [
        [ml.MEDFileFieldMultiTS(fn, fieldName, False) for fn in filesToMerge]
        for fieldName in allFields
    ]

    for iField, listOfFmts in enumerate(fmts):
        refField = listOfFmts[0]
        nbTs = len(refField)
        for iTs in range(nbTs):
            with contextlib.ExitStack() as stack:
                for iPart in range(len(listOfFmts)):
                    stack.enter_context(listOfFmts[iPart][iTs])
                logger.info(
                    f"Dealing field {refField.getName()!r} time step {listOfFmts[0][iTs].getTime()[2]}"
                )
                mcf = [
                    fmts[iTs].field(meshes[iPart])
                    for iPart, fmts in enumerate(listOfFmts)
                ]
            fagg = ml.MEDCouplingFieldDouble.MergeFields(mcf)
            if fagg.getDiscretization().getEnum() != ml.ON_NODES:
                m = fagg.getMesh().deepCopy()
                o2n = m.sortCellsInMEDFileFrmt()
                n2o = o2n.invertArrayO2N2N2O(m.getNumberOfCells())
                fagg = fagg[n2o]
            f1tsOut = ml.MEDFileField1TS()
            f1tsOut.setFieldNoProfileSBT(fagg)
            f1tsOut.writeXX(fnameOut, 0, *inpVersion)


def FuseCellsAndNodesInMEDFile(
    fnameIn, fnameOut, compType=2, eps=1e-6, logLev=logging.INFO, infoWrapNodes=None
):
    """
    This method read fnameIn MED file, perform operation of fusion and write the result back MED fnameOut file.

    Warning : fnameIn is expected to follow MED file convention. Positive fam values for nodes, negative for cells.

    Warning : fnameIn and fnameOut are expected to be separate files.

    See MEDFileUMesh.fuseNodesAndCells for doc of other params.
    """
    from distutils.version import StrictVersion
    import MEDLoader as ml

    def reduceOnCells(fIn, n2os):
        fOut = fIn[n2os[0]]
        return fOut

    def reduceOnNodes(fIn, n2os):
        fIn.setArray(fIn.getArray()[n2os[1]])
        return fIn

    inpVersion = StrictVersion(ml.MEDFileVersionOfFileStr(fnameIn)).version
    logger = getLogger(logLev)
    mm = ml.MEDFileMesh.New(fnameIn)
    mm.removeOrphanFamilies()
    mmOut, n2os = mm.fuseNodesAndCellsAdv(compType, eps, logLev, infoWrapNodes)
    allFields = ml.GetAllFieldNames(fnameIn)
    mmOut.writeXX(fnameOut, 2, *inpVersion)
    logger.info(f"Writing mesh into {fnameOut}")
    fmtss = [
        ml.MEDFileFieldMultiTS(fnameIn, fieldName, False) for fieldName in allFields
    ]
    for fmts in fmtss:
        for f1ts in fmts:
            with f1ts:
                logger.info(
                    f"Dealing field {f1ts.getName()!r} time step {f1ts.getTime()[2]}"
                )
                fIn = f1ts.field(mm)
                if fIn.getDiscretization().getEnum() == ml.ON_NODES:
                    fOut = reduceOnNodes(fIn, n2os)
                else:
                    fOut = reduceOnCells(fIn, n2os)
            f1tsOut = ml.MEDFileField1TS()
            f1tsOut.setFieldNoProfileSBT(fOut)
            f1tsOut.writeXX(fnameOut, 0, *inpVersion)


def MEDFileUMeshTetrahedrize(self, splitType, logLev=logging.INFO):
    """
    [EDF30178] : Method splitting hexa,prisms and underlying quads into resp and underlying triangles
    """
    import MEDLoader as ml

    def getLogger(level=logging.INFO):
        FORMAT = "%(levelname)s : %(asctime)s : [%(filename)s:%(funcName)s:%(lineno)s] : %(message)s"
        logging.basicConfig(format=FORMAT, level=level)
        return logging.getLogger()

    logger = getLogger(logLev)

    def HexaSpliter(splitType):
        """
        :param splitType : see MEDCouplingUMesh.simplexize
        """
        m3 = ml.MEDCouplingUMesh("", 3)
        m3.allocateCells()
        m3.insertNextCell(ml.NORM_HEXA8, list(range(8)))
        m3.simplexize(splitType)
        m3 = ml.MEDCoupling1SGTUMesh(m3)
        conn = m3.getNodalConnectivity()
        conn.rearrange(4)
        return conn.getValuesAsTuple()

    def Penta6Spliter(splitType):
        return [(3, 5, 4, 1), (1, 3, 5, 0), (0, 5, 1, 2)]

    def SplitByType(geoType, splitType):
        m = {ml.NORM_HEXA8: HexaSpliter, ml.NORM_PENTA6: Penta6Spliter}
        return m[geoType](splitType)

    def SplitMeshByType(splitType, m0st, famSt=None):
        """
        :param m0st: MEDCoupling1SGTUMesh instance to be split
        :param famSt: DataArrayInt storing input family field attached to m0st
        """
        conn = m0st.getNodalConnectivity()[:]
        conn.rearrange(m0st.getNumberOfNodesPerCell())
        geoType = m0st.getCellModelEnum()
        subTetra = SplitByType(geoType, splitType)
        famOut = None
        if famSt:
            famOut = famSt.duplicateEachTupleNTimes(len(subTetra))
        m0stTetras = ml.MEDCoupling1SGTUMesh(m0st.getName(), ml.NORM_TETRA4)
        m0stTetras.setCoords(self.getCoords())
        connTetras = ml.DataArrayInt.Meld([conn[:, elt] for elt in subTetra])
        connTetras.rearrange(1)
        m0stTetras.setNodalConnectivity(connTetras)
        return m0stTetras.buildUnstructured(), famOut

    def LocateTwoTrisForEachQuad(quads, tris):
        """
        This function locate for each quad in quads the 2 triangles among triangles into tris.

        :param quads: 4 components DataArrayInt storing nodal conn of quad4
        :param tris: 3 components DataArrayInt storing nodal conn containing division of quad4 to locate
        """
        from itertools import combinations

        quads.sortPerTuple(True)
        tris.sortPerTuple(True)
        curCompoId = ml.DataArrayInt(len(quads))
        curCompoId[:] = 0
        res = ml.DataArrayInt(len(quads) * 2)
        res[:] = -1
        for elt in combinations(range(4), 3):
            arr = ml.DataArrayInt.Aggregate([quads[:, elt], tris])
            offset = len(quads)
            c, ci = arr.findCommonTuples(offset)
            if not ci.deltaShiftIndex().isUniform(2):
                raise RuntimeError("Duplication of tris detected should never happen !")
            c.rearrange(2)
            if not c[:, 0].findIdsGreaterOrEqualTo(offset).empty():
                raise RuntimeError("Duplication of tris detected should never happen !")
            if not curCompoId[c[:, 0]].findIdsGreaterOrEqualTo(2).empty():
                raise RuntimeError(
                    "Internal Error : Quad4 is mapped into more than 2 sub cell triangles ! Something is wrong ! Presence of 3D overlapping cells ?"
                )
            res[2 * c[:, 0] + curCompoId[c[:, 0]]] = c[:, 1] - offset
            curCompoId[c[:, 0]] += 1
        if not curCompoId.isUniform(2):
            raise RuntimeError(
                "It smells very bad ! Impossible to find 2 triangles for some of quadrangles !"
            )
        res.rearrange(2)
        return res

    def deal3D(mmOut, splitType):
        """
        : return : 3D cells in self having a QUAD4 as subcell candidate of spliting
        """
        m0 = self[0]
        m0s = [ml.MEDCoupling1SGTUMesh(elt) for elt in m0.splitByType()]
        fams0 = self.getFamilyFieldAtLevel(0)
        outSubMesh = []
        outFams = []
        startCellId = 0
        for m0st in m0s:
            endCellId = startCellId + m0st.getNumberOfCells()
            famSt = fams0[startCellId:endCellId]
            geoType = m0st.getCellModelEnum()
            if geoType == ml.NORM_TETRA4:
                outSubMesh.append(m0st.buildUnstructured())
                outFams.append(famSt)
                continue
            m0StSplit, famStOut = SplitMeshByType(splitType, m0st, famSt)
            outFams.append(famStOut)
            outSubMesh.append(m0StSplit)
            startCellId = endCellId
        m0tetra = ml.MEDCouplingUMesh.MergeUMeshesOnSameCoords(outSubMesh)
        fam0tetra = ml.DataArrayInt.Aggregate(outFams)
        m0tetra.setDescription(self.getDescription())
        m0tetra.setName(self.getName())
        mmOut[0] = m0tetra
        mmOut.setFamilyFieldArr(0, fam0tetra)
        return ml.MEDCouplingUMesh.MergeUMeshesOnSameCoords(
            [
                elt.buildUnstructured()
                for elt in m0s
                if elt.getCellModelEnum() != ml.NORM_TETRA4
            ]
        )

    def deal2D(mmOut, meshContainingQuadsAsSubCells, splitType):
        m1 = self[-1]
        m1s = [ml.MEDCoupling1SGTUMesh(elt) for elt in m1.splitByType()]
        managed2DTypes = [ml.NORM_TRI3, ml.NORM_QUAD4]
        quads4 = [elt for elt in m1s if elt.getCellModelEnum() == ml.NORM_QUAD4]
        if not all(
            [elt.getCellModelEnum() in [ml.NORM_TRI3, ml.NORM_QUAD4] for elt in m1s]
        ):
            typesStr = [
                ml.MEDCouplingUMesh.GetReprOfGeometricType(elt.getCellModelEnum())
                for elt in m1s
            ]
            managedTypesStr = [
                ml.MEDCouplingUMesh.GetReprOfGeometricType(elt)
                for elt in managed2DTypes
            ]
            raise RuntimeError(
                f"Some geotype in -1 level ( {typesStr} ) are not in managed types ( {managedTypesStr} )"
            )
        if len(quads4) == 1:
            quads4 = quads4[0]
            pass
        logger.debug("Starting to deduce triangulation of quads in -1 level")
        logger.debug(
            "Starting to compute sub cells of 3D cells containing QUAD4 as subcell"
        )
        two2DCellContainingQuads, _, _, rd, rdi = (
            meshContainingQuadsAsSubCells.buildDescendingConnectivity()
        )
        tmp = ml.MEDCouplingUMesh.MergeUMeshesOnSameCoords(
            [quads4.buildUnstructured(), two2DCellContainingQuads]
        )
        offset = quads4.getNumberOfCells()
        logger.debug("Try to reduce list of 3D cells containing QUAD4 as subcell")
        cce, ccei = tmp.findCommonCells(2, offset)
        if not ccei.deltaShiftIndex().isUniform(2):
            raise RuntimeError("Case of fusable quad4 not managed")
        cce.rearrange(2)
        if not cce[:, 0].findIdsGreaterOrEqualTo(offset).empty():
            raise RuntimeError("Case of fusable quad4 not managed")
        cells3DToKeep, _ = ml.DataArrayInt.ExtractFromIndexedArrays(
            cce[:, 1] - offset, rd, rdi
        )
        cells3DToKeep.sort()
        cells3DToKeep = cells3DToKeep.buildUnique()
        threedCellsLyingOnQuads = meshContainingQuadsAsSubCells[cells3DToKeep]
        threedCellsLyingOnQuads.sortCellsInMEDFileFrmt()
        logger.debug("Start to compute the most compact list of tetras")
        allSubTetras = ml.MEDCouplingUMesh.MergeUMeshesOnSameCoords(
            [
                SplitMeshByType(splitType, ml.MEDCoupling1SGTUMesh(elt))[0]
                for elt in threedCellsLyingOnQuads.splitByType()
            ]
        )
        allSubTris = ml.MEDCoupling1SGTUMesh(
            allSubTetras.buildDescendingConnectivity()[0]
        )
        cSubTris = allSubTris.getNodalConnectivity()[:]
        cSubTris.rearrange(3)
        cQuads4 = quads4.getNodalConnectivity()[:]
        cQuads4.rearrange(4)
        logger.debug(
            "Start to find the right split of input quads to respect conformity with previous 3D splitting"
        )
        res = LocateTwoTrisForEachQuad(cQuads4, cSubTris)

        m1Out = ml.MEDCoupling1SGTUMesh(self.getName(), ml.NORM_TRI3)
        m1Out.copyTinyInfoFrom(self[0])
        m1Out.setCoords(self.getCoords())
        res.rearrange(1)
        cSubTris[res]
        connOut = cSubTris[res]
        connOut.rearrange(1)
        m1Out.setNodalConnectivity(connOut)
        m1Out = ml.MEDCouplingUMesh.MergeUMeshesOnSameCoords(
            [
                elt.buildUnstructured()
                for elt in m1s
                if elt.getCellModelEnum() == ml.NORM_TRI3
            ]
            + [m1Out.buildUnstructured()]
        )
        m1Out.copyTinyInfoFrom(self[0])
        mmOut[-1] = m1Out
        famM1 = self.getFamilyFieldAtLevel(-1)
        if famM1:
            outFams = []
            logger.debug("Start dealing families of 2D cells")
            startCellId = 0
            for m1st in m1s:
                endCellId = startCellId + m1st.getNumberOfCells()
                famSt = famM1[startCellId:endCellId]
                startCellId = endCellId
                geoType = m1st.getCellModelEnum()
                if geoType == ml.NORM_TRI3:
                    outFams.append(famSt)
                elif geoType == ml.NORM_QUAD4:
                    outFams.append(famSt.duplicateEachTupleNTimes(2))
                else:
                    raise RuntimeError("Not managed geo type !")
            mmOut.setFamilyFieldArr(-1, ml.DataArrayInt.Aggregate(outFams))

    #
    if self.getMeshDimension() != 3:
        raise RuntimeError(
            f"Expecting mesh with dimension 3 ! Dimension is {self.getMeshDimension()}"
        )
    mmOut = ml.MEDFileUMesh()
    levs = self.getNonEmptyLevels()
    logger.info("Treating 3D level")
    meshContainingQuadsAsSubCells = deal3D(mmOut, splitType)
    if -1 in levs:
        deal2D(mmOut, meshContainingQuadsAsSubCells, splitType)
    # dealing remaining levs not impacting by tetrahedrization
    for remainingLev in [elt for elt in self.getNonEmptyLevels() if elt not in [0, -1]]:
        logger.debug(f"Dealing with level {remainingLev}")
        mLev = self[remainingLev]
        mmOut[remainingLev] = mLev
        famField = self.getFamilyFieldAtLevel(remainingLev)
        if famField:
            mmOut.setFamilyFieldArr(remainingLev, famField)
    #
    mmOut.copyFamGrpMapsFrom(self)
    return mmOut


def MEDFileUMeshReduceToCells(self, level, keepCells, removeOrphanNodes=True):
    """
    Method returning a new MEDFileUMesh, restriction of self to level and keepCell cells at this level.
    This method also

    :param level: Specifies the top level of the returned MEDFileUMesh expected
    :param keepCells: A DataArrayInt specifying cell ids at level level of self
    :param removeOrphanNodes: Specifies if orphan nodes should be removed at the end

    see also MEDFileUMesh.extractPart
    """
    import MEDLoader as ml

    subLevs = [l for l in self.getNonEmptyLevels() if l <= level]
    subMeshes = [self[lev] for lev in subLevs]
    allFamilyFields = [self.getFamilyFieldAtLevel(lev) for lev in subLevs]
    allRefMesh = subMeshes[0]
    refMesh = allRefMesh[keepCells]

    mmOut = ml.MEDFileUMesh()
    # level 0
    mmOut[0] = refMesh
    mmOut.setFamilyFieldArr(0, allFamilyFields[0][keepCells])

    # subLevels
    for curLev, meshLev, famFieldLev in zip(
        subLevs[1:], subMeshes[1:], allFamilyFields[1:]
    ):
        allMeshLev, d, di, rd, rdi = allRefMesh.explodeMeshTo(curLev - level)
        a, b = allMeshLev.areCellsIncludedIn(meshLev, 2)
        if not a:
            raise RuntimeError("Error in mesh {}")
        dlev, dlevi = ml.DataArrayInt.ExtractFromIndexedArrays(keepCells, d, di)
        dlev2 = dlev.buildUniqueNotSorted()
        cellsToKeepLev = ml.DataArrayInt.BuildIntersection([dlev2, b])
        cellsToKeepLev = b.indicesOfSubPart(cellsToKeepLev)
        cellsToKeepLev.sort()
        mmOut[curLev] = meshLev[cellsToKeepLev]
        mmOut.setFamilyFieldArr(curLev, famFieldLev[cellsToKeepLev])

    allFamNodes = mmOut.getFamilyFieldAtLevel(1)
    if allFamNodes:
        mmOut.setFamilyFieldArr(1, allFamNodes[:])

    if removeOrphanNodes:
        mmOut.zipCoords()

    mmOut.copyFamGrpMapsFrom(self)
    return mmOut
