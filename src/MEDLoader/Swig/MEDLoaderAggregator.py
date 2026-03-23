#  -*- coding: utf-8 -*-
# Copyright (C) 2026  CEA, EDF
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

from dataclasses import dataclass
import logging
import abc

__iiiPfl = 0


def GetUniquePflName() -> str:
    global __iiiPfl
    ret = f"pfl_{__iiiPfl}"
    __iiiPfl += 1
    return ret


def positionLogger(level):
    FORMAT = "%(levelname)s : %(asctime)s : [%(filename)s:%(funcName)s:%(lineno)s] : %(message)s"
    logging.basicConfig(format=FORMAT, level=level)
    logging.getLogger().setLevel(level)


def getLogger():
    return logging.getLogger()


def MEDFileUMeshGetDistributionOfTypesHelper(self):
    """
    :param self: ml.MEDFileUMesh
    Example of return : [(14, 20), (18, 40), (3, 16), (4, 17), (1, 9)] means : 20 TETRA4, 40 HEXA8...
    """
    import MEDLoader as ml

    allTypes = ml.DataArrayInt(
        sum(
            [
                list(self.getDistributionOfTypes(lev))
                for lev in self.getNonEmptyLevels()
            ],
            [],
        )
    )
    allTypes.rearrange(3)
    return allTypes[:, [0, 1]].getValuesAsTuple()


class DefaultIterator(abc.ABC):
    @abc.abstractmethod
    def enterPart(self, iPart: int, f1ts):
        """
        :param f1ts: MEDFileAnyTypeField1TS
        """
        pass

    @abc.abstractmethod
    def leavingPart(self):
        pass

    @abc.abstractmethod
    def enterGeoType(self, geoTypeEnum: int):
        pass

    @abc.abstractmethod
    def hit(self, spatialDisc: int, startPos: int, endPos: int, pfl: str, locStr: str):
        pass

    @abc.abstractmethod
    def leavingGeoType(self):
        pass


class FieldAssignator(abc.ABC):
    """
    Abstract functor class behaving as a callback to call the MEDFileField1TS.setFieldProfile at the lowest level
    """

    @abc.abstractmethod
    def assign(self, f1tsToAssign, f, session, lev: int, pfl):
        """
        :param f1tsToAssign: ml.MEDFileField1TS
        :param f: ml.MEDCouplingFieldDouble
        :param session: CommonSession
        :param lev: level
        :param pfl: ml.DataArrayInt
        """
        raise RuntimeError("To be overloaded")

    def manageProfile(self, session, pfl):
        """
        :param session: CommonSession
        :param pfl: ml.DataArrayInt
        :return: ml.DataArrayInt
        """
        candPfl = ProfileHashable(pfl)
        if candPfl in session._father.pfl_manager:
            candPfl = session._father.pfl_manager.get(candPfl)
        else:
            session._father.pfl_manager[candPfl] = candPfl
        return candPfl.get()


class SimpleFieldAssignator(FieldAssignator):
    def assign(self, f1tsToAssign, f, session, lev: int, pfl):
        """
        :param f1tsToAssign: ml.MEDFileField1TS
        :param f: ml.MEDCouplingFieldDouble
        :param session: CommonSession
        :param lev: level
        :param pfl: ml.DataArrayInt
        """

        zePfl = self.manageProfile(session, pfl)
        f1tsToAssign.setFieldProfile(
            f, session._father.global_mesh, lev, zePfl
        )  # <- ze call


class FieldFuseAssignator(FieldAssignator):
    def __init__(self, mmMerged, n2os):
        """
        :param n2os: dict see MEDFileUMesh.fuseNodesAndCellsAdv
        """
        self._mm_merged = mmMerged
        self._n2os = n2os

    def assign(self, f1tsToAssign, f, session, lev: int, pfl):
        import MEDLoader as ml

        n2o = self._n2os[lev]
        pflName = pfl.getName()
        # first detect elements in pfl still existing after fuse
        pfl_captured = n2o.buildIntersection(pfl)
        # in new numbering ref locate elements in pfl still existing -> new profile
        pfl_in_new_ref = n2o.findIdForEach(pfl_captured)
        pfl_in_new_ref.setName(pflName)
        # to select tupleIds of f.getArray() to keep locate elements still existing in pfl reference
        tupleIdsToKeep = pfl.findIdForEach(pfl_captured)
        f.setArray(f.getArray()[tupleIdsToKeep])
        if lev <= 0:
            support = self._mm_merged[lev][pfl_in_new_ref]
            support.setName(self._mm_merged.getName())
            f.setMesh(support)
        else:
            f.setMesh(
                ml.MEDCouplingUMesh.Build0DMeshFromCoords(
                    self._mm_merged.getCoords()[pfl_in_new_ref]
                )
            )
        f.checkConsistencyLight()
        #
        zePfl = self.manageProfile(session, pfl_in_new_ref)
        f1tsToAssign.setFieldProfile(f, self._mm_merged, lev, zePfl)  # <- ze call


class CommonSession(abc.ABC):
    def __init__(self, father):
        self._father = father
        self._subparts = []

    @abc.abstractmethod
    def getSpationDiscretization(self):
        raise RuntimeError("Must be overloaded")

    @abc.abstractmethod
    def getGeoSupport(self, session):
        """
        :return: ml.MEDCouplingUMesh
        """
        raise RuntimeError("Must be overloaded")

    def append(
        self,
        mesh,
        geoType: int,
        posInAggMesh: int,
        arr,
        startPos: int,
        endPos: int,
        pfl,
    ):
        """
        :param mesh: ml.MEDFileMesh
        :param arr: ml.DataArrayDouble
        :param pfl: ml.DataArrayInt
        """
        self._subparts.append(
            SubPart(
                self,
                self._father._cur_part,
                mesh,
                geoType,
                posInAggMesh,
                arr,
                startPos,
                endPos,
                pfl,
            )
        )

    def build(self, f1ts, fieldAssign: FieldAssignator):
        """
        :param f1ts: ml.MEDFileField1TS
        """
        import MEDLoader as ml

        if len(self._subparts) < 1:
            raise RuntimeError("Internal Error")
        # very important transpose _subparts. From (P0, TRI3), (P0, QUAD4), (P1,TRI3), (P1, QUAD4), (P2, TRI3), (P3,QUAD4)
        # to (P0,TRI3), (P1,TRI3), (P2,TRI3), (P0,QUAD4), (P1,QUAD4), (P3,QUAD4)
        allGeoTypes = ml.AllGeometricTypes()
        allGeoTypes.append(ml.NORM_ERROR)
        dictGeoType = {gt: i for i, gt in enumerate(allGeoTypes)}
        sortedSubparts = sorted(
            self._subparts, key=lambda sbp: dictGeoType[sbp.geoType]
        )
        # end of transposition
        base = sortedSubparts[0]
        lev = base.getLevel()
        baseArray = base.constructArray()
        baseArrays = [baseArray]
        pfl = base.getProfile()
        pfls = [pfl]
        geoSupport = self.getGeoSupport(base)
        geoSupports = [geoSupport]
        getLogger().debug("Start iterating on subparts")
        for sbp in sortedSubparts[1:]:
            baseArrays.append(sbp.constructArray())
            pfls.append(sbp.getProfile())
            geoSupports.append(self.getGeoSupport(sbp))
            pass
        dftValue = 0.0
        maxNbCompo, arrWithMaxNbOfCompo = max(
            [(elt.getNumberOfComponents(), elt) for elt in baseArrays],
            key=lambda x: x[0],
        )
        minNbCompo = min([elt.getNumberOfComponents() for elt in baseArrays])
        if maxNbCompo != minNbCompo:
            for i in range(len(baseArrays)):
                if baseArrays[i].getNumberOfComponents() != maxNbCompo:
                    arrRet = baseArrays[i].changeNbOfComponents(maxNbCompo, dftValue)
                    arrRet.setInfoOnComponents(
                        arrWithMaxNbOfCompo.getInfoOnComponents()
                    )
                    baseArrays[i] = arrRet
        baseArray = baseArray.__class__.Aggregate(baseArrays)
        pfl = pfl.__class__.Aggregate(pfls)
        geoSupport = ml.MEDCouplingUMesh.MergeUMeshes(geoSupports)
        getLogger().debug("Start iterating on subparts")
        pfl.setName(GetUniquePflName())
        f = ml.MEDCouplingFieldDouble(self.getSpationDiscretization())
        f.setName(f1ts.getName())
        dt, it, zeTime = f1ts.getTime()
        f.setTime(zeTime, dt, it)
        f.setArray(baseArray)
        f.setMesh(geoSupport)
        fieldAssign.assign(f1ts, f, self, lev, pfl)  # <- ze call
        pass

    pass


class OnNodesSession(CommonSession):
    def __init__(self, father):
        super().__init__(father)

    @property
    def dim(self):
        return -1

    def getSpationDiscretization(self):
        import MEDLoader as ml

        return ml.ON_NODES

    def getGeoSupport(self, subPart):
        """
        :return: ml.MEDCouplingUMesh:
        """
        import MEDLoader as ml

        pfl = subPart.getProfileInLocalRef()
        coo = subPart.mesh.getCoords()[pfl]
        ret = ml.MEDCouplingUMesh.Build0DMeshFromCoords(coo)
        ret.setName(subPart.mesh.getName())
        return ret


class OnCellsSession(CommonSession):
    def __init__(self, father, dim: int):
        super().__init__(father)
        self._dim = dim

    @property
    def dim(self):
        return self._dim

    def getSpationDiscretization(self):
        import MEDLoader as ml

        return ml.ON_CELLS

    def getGeoSupport(self, subPart):
        """
        :return: ml.MEDCouplingUMesh
        """
        pfl = subPart.getProfileInLocalRef()
        geoSupport = subPart.mesh[self.dim - subPart.mesh.getMeshDimension()][pfl]
        return geoSupport


@dataclass(frozen=True)
class SubPart:
    """
    Represent a contribution inside a geometric type corresponding of a filePart
    """

    father: CommonSession
    cur_part: int
    # mesh : ml.MEDFileMesh
    mesh: list
    geoType: int
    posInAggMesh: int
    # arr : ml.DataArrayDouble
    arr: list
    startPos: int
    endPos: int
    # pfl may be None ml.DataArrayInt
    pfl: list

    def getLevel(self):
        import MEDLoader as ml

        return (
            ml.MEDCouplingUMesh.GetDimensionOfGeometricType(self.geoType)
            - self.mesh.getMeshDimension()
            if self.geoType != ml.NORM_ERROR
            else 1
        )

    def constructArray(self):
        return self.arr[self.startPos : self.endPos]

    def getNbOfEntitiesInPart(self):
        import MEDLoader as ml

        if self.geoType != ml.NORM_ERROR:
            return self.mesh.getNumberOfCellsWithType(self.geoType)
        else:
            return self.mesh.getNumberOfNodes()

    def getProfileInLocalRef(self):
        """
        :return: ml.DataArrayInt
        """
        import MEDLoader as ml

        ret = self.pfl
        if ret is None:
            ret = ml.DataArrayInt.New(self.getNbOfEntitiesInPart())
            ret.iota()
            pass
        else:
            ret = ret[:]
        return ret

    def getProfile(self):
        """
        :return: ml.DataArrayInt
        """
        ret = self.getProfileInLocalRef()
        ret += self.posInAggMesh
        return ret


class AssignmentSession(DefaultIterator):
    def __init__(self, listOfMeshes: list, mmagg, pflMngr: dict):
        """
        :param mmagg: ml.MEDFileUMesh
        """
        import MEDLoader as ml
        from collections import defaultdict

        # immutable accross walk.
        self._mm_agg = mmagg
        # immutable accross walk.
        self._list_of_meshes = listOfMeshes
        # mutable across walk accumalated. list of instances of CommonSession
        self._sessions = []
        # mutable across walk
        self._cur_part = -1
        self._gt = None
        # mutable across walk
        self._f1ts = None
        # mutable across walk accumalated
        self._offsets = {}
        # immutable accross walk. { dim : { geotypeEnum : offset inside current dim } } in aggregated mesh referential
        self._agg_distribution = defaultdict(dict)
        #
        self._pfl_manager = pflMngr
        #
        tmp = defaultdict(int)
        aggDistribution = MEDFileUMeshGetDistributionOfTypesHelper(
            self._mm_agg
        )  # distribution of aggregated mesh
        for gt, nbGt in aggDistribution:
            dim = ml.MEDCouplingUMesh.GetDimensionOfGeometricType(gt)
            self._agg_distribution[dim][gt] = tmp[dim]
            tmp[dim] += nbGt
        pass

    @property
    def global_mesh(self):
        """
        :return: ml.MEDFileUMesh
        """
        return self._mm_agg

    @property
    def pfl_manager(self):
        return self._pfl_manager

    def build(self, f1ts, fieldAssign: FieldAssignator):
        """
        :param f1ts: ml.MEDFileField1TS
        """

        def SessionKey(session) -> int:
            return session.dim

        for session in reversed(sorted(self._sessions, key=SessionKey)):
            session.build(f1ts, fieldAssign)
            pass

    def prepareOffsets(self):
        import MEDLoader as ml

        mesh = self._list_of_meshes[self._cur_part]
        typeDist = MEDFileUMeshGetDistributionOfTypesHelper(mesh)
        for typ, _ in typeDist:
            if typ not in self._offsets:
                self._offsets[typ] = 0
        if ml.NORM_ERROR not in self._offsets:
            self._offsets[ml.NORM_ERROR] = 0

    def updateOffsets(self):
        import MEDLoader as ml

        mesh = self._list_of_meshes[self._cur_part]
        typeDist = MEDFileUMeshGetDistributionOfTypesHelper(mesh)
        for typ, nbOfType in typeDist:
            self._offsets[typ] += nbOfType
        self._offsets[ml.NORM_ERROR] += mesh.getNumberOfNodes()

    def getSessionGivenSpatialDisc(self, spatialDisc: int) -> CommonSession:
        import MEDLoader as ml

        if spatialDisc == ml.ON_NODES:
            sess = [elt for elt in self._sessions if isinstance(elt, OnNodesSession)]
            if len(sess) > 1:
                raise RuntimeError("Internal error")
            if len(sess) == 1:
                return sess[0]
            else:
                sess = OnNodesSession(self)
                self._sessions.append(sess)
                return sess
        elif spatialDisc == ml.ON_CELLS:
            dimRequested = ml.MEDCouplingUMesh.GetDimensionOfGeometricType(self._gt)
            sess = [elt for elt in self._sessions if isinstance(elt, OnCellsSession)]
            sess = [elt for elt in sess if elt.dim == dimRequested]
            if len(sess) > 1:
                raise RuntimeError("Internal error")
            if len(sess) == 1:
                return sess[0]
            else:
                sess = OnCellsSession(self, dim=dimRequested)
                self._sessions.append(sess)
                return sess
        else:
            raise RuntimeError(f"Spatial disc {spatialDisc} not managed yet")

    def enterPart(self, iPart: int, f1ts):
        """
        :param f1ts : ml.MEDFileAnyTypeField1TS
        """
        self._cur_part = iPart
        self._f1ts = f1ts
        self.prepareOffsets()
        pass

    def enterGeoType(self, geoTypeEnum: int):
        self._gt = geoTypeEnum
        pass

    def __hitCommon(self, spatialDisc: int, pfl: str):
        sess = self.getSessionGivenSpatialDisc(spatialDisc)
        zePfl = None if pfl == "" else self._f1ts.getProfile(pfl)
        if self._gt not in self._offsets:
            raise RuntimeError("Internal error")
        return sess, zePfl

    def __hit_on_nodes(self, startPos: int, endPos: int, pfl: str):
        import MEDLoader as ml

        sess, zePfl = self.__hitCommon(ml.ON_NODES, pfl)
        sess.append(
            self._list_of_meshes[self._cur_part],
            ml.NORM_ERROR,
            self._offsets[ml.NORM_ERROR],
            self._f1ts.getUndergroundDataArray(),
            startPos,
            endPos,
            zePfl,
        )

    def __hit_on_cells(
        self, spatialDisc: int, startPos: int, endPos: int, pfl: str, locStr: str
    ):
        import MEDLoader as ml

        sess, zePfl = self.__hitCommon(spatialDisc, pfl)
        dim = ml.MEDCouplingUMesh.GetDimensionOfGeometricType(self._gt)
        if dim not in self._agg_distribution:
            raise RuntimeError("Internal error")
        sess.append(
            self._list_of_meshes[self._cur_part],
            self._gt,
            self._agg_distribution[dim][self._gt] + self._offsets[self._gt],
            self._f1ts.getUndergroundDataArray(),
            startPos,
            endPos,
            zePfl,
        )

    def hit(self, spatialDisc: int, startPos: int, endPos: int, pfl: str, locStr: str):
        import MEDLoader as ml

        if self._gt != ml.NORM_ERROR:
            self.__hit_on_cells(spatialDisc, startPos, endPos, pfl, locStr)
        else:
            self.__hit_on_nodes(startPos, endPos, pfl)

    def leavingGeoType(self):
        self._gt = None
        pass

    def leavingPart(self):
        self.updateOffsets()
        self._f1ts = None
        self._cur_part = -1
        pass


def FieldWalkerTexasRanger(listOfMEDFileField1TS: list, iterator):
    """
    Iterator pattern
    """
    for ifile, f1ts in enumerate(listOfMEDFileField1TS):
        iterator.enterPart(ifile, f1ts)
        dist = f1ts.getFieldSplitedByType()
        for perGeoType in dist:
            gt, listOfDiscs = perGeoType
            iterator.enterGeoType(gt)
            for zeSpDisc in listOfDiscs:
                spatialDisc, (startPos, endPos), pfl, locStr = zeSpDisc
                iterator.hit(spatialDisc, startPos, endPos, pfl, locStr)
            iterator.leavingGeoType()
        iterator.leavingPart()


def AggregateFieldsNoFusionOnListOfF1TS(
    listOfMEDFileField1TS: list, listOfMeshes: list, mmagg, pflMngr: dict
):
    """
    :param mmagg : ml.MEDFileUMesh
    :return: ml.MEDFileField1TS
    """
    import MEDLoader as ml

    if len(listOfMEDFileField1TS) < 1:
        raise RuntimeError("List is excepted to be of size >=1 !")
    asignment = AssignmentSession(listOfMeshes, mmagg, pflMngr)
    getLogger().debug(
        f"Start walking accross {len(listOfMEDFileField1TS)} instances of MEDFileField1TS to prepare assignation"
    )
    FieldWalkerTexasRanger(listOfMEDFileField1TS, asignment)
    getLogger().debug(
        f"End of walking accross {len(listOfMEDFileField1TS)} instances of MEDFileField1TS to prepare assignation"
    )
    ret = ml.MEDFileField1TS()
    ret.setName(listOfMEDFileField1TS[0].getName())
    ret.setTime(*(listOfMEDFileField1TS[0].getTime()))
    assign = SimpleFieldAssignator()
    getLogger().debug(f"Start feeding output MEDFileField1TS result")
    asignment.build(ret, assign)
    getLogger().debug(f"End feeding output MEDFileField1TS result")
    return ret


class ProfileHashable:
    def __init__(self, arr):
        """
        :param arr: ml.DataArrayInt
        """
        self._arr = arr

    def __eq__(self, other):
        if self._arr.getNbOfElems() != other._arr.getNbOfElems():
            return False
        return self._arr.isEqualWithoutConsideringStr(other._arr)

    def __hash__(self):
        return self._arr.getHashCode2()

    def get(self):
        return self._arr


def AggregateMEDFilesNoFusion(pat: str, fnameOut: str, logLev=logging.INFO):
    """
    This method is useful to aggregate split MED files into a single one.

    This method fuse content of MED files in pat and put the result into fnameOut. For the moment profiles are not managed. And only one mesh is supported. All MED files
    for pat are expected to have the same structure.

    Pay attention nodes and cells may be duplicated by this method. To remove cells/nodes duplication call fuseCellsAndNodes

    :param pat: pattern of MED files to be aggregated
    :param fnameOut: output file storing the result
    """
    from MEDLoaderFinalize import FindIdFromPathAndPattern
    import MEDLoader as ml
    import contextlib
    from glob import glob
    from distutils.version import StrictVersion

    positionLogger(logLev)
    filesToMerge = sorted(glob(pat), key=lambda x: FindIdFromPathAndPattern(x, pat))
    inpVersion = StrictVersion(ml.MEDFileVersionOfFileStr(filesToMerge[0])).version
    getLogger().debug(
        f"Start to load all {len(filesToMerge)} meshes in memory to perform aggregation"
    )
    meshes = [ml.MEDFileMesh.New(elt) for elt in filesToMerge]
    getLogger().debug(
        f"End of load all {len(filesToMerge)} meshes in memory to perform aggregation. Start aggregation of meshes"
    )
    mm = ml.MEDFileUMesh.Aggregate(meshes)
    getLogger().debug(f"End aggregation of meshes. Start writing back into {fnameOut}")
    mm.writeXX(fnameOut, 2, *inpVersion)
    getLogger().debug(f"End writing back into {fnameOut}")
    allFields = ml.GetAllFieldNames(filesToMerge[0])
    getLogger().debug(f"Start reading structure of fields of {len(filesToMerge)} files")
    fmts = [
        [ml.MEDFileFieldMultiTS(fn, fieldName, False) for fn in filesToMerge]
        for fieldName in allFields
    ]
    getLogger().debug(
        f"End of reading structure of fields of {len(filesToMerge)} files"
    )
    pflMngr = {}

    for iField, listOfFmts in enumerate(fmts):
        refField = listOfFmts[0]
        nbTs = len(refField)
        for iTs in range(nbTs):
            with contextlib.ExitStack() as stack:
                getLogger().debug(
                    f"Start reading {iTs}/{nbTs} time step fields of {len(filesToMerge)} files for {refField.getName()}"
                )
                for iPart in range(len(listOfFmts)):
                    stack.enter_context(listOfFmts[iPart][iTs])
                listOfF1ts = [
                    listOfFmts[iPart][iTs] for iPart in range(len(listOfFmts))
                ]
                getLogger().info(
                    f"Dealing field {refField.getName()!r} time step {listOfFmts[0][iTs].getTime()[2]}"
                )
                f1tsagg = AggregateFieldsNoFusionOnListOfF1TS(
                    listOfF1ts, meshes, mm, pflMngr
                )
                getLogger().debug(
                    f"Start writing {iTs}/{nbTs} time step field {refField.getName()!r} of {len(filesToMerge)} files for {refField.getName()}"
                )
                f1tsagg.writeXX(fnameOut, 0, *inpVersion)
                getLogger().debug(
                    f"End writing {iTs}/{nbTs} time step field {refField.getName()!r} of {len(filesToMerge)} files for {refField.getName()}"
                )


def FuseCellAndNodesField1TS_NoProfile(f1tsIn, mm, n2os):
    """
    :param f1tsIn: ml.MEDFileField1TS single time step to be merged
    :param mm: original ml.MEDFileMesh not merged
    :param mmMerged: ml.MEDFileMesh merged
    :param n2os: dict see MEDFileUMesh.fuseNodesAndCellsAdv
    :return: ml.MEDFileField1TS
    """

    def reduceOnCells(fIn, n2os):
        fOut = fIn[n2os[0]]
        return fOut

    def reduceOnNodes(fIn, n2os):
        fIn.setArray(fIn.getArray()[n2os[1]])
        return fIn

    import MEDLoader as ml

    with f1tsIn:
        fIn = f1tsIn.field(mm)
        if fIn.getDiscretization().getEnum() == ml.ON_NODES:
            fOut = reduceOnNodes(fIn, n2os)
        else:
            fOut = reduceOnCells(fIn, n2os)
    f1tsOut = ml.MEDFileField1TS()
    f1tsOut.setFieldNoProfileSBT(fOut)
    return f1tsOut


def FuseCellAndNodesField1TS_Profile(f1tsIn, mm, mmMerged, n2os, pflMngr: dict):
    import MEDLoader as ml

    asignment = AssignmentSession([mm], mm, pflMngr)
    FieldWalkerTexasRanger([f1tsIn], asignment)
    ret = ml.MEDFileField1TS()
    ret.setName(f1tsIn.getName())
    ret.setTime(*(f1tsIn.getTime()))
    assign = FieldFuseAssignator(mmMerged, n2os)
    with f1tsIn:
        asignment.build(ret, assign)
    return ret


def FuseCellAndNodesField1TS(
    f1tsIn, mm, mmMerged, n2os, pflMngr: dict, logLev=logging.INFO
):
    """
    :param f1tsIn: ml.MEDFileField1TS single time step to be merged
    :param mm: original ml.MEDFileMesh not merged
    :param mmMerged: ml.MEDFileMesh merged
    :param n2os: tuple see MEDFileUMesh.fuseNodesAndCellsAdv

    :return: ml.MEDFileField1TS
    """
    positionLogger(logLev)
    getLogger().info(
        f"Dealing field {f1tsIn.getName()!r} time step {f1tsIn.getTime()[2]}"
    )
    if len(f1tsIn.getPflsReallyUsed()) > 0:
        return FuseCellAndNodesField1TS_Profile(f1tsIn, mm, mmMerged, n2os, pflMngr)
    else:
        return FuseCellAndNodesField1TS_NoProfile(f1tsIn, mm, n2os)
