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

# EDF34807

from MEDLoader import *
import unittest
from MEDLoaderDataForTest import WriteInTmpDir
from dataclasses import dataclass
from math import modf
import logging


@dataclass(frozen=True)
class MeshInfo:
    name: str
    typeDistribution: list
    nodeOffsetInConn: int
    cooOffset: float

    @property
    def maxDim(self):
        return max(
            [
                MEDCouplingUMesh.GetDimensionOfGeometricType(elt[0])
                for elt in self.typeDistribution
            ]
        )


@dataclass(frozen=True)
class FieldInfo:
    name: str
    mi: MeshInfo
    distribution: dict
    valueOffset: float

    def attachTimeInfo(self, f1ts: MEDFileAnyTypeField1TS):
        f1ts.setTime(1, 2, 3.5)

    def putInfoOnComponents(self, arr: DataArrayDouble):
        arr.setInfoOnComponents([2 * self.name])

    def createFieldDoubleInstance(self) -> MEDCouplingFieldDouble:
        typeOfSpDisc = None
        if len([elt for elt in self.distribution if elt != NORM_ERROR]) == len(
            self.distribution
        ):
            typeOfSpDisc = ON_CELLS
        elif len([elt for elt in self.distribution if elt == NORM_ERROR]) == len(
            self.distribution
        ):
            typeOfSpDisc = ON_NODES
        else:
            raise RuntimeError("Invalid distribution of spatial disc of fields")
        ret = MEDCouplingFieldDouble(typeOfSpDisc)
        ret.setName(self.name)
        return ret


def GenerateMeshFuseReady(mi: MeshInfo) -> MEDFileUMesh:
    """
    Conn of meshes generated is something correct. Nodal conn is in a valid range regarding coordinantes.
    /important : nb of cells of lev 0 must be greater or equal to the nb of cells of lesser levels

    Coo = coo of highest level. nbNodes == nbCells of highest level. iota + decimal part. decimal part locate the process id.
    """

    def generateUMesh(
        mcCellTypeEnum: int,
        nbOfCells: int,
        totalNbOfCells: int,
        nodeOffsetInConn: int,
        cooOffset: float,
    ) -> MEDCouplingUMesh:
        ret = MEDCoupling1SGTUMesh("", mcCellTypeEnum)
        nbOfNodesPerType = MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(
            mcCellTypeEnum
        )
        arr = DataArrayInt(nbOfCells)
        arr.iota()
        arr += nodeOffsetInConn
        conn = DataArrayInt.Meld(
            [(arr + i) % totalNbOfCells for i in range(nbOfNodesPerType)]
        )
        conn.rearrange(1)
        ret.setNodalConnectivity(conn)
        #
        nbOfNodes = nbOfCells  # not a bug : for test purpose nbNodes==nbOfCells
        coo = DataArrayDouble(nbOfNodes)
        coo.iota()
        coo += cooOffset
        coo = DataArrayDouble.Meld([coo, coo + float(2**30), coo + float(2**31)])
        ret.setCoords(coo)
        return ret.buildUnstructured()

    def buildCoordinates(m: MEDCouplingUMesh, cooOffset: float) -> DataArrayDouble:
        # not a bug : for test purpose nbNodes==nbOfCells
        coo = DataArrayDouble(m.getNumberOfCells())
        coo.iota()
        coo += cooOffset
        coo = DataArrayDouble.Meld([coo, coo + float(2**30), coo + float(2**31)])
        coo.setInfoOnComponents(["XX", "YYY", "ZZZZ"])
        return coo

    maxDim = mi.maxDim
    mm = MEDFileUMesh.New()
    mm.setName(mi.name)
    CST_OFFSET_IN_CONN = 10000
    CST_OFFSET_IN_CONN_CT = 1000
    coo = None
    for dim in range(0, -maxDim - 1, -1):
        tdOnDim = [
            elt
            for elt in mi.typeDistribution
            if MEDCouplingUMesh.GetDimensionOfGeometricType(elt[0]) - maxDim == dim
        ]
        if any(tdOnDim):
            # totalNbOfCells also equal to totalNbOfNodes regarding generateUMesh and MergeUMeshes
            totalNbOfCells = sum([nb for _, nb in tdOnDim])
            ms = []
            offset = 0
            offsetNode = 0
            for iType, (mcCellTypeEnum, nbOfCells) in enumerate(tdOnDim):
                m = generateUMesh(
                    mcCellTypeEnum,
                    nbOfCells,
                    totalNbOfCells,
                    nodeOffsetInConn=offset,
                    cooOffset=mi.cooOffset,
                )
                offset += m.getNumberOfCells()
                m.shiftNodeNumbersInConn(-offsetNode)
                offsetNode += m.getNumberOfNodes()
                ms.append(m)
            #

            m = MEDCouplingUMesh.MergeUMeshes(ms)
            m.sortCellsInMEDFileFrmt()
            if coo is None:
                coo = buildCoordinates(m, mi.cooOffset)
                m.setCoords(coo)
            else:
                m.setCoords(coo)
            m.setName(mi.name)
            mm[dim] = m
        pass
    return mm


def GenerateMesh(mi: MeshInfo) -> MEDFileUMesh:
    """
    Conn = ( dim + 1 )*10000 + ( ranktypeof cell + 1 ) * 1000 + 30 * ( id inside type of cell ) + offset coordinates (due to final MergeUMeshes)
    Coo = coo of highest level. nbNodes == nbCells of highest level. iota + decimal part. decimal part locate the process id.
    """

    def generate1GTUMesh(
        mcCellTypeEnum: int, nbOfCells: int, nodeOffsetInConn: int, cooOffset: float
    ) -> MEDCoupling1SGTUMesh:
        CST = 30
        ret = MEDCoupling1SGTUMesh("", mcCellTypeEnum)
        nbOfNodesPerType = MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(
            mcCellTypeEnum
        )
        arr = DataArrayInt(nbOfCells)
        arr.iota()
        arr *= CST
        arr += nodeOffsetInConn + mi.nodeOffsetInConn
        conn = DataArrayInt.Meld([arr + i for i in range(nbOfNodesPerType)])
        conn.rearrange(1)
        ret.setNodalConnectivity(conn)
        #
        nbOfNodes = nbOfCells  # not a bug : for test purpose nbNodes==nbOfCells
        coo = DataArrayDouble(nbOfNodes)
        coo.iota()
        coo += cooOffset
        coo = DataArrayDouble.Meld([coo, coo + float(2**30), coo + float(2**31)])
        ret.setCoords(coo)
        return ret

    maxDim = mi.maxDim
    mm = MEDFileUMesh.New()
    mm.setName(mi.name)
    CST_OFFSET_IN_CONN = 10000
    CST_OFFSET_IN_CONN_CT = 1000
    coo = None
    for dim in range(0, -maxDim - 1, -1):
        tdOnDim = [
            elt
            for elt in mi.typeDistribution
            if MEDCouplingUMesh.GetDimensionOfGeometricType(elt[0]) - maxDim == dim
        ]
        if any(tdOnDim):
            ms = [
                generate1GTUMesh(
                    mcCellTypeEnum,
                    nbOfCells,
                    (maxDim + dim + 1) * CST_OFFSET_IN_CONN
                    + (iType + 1) * CST_OFFSET_IN_CONN_CT,
                    mi.cooOffset,
                )
                for iType, (mcCellTypeEnum, nbOfCells) in enumerate(tdOnDim)
            ]
            m = MEDCouplingUMesh.MergeUMeshes([elt.buildUnstructured() for elt in ms])
            m.sortCellsInMEDFileFrmt()
            if coo is None:
                coo = m.getCoords()
            else:
                m.setCoords(coo)
            m.setName(mi.name)
            mm[dim] = m
        pass
    return mm


def GenerateCellFieldHelper0(
    fi: FieldInfo, partOfFieldToKeep: list, mm: MEDFileUMesh, lev: int
):
    """
    returns :
    - list of paire containing distribution of geotype of mm at level lev (can be seen as a dict sorted like MED file does)
    - list of paire containing geotype and DataArrayInt of cells activated inside geotype in mm referential a lev level. (can be seen as a dict sorted like MED file does)
    """

    def postTreat(value, offset: int, nbOfCellsInCurGT: int):
        if value is not None:
            return value + offset
        else:
            ret = DataArrayInt(nbOfCellsInCurGT)
            ret.iota()
            ret += offset
            return ret

    curLevSplit = [MEDCoupling1SGTUMesh(elt) for elt in mm[lev].splitByType()]
    cts = [elt.getCellModelEnum() for elt in curLevSplit]
    zeMap = {k: v for k, v in fi.mi.typeDistribution}
    if [elt.getNumberOfCells() for elt in curLevSplit] != [
        zeMap[elt.getCellModelEnum()] for elt in curLevSplit
    ]:
        raise RuntimeError("Internal error....")
    ret = [(elt.getCellModelEnum(), elt.getNumberOfCells()) for elt in curLevSplit]
    offset = DataArrayInt([elt[1] for elt in ret])
    offset.computeOffsets()
    offset = offset.getValues()
    tmp = {k: v for k, v in ret}  # map key = geoType value is nbOfCells with this type
    tmp2 = {
        k: off for (k, v), off in zip(ret, offset)
    }  # map key = geoType value is offset regarding whole mesh
    ret2 = [(k, postTreat(v, tmp2[k], tmp[k])) for k, v in partOfFieldToKeep]
    return ret, ret2


__pflid = 0


def getProfileName() -> str:
    global __pflid
    ret = f"pfl_{__pflid}"
    __pflid += 1
    return ret


def GenerateFieldNoProfileCommon(
    fi: FieldInfo, arr: DataArrayDouble, mesh: MEDCouplingMesh
) -> MEDFileAnyTypeField1TS:
    fi.putInfoOnComponents(arr)
    f = fi.createFieldDoubleInstance()
    f.setArray(arr)
    f.setMesh(mesh)
    f1ts = MEDFileField1TS()
    f1ts.setFieldNoProfileSBT(f)  # <- ze call
    fi.attachTimeInfo(f1ts)
    return f1ts


def GenerateNodeFieldNoProfile(
    fi: FieldInfo, partOfFieldToKeep: list, mm: MEDFileUMesh, offsetValue: int
) -> MEDFileAnyTypeField1TS:
    if len(partOfFieldToKeep) != 1:
        raise RuntimeError("Only single part for node supported")
    arr = DataArrayDouble(mm.getNumberOfNodes())
    arr.iota()
    arr += fi.valueOffset + offsetValue
    return GenerateFieldNoProfileCommon(fi, arr, mm[0])


def GenerateCellFieldNoProfile(
    fi: FieldInfo, partOfFieldToKeep: list, mm: MEDFileUMesh, offsetValue: int, lev: int
) -> MEDFileAnyTypeField1TS:
    """
    Ex: partOfFieldToKeep = [(NORM_TRI3, None), (NORM_QUAD4, DataArrayInt())]
    """
    zeTd, _ = GenerateCellFieldHelper0(fi, partOfFieldToKeep, mm, lev)
    zeArr = []
    for iCt, (ct, nbCellsInCt) in enumerate(zeTd):
        arr = DataArrayDouble(nbCellsInCt)
        arr.iota()
        arr += fi.valueOffset + 100 * iCt + offsetValue
        zeArr.append(arr)
    arr = DataArrayDouble.Aggregate(zeArr)
    return GenerateFieldNoProfileCommon(fi, arr, mm[lev])


def GenerateFieldNoProfile(
    fi: FieldInfo, partOfFieldToKeep: list, mm: MEDFileUMesh, offsetValue: int, lev: int
) -> MEDFileAnyTypeField1TS:
    if lev == 1:
        return GenerateNodeFieldNoProfile(fi, partOfFieldToKeep, mm, offsetValue)
    else:
        return GenerateCellFieldNoProfile(fi, partOfFieldToKeep, mm, offsetValue, lev)


def GenerateFieldsProfile(
    fi: FieldInfo, partOfFieldToKeep: list, mm: MEDFileUMesh, offsetValue: int, lev: int
) -> MEDFileAnyTypeField1TS:
    if lev == 1:
        return GenerateNodeFieldsProfile(fi, partOfFieldToKeep, mm, offsetValue, lev)
    else:
        return GenerateCellFieldsProfile(fi, partOfFieldToKeep, mm, offsetValue, lev)


def GenerateFieldsProfileCommon(
    fi: FieldInfo,
    arr: DataArrayDouble,
    mm: MEDFileUMesh,
    lev: int,
    mesh: MEDCouplingMesh,
    pfl: DataArrayInt,
) -> MEDFileAnyTypeField1TS:
    fi.putInfoOnComponents(arr)
    f = fi.createFieldDoubleInstance()
    f.setName(fi.name)
    f.setArray(arr)
    mesh.setName(fi.mi.name)
    f.setMesh(mesh)
    f1ts = MEDFileField1TS()
    f1ts.setFieldProfile(f, mm, lev, pfl)  # <- ze call
    fi.attachTimeInfo(f1ts)
    return f1ts


def GenerateNodeFieldsProfile(
    fi: FieldInfo, partOfFieldToKeep: list, mm: MEDFileUMesh, offsetValue: int, lev: int
) -> MEDFileAnyTypeField1TS:
    if len(partOfFieldToKeep) != 1:
        raise RuntimeError("Only single part for node supported")
    _, pfl = partOfFieldToKeep[0]
    pfl2 = pfl.deepCopy()
    pfl2.setName(getProfileName())
    arr = pfl2.convertToDblArr()
    arr += fi.valueOffset + 5 * 1000
    mesh = MEDCouplingUMesh.Build0DMeshFromCoords(mm.getCoords()[pfl])
    return GenerateFieldsProfileCommon(fi, arr, mm, 0, mesh, pfl2)


def GenerateCellFieldsProfile(
    fi: FieldInfo, partOfFieldToKeep: list, mm: MEDFileUMesh, offsetValue: int, lev: int
) -> MEDFileAnyTypeField1TS:
    """
    Ex: partOfFieldToKeep = [(NORM_TRI3, None), (NORM_QUAD4, DataArrayInt())]
    """
    zeTd, partOfFieldToKeepPost = GenerateCellFieldHelper0(
        fi, partOfFieldToKeep, mm, lev
    )
    pfl = DataArrayInt.Aggregate([elt for gt, elt in partOfFieldToKeepPost])
    pfl.setName(getProfileName())
    zeArr = []
    for iCt, (ct, pflLoc) in enumerate(partOfFieldToKeepPost):
        arr = pflLoc.convertToDblArr()
        arr += fi.valueOffset + 100 * iCt + offsetValue
        zeArr.append(arr)
    arr = DataArrayDouble.Aggregate(zeArr)
    return GenerateFieldsProfileCommon(fi, arr, mm, lev, mm[lev], pfl)


def IsProfile(partOfFieldToKeep: list, mm: MEDFileUMesh, lev: int) -> bool:
    if len([part for ct, part in partOfFieldToKeep if part is not None]) > 0:
        return True  # presence at least of one part with profile -> profile
    # at this point all gt have None as profile.
    # first if we are in NODES mode -> no profile.
    if lev == 1:
        return False
    # detect if all mm[lev] geotypes are fetched by partOfFieldToKeep. If yes no profile. If not : profile
    curLevSplit = [MEDCoupling1SGTUMesh(elt) for elt in mm[lev].splitByType()]
    cts = [elt.getCellModelEnum() for elt in curLevSplit]
    return set(cts) != set([ct for ct, _ in partOfFieldToKeep])


def GenerateField(
    fi: FieldInfo, partOfFieldToKeep: list, mm: MEDFileUMesh, offsetValue: int, lev: int
) -> MEDFileAnyTypeField1TS:
    if IsProfile(partOfFieldToKeep, mm, lev):
        return GenerateFieldsProfile(fi, partOfFieldToKeep, mm, offsetValue, lev)
    else:
        return GenerateFieldNoProfile(fi, partOfFieldToKeep, mm, offsetValue, lev)


def GenerateFile(fileName: str, fi: FieldInfo, generator):
    """
    :param generator: is a py func taking MeshInfo and returning MEDFileUMesh
    Fields values for profiles = 10^4 * timeId + (dim+1) * 10^3  + geotypeid * 10^2 + cellId in geotype in profile
    Fields values for non profiles = 10^4 * timeId + (dim+1) * 10^3  + geotypeid * 10^2 + cellId in geotype
    """
    mm = generator(fi.mi)
    levs = mm.getNonEmptyLevels()
    maxDim = fi.mi.maxDim
    ret = []
    for lev in levs:
        partOfFieldToKeep = [
            (ct, part)
            for ct, part in fi.distribution.items()
            if MEDCouplingUMesh.GetDimensionOfGeometricType(ct) == lev + fi.mi.maxDim
        ]
        if any(partOfFieldToKeep):
            ret.append(
                GenerateField(fi, partOfFieldToKeep, mm, (maxDim + lev + 1) * 1000, lev)
            )
            pass
        pass
    # field on nodes
    partOfFieldToKeep = [
        (ct, part) for ct, part in fi.distribution.items() if ct == NORM_ERROR
    ]
    if any(partOfFieldToKeep):
        ret.append(GenerateField(fi, partOfFieldToKeep, mm, 5 * 1000, 1))
    mm.write(fileName, 2)
    for f1ts in ret:
        f1ts.write(fileName, 0)


class MEDLoaderAggregatorTest(unittest.TestCase):
    @WriteInTmpDir
    def testAggregation0(self):
        """
        Case of aggregation of field on 2 procs mixing non profile part and profile part
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_TRI3 : None, NORM_QUAD4 : None }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMesh )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,6), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_TRI3 : None, NORM_HEXA8 : DataArrayInt([0,1,3,8]) }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMesh )

        f1ts_0 = MEDFileField1TS(file0_name)
        f1ts_0.getUndergroundDataArrayExt()
        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.WARNING)
        mmagg = MEDFileMesh.New( merge_name )
        f1ts = MEDFileField1TS.New( merge_name )
        #
        self.assertTrue( mmagg.getNumberOfNodes() == 60 )
        self.assertTrue( (mmagg.getCoords()[:,1] - mmagg.getCoords()[:,0]).isUniform(2**30,1e-200) )
        self.assertTrue( (mmagg.getCoords()[:,2] - mmagg.getCoords()[:,0]).isUniform(2**31,1e-200) )
        self.assertTrue( mmagg.getNonEmptyLevels() == (0,-1,-2) )
        self.assertTrue( DataArrayDouble([modf(elt)[0] for elt in (mmagg.getCoords()[:50,0]).getValues()]).isUniform(0.5,1e-200) )
        self.assertTrue( DataArrayDouble([modf(elt)[0] for elt in (mmagg.getCoords()[50:60,0]).getValues()]).isUniform(0.25,1e-200) )
        conn0 = DataArrayInt( [14, 141000, 141001, 141002, 141003, 14, 141030, 141031, 141032, 141033, 14, 141060, 141061, 141062, 141063, 14, 141090, 141091, 141092, 141093, 14, 141120, 141121, 141122, 141123, 14, 141150, 141151, 141152, 141153, 14, 141180, 141181, 141182, 141183, 14, 141210, 141211, 141212, 141213, 14, 141240, 141241, 141242, 141243, 14, 141270, 141271, 141272, 141273, 14, 141300, 141301, 141302, 141303, 14, 141330, 141331, 141332, 141333, 14, 141360, 141361, 141362, 141363, 14, 141390, 141391, 141392, 141393, 14, 141420, 141421, 141422, 141423, 14, 141450, 141451, 141452, 141453, 14, 141480, 141481, 141482, 141483, 14, 141510, 141511, 141512, 141513, 14, 141540, 141541, 141542, 141543, 14, 141570, 141571, 141572, 141573, 18, 142020, 142021, 142022, 142023, 142024, 142025, 142026, 142027, 18, 142050, 142051, 142052, 142053, 142054, 142055, 142056, 142057, 18, 142080, 142081, 142082, 142083, 142084, 142085, 142086, 142087, 18, 142110, 142111, 142112, 142113, 142114, 142115, 142116, 142117, 18, 142140, 142141, 142142, 142143, 142144, 142145, 142146, 142147, 18, 142170, 142171, 142172, 142173, 142174, 142175, 142176, 142177, 18, 142200, 142201, 142202, 142203, 142204, 142205, 142206, 142207, 18, 142230, 142231, 142232, 142233, 142234, 142235, 142236, 142237, 18, 142260, 142261, 142262, 142263, 142264, 142265, 142266, 142267, 18, 142290, 142291, 142292, 142293, 142294, 142295, 142296, 142297, 18, 142320, 142321, 142322, 142323, 142324, 142325, 142326, 142327, 18, 142350, 142351, 142352, 142353, 142354, 142355, 142356, 142357, 18, 142380, 142381, 142382, 142383, 142384, 142385, 142386, 142387, 18, 142410, 142411, 142412, 142413, 142414, 142415, 142416, 142417, 18, 142440, 142441, 142442, 142443, 142444, 142445, 142446, 142447, 18, 142470, 142471, 142472, 142473, 142474, 142475, 142476, 142477, 18, 142500, 142501, 142502, 142503, 142504, 142505, 142506, 142507, 18, 142530, 142531, 142532, 142533, 142534, 142535, 142536, 142537, 18, 142560, 142561, 142562, 142563, 142564, 142565, 142566, 142567, 18, 142590, 142591, 142592, 142593, 142594, 142595, 142596, 142597, 18, 142620, 142621, 142622, 142623, 142624, 142625, 142626, 142627, 18, 142650, 142651, 142652, 142653, 142654, 142655, 142656, 142657, 18, 142680, 142681, 142682, 142683, 142684, 142685, 142686, 142687, 18, 142710, 142711, 142712, 142713, 142714, 142715, 142716, 142717, 18, 142740, 142741, 142742, 142743, 142744, 142745, 142746, 142747, 18, 142770, 142771, 142772, 142773, 142774, 142775, 142776, 142777, 18, 142800, 142801, 142802, 142803, 142804, 142805, 142806, 142807, 18, 142830, 142831, 142832, 142833, 142834, 142835, 142836, 142837, 18, 142860, 142861, 142862, 142863, 142864, 142865, 142866, 142867, 18, 142890, 142891, 142892, 142893, 142894, 142895, 142896, 142897, 18, 241050, 241051, 241052, 241053, 241054, 241055, 241056, 241057, 18, 241080, 241081, 241082, 241083, 241084, 241085, 241086, 241087, 18, 241110, 241111, 241112, 241113, 241114, 241115, 241116, 241117, 18, 241140, 241141, 241142, 241143, 241144, 241145, 241146, 241147, 18, 241170, 241171, 241172, 241173, 241174, 241175, 241176, 241177, 18, 241200, 241201, 241202, 241203, 241204, 241205, 241206, 241207, 18, 241230, 241231, 241232, 241233, 241234, 241235, 241236, 241237, 18, 241260, 241261, 241262, 241263, 241264, 241265, 241266, 241267, 18, 241290, 241291, 241292, 241293, 241294, 241295, 241296, 241297, 18, 241320, 241321, 241322, 241323, 241324, 241325, 241326, 241327] )
        self.assertTrue( mmagg[0].getNodalConnectivity().isEqual( conn0 ) )
        conn0I = DataArrayInt( [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 109, 118, 127, 136, 145, 154, 163, 172, 181, 190, 199, 208, 217, 226, 235, 244, 253, 262, 271, 280, 289, 298, 307, 316, 325, 334, 343, 352, 361, 370, 379, 388, 397, 406, 415, 424, 433, 442, 451, 460 ] )
        self.assertTrue( mmagg[0].getNodalConnectivityIndex().isEqual( conn0I ) )
        #
        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        self.assertTrue( f1ts.getUndergroundDataArray().getInfoOnComponents()==['zeFieldzeField'] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 3 )
        #
        tri3Part = [part for part in fieldSpectrum if part[0] == NORM_TRI3]
        self.assertTrue( len(tri3Part) == 1 )# check single TRI3 section if not smell bad
        tri3Part = tri3Part[0]
        tri3Part = tri3Part[1]
        self.assertTrue( len(tri3Part)==1 )# single spatial disc for TRI3
        tri3Part = tri3Part[0]
        spatialDisc, (start,endd), pfl, loc = tri3Part
        self.assertTrue( spatialDisc == ON_CELLS )
        self.assertTrue( pfl == "" ) # no profile because all TRI3 in aggregated are fetched
        self.assertTrue( loc == "" )
        arr = f1ts.getUndergroundDataArray()[start:endd]
        self.assertTrue(arr.isEqualWithoutConsideringStr(DataArrayDouble([13000, 13001, 13002, 13003, 13004, 13005, 13006, 13007, 13008, 13009, 23000, 23001, 23002, 23003, 23004, 23005]),1e-20))
        #
        quad4Part = [part for part in fieldSpectrum if part[0] == NORM_QUAD4]
        self.assertTrue( len(quad4Part) == 1 )# check single QUAD4 section if not smell bad
        quad4Part = quad4Part[0]
        quad4Part = quad4Part[1]
        self.assertTrue( len(quad4Part)==1 )# single spatial disc for QUAD4
        quad4Part = quad4Part[0]
        spatialDisc, (start,endd), pfl, loc = quad4Part
        self.assertTrue( spatialDisc == ON_CELLS )
        self.assertTrue( pfl != "" ) # profile because all QUAD4 in aggregated are not fetched
        self.assertTrue( f1ts.getProfile(pfl).isEqualWithoutConsideringStr( DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11]) ) )
        self.assertTrue( loc == "" )
        arr = f1ts.getUndergroundDataArray()[start:endd]
        self.assertTrue(arr.isEqualWithoutConsideringStr(DataArrayDouble([13100, 13101, 13102, 13103, 13104, 13105, 13106, 13107, 13108, 13109, 13110, 13111 ]),1e-20))
        #
        hexa8Part = [part for part in fieldSpectrum if part[0] == NORM_HEXA8]
        self.assertTrue( len(hexa8Part) == 1 )# check single HEXA8 section if not smell bad
        hexa8Part = hexa8Part[0]
        hexa8Part = hexa8Part[1]
        self.assertTrue( len(hexa8Part)==1 )# single spatial disc for HEXA8
        hexa8Part = hexa8Part[0]
        spatialDisc, (start,endd), pfl, loc = hexa8Part
        self.assertTrue( spatialDisc == ON_CELLS )
        self.assertTrue( f1ts.getProfile(pfl).isEqualWithoutConsideringStr( DataArrayInt([30,31,33,38]) ) )
        arr = f1ts.getUndergroundDataArray()[start:endd]
        self.assertTrue( arr.isEqualWithoutConsideringStr(DataArrayDouble([24000, 24001, 24003, 24008]), 1e-20 ) )
        # fmt: on

    @WriteInTmpDir
    def testAggregation1(self):
        """
        Case of aggregation of field on 2 procs mixing non profile part and profile part
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_TRI3 : None, NORM_QUAD4 : None }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMesh )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,6), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_TRI3 : DataArrayInt([2,4,5]), NORM_HEXA8 : DataArrayInt([0,1,3,8]) }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMesh )

        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.WARNING)
        f1ts = MEDFileField1TS.New( merge_name )

        #
        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        self.assertTrue( f1ts.getUndergroundDataArray().getInfoOnComponents()==['zeFieldzeField'] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 3 )
        #
        tri3Part = [part for part in fieldSpectrum if part[0] == NORM_TRI3]
        self.assertTrue( len(tri3Part) == 1 )# check single TRI3 section if not smell bad
        tri3Part = tri3Part[0]
        tri3Part = tri3Part[1]
        self.assertTrue( len(tri3Part)==1 )# single spatial disc for TRI3
        tri3Part = tri3Part[0]
        spatialDisc, (start,endd), pfl, loc = tri3Part
        self.assertTrue( spatialDisc == ON_CELLS )
        self.assertTrue( pfl != "" ) # no profile because all TRI3 in aggregated are fetched
        self.assertTrue( loc == "" )
        self.assertTrue( f1ts.getProfile(pfl).isEqualWithoutConsideringStr( DataArrayInt([0,1,2,3,4,5,6,7,8,9,12,14,15]) ) )
        arr = f1ts.getUndergroundDataArray()[start:endd]
        self.assertTrue(arr.isEqualWithoutConsideringStr(DataArrayDouble([13000, 13001, 13002, 13003, 13004, 13005, 13006, 13007, 13008, 13009, 23002, 23004, 23005]),1e-20))
        #
        quad4Part = [part for part in fieldSpectrum if part[0] == NORM_QUAD4]
        self.assertTrue( len(quad4Part) == 1 )# check single QUAD4 section if not smell bad
        quad4Part = quad4Part[0]
        quad4Part = quad4Part[1]
        self.assertTrue( len(quad4Part)==1 )# single spatial disc for QUAD4
        quad4Part = quad4Part[0]
        spatialDisc, (start,endd), pfl, loc = quad4Part
        self.assertTrue( spatialDisc == ON_CELLS )
        self.assertTrue( pfl != "" ) # profile because all QUAD4 in aggregated are not fetched
        self.assertTrue( f1ts.getProfile(pfl).isEqualWithoutConsideringStr( DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11]) ) )
        self.assertTrue( loc == "" )
        arr = f1ts.getUndergroundDataArray()[start:endd]
        self.assertTrue(arr.isEqualWithoutConsideringStr(DataArrayDouble([13100, 13101, 13102, 13103, 13104, 13105, 13106, 13107, 13108, 13109, 13110, 13111 ]),1e-20))
        #
        hexa8Part = [part for part in fieldSpectrum if part[0] == NORM_HEXA8]
        self.assertTrue( len(hexa8Part) == 1 )# check single HEXA8 section if not smell bad
        hexa8Part = hexa8Part[0]
        hexa8Part = hexa8Part[1]
        self.assertTrue( len(hexa8Part)==1 )# single spatial disc for HEXA8
        hexa8Part = hexa8Part[0]
        spatialDisc, (start,endd), pfl, loc = hexa8Part
        self.assertTrue( spatialDisc == ON_CELLS )
        self.assertTrue( f1ts.getProfile(pfl).isEqualWithoutConsideringStr( DataArrayInt([30,31,33,38]) ) )
        arr = f1ts.getUndergroundDataArray()[start:endd]
        self.assertTrue( arr.isEqualWithoutConsideringStr(DataArrayDouble([24000, 24001, 24003, 24008]), 1e-20 ) )
        # fmt: on

    @WriteInTmpDir
    def testAggregation2(self):
        """
        Case of aggregation of field on nodes 2 procs mixing non profile part and profile part
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_ERROR : None }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMesh )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,6), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_ERROR : None }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMesh )

        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.WARNING)
        f1ts = MEDFileField1TS.New( merge_name )
        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        arr = f1ts.getUndergroundDataArray()
        self.assertTrue( arr.getInfoOnComponents()==['zeFieldzeField'] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 1 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_ERROR )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_NODES )
        self.assertTrue( pfl == "" ) # profile because all nodes in aggregated
        self.assertTrue( (start,endd) == (0,60) )
        self.assertTrue( arr.isEqualWithoutConsideringStr( DataArrayDouble( [15000, 15001, 15002, 15003, 15004, 15005, 15006, 15007, 15008, 15009, 15010, 15011, 15012, 15013, 15014, 15015, 15016, 15017, 15018, 15019, 15020, 15021, 15022, 15023, 15024, 15025, 15026, 15027, 15028, 15029, 15030, 15031, 15032, 15033, 15034, 15035, 15036, 15037, 15038, 15039, 15040, 15041, 15042, 15043, 15044, 15045, 15046, 15047, 15048, 15049, 25000, 25001, 25002, 25003, 25004, 25005, 25006, 25007, 25008, 25009]) , 1e-200) )
        # fmt: on

    @WriteInTmpDir
    def testAggregation3(self):
        """
        Case of aggregation of field on nodes 2 procs mixing non profile part and profile part
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_ERROR : DataArrayInt([0,2,4,7,14,21,28]) }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMesh )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,6), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_ERROR : None }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMesh )

        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.WARNING)
        f1ts = MEDFileField1TS.New( merge_name )
        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        arr = f1ts.getUndergroundDataArray()
        self.assertTrue( arr.getInfoOnComponents()==['zeFieldzeField'] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 1 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_ERROR )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_NODES )
        self.assertTrue( pfl != "" )
        self.assertTrue( f1ts.getProfile( pfl ).isEqualWithoutConsideringStr( DataArrayInt( [0, 2, 4, 7, 14, 21, 28, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59] ) ) )
        self.assertTrue( (start,endd) == (0,17) )
        self.assertTrue( arr.isEqualWithoutConsideringStr( DataArrayDouble( [15000, 15002, 15004, 15007, 15014, 15021, 15028, 25000, 25001, 25002, 25003, 25004, 25005, 25006, 25007, 25008, 25009]) , 1e-200) )
        # fmt: on

    @WriteInTmpDir
    def testAggregation4(self):
        """
        Case of aggregation of field on nodes 2 procs mixing non profile part and profile part
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_ERROR : None }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMesh )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,6), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_ERROR : DataArrayInt([2,4,7]) }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMesh )

        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.WARNING)
        f1ts = MEDFileField1TS.New( merge_name )

        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        arr = f1ts.getUndergroundDataArray()
        self.assertTrue( arr.getInfoOnComponents()==['zeFieldzeField'] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 1 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_ERROR )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_NODES )
        self.assertTrue( pfl != "" )
        self.assertTrue( f1ts.getProfile( pfl ).isEqualWithoutConsideringStr( DataArrayInt( [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 52, 54, 57] ) ) )
        self.assertTrue( (start,endd) == (0,53) )
        self.assertTrue( arr.isEqualWithoutConsideringStr( DataArrayDouble( [15000, 15001, 15002, 15003, 15004, 15005, 15006, 15007, 15008, 15009, 15010, 15011, 15012, 15013, 15014, 15015, 15016, 15017, 15018, 15019, 15020, 15021, 15022, 15023, 15024, 15025, 15026, 15027, 15028, 15029, 15030, 15031, 15032, 15033, 15034, 15035, 15036, 15037, 15038, 15039, 15040, 15041, 15042, 15043, 15044, 15045, 15046, 15047, 15048, 15049, 25002, 25004, 25007] ) , 1e-200) )
        # fmt: on

    @WriteInTmpDir
    def testAggregation5(self):
        """
        Case of aggregation and fusion based on fields on profiles on NODES
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        merge_name_no_dup = "merge_no_dup.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_ERROR : None }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMeshFuseReady )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,4), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_ERROR : DataArrayInt([1,4,5,6,8]) }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMeshFuseReady )

        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.WARNING)

        c = DataArrayInt([10,50,11,51,31,53,37,56,57,59])
        ci = DataArrayInt.Range(0,12,2)
        # 60 == 50 (#nodes Part0) + 10 (#nodes Part1)
        o2nNodes, newNbNodes = DataArrayInt.ConvertIndexArrayToO2N(60,c,ci)
        n2oNodes = o2nNodes.invertArrayO2N2N2O(newNbNodes)
        self.assertTrue( MEDFileMesh.New( merge_name ).getCoords().getInfoOnComponents() == ["XX","YYY","ZZZZ"])
        FuseCellsAndNodesInMEDFile( merge_name, merge_name_no_dup, infoWrapNodes = ( c, ci, o2nNodes, n2oNodes), logLev = logging.WARNING )
        f1ts = MEDFileField1TS.New( merge_name_no_dup )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 1 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_ERROR )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_NODES )
        # 55 - 4. There 4 couples among 5 impacting the field
        self.assertTrue( start == 0 and endd == 51 )
        self.assertTrue( f1ts.getProfile( pfl ).isEqualWithoutConsideringStr( DataArrayInt( [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 11, 51, 52, 37, 54] ) ) )
        valuesRef = DataArrayDouble( [15000, 15001, 15002, 15003, 15004, 15005, 15006, 15007, 15008, 15009, 15012, 15013, 15014, 15015, 15016, 15017, 15018, 15019, 15020, 15021, 15022, 15023, 15024, 15025, 15026, 15027, 15028, 15029, 15030, 15032, 15033, 15034, 15035, 15036, 15038, 15039, 15040, 15041, 15042, 15043, 15044, 15045, 15046, 15047, 15048, 15049, 25001, 25004, 25005, 25006, 25008] )
        self.assertTrue( f1ts.getUndergroundDataArray().isEqualWithoutConsideringStr( valuesRef, 1e-200) )
        mm_merged = MEDFileMesh.New( merge_name_no_dup )
        self.assertTrue( mm_merged.getCoords().getInfoOnComponents() == ["XX","YYY","ZZZZ"])
        c_ref = DataArrayInt( [14, 0, 1, 2, 3, 14, 1, 2, 3, 4, 14, 2, 3, 4, 5, 14, 3, 4, 5, 6, 14, 4, 5, 6, 7, 14, 5, 6, 7, 8, 14, 6, 7, 8, 9, 14, 7, 8, 9, 10, 14, 8, 9, 10, 11, 14, 9, 10, 11, 12, 14, 10, 11, 12, 13, 14, 11, 12, 13, 14, 14, 12, 13, 14, 15, 14, 13, 14, 15, 16, 14, 14, 15, 16, 17, 14, 15, 16, 17, 18, 14, 16, 17, 18, 19, 14, 17, 18, 19, 20, 14, 18, 19, 20, 21, 14, 19, 20, 21, 22, 18, 20, 21, 22, 23, 24, 25, 26, 27, 18, 21, 22, 23, 24, 25, 26, 27, 28, 18, 22, 23, 24, 25, 26, 27, 28, 29, 18, 23, 24, 25, 26, 27, 28, 29, 30, 18, 24, 25, 26, 27, 28, 29, 30, 31, 18, 25, 26, 27, 28, 29, 30, 31, 32, 18, 26, 27, 28, 29, 30, 31, 32, 33, 18, 27, 28, 29, 30, 31, 32, 33, 34, 18, 28, 29, 30, 31, 32, 33, 34, 35, 18, 29, 30, 31, 32, 33, 34, 35, 36, 18, 30, 31, 32, 33, 34, 35, 36, 37, 18, 31, 32, 33, 34, 35, 36, 37, 38, 18, 32, 33, 34, 35, 36, 37, 38, 39, 18, 33, 34, 35, 36, 37, 38, 39, 40, 18, 34, 35, 36, 37, 38, 39, 40, 41, 18, 35, 36, 37, 38, 39, 40, 41, 42, 18, 36, 37, 38, 39, 40, 41, 42, 43, 18, 37, 38, 39, 40, 41, 42, 43, 44, 18, 38, 39, 40, 41, 42, 43, 44, 45, 18, 39, 40, 41, 42, 43, 44, 45, 46, 18, 40, 41, 42, 43, 44, 45, 46, 47, 18, 41, 42, 43, 44, 45, 46, 47, 48, 18, 42, 43, 44, 45, 46, 47, 48, 49, 18, 43, 44, 45, 46, 47, 48, 49, 0, 18, 44, 45, 46, 47, 48, 49, 0, 1, 18, 45, 46, 47, 48, 49, 0, 1, 2, 18, 46, 47, 48, 49, 0, 1, 2, 3, 18, 47, 48, 49, 0, 1, 2, 3, 4, 18, 48, 49, 0, 1, 2, 3, 4, 5, 18, 49, 0, 1, 2, 3, 4, 5, 6, 18, 53, 10, 11, 50, 31, 51, 52, 37, 18, 11, 50, 31, 51, 52, 37, 53, 54, 18, 50, 31, 51, 52, 37, 53, 54, 53, 18, 31, 51, 52, 37, 53, 54, 53, 10, 18, 51, 52, 37, 53, 54, 53, 10, 11, 18, 52, 37, 53, 54, 53, 10, 11, 50, 18, 37, 53, 54, 53, 10, 11, 50, 31, 18, 53, 54, 53, 10, 11, 50, 31, 51, 18, 54, 53, 10, 11, 50, 31, 51, 52] )
        ci_ref = DataArrayInt( [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 109, 118, 127, 136, 145, 154, 163, 172, 181, 190, 199, 208, 217, 226, 235, 244, 253, 262, 271, 280, 289, 298, 307, 316, 325, 334, 343, 352, 361, 370, 379, 388, 397, 406, 415, 424, 433, 442, 451] )
        # 60 cells & 60 nodes -> 59 cells & 55 nodes
        coo_ref = DataArrayDouble( [(0.5, 1073741824.5, 2147483648.5), (1.5, 1073741825.5, 2147483649.5), (2.5, 1073741826.5, 2147483650.5), (3.5, 1073741827.5, 2147483651.5), (4.5, 1073741828.5, 2147483652.5), (5.5, 1073741829.5, 2147483653.5), (6.5, 1073741830.5, 2147483654.5), (7.5, 1073741831.5, 2147483655.5), (8.5, 1073741832.5, 2147483656.5), (9.5, 1073741833.5, 2147483657.5), (0.25, 1073741824.25, 2147483648.25), (1.25, 1073741825.25, 2147483649.25), (12.5, 1073741836.5, 2147483660.5), (13.5, 1073741837.5, 2147483661.5), (14.5, 1073741838.5, 2147483662.5), (15.5, 1073741839.5, 2147483663.5), (16.5, 1073741840.5, 2147483664.5), (17.5, 1073741841.5, 2147483665.5), (18.5, 1073741842.5, 2147483666.5), (19.5, 1073741843.5, 2147483667.5), (20.5, 1073741844.5, 2147483668.5), (21.5, 1073741845.5, 2147483669.5), (22.5, 1073741846.5, 2147483670.5), (23.5, 1073741847.5, 2147483671.5), (24.5, 1073741848.5, 2147483672.5), (25.5, 1073741849.5, 2147483673.5), (26.5, 1073741850.5, 2147483674.5), (27.5, 1073741851.5, 2147483675.5), (28.5, 1073741852.5, 2147483676.5), (29.5, 1073741853.5, 2147483677.5), (30.5, 1073741854.5, 2147483678.5), (3.25, 1073741827.25, 2147483651.25), (32.5, 1073741856.5, 2147483680.5), (33.5, 1073741857.5, 2147483681.5), (34.5, 1073741858.5, 2147483682.5), (35.5, 1073741859.5, 2147483683.5), (36.5, 1073741860.5, 2147483684.5), (6.25, 1073741830.25, 2147483654.25), (38.5, 1073741862.5, 2147483686.5), (39.5, 1073741863.5, 2147483687.5), (40.5, 1073741864.5, 2147483688.5), (41.5, 1073741865.5, 2147483689.5), (42.5, 1073741866.5, 2147483690.5), (43.5, 1073741867.5, 2147483691.5), (44.5, 1073741868.5, 2147483692.5), (45.5, 1073741869.5, 2147483693.5), (46.5, 1073741870.5, 2147483694.5), (47.5, 1073741871.5, 2147483695.5), (48.5, 1073741872.5, 2147483696.5), (49.5, 1073741873.5, 2147483697.5), (2.25, 1073741826.25, 2147483650.25), (4.25, 1073741828.25, 2147483652.25), (5.25, 1073741829.25, 2147483653.25), (9.25, 1073741833.25, 2147483657.25), (8.25, 1073741832.25, 2147483656.25)] )
        self.assertTrue( mm_merged[0].getNodalConnectivity().isEqual( c_ref ) )
        self.assertTrue( mm_merged[0].getNodalConnectivityIndex().isEqual( ci_ref ) )
        self.assertTrue( mm_merged[0].getCoords().isEqualWithoutConsideringStr( coo_ref, 1e-200 ) )
        # fmt: on

    @WriteInTmpDir
    def testAggregation6(self):
        """
        Case of aggregation and fusion based on fields on profiles on CELLS
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        merge_name_no_dup = "merge_no_dup.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_TETRA4 : None, NORM_HEXA8 : DataArrayInt([7,8,9,27]) }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMeshFuseReady )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,4), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_HEXA8 : DataArrayInt([0,5,8,9]) }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMeshFuseReady )

        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.WARNING)

        c = DataArrayInt([10,50,11,51,31,53,37,56,57,59])
        ci = DataArrayInt.Range(0,12,2)
        # 60 == 50 (#nodes Part0) + 10 (#nodes Part1)
        o2nNodes, newNbNodes = DataArrayInt.ConvertIndexArrayToO2N(60,c,ci)
        n2oNodes = o2nNodes.invertArrayO2N2N2O(newNbNodes)
        self.assertTrue( MEDFileMesh.New( merge_name ).getCoords().getInfoOnComponents() == ["XX","YYY","ZZZZ"])
        FuseCellsAndNodesInMEDFile( merge_name, merge_name_no_dup, infoWrapNodes = ( c, ci, o2nNodes, n2oNodes), logLev = logging.WARNING )
        f1ts = MEDFileField1TS.New( merge_name_no_dup )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 2 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_TETRA4 )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_CELLS )
        # No profile because all TETRA4 are fetched
        self.assertTrue( start == 0 and endd == 20 )
        self.assertTrue( pfl == "" )
        #
        gt, allGeoDistOnGt = fieldSpectrum[1]
        self.assertTrue( gt == NORM_HEXA8 )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_CELLS )
        # Profile had size == 8 before fusion. 7 after (cell 50 and 59 have merged linked to fusion of nodes on part1)
        self.assertTrue( start == 20 and endd == 27 )
        self.assertTrue( pfl != "" )
        self.assertTrue( f1ts.getProfile( pfl ).isEqualWithoutConsideringStr( DataArrayInt( [7, 8, 9, 27, 35, 38, 30] ) ) )
        values_ref = DataArrayDouble( [14000, 14001, 14002, 14003, 14004, 14005, 14006, 14007, 14008, 14009, 14010, 14011, 14012, 14013, 14014, 14015, 14016, 14017, 14018, 14019, 14127, 14128, 14129, 14147, 24005, 24008, 24009] )
        self.assertTrue( f1ts.getUndergroundDataArray().isEqualWithoutConsideringStr( values_ref, 1e-200 ) )
        self.assertTrue( f1ts.getUndergroundDataArray().getInfoOnComponents() == ['zeFieldzeField'] )
        # fmt: on

    @WriteInTmpDir
    def testggregation7(self):
        """
        Like testAggregation3 but with number of components differents.
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_ERROR : DataArrayInt([0,2,4,7,14,21,28]) }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMesh )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,6), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_ERROR : None }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMesh )
        # Aim of the test : Patching number of components of field. And to be vicious put 2 components for proc 1 and proc0 stays with 1 component.
        mfd = MEDFileData( file1_name )
        f1ts = mfd.getFields()[zeFieldName][0]
        arr = f1ts.getUndergroundDataArray()
        arr2 = arr.changeNbOfComponents(2,0.) ; arr2[:,1] = 100000
        arr2.setInfoOnComponents( ["aa","bbb"] )
        arr.deepCopyFrom( arr2 )
        mfd.getFields()[zeFieldName].setInfo( arr2.getInfoOnComponents() )
        mfd.write( file1_name, 2 ) ; del mfd
        #

        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.ERROR)
        f1ts = MEDFileField1TS.New( merge_name )
        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        arr = f1ts.getUndergroundDataArray()
        self.assertTrue( arr.getInfoOnComponents()==["aa","bbb"] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 1 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_ERROR )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_NODES )
        self.assertTrue( pfl != "" )
        self.assertTrue( f1ts.getProfile( pfl ).isEqualWithoutConsideringStr( DataArrayInt( [0, 2, 4, 7, 14, 21, 28, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59] ) ) )
        self.assertTrue( (start,endd) == (0,17) )
        self.assertTrue( arr.isEqualWithoutConsideringStr( DataArrayDouble( [(15000, 0), (15002, 0), (15004, 0), (15007, 0), (15014, 0), (15021, 0), (15028, 0), (25000, 100000), (25001, 100000), (25002, 100000), (25003, 100000), (25004, 100000), (25005, 100000), (25006, 100000), (25007, 100000), (25008, 100000), (25009, 100000)]) , 1e-200) )

        # now inversely proc0 has 2 components and proc1 has 1 component
        mfd1 = MEDFileData( file1_name )
        f1ts = mfd1.getFields()[zeFieldName][0]
        arr = f1ts.getUndergroundDataArray()
        arr2 = arr.changeNbOfComponents(1,0.)
        arr2.setInfoOnComponents( ["proc1_compo0"] )
        arr.deepCopyFrom( arr2 )
        mfd1.getFields()[zeFieldName].setInfo( arr2.getInfoOnComponents() )
        mfd1.write( file1_name, 2 ) ; del mfd1
        mfd0 = MEDFileData( file0_name )
        f1ts = mfd0.getFields()[zeFieldName][0]
        arr = f1ts.getUndergroundDataArray()
        arr2 = arr.changeNbOfComponents(2,0.) ; arr2[:,1] = 300000
        arr2.setInfoOnComponents( ["proc0_compo0", "proc0_compo1"] )
        arr.deepCopyFrom( arr2 )
        mfd0.getFields()[zeFieldName].setInfo( arr2.getInfoOnComponents() )
        mfd0.write( file0_name, 2 )

        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.ERROR)
        f1ts = MEDFileField1TS.New( merge_name )
        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        arr = f1ts.getUndergroundDataArray()
        self.assertTrue( arr.getInfoOnComponents()==["proc0_compo0","proc0_compo1"] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 1 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_ERROR )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_NODES )
        self.assertTrue( pfl != "" )
        self.assertTrue( f1ts.getProfile( pfl ).isEqualWithoutConsideringStr( DataArrayInt( [0, 2, 4, 7, 14, 21, 28, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59] ) ) )
        self.assertTrue( (start,endd) == (0,17) )
        self.assertTrue( arr.isEqualWithoutConsideringStr( DataArrayDouble([(15000, 300000), (15002, 300000), (15004, 300000), (15007, 300000), (15014, 300000), (15021, 300000), (15028, 300000), (25000, 0), (25001, 0), (25002, 0), (25003, 0), (25004, 0), (25005, 0), (25006, 0), (25007, 0), (25008, 0), (25009, 0)]) , 1e-200) )
        # fmt: on

    @WriteInTmpDir
    def testggregation8(self):
        """
        Like testAggregation3 but with number of levels different for rank1 (richer) and rank0
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_ERROR : DataArrayInt([0,2,4,7,14,21,28]) }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMesh )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,6), (NORM_QUAD4,5), (NORM_HEXA8,10)] # interesting point is here. Presence of -2 level for proc1. -2 is not present pour proc0.
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_ERROR : None }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMesh )
        #
        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.ERROR)
        #
        mm = MEDFileMesh.New( merge_name )
        self.assertEqual( mm.getNumberOfNodes(), 60 )
        self.assertEqual( mm.getNumberOfCellsAtLevel(0), 60 )
        self.assertEqual( mm.getNumberOfCellsAtLevel(-1), 33 )
        self.assertEqual( mm.getNumberOfCellsAtLevel(-2), 5 )
        toTest = MEDCoupling1SGTUMesh( mm[-2] )
        self.assertEqual( toTest.getCellModelEnum(), NORM_SEG2 )
        self.assertTrue( toTest.getNodalConnectivity().isEqual( DataArrayInt([221000, 221001, 221030, 221031, 221060, 221061, 221090, 221091, 221120, 221121]) ) )
        #
        f1ts = MEDFileField1TS.New( merge_name )
        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        arr = f1ts.getUndergroundDataArray()
        self.assertTrue( arr.getInfoOnComponents()==['zeFieldzeField'] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 1 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_ERROR )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_NODES )
        self.assertTrue( pfl != "" )
        self.assertTrue( f1ts.getProfile( pfl ).isEqualWithoutConsideringStr( DataArrayInt( [0, 2, 4, 7, 14, 21, 28, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59] ) ) )
        self.assertTrue( (start,endd) == (0,17) )
        self.assertTrue( arr.isEqualWithoutConsideringStr( DataArrayDouble( [15000, 15002, 15004, 15007, 15014, 15021, 15028, 25000, 25001, 25002, 25003, 25004, 25005, 25006, 25007, 25008, 25009]) , 1e-200) )
        # fmt: on

    @WriteInTmpDir
    def testAggregation9(self):
        """
        Like testAggregation3 but here profile frugality
        """
        # fmt: off
        file0_name = "field0.med"
        file1_name = "field1.med"
        merge_name = "merge.med"
        #
        zeFieldName = "zeField"
        td0 = [(NORM_SEG2,4), (NORM_TRI3,10), (NORM_QUAD4,12), (NORM_TETRA4,20), (NORM_HEXA8,30)]
        mi0 = MeshInfo( name = "mesh0", typeDistribution = td0, nodeOffsetInConn = 100000, cooOffset = 0.5 )
        fd0 = { NORM_ERROR : DataArrayInt([0,2,4,7,14,21,28]) }
        fi0 = FieldInfo( name = zeFieldName, mi = mi0, distribution = fd0, valueOffset = 10000. )
        GenerateFile( file0_name, fi0, GenerateMesh )

        td1 = [(NORM_SEG2,5),(NORM_TRI3,6), (NORM_QUAD4,5), (NORM_HEXA8,10)]
        mi1 = MeshInfo( name = "mesh0", typeDistribution = td1, nodeOffsetInConn = 200000, cooOffset = 0.25 )
        fd1 = { NORM_ERROR : None }
        fi1 = FieldInfo( name = zeFieldName, mi = mi1, distribution = fd1, valueOffset = 20000. )
        GenerateFile( file1_name, fi1, GenerateMesh )

        # create an additional time step
        newTime = ( 3,4, 7.6 )
        for fname in [file0_name, file1_name]:
            f1ts = MEDFileField1TS( fname )
            f1ts.setTime( *newTime )
            f1ts.write( fname, 0 )
        #
        AggregateMEDFilesNoFusion("field*.med",merge_name, logLev = logging.WARNING)
        fs = MEDFileFields( merge_name )
        self.assertEqual( len( fs.getPfls() ), 1 ) # <- the aim of test is here. If equal 2 means that profile management failed
        #
        f1ts = MEDFileField1TS.New( merge_name )
        self.assertTrue( f1ts.getTime() == [1,2,3.5])
        arr = f1ts.getUndergroundDataArray()
        self.assertTrue( arr.getInfoOnComponents()==['zeFieldzeField'] )
        fieldSpectrum = f1ts.getFieldSplitedByType()
        self.assertTrue( len(fieldSpectrum) == 1 )
        gt, allGeoDistOnGt = fieldSpectrum[0]
        self.assertTrue( gt == NORM_ERROR )
        self.assertTrue( len( allGeoDistOnGt ) == 1 )
        spatialDisc, (start,endd), pfl, loc = allGeoDistOnGt[0]
        self.assertTrue( spatialDisc == ON_NODES )
        self.assertTrue( pfl != "" )
        self.assertTrue( f1ts.getProfile( pfl ).isEqualWithoutConsideringStr( DataArrayInt( [0, 2, 4, 7, 14, 21, 28, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59] ) ) )
        self.assertTrue( (start,endd) == (0,17) )
        self.assertTrue( arr.isEqualWithoutConsideringStr( DataArrayDouble( [15000, 15002, 15004, 15007, 15014, 15021, 15028, 25000, 25001, 25002, 25003, 25004, 25005, 25006, 25007, 25008, 25009]) , 1e-200) )
        # fmt: on

    pass


if __name__ == "__main__":
    unittest.main()
