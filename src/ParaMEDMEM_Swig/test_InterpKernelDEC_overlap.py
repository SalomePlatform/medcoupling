#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2026  CEA, EDF
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

from medcoupling import *
from ParaMEDMEMTestTools import WriteInTmpDir
import sys, os
import unittest
import math
from mpi4py import MPI


class ParaMEDMEM_IKo_DEC_Test1(unittest.TestCase):
    """Small case for easy debug"""

    def assertListAlmostEqual(self, list1, list2, places=7):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, places=places)

    def generateFullSource(self):
        """The complete source mesh: 4 squares each divided in 2 diagonaly (so 8 cells in total)"""
        m1 = MEDCouplingCMesh("src_mesh")
        dx = DataArrayDouble([0, 2, 4, 6])
        dy = DataArrayDouble([0, 2])
        m1.setCoords(dx, dy)
        msh = m1.buildUnstructured()
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        fld.setMesh(msh)
        fld.setName("source_F")
        da = DataArrayDouble(msh.getNumberOfCells())
        da.iota()
        da *= 2
        fld.setArray(da)
        numGl = DataArrayInt([5, 158, 1])
        return msh, fld, numGl

    def generateFullTarget(self):
        """The complete target mesh: 2 squares"""
        m1 = MEDCouplingCMesh("tgt_msh")
        dx = DataArrayDouble([0, 3, 6])
        dy = DataArrayDouble([0, 2])
        m1.setCoords(dx, dy)
        msh = m1.buildUnstructured()
        # Receiving side must prepare an empty field that will be filled by DEC:
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        da = DataArrayDouble(msh.getNumberOfCells())
        fld.setArray(da)
        fld.setName("tgt_F")
        fld.setMesh(msh)
        numGl = DataArrayInt([12, 27])
        return msh, fld, numGl

    #
    # Below, the two functions emulating the set up of a piece of the source and target mesh
    # on each proc. Obviously in real world problems, this comes from your code and is certainly
    # not computed by cuting again from scratch the full-size mesh!!
    #
    def getPartialSource(self, rank):
        """Will return an empty mesh piece for rank=2 and 3"""
        msh, f, numGl = self.generateFullSource()
        if rank == 0:
            sub_m, sub_f, sub_num = msh[0:2], f[0:2], numGl[0:2]
        elif rank == 1:
            sub_m, sub_f, sub_num = msh[1:3], f[1:3], numGl[1:3]
        sub_m.zipCoords()
        return sub_m, sub_f, sub_num

    def getPartialTarget(self, rank):
        """One square for each rank"""
        msh, f, numGl = self.generateFullTarget()
        if rank == 2:
            sub_m, sub_f, sub_num = msh[[0, 1]], f[[0, 1]], numGl[[0, 1]]
        elif rank == 3:
            sub_m, sub_f, sub_num = msh[[1, 0]], f[[1, 0]], numGl[[1, 0]]
        sub_m.zipCoords()

        return sub_m, sub_f, sub_num

    @WriteInTmpDir
    def testInterpKernelDEC_overlap_2D_py_1(self):
        """This test illustrates a basic use of the InterpKernelDECWithOverlap.
        Look at the C++ documentation of the class for more informations.
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 4:
            # print("Should be run on 4 procs!")
            return

        # Define two processor groups
        nproc_source = 2
        procs_source = list(range(nproc_source))
        procs_target = list(range(nproc_source, size))

        interface = CommInterface()
        source_group = MPIProcessorGroup(interface, procs_source)
        target_group = MPIProcessorGroup(interface, procs_target)
        idec = InterpKernelDECWithOverlap(source_group, target_group)

        #
        # OK, let's go DEC !!
        #
        if source_group.containsMyRank():
            _, fieldS, numS = self.getPartialSource(rank)
            fieldS.setNature(IntensiveMaximum)  # The only policy supported for now ...
            WriteField("./source_field_part_%d.med" % rank, fieldS, True)
            idec.attachLocalField(fieldS, numS)
            idec.synchronize()
            idec.sendData()

            # send second array with same mesh
            arr = fieldS.getArray()
            arr *= 2

            idec.attachLocalField(fieldS, numS)
            idec.sendData()

        if target_group.containsMyRank():
            mshT, fieldT, numT = self.getPartialTarget(rank)
            fieldT.setNature(IntensiveMaximum)
            WriteUMesh("./target_mesh_part_%d.med" % rank, mshT, True)
            idec.attachLocalField(fieldT, numT)
            idec.synchronize()
            idec.recvData()
            WriteField("./target_field_part_%d.med" % rank, fieldT, True)
            # Now the actual checks:
            if rank == 2:
                self.assertListAlmostEqual(
                    fieldT.getArray().getValues(), [2.0 / 3.0, 10.0 / 3.0], places=7
                )
            elif rank == 3:
                self.assertListAlmostEqual(
                    fieldT.getArray().getValues(), [10.0 / 3.0, 2.0 / 3.0], places=7
                )

            idec.attachLocalField(fieldT, numT)
            idec.recvData()
            # Now the actual checks:
            if rank == 2:
                self.assertListAlmostEqual(
                    fieldT.getArray().getValues(), [4.0 / 3.0, 20.0 / 3.0], places=7
                )
            elif rank == 3:
                self.assertListAlmostEqual(
                    fieldT.getArray().getValues(), [20.0 / 3.0, 4.0 / 3.0], places=7
                )

        # Release DEC (this involves MPI exchanges -- notably the release of the communicator -- so better be done before MPI.Finalize()
        idec.release()
        source_group.release()
        target_group.release()
        MPI.COMM_WORLD.Barrier()

    @WriteInTmpDir
    def testInterpKernelDEC_overlap_2D_py_1_ref(self):
        """This test illustrates a basic use of the InterpKernelDECWithOverlap.
        Look at the C++ documentation of the class for more informations.
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 2:
            # print("Should be run on 2 procs!")
            return

        # Define two processor groups
        procs_source = (0,)
        procs_target = (1,)

        interface = CommInterface()
        source_group = MPIProcessorGroup(interface, procs_source)
        target_group = MPIProcessorGroup(interface, procs_target)
        idec = InterpKernelDECWithOverlap(source_group, target_group)

        #
        # OK, let's go DEC !!
        #
        if source_group.containsMyRank():
            _, fieldS, _ = self.generateFullSource()
            fieldS.setNature(IntensiveMaximum)  # The only policy supported for now ...
            WriteField(
                "./source_field_part_seq.med",
                fieldS,
                True,
            )
            idec.attachLocalField(fieldS)
            idec.synchronize()
            idec.sendData()

            # send second array with same mesh
            arr = fieldS.getArray()
            arr *= 2

            idec.sendData()

        if target_group.containsMyRank():
            mshT, fieldT, _ = self.generateFullTarget()
            fieldT.setNature(IntensiveMaximum)
            WriteUMesh("./target_mesh_part_seq.med", mshT, True)
            idec.attachLocalField(fieldT)
            idec.synchronize()
            idec.recvData()
            WriteField("./target_field_part_seq.med", fieldT, True)
            # Now the actual checks:
            self.assertListAlmostEqual(
                fieldT.getArray().getValues(), [2.0 / 3.0, 10.0 / 3.0], places=7
            )

            idec.recvData()
            self.assertListAlmostEqual(
                fieldT.getArray().getValues(), [4.0 / 3.0, 20.0 / 3.0], places=7
            )

        # Release DEC (this involves MPI exchanges -- notably the release of the communicator -- so better be done before MPI.Finalize()
        idec.release()
        source_group.release()
        target_group.release()
        MPI.COMM_WORLD.Barrier()


class ParaMEDMEM_IKo_DEC_Test2(unittest.TestCase):
    """Generic case"""

    def checkValues(self, ref, field, cells, places=7):
        values = field.getArray().getValues()
        nbCmp = field.getNumberOfComponents()
        for i in range(len(cells)):
            c_id = cells[i]
            for j in range(nbCmp):
                self.assertAlmostEqual(
                    ref[c_id * nbCmp + j], values[i * nbCmp + j], places=places
                )

    def generateFullSource(self):
        """The complete source mesh: 3x3 square mesh"""
        m1 = MEDCouplingCMesh("src_mesh")
        dx = DataArrayDouble([0, 2, 4, 6])
        dy = DataArrayDouble([0, 3, 9, 12])
        m1.setCoords(dx, dy)
        msh = m1.buildUnstructured()
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        fld.setMesh(msh)
        fld.setName("source_F")
        da = DataArrayDouble(msh.getNumberOfCells(), 2)
        da.rearrange(1)
        da.iota()
        da.rearrange(2)
        da[:, 0] *= 2
        da[:, 1] *= 3
        fld.setArray(da)
        numGl = DataArrayInt([5, 158, 1, 6, 7, 8, 2, 257, 0])
        return msh, fld, numGl

    def generateFullTarget(self):
        """The complete target mesh: 2 squares"""
        m1 = MEDCouplingCMesh("ta_msh")
        dx = DataArrayDouble([0, 3, 5, 6])
        dy = DataArrayDouble([0, 2, 7, 12])
        m1.setCoords(dx, dy)
        msh = m1.buildUnstructured()
        # Receiving side must prepare an empty field that will be filled by DEC:
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        da = DataArrayDouble(msh.getNumberOfCells(), 2)
        fld.setArray(da)
        fld.setName("target_F")
        fld.setMesh(msh)
        numGl = DataArrayInt([12, 27, 1, 14, 15875624125, 4, 98, 99, 97])
        return msh, fld, numGl

    #
    # Below, the two functions emulating the set up of a piece of the source and target mesh
    # on each proc. Obviously in real world problems, this comes from your code and is certainly
    # not computed by cuting again from scratch the full-size mesh!!
    #
    def getPartialSource(self, cells):
        # print("src: ", cells)
        msh, f, numGl = self.generateFullSource()
        sub_m, sub_f, sub_num = msh[cells], f[cells], numGl[cells]
        sub_m.zipCoords()
        return sub_m, sub_f, sub_num

    def getPartialTarget(self, cells):
        # print("trg: ", cells)
        msh, f, numGl = self.generateFullTarget()
        sub_m, sub_f, sub_num = msh[cells], f[cells], numGl[cells]
        sub_m.zipCoords()
        return sub_m, sub_f, sub_num

    @WriteInTmpDir
    def testInterpKernelDEC_overlap_2D_py_2(self):
        """This test illustrates a basic use of the InterpKernelDECWithOverlap.
        Look at the C++ documentation of the class for more informations.
        """

        cases = [
            # reference case - 1x1 mpi without overlapping
            [[0], [list(range(9))], False, [1], [list(range(9))], False],
            # 1x1 mpi without overlapping
            [[0], [list(range(9))], True, [1], [list(range(9))], True],
            # 2x3 MPI without overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 4], [5, 6, 7, 8]],
                False,
                [1, 3, 4],
                [[0, 2, 4], [1, 3], [5, 6, 7, 8]],
                False,
            ],
            # 2x3 MPI without overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 4], [5, 6, 7, 8]],
                False,
                [1, 3, 4],
                [[0, 2, 4], [1, 3], [5, 6, 7, 8]],
                True,
            ],
            # 2x3 MPI without overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 4], [5, 6, 7, 8]],
                True,
                [1, 3, 4],
                [[0, 2, 4], [1, 3], [5, 6, 7, 8]],
                False,
            ],
            # 2x3 MPI without overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 4], [5, 6, 7, 8]],
                True,
                [1, 3, 4],
                [[0, 2, 4], [1, 3], [5, 6, 7, 8]],
                True,
            ],
            # 2x3 MPI wih overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 5, 6], [5, 6, 7, 8, 0, 4, 2]],
                True,
                [1, 3, 4],
                [[0, 2, 4, 5, 6], [1, 3], [5, 6, 7, 8, 2, 0]],
                True,
            ],
            # 2x3 MPI without/with overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 4], [5, 6, 7, 8]],
                True,
                [1, 3, 4],
                [[0, 2, 4, 5, 6], [1, 3], [5, 6, 7, 8, 2, 0]],
                True,
            ],
            # 2x3 MPI without/with overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 4], [5, 6, 7, 8]],
                False,
                [1, 3, 4],
                [[0, 2, 4, 5, 6], [1, 3], [5, 6, 7, 8, 2, 0]],
                True,
            ],
            # 2x3 MPI with/without overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 5, 6], [5, 6, 7, 8, 0, 4, 2]],
                True,
                [1, 3, 4],
                [[0, 2, 4], [1, 3], [5, 6, 7, 8]],
                True,
            ],
            # 2x3 MPI with/without overlapping,
            [
                [0, 2],
                [[0, 1, 2, 3, 5, 6], [5, 6, 7, 8, 0, 4, 2]],
                True,
                [1, 3, 4],
                [[0, 2, 4], [1, 3], [5, 6, 7, 8]],
                False,
            ],
            # 2x3 MPI with full-overlapping,
            [
                [0, 2],
                [list(range(9)), list(range(9))],
                True,
                [1, 3, 4],
                [list(range(9)), list(range(9)), list(range(9))],
                True,
            ],
            # 2x3 MPI wih overlapping,
            [
                [1, 2, 4],
                [[0, 2, 4, 5, 6], [1, 3, 5], [5, 6, 7, 8, 2, 0]],
                True,
                [0, 3],
                [[0, 1, 2, 3, 5, 6], [5, 6, 7, 8, 0, 4, 2]],
                True,
            ],
        ]

        for case in cases:
            self.base(*case)

    @WriteInTmpDir
    def base(
        self,
        procs_source,
        cells_source,
        use_overlap_src,
        procs_target,
        cells_target,
        use_overlap_trg,
    ):
        """This test illustrates a basic use of the InterpKernelDECWithOverlap.
        Look at the C++ documentation of the class for more informations.
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank

        if size != len(procs_source) + len(procs_target):
            # print(f"Should be run on {len(procs_source) + len(procs_target)} procs!")
            return

        interface = CommInterface()
        source_group = MPIProcessorGroup(interface, procs_source)
        target_group = MPIProcessorGroup(interface, procs_target)
        idec = InterpKernelDECWithOverlap(source_group, target_group)

        # reference solution without overlapping on 1 mpi
        _, full_src, _ = self.generateFullSource()
        ref_src = full_src.getArray().getValues()
        ref_trg = [
            1.3333333333333333,
            5.0,
            6.0,
            12.0,
            8.0,
            15.0,
            10.933333333333334,
            19.4,
            15.600000000000001,
            26.4,
            17.6,
            29.4,
            20.53333333333333,
            33.8,
            25.200000000000003,
            40.8,
            27.2,
            43.8,
        ]

        #
        # OK, let's go DEC !!
        #
        if source_group.containsMyRank():
            local_rank = source_group.myRank()
            _, fieldS, numS = self.getPartialSource(cells_source[local_rank])
            fieldS.setNature(IntensiveMaximum)
            WriteField("./source_field_part_%d.med" % rank, fieldS, True)

            if use_overlap_src:
                idec.attachLocalField(fieldS, numS)
            else:
                idec.attachLocalField(fieldS)

            # first send
            idec.synchronize()
            idec.sendData()

            self.checkValues(ref_src, fieldS, cells_source[local_rank])

            # FIXME
            # first recv
            # arr = fieldS.getArray()
            # arr *= 0.0
            # self.checkValues([0.0] * len(ref_src), fieldS, cells_source[local_rank])

            # idec.recvData()
            # self.checkValues(ref_src, fieldS, cells_source[local_rank])

            # second send - modify array but not local field
            arr = fieldS.getArray()
            arr *= 2.0
            self.checkValues(
                [2.0 * v for v in ref_src], fieldS, cells_source[local_rank]
            )
            idec.sendData()

        if target_group.containsMyRank():
            local_rank = target_group.myRank()
            mshT, fieldT, numT = self.getPartialTarget(cells_target[local_rank])
            fieldT.setNature(IntensiveMaximum)
            WriteUMesh("./target_mesh_part_%d.med" % rank, mshT, True)

            if use_overlap_trg:
                idec.attachLocalField(fieldT, numT)
            else:
                idec.attachLocalField(fieldT)
            idec.synchronize()

            # first recv
            idec.recvData()
            WriteField("./target_field_part_%d.med" % rank, fieldT, True)

            self.checkValues(ref_trg, fieldT, cells_target[local_rank])

            # first send
            # FIXME
            # idec.sendData()
            # self.checkValues(ref_trg, fieldT, cells_target[local_rank])

            # second recv
            idec.recvData()
            self.checkValues(
                [2.0 * v for v in ref_trg], fieldT, cells_target[local_rank]
            )

        # Release DEC (this involves MPI exchanges -- notably the release of the communicator -- so better be done before MPI.Finalize()
        idec.release()
        source_group.release()
        target_group.release()
        MPI.COMM_WORLD.Barrier()


if __name__ == "__main__":
    unittest.main()
    MPI.Finalize()
