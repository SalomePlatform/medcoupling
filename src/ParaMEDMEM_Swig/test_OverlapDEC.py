#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2025  CEA, EDF
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

class ParaMEDMEM_O_DEC_Tests(unittest.TestCase):
    """ This test illustrates a basic use of the OverlapDEC and shows notably that not all 
    processors must possess a piece of the source and/or target mesh. 
    Look at the C++ documentation of the class for more informations.
    In this case, the source mesh is only stored on 2 procs, whereas the target is on 4.
    Since only a single group of processor is defined in the setup, the 2 idle procs on the source side are just providing an empty mesh,
    thus indicating that they don't participate in the source definition. 
    
    Main method is testOverlapDEC_2D_py_1()
    """

    def generateFullSource(self, nb_compo=1):
        """ The complete source mesh: 4 squares each divided in 2 diagonaly (so 8 cells in total) """
        msh  = self.generateFullTarget()
        msh.simplexize(0)
        msh.setName("src_mesh")
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        fld.setMesh(msh); fld.setName("source_F")
        nc = msh.getNumberOfCells()
        da = DataArrayDouble(nc*nb_compo)
        da.rearrange(nb_compo)
        for c in range(nb_compo):
          da[:, c] = list(range(nc))
        da *= 2     # To compensate for the later division of the volume by 2 betw target and source cells.
        fld.setArray(da)
        return msh, fld

    def generateFullTarget(self):
        """ The complete target mesh: 4 squares """
        m1 = MEDCouplingCMesh("tgt_msh")
        da = DataArrayDouble([0,1,2])
        m1.setCoords(da, da)
        msh = m1.buildUnstructured()
        return msh

    #
    # Below, the two functions emulating the set up of a piece of the source and target mesh
    # on each proc. Obviously in real world problems, this comes from your code and is certainly
    # not computed by cutting again from scratch the full-size mesh!!
    #
    def getPartialSource(self, rank, nb_compo=1):
        """ Will return an empty mesh piece for rank=2 and 3 """
        msh, f = self.generateFullSource(nb_compo)
        if rank in [2,3]:
            sub_m, sub_f = msh[[]], f[[]]  # Little trick to select nothing in the mesh, thus producing an empty mesh
        elif rank == 0:
            sub_m, sub_f = msh[0:4], f[0:4]
        elif rank == 1:
            sub_m, sub_f = msh[4:8], f[4:8]
        sub_m.zipCoords()
        return sub_m, sub_f

    def getPartialTarget(self, rank, nb_compo=1):
        """ One square for each rank """
        msh = self.generateFullTarget()
        sub_m = msh[rank]
        sub_m.zipCoords()
        # Receiving side must prepare an empty field that will be filled by DEC:
        fld = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        da = DataArrayDouble(sub_m.getNumberOfCells(), nb_compo)
        fld.setArray(da)
        fld.setName("tgt_F")
        fld.setMesh(sub_m)
        return sub_m, fld

    def testOverlapDEC_ctor(self):
        """ Test the various Python ctors """
        size = MPI.COMM_WORLD.size
        if size != 4:
            print("Should be run on 4 procs!")
            return
        # Define processor group
        proc_group = list(range(size))
        # With 2 iterables:
        o1 = OverlapDEC.New(proc_group)
        # Should also work directly:
        o2 = OverlapDEC(proc_group)
        # With an iterable and a custom comm:
        o3 = OverlapDEC.New(proc_group, MPI.COMM_WORLD)
        # Also work directly with the **hack** on the comm:
        o4 = OverlapDEC(proc_group, MPI._addressof(MPI.COMM_WORLD))
        self.assertRaises(Exception, OverlapDEC, proc_group, MPI.COMM_WORLD)
        o4.release(); o3.release(); o2.release(); o1.release()

    @WriteInTmpDir
    def testOverlapDEC_2D_py_1(self):
        """ The main method of the test. """
        ncompo = 4   # Dummy field with 4 components
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 4:
            raise RuntimeError("Should be run on 4 procs!")

        for algo in range(3):
            # Define (single) processor group - note the difference with InterpKernelDEC which needs two groups.
            proc_group = list(range(size))   # No need for ProcessorGroup object here.
            odec = OverlapDEC(proc_group)
            odec.setWorkSharingAlgo(algo)    # Default is 1 - put here to test various different configurations of send/receive patterns

            # Write out full size meshes/fields for inspection
            if rank == 0:
                _, fld = self.generateFullSource(ncompo)
                mshT = self.generateFullTarget()
                WriteField("./source_field_FULL.med", fld, True)
                WriteUMesh("./target_mesh_FULL.med", mshT, True)

            MPI.COMM_WORLD.Barrier()  # really necessary??

            #
            # OK, let's go DEC !!
            #
            _, fieldS = self.getPartialSource(rank, ncompo)
            fieldS.setNature(IntensiveMaximum)   # The only policy supported for now ...
            mshT, fieldT = self.getPartialTarget(rank, ncompo)
            fieldT.setNature(IntensiveMaximum)
            if rank not in [2,3]:
                WriteField("./source_field_part_%d.med" % rank, fieldS, True)
            WriteUMesh("./target_mesh_part_%d.med" % rank, mshT, True)

            odec.attachSourceLocalField(fieldS)
            odec.attachTargetLocalField(fieldT)
            odec.synchronize()
            odec.sendRecvData()

            # Now the actual checks:
            ref_vals = [1.0, 5.0, 9.0, 13.0]
            self.assertEqual(fieldT.getArray().getValues(), [ref_vals[rank]]*ncompo)

            # Release DEC (this involves MPI exchanges -- notably the release of the communicator -- so better be done before MPI.Finalize()
            odec.release()

            MPI.COMM_WORLD.Barrier()

    @WriteInTmpDir
    def testOverlapDEC_2D_py_2(self):
        """ Test on a question that often comes back: should I re-synchronize() / re-attach() each time that
        I want to send a new field? 
        Basic answer: 
          - you do not have to re-synchronize, but you can re-attach a new field, as long as it has the same support.
        WARNING: this differs in InterpKernelDEC ...
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 4:
            raise RuntimeError("Should be run on 4 procs!")

        # Define (single) processor group - note the difference with InterpKernelDEC which needs two groups.
        proc_group = list(range(size))   # No need for ProcessorGroup object here.
        odec = OverlapDEC(proc_group)

        MPI.COMM_WORLD.Barrier()  # really necessary??

        #
        # OK, let's go DEC !!
        #
        _, fieldS = self.getPartialSource(rank)
        fieldS.setNature(IntensiveMaximum)   # The only policy supported for now ...
        mshT, fieldT = self.getPartialTarget(rank)
        fieldT.setNature(IntensiveMaximum)

        mul = 1
        for t in range(3):  # Emulating a time loop ...
            if t == 0:
                odec.attachSourceLocalField(fieldS)   # only once!
                odec.attachTargetLocalField(fieldT)   # only once!
                odec.synchronize()                    # only once!
            else:
                das = fieldS.getArray()               # but we can still hack the underlying field values ...
                das *= 2
                mul *= 2
            odec.sendRecvData()                       # each time!

            # Now the actual checks:
            ref_vals = [1.0, 5.0, 9.0, 13.0]
            self.assertEqual(fieldT.getArray().getValues(), [ref_vals[rank]*mul])

        # Release DEC (this involves MPI exchanges -- notably the release of the communicator -- so better be done before MPI.Finalize()
        odec.release()

        MPI.COMM_WORLD.Barrier()

if __name__ == "__main__":
    unittest.main()
    MPI.Finalize()
    # tt = ParaMEDMEM_O_DEC_Tests()
    # tt.testOverlapDEC_2D_py_1()
