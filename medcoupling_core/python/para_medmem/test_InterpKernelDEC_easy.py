#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2023  CEA, EDF
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
#from ParaMEDMEMTestTools import WriteInTmpDir
import sys, os
import unittest
import math
from mpi4py import MPI

def ranksByGroup(groupString, jobPerWorldRank):
    ranks=[]
    for key,value in jobPerWorldRank.items():
        if (groupString == value ):
            ranks.append(key)
    return ranks


class ParaMEDMEM_IK_DEC_Tests(unittest.TestCase):
    def test_InterpKernelDEC_easy_comm_creation(self):
        """
        [EDF26706] :
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 5:
            print("Should be run on 5 procs!")
            return
        jobPerWorldRank = {0:"A",1:"B",2:"B",3:"C",4:"C"}
        interface = CommInterface()
        group = ByStringMPIProcessorGroup(interface, jobPerWorldRank[rank])
        decBC = InterpKernelDEC(group,"B<->C")
        decAC = InterpKernelDEC(group,"A<->C")
        eval("Easy_comm_creation_{}".format(rank))(decBC,decAC)
        #
        MPI.COMM_WORLD.Barrier() 
    
    def test_InterpKernelDEC_easy_comm_creation_2(self):
        """
        [EDF26706] :
        """
        size = MPI.COMM_WORLD.size
        rank = MPI.COMM_WORLD.rank
        if size != 5:
            print("Should be run on 5 procs!")
            return
        jobPerWorldRank = {0:"A",1:"B",2:"B",3:"C",4:"C"}
        interface = CommInterface()
        group = ByStringMPIProcessorGroup(interface, jobPerWorldRank[rank])
        decBC = InterpKernelDEC(group,"B","C")
        decAC = InterpKernelDEC(group,"A","C")
        eval("Easy_comm_creation_{}".format(rank))(decBC,decAC)
        #
        MPI.COMM_WORLD.Barrier()    

def Easy_comm_creation_0(decBC,decAC):
    """ Proc 0 of A"""
    m = MEDCouplingCMesh() ; m.setCoords(DataArrayDouble([0,1]),DataArrayDouble([0,1])) ; m = m.buildUnstructured()
    field = MEDCouplingFieldDouble(ON_CELLS)
    field.setNature(IntensiveMaximum)
    field.setMesh( m )
    field.setArray( DataArrayDouble([1.2]))
    decAC.attachLocalField( field )
    decAC.synchronize()
    decAC.sendData()
    pass

def Easy_comm_creation_1(decBC,decAC):
    """ Proc 0 of B"""
    m = MEDCouplingCMesh() ; m.setCoords(DataArrayDouble([2,3]),DataArrayDouble([1,2])) ; m = m.buildUnstructured()
    field = MEDCouplingFieldDouble(ON_CELLS)
    field.setNature(IntensiveMaximum)
    field.setMesh( m )
    field.setArray( DataArrayDouble([2.3]))
    decBC.attachLocalField( field )
    decBC.synchronize()
    decBC.sendData()
    pass

def Easy_comm_creation_2(decBC,decAC):
    """ Proc 1 of B"""
    m = MEDCouplingCMesh() ; m.setCoords(DataArrayDouble([3,4]),DataArrayDouble([1,2])) ; m = m.buildUnstructured()
    field = MEDCouplingFieldDouble(ON_CELLS)
    field.setNature(IntensiveMaximum)
    field.setMesh( m )
    field.setArray( DataArrayDouble([3.3]))
    decBC.attachLocalField( field )
    decBC.synchronize()
    decBC.sendData()
    pass

def Easy_comm_creation_3(decBC,decAC):
    """ Proc 0 of C"""
    m = MEDCouplingCMesh() ; m.setCoords(DataArrayDouble([0.5,3.5]),DataArrayDouble([0,1.5])) ; m = m.buildUnstructured()
    field = MEDCouplingFieldDouble(ON_CELLS)
    field.setNature(IntensiveMaximum)
    field.setMesh( m )
    field.setArray( DataArrayDouble([0.]))
    decBC.attachLocalField( field )
    decAC.attachLocalField( field )
    decBC.synchronize()
    decAC.synchronize()
    decBC.recvData()
    print(field.getArray().getValues())
    decAC.recvData()
    print(field.getArray().getValues())
    pass

def Easy_comm_creation_4(decBC,decAC):
    """ Proc 1 of C"""
    m = MEDCouplingCMesh() ; m.setCoords(DataArrayDouble([0.7,3.5]),DataArrayDouble([0,1.5])) ; m = m.buildUnstructured()
    field = MEDCouplingFieldDouble(ON_CELLS)
    field.setNature(IntensiveMaximum)
    field.setMesh( m )
    field.setArray( DataArrayDouble([0.]))
    decBC.attachLocalField( field )
    decAC.attachLocalField( field )
    decBC.synchronize()
    decAC.synchronize()
    decBC.recvData()
    print(field.getArray().getValues())
    decAC.recvData()
    print(field.getArray().getValues())
    pass

if __name__ == "__main__":
    unittest.main()
    MPI.Finalize()

