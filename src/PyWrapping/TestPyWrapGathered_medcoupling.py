#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2017  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (EDF R&D)

from medcoupling import *
import unittest

class FileCreator(object):
    def __init__(self,tester,fname):
        self._tester=tester
        self._fname=fname
        pass
        
    def fileName(self):
        return self._fname
    
    def __enter__(self):
        import os
        if os.path.exists(self._fname):
            os.remove(self._fname)
            pass
        return self
    
    def __exit__(self, type, value, traceback):
        import os
        if not os.path.exists(self._fname):
            self._tester.assertTrue(False)
            pass
        else:
            os.remove(self._fname)
        pass
        
class medcouplingTest(unittest.TestCase):
    
    @unittest.skipUnless(HasMEDFileExt(),"Requires link to MED file")
    def test1(self):
        import sys
        fname="mctest1.med"
        arr=DataArrayDouble(10) ; arr.iota()
        m=MEDCouplingCMesh()
        m.setCoords(arr,arr)
        m.setName("mesh")
        with FileCreator(self,fname) as fc:
            m.write(fc.fileName())
        m=m.buildUnstructured()
        with FileCreator(self,fname) as fc:
            m.write(fc.fileName())
        f=MEDCouplingFieldDouble(ON_NODES) ; f.setMesh(m) ; f.setArray(m.getCoords()) ; f.setName("field")
        with FileCreator(self,fname) as fc:
            f.write(fc.fileName())
        f=MEDCouplingFieldFloat(ON_NODES) ; f.setMesh(m)
        d=DataArrayFloat(m.getNumberOfNodes()) ; d.iota()
        f.setArray(d) ; f.setName("field1")
        with FileCreator(self,fname) as fc:
            f.write(fc.fileName())
        pass

    @unittest.skipUnless(HasRenumberExt(),"Requires Boost or Metis to activate Renumberer")
    def test2(self):
        arr=DataArrayDouble(10) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured() ; m.setName("mesh")
        #
        renf=RenumberingFactory("Boost")
        neigh,neighi=m.computeNeighborsOfCells()
        n2o,o2n=renf.renumber(neigh,neighi)
        mRenum=m[n2o]
        pass

    @unittest.skipUnless(HasPartitionerExt(),"Requires Partitioner activation")
    def test3(self):
        for alg in MEDPartitioner.AvailableAlgorithms():
            st="Graph.%s"%alg.upper()
            print(st)
            self.partitionerTesterHelper(eval(st))
            pass
        pass
    
    @unittest.skipUnless(HasParallelInterpolatorExt(),"Requires // interpolator activated")
    def test4(self):
        interface=CommInterface()
        pass

    def partitionerTesterHelper(self,algoSelected):
        arr=DataArrayDouble(10) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured() ; m.setName("mesh")
        a,b=m.computeNeighborsOfCells()
        sk=MEDCouplingSkyLineArray(b,a)
        g=MEDPartitioner.Graph(sk,algoSelected)
        g.partGraph(4)
        procIdOnCells=g.getPartition().getValuesArray()
        m0=m[procIdOnCells.findIdsEqual(0)] ; m0.setName("m0")
        pass
    
    pass

if __name__ == "__main__":
    if HasParallelInterpolatorExt():
        try:
            from mpi4py import MPI # if not imported test3 may failed due to MPI call of partitioner algorithms.
        except:
            pass
        pass
    unittest.main()
    pass

