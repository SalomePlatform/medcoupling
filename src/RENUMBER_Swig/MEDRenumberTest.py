# Copyright (C) 2012-2016  CEA/DEN, EDF R&D
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

from MEDRenumber import *
import unittest

class MEDRenumberTest(unittest.TestCase):

    @unittest.skipUnless("BOOST" in RenumberAvailableMethods(),"requires BOOST prerequisite !")
    def test1(self):
        from MEDCoupling import MEDCouplingCMesh
        ren=RenumberingFactory("BOOST")
        arr=DataArrayDouble(10) ; arr.iota()
        c=MEDCouplingCMesh() ; c.setCoords(arr,arr)
        m=c.buildUnstructured()
        a,b=m.computeNeighborsOfCells()
        n2o,o2n=ren.renumber(a,b)
        self.assertTrue(o2n.isEqual(DataArrayInt([0,2,5,9,14,20,27,35,44,1,4,8,13,19,26,34,43,52,3,7,12,18,25,33,42,51,59,6,11,17,24,32,41,50,58,65,10,16,23,31,40,49,57,64,70,15,22,30,39,48,56,63,69,74,21,29,38,47,55,62,68,73,77,28,37,46,54,61,67,72,76,79,36,45,53,60,66,71,75,78,80])))
        pass

    def setUp(self):
        pass
    pass

unittest.main()
