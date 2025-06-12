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

import unittest


class ICoCoICoCoMEDDoubleFieldTest(unittest.TestCase):
    def generate_fields_double(self):
        """Dummy MCFieldDouble"""
        msh = MEDCouplingCMesh("toto_mesh")
        msh.setCoords(DataArrayDouble([0.0, 1.0, 2.0]))
        msh = msh.buildUnstructured()
        f1 = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        f1.setMesh(msh)
        f1.setName("toto")
        f1.setArray(DataArrayDouble([0.0, 1.0, 2.0, 3.0]))

        f2 = f1.deepCopy()
        da = f2.getArray()
        da += 3.5
        da2 = f2.getArray()
        return f1, f2

    def generate_fields_int(self):
        """Dummy MCFieldDouble"""
        msh = MEDCouplingCMesh("toto_mesh")
        msh.setCoords(DataArrayDouble([0.0, 1.0, 2.0]))
        msh = msh.buildUnstructured()
        f1 = MEDCouplingFieldInt32(ON_CELLS, ONE_TIME)
        f1.setMesh(msh)
        f1.setName("toto")
        f1.setArray(DataArrayInt32([0, 1, 2, 3]))

        f2 = f1.deepCopy()
        da = f2.getArray()
        da += 3
        da2 = f2.getArray()
        return f1, f2

    def test1(self):
        lst_typ = [ICoCoMEDDoubleField, ICoCoMEDIntField]
        f1d, f2d = self.generate_fields_double()
        f1i, f2i = self.generate_fields_int()
        fld1_lst = [f1d, f1i]
        fld2_lst = [f2d, f2i]
        for ICoCoMED_T_Field, f1, f2 in zip(lst_typ, fld1_lst, fld2_lst):
            mf = ICoCoMED_T_Field()
            mf.setName("titi")
            self.assertEqual(mf.getName(), "titi")
            mfd = mf.getMCField()
            self.assertTrue(mfd is None)
            mf.setMCField(f1)
            f11 = mf.getMCField()
            self.assertEqual(
                f1.getHiddenCppPointer(), f11.getHiddenCppPointer()
            )  # strictly the same
            self.assertEqual(mf.getName(), "toto")  # name is taken from MC object
            mf.setMCField(f2)
            f22 = mf.getMCField()
            self.assertEqual(
                f2.getHiddenCppPointer(), f22.getHiddenCppPointer()
            )  # strictly the same

            mf = ICoCoMED_T_Field(f1)  # ctor with MC object
            mfd = mf.getMCField()
            self.assertEqual(
                mfd.getHiddenCppPointer(), f1.getHiddenCppPointer()
            )  # strictly the same
            self.assertEqual(mf.getName(), "toto")  # name is taken from MC object

            mf.setMCField(None)
            mfd = mf.getMCField()
            self.assertTrue(mfd is None)
            self.assertEqual(mf.getName(), "")  # name is reset

            mf.setMCField(f2)
            f22 = mf.getMCField()
            self.assertEqual(
                f2.getHiddenCppPointer(), f22.getHiddenCppPointer()
            )  # strictly the same

            mf.setName("aa")
            mf2 = ICoCoMED_T_Field(mf)  # copy ctor
            f22 = mf2.getMCField()
            self.assertEqual(
                f2.getHiddenCppPointer(), f22.getHiddenCppPointer()
            )  # strictly the same
            self.assertEqual(mf2.getName(), "aa")

            mf2 = mf  # assignement op
            f22 = mf2.getMCField()
            self.assertEqual(
                f2.getHiddenCppPointer(), f22.getHiddenCppPointer()
            )  # strictly the same
            self.assertEqual(mf2.getName(), "aa")

            mf2.setMCField(None)
            self.assertEqual(mf2.getName(), "")


if __name__ == "__main__":
    unittest.main()
