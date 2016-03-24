#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

from MEDLoader import *
import os
import sys
import unittest

class RenumberingTest(unittest.TestCase):
    def testBoost2D(self):
        filename="Test2D.med"
        meshname="Mesh_1"
        method="BOOST"
        string_to_execute=self.dir_renumber+" "+self.dir_mesh+"/"+filename+" "+meshname+" "+method+" "+self.dir_mesh+"/out_"+filename
        os.system(string_to_execute)
        mm=MEDFileMesh.New(self.dir_mesh+"/out_"+filename,meshname)
        m=mm.getMeshAtLevel(0)
        ff=MEDFileField1TS(self.dir_mesh+"/out_"+filename,"Test field")
        field_ini=DataArrayDouble([(2,3),(12,13),(14,15),(4,5),(6,7),(8,9),(16,17),(0,1),(10,11)])
        ff.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
        f=ff.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
        field=f.getArray().isEqual(field_ini,1e-15)
        connectivite=[4,1,5,12,10,4,10,12,13,11,4,5,4,14,12,4,11,13,9,3,4,12,14,15,13,4,4,0,6,14,4,13,15,8,9,4,14,6,7,15,4,15,7,2,8]
        connectivite_index=[0,5,10,15,20,25,30,35,40,45]
        Boost2D=m.getNodalConnectivity().getValues()==connectivite and m.getNodalConnectivityIndex().getValues()==connectivite_index and field
        self.assertTrue(Boost2D)
        os.remove(self.dir_mesh+"/out_"+filename)
        pass

    def tessMetis2D(self):#not activated yet
        filename="Test2D.med"
        meshname="Mesh_1"
        method="METIS"
        string_to_execute=self.dir_renumber+" "+self.dir_mesh+"/"+filename+" "+meshname+" "+method+" "+self.dir_mesh+"/out_"+filename
        os.system(string_to_execute)
        m = MESH(MED_DRIVER,self.dir_mesh+"/out_"+filename,meshname)
        nbcell2dmetis=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
        connectivite=[12,14,10,4,2,6,13,11,11,13,14,12,16,8,3,9,5,1,7,15,15,7,8,16,14,16,9,10,6,5,15,13,13,15,16,14]
        connectivite_index=[1,5,9,13,17,21,25,29,33,37]
        conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_QUAD4)
        conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
        conn2dmetis=(list(conn)==connectivite)
        conn_index2dmetis=(list(conn_index)==connectivite_index)
        Metis2D=conn2dmetis and conn_index2dmetis and (nbcell2dmetis==9)
        self.assertTrue(Metis2D)
        os.remove(self.dir_mesh+"/out_"+filename)
        pass

    def testBoost2DPolygon(self):
        filename="Test2Dpoly.med"
        meshname="Mesh_1"
        method="BOOST"
        string_to_execute=self.dir_renumber+" "+self.dir_mesh+"/"+filename+" "+meshname+" "+method+" "+self.dir_mesh+"/out_"+filename
        os.system(string_to_execute)
        mm=MEDFileMesh.New(self.dir_mesh+"/out_"+filename,meshname)
        m=mm.getMeshAtLevel(0)
        nbcell2dpolyboost=m.getNumberOfCells()
        connectivite=[5,1,4,8,9,5,10,9,8,11,5,4,5,7,8,5,3,10,11,15,5,11,8,7,12,5,5,0,6,7,5,15,11,12,14,5,12,7,6,13,5,14,12,13,2]
        connectivite_index=[0,5,10,15,20,25,30,35,40,45]
        conn=m.getNodalConnectivity().getValues()
        conn_index=m.getNodalConnectivityIndex().getValues()
        conn2dpolyboost=(list(conn)==connectivite)
        conn_index2dpolyboost=(list(conn_index)==connectivite_index)
        PolyBoost2D=conn2dpolyboost and conn_index2dpolyboost and (nbcell2dpolyboost==9)
        self.assertTrue(PolyBoost2D)
        os.remove(self.dir_mesh+"/out_"+filename)
        pass

    def tessMetis2DPolygon(self):#not activated yet
        filename="Test2Dpoly.med"
        meshname="Mesh_1"
        method="METIS"
        string_to_execute=self.dir_renumber+" "+self.dir_mesh+"/"+filename+" "+meshname+" "+method+" "+self.dir_mesh+"/out_"+filename
        os.system(string_to_execute)
        m = MESH(MED_DRIVER,self.dir_mesh+"/out_"+filename,meshname)
        nbcell2dpolymetis=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
        connectivite=[6,1,7,8,2,5,9,10,5,6,8,9,15,13,14,3,4,11,12,16,16,12,13,15,11,10,9,12,12,9,8,13,13,8,7,14]
        connectivite_index=[1,5,9,13,17,21,25,29,33,37]
        conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_POLYGON)
        conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL)
        conn2dpolymetis=(list(conn)==connectivite)
        conn_index2dpolymetis=(list(conn_index)==connectivite_index)
        PolyMetis2D=conn2dpolymetis and conn_index2dpolymetis and (nbcell2dpolymetis==9)
        self.assertTrue(PolyMetis2D)
        os.remove(self.dir_mesh+"/out_"+filename)
        pass

    def testBoost3D(self):
        filename="Test3D.med"
        meshname="Mesh_1"
        method="BOOST"
        string_to_execute=self.dir_renumber+" "+self.dir_mesh+"/"+filename+" "+meshname+" "+method+" "+self.dir_mesh+"/out_"+filename
        os.system(string_to_execute)
        mm=MEDFileMesh.New(self.dir_mesh+"/out_"+filename,meshname)
        m=mm.getMeshAtLevel(0)
        nbcell3dboost=m.getNumberOfCells()
        connectivite=[18,22,12,4,17,26,21,13,25,18,16,5,12,22,24,15,21,26,18,26,21,13,25,23,14,6,19,18,8,22,17,0,20,26,25,9,18,24,15,21,26,18,7,14,23,18,1,16,22,8,11,24,26,20,18,20,26,25,9,10,23,19,2,18,11,24,26,20,3,18,23,10]
        connectivite_index=[0,9,18,27,36,45,54,63,72]
        conn=m.getNodalConnectivity().getValues()
        conn_index=m.getNodalConnectivityIndex().getValues()
        conn3dboost=(list(conn)==connectivite)
        conn_index3dboost=(list(conn_index)==connectivite_index)
        Boost3D=conn3dboost and conn_index3dboost and (nbcell3dboost==8)
        self.assertTrue(Boost3D)
        os.remove(self.dir_mesh+"/out_"+filename)
        pass

    def tessMetis3D(self):#not activated yet
        filename="Test3D.med"
        meshname="Mesh_1"
        method="METIS"
        string_to_execute=self.dir_renumber+" "+self.dir_mesh+"/"+filename+" "+meshname+" "+method+" "+self.dir_mesh+"/out_"+filename
        os.system(string_to_execute)
        m = MESH(MED_DRIVER,self.dir_mesh+"/out_"+filename,meshname)
        nbcell3dmetis=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
        connectivite=[12,25,27,21,4,19,24,11,27,22,14,26,24,15,7,20,17,6,13,23,25,16,22,27,9,23,18,1,21,27,26,10,23,13,5,18,27,22,14,26,25,16,22,27,19,8,15,24,2,17,23,9,12,25,27,21,21,27,26,10,11,24,20,3]
        connectivite_index=[1,9,17,25,33,41,49,57,65]
        conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_HEXA8)
        conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
        conn3dmetis=(list(conn)==connectivite)
        conn_index3dmetis=(list(conn_index)==connectivite_index)
        Metis3D=conn3dmetis&conn_index3dmetis&(nbcell3dmetis==8)
        self.assertTrue(Metis3D)
        os.remove(self.dir_mesh+"/out_"+filename)
        pass

    def testBoost3DPoly(self):
        filename="Test3Dpoly.med"
        meshname="Mesh_1"
        method="BOOST"
        string_to_execute=self.dir_renumber+" "+self.dir_mesh+"/"+filename+" "+meshname+" "+method+" "+self.dir_mesh+"/out_"+filename
        os.system(string_to_execute)
        mm=MEDFileMesh.New(self.dir_mesh+"/out_"+filename,meshname)
        m=mm.getMeshAtLevel(0)
        nbcell3dpolyboost=m.getNumberOfCells()
        connectivite=[31,22,12,4,17,-1,26,25,13,21,-1,22,26,21,12,-1,12,21,13,4,-1,4,13,25,17,-1,17,25,26,22,31,16,5,12,22,-1,24,26,21,15,-1,16,24,15,5,-1,5,15,21,12,-1,12,21,26,22,-1,22,26,24,16,31,26,21,13,25,-1,23,19,6,14,-1,26,23,14,21,-1,21,14,6,13,-1,13,6,19,25,-1,25,19,23,26,31,8,22,17,0,-1,20,9,25,26,-1,8,20,26,22,-1,22,26,25,17,-1,17,25,9,0,-1,0,9,20,8,31,24,15,21,26,-1,18,23,14,7,-1,24,18,7,15,-1,15,7,14,21,-1,21,14,23,26,-1,26,23,18,24,31,1,16,22,8,-1,11,20,26,24,-1,1,11,24,16,-1,16,24,26,22,-1,22,26,20,8,-1,8,20,11,1,31,20,26,25,9,-1,10,2,19,23,-1,20,10,23,26,-1,26,23,19,25,-1,25,19,2,9,-1,9,2,10,20,31,11,24,26,20,-1,3,10,23,18,-1,11,3,18,24,-1,24,18,23,26,-1,26,23,10,20,-1,20,10,3,11]
        connectivite_index=[0,30,60,90,120,150,180,210,240]
        conn=m.getNodalConnectivity().getValues()
        conn_index=m.getNodalConnectivityIndex().getValues()
        conn3dpolyboost=(connectivite==list(conn))
        conn_index3dpolyboost=(connectivite_index==list(conn_index))
        PolyBoost3D=(conn3dpolyboost and conn_index3dpolyboost and (nbcell3dpolyboost==8))
        self.assertTrue(PolyBoost3D)
        os.remove(self.dir_mesh+"/out_"+filename)
        pass

    def tessBoost3DPoly(self):#not activated yet
        filename="Test3Dpoly.med"
        meshname="Mesh_1"
        method="METIS"
        string_to_execute=self.dir_renumber+" "+self.dir_mesh+"/"+filename+" "+meshname+" "+method+" "+self.dir_mesh+"/out_"+filename
        os.system(string_to_execute)
        m = MESH(MED_DRIVER,self.dir_mesh+"/out_"+filename,meshname)
        nbcell3dpolymetis=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
        connectivite=[12,25,27,21,-1,4,11,24,19,-1,12,4,19,25,-1,25,19,24,27,-1,27,24,11,21,-1,21,11,4,12,
                      27,22,14,26,-1,24,20,7,15,-1,27,24,15,22,-1,22,15,7,14,-1,14,7,20,26,-1,26,20,24,27,
                      17,6,13,23,-1,25,27,22,16,-1,17,25,16,6,-1,6,16,22,13,-1,13,22,27,23,-1,23,27,25,17,
                      9,23,18,1,-1,21,10,26,27,-1,9,21,27,23,-1,23,27,26,18,-1,18,26,10,1,-1,1,10,21,9,
                      23,13,5,18,-1,27,26,14,22,-1,23,27,22,13,-1,13,22,14,5,-1,5,14,26,18,-1,18,26,27,23,
                      25,16,22,27,-1,19,24,15,8,-1,25,19,8,16,-1,16,8,15,22,-1,22,15,24,27,-1,27,24,19,25,
                      2,17,23,9,-1,12,21,27,25,-1,2,12,25,17,-1,17,25,27,23,-1,23,27,21,9,-1,9,21,12,2,
                      21,27,26,10,-1,11,3,20,24,-1,21,11,24,27,-1,27,24,20,26,-1,26,20,3,10,-1,10,3,11,21]
        connectivite_index=[1, 30, 59, 88, 117, 146, 175, 204, 233]
        conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_POLYHEDRA)
        conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
        conn3dpolymetis=(list(conn)==connectivite)
        conn_index3dpolymetis=(list(conn_index)==connectivite_index)
        PolyMetis3D=(conn3dpolymetis and conn_index3dpolymetis and (nbcell3dpolymetis==8))
        self.assertTrue(PolyMetis3D)
        os.remove(self.dir_mesh+"/out_"+filename)
        pass

    def setUp(self):
        med_root_dir=os.getenv("MEDCOUPLING_ROOT_DIR")
        self.dir_renumber=os.path.join(med_root_dir, "bin/renumber")
        self.dir_mesh=os.path.join(med_root_dir, "share","resources","med")
        pass
    pass

unittest.main()
