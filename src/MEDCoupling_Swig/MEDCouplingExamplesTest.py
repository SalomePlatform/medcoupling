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

from MEDCoupling import *
import unittest
from math import pi, sqrt

class MEDCouplingBasicsTest(unittest.TestCase):

    def testExample_MEDCouplingFieldDouble_WriteVTK(self):
        #! [PySnippet_MEDCouplingFieldDouble_WriteVTK_1]
        # mesh
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr) # mesh becomes a 2D one

        # 3 fields (lying on the same mesh!)
        field1 = mesh.getMeasureField( True )
        field2 = mesh.buildOrthogonalField()
        field3 = mesh.fillFromAnalytic( ON_CELLS, 2, "IVec * x + JVec * y" )
        field2.setName( "Normal" ) # name is necessary!
        field3.setName( "Barycenter" ) # name is necessary!

        # WriteVTK
        fileName = "testExample_MEDCouplingFieldDouble_WriteVTK"
        fs = [ field1, field2, field3 ] # field series
        writtenFileName=MEDCouplingFieldDouble.WriteVTK( fileName, fs )
        print("The file name with correct extension is : %s"%(writtenFileName))
        #! [PySnippet_MEDCouplingFieldDouble_WriteVTK_1]
        import os
        os.remove( writtenFileName )

        return

    def testExample_MEDCouplingFieldDouble_MaxFields(self):
        #! [PySnippet_MEDCouplingFieldDouble_MaxFields_1]
        vals1   = [0.,2., 4.,6.] # for field 1
        vals2   = [2.,0., 6.,4.] # for field 2
        valsMax = [2.,2., 6.,6.] # expected max field
        valsMin = [0.,0., 4.,4.] # expected min field

        # field 1
        valsArr1=DataArrayDouble(vals1,2,2) # 2 tuples per 2 components
        field1 = MEDCouplingFieldDouble( ON_NODES )
        field1.setArray( valsArr1 )

        # field 2
        valsArr2=DataArrayDouble(vals2,2,2) # 2 tuples per 2 components
        field2 = MEDCouplingFieldDouble( ON_NODES )
        field2.setArray( valsArr2 )

        # max field 
        fieldMax = MEDCouplingFieldDouble.MaxFields( field1, field2 )
        self.assertTrue( fieldMax.getArray().getValues() == valsMax )

        # min field 
        fieldMin = MEDCouplingFieldDouble.MinFields( field1, field2 )
        self.assertTrue( fieldMin.getArray().getValues() == valsMin )
        #! [PySnippet_MEDCouplingFieldDouble_MaxFields_1]

    def testExample_MEDCouplingFieldDouble_MergeFields(self):
        #! [PySnippet_MEDCouplingFieldDouble_MergeFields_1]
        # mesh 1
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh1=MEDCouplingCMesh()
        mesh1.setCoords(coordsArr)
        # field 1
        field1 = mesh1.fillFromAnalytic( ON_CELLS, 1, "x")

        # mesh 2 and field 2
        field2 = field1.cloneWithMesh( True )
        vec = [5.]
        field2.getMesh().translate(vec) # translate mesh2
        field2.applyFunc("x + 5") # "translate" field2

        # concatenate field1 and field2
        field3 = MEDCouplingFieldDouble.MergeFields( field1, field2 )
        field4 = MEDCouplingFieldDouble.MergeFields( [ field1, field2] )
        #! [PySnippet_MEDCouplingFieldDouble_MergeFields_1]
        return

    def testExample_MEDCouplingFieldDouble_substractInPlaceDM(self):
        #! [PySnippet_MEDCouplingFieldDouble_substractInPlaceDM_1]
        coords1=[0.,1.,2.,3.]
        coords2=[2.,1.,0.,3.] #0 <==> #2
        # mesh 1
        mesh1=MEDCouplingUMesh()
        coordsArr=DataArrayDouble(coords1, 4, 1)
        mesh1.setCoords(coordsArr)
        mesh1.setMeshDimension(0)
        mesh1.allocateCells(0)
        mesh1.finishInsertingCells()
        # mesh 2
        mesh2=mesh1.deepCopy()
        mesh2.getCoords().setValues(coords2, 4, 1)
        #! [PySnippet_MEDCouplingFieldDouble_substractInPlaceDM_1]
        #! [PySnippet_MEDCouplingFieldDouble_substractInPlaceDM_2]
        field1 = mesh1.fillFromAnalytic(ON_NODES,1,"x") # field1 values == coords1
        field2 = mesh2.fillFromAnalytic(ON_NODES,1,"x") # field2 values == coords2
        levOfCheck = 10 # nodes can be permuted
        field1.substractInPlaceDM( field2, levOfCheck, 1e-13, 0 ) # values #0 and #2 must swap
        #! [PySnippet_MEDCouplingFieldDouble_substractInPlaceDM_2]
        #! [PySnippet_MEDCouplingFieldDouble_substractInPlaceDM_3]
        field2.applyFunc( 1, 0.0 ) # all field2 values == 0.0
        self.assertTrue( field1.isEqual( field2, 1e-13, 1e-13 )) # field1 == field2 == 0.0
        #! [PySnippet_MEDCouplingFieldDouble_substractInPlaceDM_3]
        return

    def testExample_MEDCouplingFieldDouble_changeUnderlyingMesh(self):
        #! [PySnippet_MEDCouplingFieldDouble_changeUnderlyingMesh_1]
        coords1=[0.,1.,2.,3.]
        coords2=[2.,1.,0.,3.] #0 <==> #2
        # mesh 1
        mesh1=MEDCouplingUMesh()
        coordsArr=DataArrayDouble(coords1, 4, 1)
        mesh1.setCoords(coordsArr)
        mesh1.setMeshDimension(0)
        mesh1.allocateCells(0)
        mesh1.finishInsertingCells()
        # mesh 2
        mesh2=mesh1.deepCopy()
        mesh2.getCoords().setValues(coords2, 4, 1)
        #! [PySnippet_MEDCouplingFieldDouble_changeUnderlyingMesh_1]
        #! [PySnippet_MEDCouplingFieldDouble_changeUnderlyingMesh_2]
        field = mesh1.fillFromAnalytic(ON_NODES,1,"x") # field values == coords1
        levOfCheck = 10 # nodes can be permuted
        field.changeUnderlyingMesh( mesh2, levOfCheck, 1e-13, 0 ) # values #0 and #2 must swap
        self.assertTrue( field.getArray().getValues() == coords2 )
        #! [PySnippet_MEDCouplingFieldDouble_changeUnderlyingMesh_2]
        return

    def testExample_MEDCouplingFieldDouble_applyFunc_same_nb_comp(self):
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_same_nb_comp_1]
        v = [1.,2., 3.,4.]
        array = DataArrayDouble( v, 2, 2 ) # 2 tuples per 2 components
        field = MEDCouplingFieldDouble( ON_CELLS )
        field.setArray( array )
        func = "IVec * v + JVec * w*w + 10"
        field.applyFunc( 2, func )
        self.assertTrue( field.getNumberOfComponents() == 2 ) # 2 components remains
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_same_nb_comp_1]
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_same_nb_comp_2]
        v2 = field.getArray().getValues()
        self.assertAlmostEqual( v2[0], 10 + v[0], 13 )      # "10 + IVec * v"
        self.assertAlmostEqual( v2[1], 10 + v[1]*v[1], 13 ) # "10 + JVec * v*v"
        self.assertAlmostEqual( v2[2], 10 + v[2], 13 )      # "10 + IVec * v"
        self.assertAlmostEqual( v2[3], 10 + v[3]*v[3], 13 ) # "10 + JVec * v*v"
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_same_nb_comp_2]
        return

    def testExample_MEDCouplingFieldDouble_applyFunc3(self):
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc3_1]
        # create a 2D vector field
        values = [1.,1., 2.,1.]
        array = DataArrayDouble( values, 2, 2 ) # 2 tuples per 2 components
        field = MEDCouplingFieldDouble( ON_CELLS )
        field.setArray( array )
        # transform the field to a 3D vector field
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        varNames=["a","b"] # names used to refer to X and Y components
        field.applyFuncNamedCompo( 3, varNames, func ) # require 3 components 
        self.assertTrue( field.getNumberOfComponents() == 3 ) # 3 components as required
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc3_1]
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc3_2]
        vec1 = field.getArray().getTuple(1) # vector #1
        a,b = values[2], values[3] # initial components of the vector #1
        self.assertAlmostEqual( vec1[0], 10 + b, 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vec1[1], 10 + a, 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vec1[2], 10 + sqrt(a*a+b*b), 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc3_2]
        return

    def testExample_MEDCouplingFieldDouble_applyFunc2(self):
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc2_1]
        # create a 2D vector field
        values = [1.,1., 2.,1.]
        array = DataArrayDouble( values, 2, 2 ) # 2 tuples per 2 components
        array.setInfoOnComponent(0,"a") # name used to refer to X component within a function
        array.setInfoOnComponent(1,"b") # name used to refer to Y component within a function
        field = MEDCouplingFieldDouble( ON_CELLS )
        field.setArray( array )
        # transform the field to a 3D vector field
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        field.applyFuncCompo( 3, func ) # require 3 components 
        self.assertTrue( field.getNumberOfComponents() == 3 ) # 3 components as required
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc2_1]
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc2_2]
        vec1 = field.getArray().getTuple(1) # vector #1
        a,b = values[2], values[3] # initial components of the vector #1
        self.assertAlmostEqual( vec1[0], 10 + b, 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vec1[1], 10 + a, 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vec1[2], 10 + sqrt(a*a+b*b), 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc2_2]
        return

    def testExample_MEDCouplingFieldDouble_applyFunc(self):
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_1]
        # create a 2D vector field
        values = [1.,1., 2.,1.]
        array = DataArrayDouble( values, 2, 2 ) # 2 tuples per 2 components
        field = MEDCouplingFieldDouble( ON_CELLS )
        field.setArray( array )
        # transform the field to a 3D vector field
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        field.applyFunc( 3, func ) # require 3 components 
        self.assertTrue( field.getNumberOfComponents() == 3 ) # 3 components as required
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_1]
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_2]
        vec1 = field.getArray().getTuple(1) # vector #1
        a,b = values[2], values[3] # initial components of the vector #1
        self.assertAlmostEqual( vec1[0], 10 + b, 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vec1[1], 10 + a, 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vec1[2], 10 + sqrt(a*a+b*b), 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_2]
        return

    def testExample_MEDCouplingFieldDouble_applyFunc_val(self):
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_val_1]
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr)
        field = MEDCouplingFieldDouble( ON_CELLS )
        field.setMesh( mesh )
        field.fillFromAnalytic(2,"IVec * x + JVec * y") # 2 components
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_val_1]
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_val_2]
        newValue = 7.
        field.applyFunc( 3, newValue ) # 3 components are required
        self.assertTrue( field.getIJ(1,0) == newValue ) # a value is as expected
        self.assertTrue( field.getNumberOfComponents() == 3 )
        self.assertTrue( field.getNumberOfTuples() == mesh.getNumberOfCells() )
        #! [PySnippet_MEDCouplingFieldDouble_applyFunc_val_2]
        return

    def testExample_MEDCouplingFieldDouble_fillFromAnalytic3(self):
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic3_1]
        coords = [0.,2.,4.,6.] #  6. is not used
        x=DataArrayDouble(coords[:3],3,1)
        y=DataArrayDouble(coords[:2],2,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(x,y)
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic3_1]
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic3_2]
        field = MEDCouplingFieldDouble( ON_CELLS )
        field.setMesh( mesh )
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        varNames=["a","b"] # names used to refer to X and Y coord components
        field.fillFromAnalyticNamedCompo(3,varNames,func)
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic3_2]
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic3_3]
        vals1 = field.getArray().getTuple(1) # values of the cell #1
        assert len( vals1 ) == 3 # 3 components in the field
        #
        bc = mesh.computeCellCenterOfMass() # func is applied to barycenters of cells
        bc1 = bc.getTuple(1) # coordinates of the second point
        #
        dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] ) # "sqrt( a*a + b*b )"
        self.assertAlmostEqual( vals1[0], 10 + bc1[1], 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vals1[1], 10 + bc1[0], 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vals1[2], 10 + dist  , 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic3_3]
        return

    def testExample_MEDCouplingFieldDouble_fillFromAnalytic2(self):
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic2_1]
        coords = [0.,2.,4.]
        x=DataArrayDouble(coords[:3],3,1)
        y=DataArrayDouble(coords[:2],2,1)
        x.setInfoOnComponent(0,"a") # name used to refer to X coordinate within a function
        y.setInfoOnComponent(0,"b") # name used to refer to Y coordinate within a function
        mesh=MEDCouplingCMesh()
        mesh.setCoords(x,y)
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic2_1]
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic2_2]
        field = MEDCouplingFieldDouble( ON_CELLS )
        field.setMesh( mesh )
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        field.fillFromAnalyticCompo(3,func)
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic2_2]
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic2_3]
        vals1 = field.getArray().getTuple(1) # values of the cell #1
        assert len( vals1 ) == 3 # 3 components in the field
        #
        bc = mesh.computeCellCenterOfMass() # func is applied to barycenters of cells
        bc1 = bc.getTuple(1) # coordinates of the second point
        #
        dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] ) # "sqrt( a*a + b*b )"
        self.assertAlmostEqual( vals1[0], 10 + bc1[1], 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vals1[1], 10 + bc1[0], 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vals1[2], 10 + dist  , 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic2_3]
        return

    def testExample_MEDCouplingFieldDouble_fillFromAnalytic(self):
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic_1]
        coords = [0.,2.,4.]
        x=DataArrayDouble(coords[:3],3,1)
        y=DataArrayDouble(coords[:2],2,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(x,y)
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic_1]
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic_2]
        field = MEDCouplingFieldDouble( ON_CELLS )
        field.setMesh( mesh )
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        field.fillFromAnalytic(3,func)
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic_2]
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic_3]
        vals1 = field.getArray().getTuple(1) # values of the cell #1
        assert len( vals1 ) == 3 # 3 components in the field
        #
        bc = mesh.computeCellCenterOfMass() # func is applied to barycenters of cells
        bc1 = bc.getTuple(1) # coordinates of the second point
        #
        dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] ) # "sqrt( a*a + b*b )"
        self.assertAlmostEqual( vals1[0], 10 + bc1[1], 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vals1[1], 10 + bc1[0], 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vals1[2], 10 + dist  , 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingFieldDouble_fillFromAnalytic_3]
        return

    def testExample_MEDCouplingFieldDouble_getValueOn_time(self):
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_time_1]
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr)
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_time_1]
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_time_2]
        field = MEDCouplingFieldDouble( ON_CELLS, LINEAR_TIME )
        field.setMesh( mesh )
        field.fillFromAnalytic(1,"10") # all values == 10.
        field.setEndArray( field.getArray() + field.getArray() ) # all values == 20.
        time1, time2 = 1.1, 22.
        field.setStartTime( time1, 0, 0 )
        field.setEndTime  ( time2, 0, 0 )
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_time_2]
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_time_3]
        pos = [ 1., 1. ] # we are in 2D space
        value = field.getValueOn( pos, 0.5*( time1 + time2 ))
        self.assertTrue( value[0] == 0.5*( 10. + 20.))
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_time_3]
        return

    def testExample_MEDCouplingFieldDouble_getValueOnMulti(self):
        #! [PySnippet_MEDCouplingFieldDouble_getValueOnMulti_1]
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr)
        field = mesh.fillFromAnalytic(ON_CELLS,1,"x+y")
        #! [PySnippet_MEDCouplingFieldDouble_getValueOnMulti_1]
        #! [PySnippet_MEDCouplingFieldDouble_getValueOnMulti_2]
        bc = mesh.computeCellCenterOfMass() # field values are located at cell barycenters
        valArray = field.getValueOnMulti( bc )
        self.assertTrue( valArray.isEqual( field.getArray(), 1e-13 ))
        #! [PySnippet_MEDCouplingFieldDouble_getValueOnMulti_2]
        return

    def testExample_MEDCouplingFieldDouble_getValueOn(self):
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_1]
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr)
        field = mesh.fillFromAnalytic(ON_CELLS,1,"x+y")
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_1]
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_2]
        bc = mesh.computeCellCenterOfMass() # field values are located at cell barycenters
        vals = [] # array to collect values returned by getValueOn()
        for i,tupl in enumerate( bc ):
            vals.extend( field.getValueOn( tupl ) )
        self.assertTrue( vals == field.getArray().getValues() )
        #! [PySnippet_MEDCouplingFieldDouble_getValueOn_2]
        return

    def testExample_MEDCouplingFieldDouble_getValueOnPos(self):
        #! [PySnippet_MEDCouplingFieldDouble_getValueOnPos_1]
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr)
        field = mesh.fillFromAnalytic(ON_CELLS,1,"x+y")
        #! [PySnippet_MEDCouplingFieldDouble_getValueOnPos_1]
        #! [PySnippet_MEDCouplingFieldDouble_getValueOnPos_2]
        val11 = field.getValueOnPos( 1,1,-1)
        bc = mesh.computeCellCenterOfMass() # field values are located at cell barycenters
        self.assertTrue( val11[0] == bc[3,0] + bc[3,1] )
        #! [PySnippet_MEDCouplingFieldDouble_getValueOnPos_2]
        return

    def testExample_MEDCouplingFieldDouble_renumberNodes(self):
        #! [PySnippet_MEDCouplingFieldDouble_renumberNodes_1]
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr)
        mesh=mesh.buildUnstructured()
        #! [PySnippet_MEDCouplingFieldDouble_renumberNodes_1]
        #! [PySnippet_MEDCouplingFieldDouble_renumberNodes_2]
        field = mesh.fillFromAnalytic(ON_NODES,2,"IVec*x+JVec*y")
        values = field.getArray()
        nodeCoords = mesh.getCoords()
        self.assertTrue( values.isEqualWithoutConsideringStr( nodeCoords, 1e-13 ))
        #! [PySnippet_MEDCouplingFieldDouble_renumberNodes_2]
        #! [PySnippet_MEDCouplingFieldDouble_renumberNodes_3]
        renumber = [8, 7, 6, 5, 4, 3, 2, 1, 0]
        field.renumberNodes(renumber,False)
        mesh2 = field.getMesh() # field now refers to another mesh
        values = field.getArray()
        nodeCoords = mesh2.getCoords()
        self.assertTrue( values.isEqualWithoutConsideringStr( nodeCoords, 1e-13 ))
        #! [PySnippet_MEDCouplingFieldDouble_renumberNodes_3]
        return


    def testExample_MEDCouplingFieldDouble_renumberCells(self):
        #! [PySnippet_MEDCouplingFieldDouble_renumberCells_1]
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr)
        mesh=mesh.buildUnstructured()
        #! [PySnippet_MEDCouplingFieldDouble_renumberCells_1]
        #! [PySnippet_MEDCouplingFieldDouble_renumberCells_2]
        field = mesh.fillFromAnalytic(ON_CELLS,2,"IVec*x+JVec*y")
        values = field.getArray()
        bc = mesh.computeCellCenterOfMass()
        self.assertTrue( values.isEqualWithoutConsideringStr( bc, 1e-13 ))
        #! [PySnippet_MEDCouplingFieldDouble_renumberCells_2]
        #! [PySnippet_MEDCouplingFieldDouble_renumberCells_3]
        renumber = [ 3, 2, 1, 0 ]
        field.renumberCells(renumber,False)
        mesh2 = field.getMesh() # field now refers to another mesh
        values = field.getArray()
        bc = mesh2.computeCellCenterOfMass()
        self.assertTrue( values.isEqualWithoutConsideringStr( bc, 1e-13 ))
        #! [PySnippet_MEDCouplingFieldDouble_renumberCells_3]
        return

    def testExample_MEDCouplingFieldDouble_buildNewTimeReprFromThis(self):
        #! [PySnippet_MEDCouplingFieldDouble_buildNewTimeReprFromThis_1]
        coords = [0.,2.,4.]
        coordsArr=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(coordsArr,coordsArr)
        field1 = mesh.fillFromAnalytic(ON_NODES,1,"x+y")
        self.assertTrue( field1.getTimeDiscretization() == ONE_TIME )
        #! [PySnippet_MEDCouplingFieldDouble_buildNewTimeReprFromThis_1]
        #! [PySnippet_MEDCouplingFieldDouble_buildNewTimeReprFromThis_2]
        field2 = field1.buildNewTimeReprFromThis(NO_TIME,False)
        self.assertTrue( field2.getTimeDiscretization() == NO_TIME )
        #! [PySnippet_MEDCouplingFieldDouble_buildNewTimeReprFromThis_2]
        return

    def testExample_MEDCouplingMesh_fillFromAnalytic3(self):
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic3_1]
        coords = [0.,2.,4.,6.] #  6. is not used
        x=DataArrayDouble(coords[:3],3,1)
        y=DataArrayDouble(coords[:2],2,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(x,y)
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic3_1]
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic3_2]
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        varNames=["a","b"] # names used to refer to X and Y coord components
        field=mesh.fillFromAnalyticNamedCompo(ON_CELLS,3,varNames,func)
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic3_2]
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic3_3]
        vals1 = field.getArray().getTuple(1) # values of the cell #1
        assert len( vals1 ) == 3 # 3 components in the field
        #
        bc = mesh.computeCellCenterOfMass() # func is applied to barycenters of cells
        bc1 = bc.getTuple(1) # coordinates of the second point
        #
        dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] ) # "sqrt( a*a + b*b )"
        self.assertAlmostEqual( vals1[0], 10 + bc1[1], 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vals1[1], 10 + bc1[0], 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vals1[2], 10 + dist  , 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic3_3]
        return

    def testExample_MEDCouplingMesh_fillFromAnalytic2(self):
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic2_1]
        coords = [0.,2.,4.,6.] #  6. is not used
        x=DataArrayDouble(coords[:3],3,1)
        y=DataArrayDouble(coords[:2],2,1)
        x.setInfoOnComponent(0,"a") # name used to refer to X coordinate within a function
        y.setInfoOnComponent(0,"b") # name used to refer to Y coordinate within a function
        mesh=MEDCouplingCMesh()
        mesh.setCoords(x,y)
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic2_1]
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic2_2]
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        field=mesh.fillFromAnalyticCompo(ON_CELLS,3,func)
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic2_2]
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic2_3]
        vals1 = field.getArray().getTuple(1) # values of the cell #1
        assert len( vals1 ) == 3 # 3 components in the field
        #
        bc = mesh.computeCellCenterOfMass() # func is applied to barycenters of cells
        bc1 = bc.getTuple(1) # coordinates of the second point
        #
        dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] ) # "sqrt( a*a + b*b )"
        self.assertAlmostEqual( vals1[0], 10 + bc1[1], 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vals1[1], 10 + bc1[0], 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vals1[2], 10 + dist  , 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic2_3]
        return

    def testExample_MEDCouplingMesh_fillFromAnalytic(self):
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic_1]
        coords = [0.,2.,4.,6.] #  6. is not used
        x=DataArrayDouble(coords[:3],3,1)
        y=DataArrayDouble(coords[:2],2,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoords(x,y)
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic_1]
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic_2]
        func = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10"
        field=mesh.fillFromAnalytic(ON_CELLS,3,func)
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic_2]
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic_3]
        vals1 = field.getArray().getTuple(1) # values of the cell #1
        assert len( vals1 ) == 3 # 3 components in the field
        #
        bc = mesh.computeCellCenterOfMass() # func is applied to barycenters of cells
        bc1 = bc.getTuple(1) # coordinates of the second point
        #
        dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] ) # "sqrt( a*a + b*b )"
        self.assertAlmostEqual( vals1[0], 10 + bc1[1], 13 ) # "10 + IVec * b"
        self.assertAlmostEqual( vals1[1], 10 + bc1[0], 13 ) # "10 + JVec * a"
        self.assertAlmostEqual( vals1[2], 10 + dist  , 13 ) # "10 + KVec * sqrt( a*a + b*b )"
        #! [PySnippet_MEDCouplingMesh_fillFromAnalytic_3]
        return

    def testExample_MEDCouplingCMesh_getCoordsAt(self):
        #! [PySnippet_MEDCouplingCMesh_getCoordsAt_1]
        coords = [1.,2.,4.]
        x=DataArrayDouble(coords,3,1)
        mesh=MEDCouplingCMesh()
        mesh.setCoordsAt(0,x)
        x2=mesh.getCoordsAt(0)
        assert coords == x2.getValues()
        #! [PySnippet_MEDCouplingCMesh_getCoordsAt_1]
        return

    def testExample_MEDCouplingUMesh_areCellsIncludedIn(self):
        #! [PySnippet_MEDCouplingUMesh_areCellsIncludedIn_1]
        mesh1=MEDCouplingUMesh()
        mesh1.setMeshDimension(2)
        mesh1.allocateCells(5)
        conn=[0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh1.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh1.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh1.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh1.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh1.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh1.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh1.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_areCellsIncludedIn_1]
        #! [PySnippet_MEDCouplingUMesh_areCellsIncludedIn_2]
        cells2 = [ 4,2,0 ]
        mesh2 = mesh1.buildPartOfMySelf(cells2, True ) # even cells selected
        #! [PySnippet_MEDCouplingUMesh_areCellsIncludedIn_2]
        #! [PySnippet_MEDCouplingUMesh_areCellsIncludedIn_3]
        compType = 0 # the strongest policy
        isOk, corr2to1 = mesh1.areCellsIncludedIn( mesh2, compType )
        assert isOk # a larger mesh1 includes a smaller mesh2
        assert corr2to1.getValues() == cells2
        #! [PySnippet_MEDCouplingUMesh_areCellsIncludedIn_3]
        #! [PySnippet_MEDCouplingUMesh_areCellsIncludedIn_4]
        isOk, corr1to2 = mesh2.areCellsIncludedIn( mesh1, compType )
        assert not isOk # the smaller mesh2 does NOT include the larger mesh1
        assert corr1to2.getValues() == [2, 3, 1, 4, 0]
        #! [PySnippet_MEDCouplingUMesh_areCellsIncludedIn_4]

    def testExample_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells(self):
        #! [PySnippet_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells_1]
        # 2D coordinates of 5 base nodes
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2]
        coordsArr=DataArrayDouble(coords,5,2)
        # coordinates of 5 top nodes
        coordsArr2 = coordsArr.deepCopy()
        # 3D coordinates of base + top nodes
        coordsArr  = coordsArr.changeNbOfComponents( 3, 0 )
        coordsArr2 = coordsArr2.changeNbOfComponents( 3, 1 )
        coordsArr = DataArrayDouble.Aggregate([coordsArr,coordsArr2])
        # mesh
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        mesh.setMeshDimension(3)
        mesh.allocateCells(2)
        # connectivity of reversed HEXA8 and PENTA6
        conn=[0,1,4,3, 5,6,9,8, 1,2,4, 6,7,9]
        mesh.insertNextCell(NORM_HEXA8, 8,conn[0:0+8])
        mesh.insertNextCell(NORM_PENTA6,6,conn[8:8+6])
        mesh.finishInsertingCells()
        #! [PySnippet_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells_1]
        #! [PySnippet_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells_2]
        fixedCells = mesh.findAndCorrectBadOriented3DExtrudedCells()
        assert len( fixedCells ) == 2 # 2 cells fixed
        fixedCells = mesh.findAndCorrectBadOriented3DExtrudedCells()
        assert len( fixedCells ) == 0 # no bad cells
        #! [PySnippet_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells_2]
        return

    def testExample_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented(self):
        #! [PySnippet_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented_1]
        # 2D coordinates of 5 base nodes
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2]
        coordsArr=DataArrayDouble(coords,5,2)
        # coordinates of 5 top nodes
        coordsArr2 = coordsArr.deepCopy()
        # 3D coordinates of base + top nodes
        coordsArr  = coordsArr.changeNbOfComponents( 3, 0 )
        coordsArr2 = coordsArr2.changeNbOfComponents( 3, 1 )
        coordsArr = DataArrayDouble.Aggregate([coordsArr,coordsArr2])
        # mesh
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        mesh.setMeshDimension(3)
        mesh.allocateCells(2)
        # connectivity of a HEXA8 + a reversed PENTA6
        conn=[0,3,4,1, 5,8,9,6, 1,2,4, 6,7,9]
        mesh.insertNextCell(NORM_POLYHED, 8,conn[0:0+8]) # "extruded" polyhedron
        mesh.insertNextCell(NORM_POLYHED,6,conn[8:8+6])
        mesh.finishInsertingCells()
        # fix connectivity of NORM_POLYHED's
        mesh.convertExtrudedPolyhedra()
        #! [PySnippet_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented_1]
        #! [PySnippet_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented_2]
        badCells = mesh.arePolyhedronsNotCorrectlyOriented()
        assert len( badCells ) == 1 # one polyhedron is KO
        # fix invalid rolyherdons
        mesh.orientCorrectlyPolyhedrons()
        # re-check the orientation
        badCells = mesh.arePolyhedronsNotCorrectlyOriented()
        assert len( badCells ) == 0 # connectivity is OK
        #! [PySnippet_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented_2]
        return

    def testExample_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented(self):
        #! [PySnippet_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        mesh.changeSpaceDimension(3)
        #! [PySnippet_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented_1]
        #! [PySnippet_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented_2]
        vec = [0.,0.,-1.]
        badCellIds=mesh.are2DCellsNotCorrectlyOriented( vec, False )
        assert len( badCellIds ) == 1 # one cell is reversed
        # fix orientation
        mesh.orientCorrectly2DCells( vec, False )
        # re-check orientation
        badCellIds=mesh.are2DCellsNotCorrectlyOriented( vec, False )
        assert len( badCellIds ) == 0 # the orientation is OK
        #! [PySnippet_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented_2]
        return

    def testExample_MEDCouplingUMesh_getCellsContainingPoints(self):
        #! [PySnippet_MEDCouplingUMesh_getCellsContainingPoints_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_getCellsContainingPoints_1]
        #! [PySnippet_MEDCouplingUMesh_getCellsContainingPoints_2]
        pos = [ 10., 10,              # point out of the mesh
                0.3, 0.3,             # point located somewhere inside the mesh
                coords[2], coords[3]] # point at the node #1
        eps = 1e-4 # ball radius
        cells,cellsIndex=mesh.getCellsContainingPoints( pos, 3, eps )
        assert cells.getValues() == [4, 0, 1]
        assert cellsIndex.getValues() == [0, 0, 1, 3]
        #! [PySnippet_MEDCouplingUMesh_getCellsContainingPoints_2]
        return


    def testExample_MEDCouplingUMesh_getCellsContainingPoint(self):
        #! [PySnippet_MEDCouplingUMesh_getCellsContainingPoint_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])  
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])  
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10]) 
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14])
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18])
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_getCellsContainingPoint_1]
        #! [PySnippet_MEDCouplingUMesh_getCellsContainingPoint_2]
        pos4 = coords[ 4*2 : ] # coordinates of the node #4
        eps = 1e-4 # ball radius
        pos = [ pos4[0]+eps, pos4[1]-eps ] # ball center
        cellIds=mesh.getCellsContainingPoint( pos, eps )
        assert len( cellIds ) == mesh.getNumberOfCells()
        #! [PySnippet_MEDCouplingUMesh_getCellsContainingPoint_2]
        return


    def testExample_MEDCouplingUMesh_buildPartOrthogonalField(self):
        #! [PySnippet_MEDCouplingUMesh_buildPartOrthogonalField_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_buildPartOrthogonalField_1]
        #! [PySnippet_MEDCouplingUMesh_buildPartOrthogonalField_2]
        part = DataArrayInt([1,2,3,4],4,1) # cell #0 is omitted
        vecField=mesh.buildPartOrthogonalField( part )
        vecArr = vecField.getArray()
        assert len( vecArr ) == len( part )
        assert vecArr.getNumberOfComponents() == 3
        #! [PySnippet_MEDCouplingUMesh_buildPartOrthogonalField_2]
        return

    def testExample_MEDCouplingUMesh_getPartMeasureField(self):
        #! [PySnippet_MEDCouplingUMesh_getPartMeasureField_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_getPartMeasureField_1]
        #! [PySnippet_MEDCouplingUMesh_getPartMeasureField_2]
        isAbs = True
        part = DataArrayInt([1,2,3,4],4,1) # cell #0 is omitted
        areaArr=mesh.getPartMeasureField( isAbs, part )
        assert areaArr[0] > 0 # orientation ignored
        areaArr=mesh.getPartMeasureField( not isAbs, part )
        assert areaArr[0] < 0 # orientation considered
        assert len( areaArr ) == len( part )
        #! [PySnippet_MEDCouplingUMesh_getPartMeasureField_2]
        #! [PySnippet_MEDCouplingUMesh_getPartMeasureField_3]
        part = DataArrayInt([1,2,3,4],4,1) # cell #0 is omitted
        baryCenters = mesh.getPartBarycenterAndOwner( part )
        assert len( baryCenters ) == len( part )
        assert baryCenters.getNumberOfComponents() == mesh.getSpaceDimension()
        #! [PySnippet_MEDCouplingUMesh_getPartMeasureField_3]
        return

    def testExample_MEDCouplingUMesh_getCellsInBoundingBox(self):
        #! [PySnippet_MEDCouplingUMesh_getCellsInBoundingBox_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        coords=[0.,0., 0.,1., 1.,1]
        coordsArr=DataArrayDouble(coords,3,2)
        mesh.setCoords(coordsArr)
        mesh.allocateCells(1)
        conn=[0,1,2]
        mesh.insertNextCell(NORM_TRI3,3,conn)
        mesh.finishInsertingCells()
        #! [PySnippet_MEDCouplingUMesh_getCellsInBoundingBox_1]
        #! [PySnippet_MEDCouplingUMesh_getCellsInBoundingBox_2]
        bbox = [1., 1., 1.001,1.001] # xMin, xMax, yMin, yMax
        cellsInBox = mesh.getCellsInBoundingBox( bbox, 0.0 )
        assert cellsInBox.getValues() == []
        cellsInBox = mesh.getCellsInBoundingBox( bbox, 0.1 )
        assert cellsInBox.getValues() == [0]
        #! [PySnippet_MEDCouplingUMesh_getCellsInBoundingBox_2]


    def testExample_MEDCouplingUMesh_renumberNodesInConn(self):
        #! [PySnippet_MEDCouplingUMesh_renumberNodesInConn_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(1)
        conn=[4,3,2,1]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])
        mesh.finishInsertingCells()
        #! [PySnippet_MEDCouplingUMesh_renumberNodesInConn_1]
        #! [PySnippet_MEDCouplingUMesh_renumberNodesInConn_2]
        old2newIds = [-1,3,2,1,0]
        mesh.renumberNodesInConn( old2newIds )
        nodes0 = mesh.getNodeIdsOfCell( 0 )
        assert nodes0 == [0,1,2,3]
        #! [PySnippet_MEDCouplingUMesh_renumberNodesInConn_2]
        return


    def testExample_MEDCouplingUMesh_renumberNodes(self):
        #! [PySnippet_MEDCouplingUMesh_renumberNodes_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.3]
        coordsArr=DataArrayDouble(coords,4,2)
        mesh.setCoords(coordsArr)
        mesh.allocateCells(0)
        mesh.finishInsertingCells()
        #! [PySnippet_MEDCouplingUMesh_renumberNodes_1]
        #! [PySnippet_MEDCouplingUMesh_renumberNodes_2]
        mesh.renumberNodes([ 2,1,0,-1 ], 3)
        coordsArr = mesh.getCoords() # get a shorten array
        assert coordsArr.getValues() == [0.7,-0.3, 0.2,-0.3, -0.3,-0.3]
        #! [PySnippet_MEDCouplingUMesh_renumberNodes_2]
        #! [PySnippet_MEDCouplingUMesh_renumberNodes_3]
        coordsArr.setValues(coords,4,2) # restore old nodes
        mesh.renumberNodesCenter([ 2,1,0,2 ], 3)
        coordsArr = mesh.getCoords() # get a shorten array
        assert coordsArr.getValues() == [0.7,-0.3, 0.2,-0.3, -0.3,0.0]
        #! [PySnippet_MEDCouplingUMesh_renumberNodes_3]
        return

    def testExample_MEDCouplingUMesh_findBoundaryNodes(self):
        #! [PySnippet_MEDCouplingUMesh_findBoundaryNodes_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])  
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])  
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10]) 
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14])
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18])
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_findBoundaryNodes_1]
        #! [PySnippet_MEDCouplingUMesh_findBoundaryNodes_2]
        nodeIdsArr=mesh.findBoundaryNodes()
        assert nodeIdsArr.getNumberOfTuples() == mesh.getNumberOfNodes() - 1 
        #! [PySnippet_MEDCouplingUMesh_findBoundaryNodes_2]
        return

    def testExample_MEDCouplingUMesh_buildBoundaryMesh(self):
        #! [PySnippet_MEDCouplingUMesh_buildBoundaryMesh_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])  
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])  
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10]) 
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14])
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18])
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_buildBoundaryMesh_1]
        #! [PySnippet_MEDCouplingUMesh_buildBoundaryMesh_2]
        mesh1=mesh.buildBoundaryMesh(True)
        mesh2=mesh.buildBoundaryMesh(False)
        assert coordsArr.isEqual( mesh1.getCoords(), 1e-13 )  # same nodes
        assert not coordsArr.isEqual( mesh2.getCoords(), 1e-13 ) # different nodes
        #! [PySnippet_MEDCouplingUMesh_buildBoundaryMesh_2]
        return

    def testExample_MEDCouplingUMesh_buildFacePartOfMySelfNode(self):
        #! [PySnippet_MEDCouplingUMesh_buildFacePartOfMySelfNode_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_buildFacePartOfMySelfNode_1]
        #! [PySnippet_MEDCouplingUMesh_buildFacePartOfMySelfNode_2]
        nodeIds = mesh.getNodeIdsOfCell( 0 )
        allNodes = True
        mesh1 = mesh.buildFacePartOfMySelfNode( nodeIds, allNodes )
        assert mesh1.getNumberOfCells() == 4 # 4 segments bounding QUAD4 #0 only
        mesh2 = mesh.buildFacePartOfMySelfNode( nodeIds, not allNodes )
        assert mesh2.getNumberOfCells() >  4 # more segments added
        #! [PySnippet_MEDCouplingUMesh_buildFacePartOfMySelfNode_2]
        return


    def testExample_MEDCouplingUMesh_buildPartOfMySelfNode(self):
        #! [PySnippet_MEDCouplingUMesh_buildPartOfMySelfNode_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_buildPartOfMySelfNode_1]
        #! [PySnippet_MEDCouplingUMesh_buildPartOfMySelfNode_2]
        nodeIds = mesh.getNodeIdsOfCell( 0 )
        allNodes = True
        mesh1 = mesh.buildPartOfMySelfNode( nodeIds, allNodes )
        mesh2 = mesh.buildPartOfMySelfNode( nodeIds, not allNodes )
        assert mesh1.getNumberOfCells() == 1 # cell #0 is found only
        assert mesh2.getNumberOfCells() == mesh.getNumberOfCells() # all cells are found
        #! [PySnippet_MEDCouplingUMesh_buildPartOfMySelfNode_2]
        return


    def testExample_MEDCouplingUMesh_getCellIdsLyingOnNodes(self):
        #! [PySnippet_MEDCouplingUMesh_getCellIdsLyingOnNodes_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_getCellIdsLyingOnNodes_1]
        #! [PySnippet_MEDCouplingUMesh_getCellIdsLyingOnNodes_2]
        nodeIds = mesh.getNodeIdsOfCell( 0 )
        allNodes = True
        cellIdsArr1 = mesh.getCellIdsLyingOnNodes( nodeIds, allNodes )
        cellIdsArr2 = mesh.getCellIdsLyingOnNodes( nodeIds, not allNodes )
        assert cellIdsArr1.getNumberOfTuples() == 1 # cell #0 is found only
        assert cellIdsArr2.getNumberOfTuples() == mesh.getNumberOfCells() # all cells are found
        #! [PySnippet_MEDCouplingUMesh_getCellIdsLyingOnNodes_2]
        return


    def testExample_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds(self):
        #! [PySnippet_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])  
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])  
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10]) 
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14])
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18])
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds_1]
        #! [PySnippet_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds_2]
        cellIds = [1,2]
        nodeIds =  mesh.getNodeIdsOfCell( cellIds[0] )
        nodeIds += mesh.getNodeIdsOfCell( cellIds[1] )
        cellIdsArr = mesh.getCellIdsFullyIncludedInNodeIds( nodeIds )
        assert cellIdsArr.getValues() == cellIds
        #! [PySnippet_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds_2]
        return


    def testExample_MEDCouplingUMesh_buildPartOfMySelf(self):
        #! [PySnippet_MEDCouplingUMesh_buildPartOfMySelf_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_buildPartOfMySelf_1]
        #! [PySnippet_MEDCouplingUMesh_buildPartOfMySelf_2]
        cellIds=[1,2]
        mesh2=mesh.buildPartOfMySelf(cellIds, True)
        mesh3=mesh.buildPartOfMySelf(cellIds, False)
        coordsArr2 = mesh2.getCoords()
        assert coordsArr.isEqual( coordsArr2, 1e-13 )  # same nodes
        coordsArr3 = mesh3.getCoords()
        assert not coordsArr.isEqual( coordsArr3, 1e-13 ) # different nodes
        assert mesh2.getNodeIdsOfCell(0) == mesh.getNodeIdsOfCell( cellIds[0]) # cell #1 was copied
        assert mesh2.getNodeIdsOfCell(1) == mesh.getNodeIdsOfCell( cellIds[1]) # cell #2 was copied
        #! [PySnippet_MEDCouplingUMesh_buildPartOfMySelf_2]
        return

    def testExample_MEDCouplingUMesh_mergeNodes(self):
        #! [PySnippet_MEDCouplingUMesh_mergeNodes_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4]) 
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7]) 
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])
        mesh.finishInsertingCells()
        coords=[0.3,-0.301, # 0
                0.2,-0.3,   # 1
                0.3,-0.302, # 2 ~~ 0
                1.1,0.0,    # 3
                1.1,0.0,    # 4 == 3
                0.3,-0.303]# 5 ~~ 0
        coordsArr=DataArrayDouble(coords,6,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_mergeNodes_1]
        #! [PySnippet_MEDCouplingUMesh_mergeNodes_2]
        arr,areNodesMerged,newNbOfNodes=mesh.mergeNodes(0.004)
        assert arr.getValues() == [0, 1, 0, 2, 2, 0]
        assert areNodesMerged
        assert newNbOfNodes == 3
        #! [PySnippet_MEDCouplingUMesh_mergeNodes_2]
        #! [PySnippet_MEDCouplingUMesh_mergeNodes_3]
        baryCoords2 = coords[2*2:] # initial coordinates of node #2
        coordsArr = mesh.getCoords() # retrieve a new shorten coord array
        self.assertNotAlmostEqual( baryCoords2[1], coordsArr.getIJ(0,1), 13 ) # Y of node #0 differs from that of baryCoords2
        # restore coordinates
        coordsArr = DataArrayDouble(coords,6,2)
        mesh.setCoords(coordsArr)
        # call mergeNodesCenter()
        mesh.mergeNodesCenter(0.004)
        coordsArr = mesh.getCoords() # retrieve a new shorten coord array
        self.assertAlmostEqual( baryCoords2[1], coordsArr.getIJ(0,1), 13 ) # Y of node #0 equals to that of baryCoords2
        #! [PySnippet_MEDCouplingUMesh_mergeNodes_3]
        return

    def testExample_MEDCouplingUMesh_zipConnectivityTraducer(self):
        #! [PySnippet_MEDCouplingUMesh_zipConnectivityTraducer_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])           # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])           # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])           # 2 == 1
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])           # 3 == 0
        mesh.insertNextCell(NORM_QUAD4,4,conn[2:4]+conn[0:2]) # 4 ~~ 0
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_zipConnectivityTraducer_1]
        #! [PySnippet_MEDCouplingUMesh_zipConnectivityTraducer_2]
        oldNbCells = mesh.getNumberOfCells()
        arr = mesh.zipConnectivityTraducer(0)
        assert mesh.getNumberOfCells() == oldNbCells-2
        assert arr.getValues() == [0, 1, 1, 0, 2]
        #! [PySnippet_MEDCouplingUMesh_zipConnectivityTraducer_2]
        return

    def testExample_MEDCouplingUMesh_zipCoordsTraducer(self):
        #! [PySnippet_MEDCouplingUMesh_zipCoordsTraducer_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])  
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])  
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10]) 
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14])
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18])
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_zipCoordsTraducer_1]
        #! [PySnippet_MEDCouplingUMesh_zipCoordsTraducer_2]
        cellIds=[1,2]
        mesh2=mesh.buildPartOfMySelf(cellIds,True)
        arr=mesh2.zipCoordsTraducer()
        assert mesh2.getNumberOfNodes() == 4 # nb of nodes decreased
        assert arr.getValues() == [-1,0,1,-1,2,3,-1,-1,-1] # -1 for unused nodes
        #! [PySnippet_MEDCouplingUMesh_zipCoordsTraducer_2]
        return

    def testExample_MEDCouplingUMesh_getNodeIdsInUse(self):
        #! [PySnippet_MEDCouplingUMesh_getNodeIdsInUse_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])  
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])  
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10]) 
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14])
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18])
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_getNodeIdsInUse_1]
        #! [PySnippet_MEDCouplingUMesh_getNodeIdsInUse_2]
        cellIds=[1,2]
        mesh2=mesh.buildPartOfMySelf(cellIds,True)
        arr,newNbOfNodes=mesh2.getNodeIdsInUse()
        assert arr.getValues() == [-1,0,1,-1,2,3,-1,-1,-1]
        #! [PySnippet_MEDCouplingUMesh_getNodeIdsInUse_2]
        #! [PySnippet_MEDCouplingUMesh_getNodeIdsInUse_3]
        arr2=arr.invertArrayO2N2N2O(newNbOfNodes)
        assert arr2.getValues() == [1,2,4,5]
        #! [PySnippet_MEDCouplingUMesh_getNodeIdsInUse_3]
        return

    def testExample_MEDCouplingUMesh_convertToPolyTypes(self):
        #! [PySnippet_MEDCouplingUMesh_convertToPolyTypes_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_convertToPolyTypes_1]
        #! [PySnippet_MEDCouplingUMesh_convertToPolyTypes_2]
        cells=[1,3]
        mesh.convertToPolyTypes(cells)
        assert mesh.getTypeOfCell(0) == NORM_QUAD4
        assert mesh.getTypeOfCell(1) == NORM_POLYGON, mesh.getTypeOfCell(1)
        assert mesh.getTypeOfCell(2) == NORM_TRI3
        assert mesh.getTypeOfCell(3) == NORM_POLYGON
        #! [PySnippet_MEDCouplingUMesh_convertToPolyTypes_2]
        return

    def testExample_MEDCouplingUMesh_buildDescendingConnectivity2(self):
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity2_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity2_1]
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity2_2]
        mesh2,desc,descIndx,revDesc,revDescIndx=mesh.buildDescendingConnectivity2()
        assert desc.getValues()        == [1,2,3,4,-3,5,6, 7,8,-5,9,10,-2,11, 12,13,-7,-10]
        assert descIndx.getValues()    == [0,4,7,10,14,18]
        assert revDesc.getValues()     == [0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4]
        assert revDescIndx.getValues() == [0,1,3,5,6,8,9,11,12,13,15,16,17,18]
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity2_2]
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity2_3]
        assert mesh2.getNodeIdsOfCell( 3-1 ) == [4, 1]  # cell #3 in FORTRAN mode
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity2_3]
        return

    def testExample_MEDCouplingUMesh_buildDescendingConnectivity(self):
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity_1]
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity_2]
        mesh2,desc,descIndx,revDesc,revDescIndx=mesh.buildDescendingConnectivity()
        assert desc.getValues()        == [0,1,2,3, 2,4,5, 6,7,4, 8,9,1,10, 11,12,6,9]
        assert descIndx.getValues()    == [0,4,7,10,14,18]
        assert revDesc.getValues()     == [0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4]
        assert revDescIndx.getValues() == [0,1,3,5,6,8,9,11,12,13,15,16,17,18]
        #! [PySnippet_MEDCouplingUMesh_buildDescendingConnectivity_2]
        return

    def testExample_MEDCouplingUMesh_getReverseNodalConnectivity(self):
        #! [PySnippet_MEDCouplingUMesh_getReverseNodalConnectivity_1]
        mesh=MEDCouplingUMesh()
        mesh.setMeshDimension(2)
        mesh.allocateCells(5)
        conn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        mesh.insertNextCell(NORM_QUAD4,4,conn[0:4])   # 0
        mesh.insertNextCell(NORM_TRI3,3, conn[4:7])   # 1
        mesh.insertNextCell(NORM_TRI3,3, conn[7:10])  # 2
        mesh.insertNextCell(NORM_QUAD4,4,conn[10:14]) # 3
        mesh.insertNextCell(NORM_QUAD4,4,conn[14:18]) # 4
        mesh.finishInsertingCells()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        coordsArr=DataArrayDouble(coords,9,2)
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingUMesh_getReverseNodalConnectivity_1]
        #! [PySnippet_MEDCouplingUMesh_getReverseNodalConnectivity_2]
        revNodal,revNodalIndx=mesh.getReverseNodalConnectivity()
        assert revNodal.getValues()     == [0,0,1,1,2,0,3,0,1,2,3,4,2,4,3,3,4,4]
        assert revNodalIndx.getValues() == [0,1,3,5,7,12,14,15,17,18]
        #! [PySnippet_MEDCouplingUMesh_getReverseNodalConnectivity_2]
        return

    def testExample_MEDCouplingUMesh_checkDeepEquivalWith(self):
        #! [PySnippet_MEDCouplingUMesh_checkDeepEquivalWith_1]
        # mesh 1
        mesh1=MEDCouplingUMesh()
        mesh1.setMeshDimension(2)
        coords=[0.0,0.0, #0
                1.0,0.0, #1
                1.0,1.0, #2
                0.0,1.0] #3
        coordsArr=DataArrayDouble(coords,4,2)
        mesh1.setCoords(coordsArr)
        mesh1.allocateCells(2)
        mesh1.insertNextCell(NORM_TRI3,3,[0,1,2]) #0
        mesh1.insertNextCell(NORM_TRI3,3,[1,2,3]) #1
        mesh1.finishInsertingCells()
        # mesh 2
        mesh2=MEDCouplingUMesh()
        mesh2.setMeshDimension(2)
        coords=[0.0,1.0,    #0 = #3
                0.0,0.0,    #1 = #0
                1.0,0.0,    #2 = #1
                1.0,1.001]  #3 ~ #2
        coordsArr2=DataArrayDouble(coords,4,2)
        mesh2.setCoords(coordsArr2)
        mesh2.allocateCells(2)
        mesh2.insertNextCell(NORM_TRI3,3,[2,3,0]) #0 = #1
        mesh2.insertNextCell(NORM_TRI3,3,[3,1,2]) #1 ~ #0
        mesh2.finishInsertingCells()
        #! [PySnippet_MEDCouplingUMesh_checkDeepEquivalWith_1]
        #! [PySnippet_MEDCouplingUMesh_checkDeepEquivalWith_2]
        cellCompPol = 1 # "permuted same orientation" - policy of medium severity
        cOld2New, nOld2New = mesh1.checkDeepEquivalWith( mesh2, cellCompPol, 0.002 )
        assert nOld2New.getValues() == [3, 0, 1, 2]
        assert cOld2New.getValues() == [1, 0]
        #! [PySnippet_MEDCouplingUMesh_checkDeepEquivalWith_2]
        #! [PySnippet_MEDCouplingUMesh_checkDeepEquivalWith_3]
        self.assertRaises( InterpKernelException, mesh1.checkDeepEquivalOnSameNodesWith, mesh2, cellCompPol, 0.002)
        mesh2.setCoords(coordsArr) # make meshes share the same coordinates array
        mesh2.allocateCells(2)
        mesh2.insertNextCell(NORM_TRI3,3,[1,2,3]) #0 = #1
        mesh2.insertNextCell(NORM_TRI3,3,[1,0,2]) #1 ~ #0
        mesh2.finishInsertingCells()
        cellCompPol = 2 # the weakest policy
        mesh1.checkDeepEquivalOnSameNodesWith( mesh2, cellCompPol, 0 )
        #! [PySnippet_MEDCouplingUMesh_checkDeepEquivalWith_3]
        return

    def testExample_MEDCouplingPointSet_scale(self):
        #! [PySnippet_MEDCouplingPointSet_scale_1]
        coords=[0.0,0.0, 1.0,0.0, 1.0,1.0, 0.0,1.0] # 2D coordinates of 4 nodes
        coordsArr=DataArrayDouble(coords,4,2)
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        initCoords = coordsArr.deepCopy()
        #! [PySnippet_MEDCouplingPointSet_scale_1]
        #! [PySnippet_MEDCouplingPointSet_scale_2]
        center = [0.,0.]
        factor = 2.
        mesh.scale(center,factor)
        #! [PySnippet_MEDCouplingPointSet_scale_2]
        #! [PySnippet_MEDCouplingPointSet_scale_3]
        coords2 = mesh.getCoords()
        assert coords2.isEqualWithoutConsideringStr( initCoords, 1.0 )
        assert not coords2.isEqualWithoutConsideringStr( initCoords, 0.9 )
        #! [PySnippet_MEDCouplingPointSet_scale_3]
        return

    def testExample_MEDCouplingPointSet_translate(self):
        #! [PySnippet_MEDCouplingPointSet_translate_1]
        coords=[0.0,0.0, 1.0,0.0, 1.0,1.0, 0.0,1.0] # 2D coordinates of 4 nodes
        coordsArr=DataArrayDouble(coords,4,2)
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        initCoords = coordsArr.deepCopy()
        #! [PySnippet_MEDCouplingPointSet_translate_1]
        #! [PySnippet_MEDCouplingPointSet_translate_2]
        vector = [1.,1.]
        mesh.translate(vector)
        #! [PySnippet_MEDCouplingPointSet_translate_2]
        #! [PySnippet_MEDCouplingPointSet_translate_3]
        coords2 = mesh.getCoords()
        assert coords2.isEqualWithoutConsideringStr( initCoords, 1 )
        assert not coords2.isEqualWithoutConsideringStr( initCoords, 0.9 )
        #! [PySnippet_MEDCouplingPointSet_translate_3]
        return

    def testExample_MEDCouplingPointSet_rotate(self):
        #! [PySnippet_MEDCouplingPointSet_rotate_1]
        coords=[0.0,0.0, 0.1,0.0, 0.1,0.1, 0.0,0.1] # 2D coordinates of 4 nodes
        coordsArr=DataArrayDouble(coords,4,2)
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingPointSet_rotate_1]
        #! [PySnippet_MEDCouplingPointSet_rotate_2]
        center = [0.,0.]
        mesh.rotate(center,-pi/2)
        #! [PySnippet_MEDCouplingPointSet_rotate_2]
        #! [PySnippet_MEDCouplingPointSet_rotate_3]
        mesh.changeSpaceDimension(3)
        center = [0.,0.,0.]
        vector = [0.,0.,1.]
        mesh.rotate(center,vector,pi/2)
        #! [PySnippet_MEDCouplingPointSet_rotate_3]
        #! [PySnippet_MEDCouplingPointSet_rotate_4]
        mesh.changeSpaceDimension(2)
        coords2 = mesh.getCoords()
        for i,c in enumerate( coords ):
            self.assertAlmostEqual( c, coords2.getIJ(0,i), 13 )
        #! [PySnippet_MEDCouplingPointSet_rotate_4]
        return

    def testExample_MEDCouplingPointSet_getBoundingBox(self):
        #! [PySnippet_MEDCouplingPointSet_getBoundingBox_1]
        cc=[0.0, 0.1, 0.2, # 3D coordinates of 2 nodes
            2.0, 2.1, 2.2]
        coordsArr=DataArrayDouble(cc,2,3)
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingPointSet_getBoundingBox_1]
        #! [PySnippet_MEDCouplingPointSet_getBoundingBox_2]
        bbox=mesh.getBoundingBox()
        assert bbox == [( cc[0], cc[3] ), # NOTE: list of 3 tuples is retirned!
                        ( cc[1], cc[4] ),
                        ( cc[2], cc[5] )]
        #! [PySnippet_MEDCouplingPointSet_getBoundingBox_2]

    def testExample_MEDCouplingPointSet_getNodeIdsNearPoint(self):
        #! [PySnippet_MEDCouplingPointSet_getNodeIdsNearPoint_1]
        # 2D coordinates of 5 nodes
        coords=[0.3,-0.301, # 0
                0.2,-0.3,   # 1
                0.3,-0.302, # 2
                1.1,0.0,    # 3
                0.3,-0.30299999999999]# 4
        coordsArr=DataArrayDouble(coords,5,2)
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingPointSet_getNodeIdsNearPoint_1]
        #! [PySnippet_MEDCouplingPointSet_getNodeIdsNearPoint_2]
        point=[0.3, -0.3]   # point close to nodes #0, #2 and #4
        ids=mesh.getNodeIdsNearPoint(point,0.003)
        assert ids.getValues() == [0,2,4]
        #! [PySnippet_MEDCouplingPointSet_getNodeIdsNearPoint_2]
        return

    def testExample_MEDCouplingPointSet_getNodeIdsNearPoints(self):
        #! [PySnippet_MEDCouplingPointSet_getNodeIdsNearPoints_1]
        # 2D coordinates of 7 nodes
        coords=[0.3,-0.301, # 0
                0.2,-0.3,   # 1
                0.3,-0.302, # 2
                1.1,0.0,    # 3
                1.1,0.0,    # 4
                1.1,0.002,  # 5
                0.3,-0.303]# 6
        coordsArr=DataArrayDouble(coords,7,2)
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingPointSet_getNodeIdsNearPoints_1]
        #! [PySnippet_MEDCouplingPointSet_getNodeIdsNearPoints_2]
        points=[0.2,-0.301,   # ~ node #1
                0.0, 0.0,
                1.1, 0.002]   # ~ nodes #3, #4 and #5
        ids,idsIndex=mesh.getNodeIdsNearPoints(points,3,0.003)
        assert ids.getValues() == [1, 3, 4, 5]
        assert idsIndex.getValues() == [0, 1, 1, 4]
        #! [PySnippet_MEDCouplingPointSet_getNodeIdsNearPoints_2]
        return

    def testExample_MEDCouplingPointSet_findCommonNodes(self):
        #! [PySnippet_MEDCouplingPointSet_findCommonNodes_1]
        coords=[0.3,-0.301, # 0
                0.2,-0.3,   # 1
                0.3,-0.302, # 2
                1.1,0.0,    # 3
                1.1,0.0,    # 4
                0.3,-0.303]# 5
        coordsArr=DataArrayDouble(coords,6,2)
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
        #! [PySnippet_MEDCouplingPointSet_findCommonNodes_1]
        #! [PySnippet_MEDCouplingPointSet_findCommonNodes_2]
        comm,commI=mesh.findCommonNodes(1e-13)
        assert comm.getValues() == [3,4]
        comm,commI=mesh.findCommonNodes(0.004)
        assert comm.getValues() == [0,2,5,3,4]
        #! [PySnippet_MEDCouplingPointSet_findCommonNodes_2]
        return

    def testExample_MEDCouplingPointSet_getCoordinatesOfNode(self):
        #! [PySnippet_MEDCouplingPointSet_getCoordinatesOfNode_1]
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3]
        coordsArr=DataArrayDouble(coords,3,2)
        mesh=MEDCouplingUMesh()
        mesh.setCoords(coordsArr)
#! [PySnippet_MEDCouplingPointSet_getCoordinatesOfNode_1]
#! [PySnippet_MEDCouplingPointSet_getCoordinatesOfNode_2]
        nodeCoords=mesh.getCoordinatesOfNode(1)
        self.assertAlmostEqual(0.2, nodeCoords[0],13)
        self.assertAlmostEqual(-0.3,nodeCoords[1],13)
#! [PySnippet_MEDCouplingPointSet_getCoordinatesOfNode_2]
        return

    def testExample_DataArrayInt_getTuple(self):
#! [Snippet_DataArrayInt_getTuple_1]
        dv=DataArrayInt()
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
        assert dv.getTuple( 1 ) == [9,10]
#! [Snippet_DataArrayInt_getTuple_1]
#! [Snippet_DataArrayInt_getTuple_2]
        for tpl in dv:
            print(tpl)
#! [Snippet_DataArrayInt_getTuple_2]
        return

    def testExample_DataArrayInt_buildPermutationArr(self):
#! [PySnippet_DataArrayInt_buildPermutationArr_1]
        a=DataArrayInt()
        a.setValues([4,5,6,7,8],5,1)
        b=DataArrayInt()
        b.setValues([5,4,8,6,7],5,1)
        c=a.buildPermutationArr(b)
#! [PySnippet_DataArrayInt_buildPermutationArr_1]
        self.assertEqual([1,0,4,2,3],c.getValues())
        return

    def testExample_DataArrayInt_invertArrayO2N2N2O(self):
#! [PySnippet_DataArrayInt_invertArrayO2N2N2O_1]
        arr1=[2,0,4,1,5,3]
        da=DataArrayInt()
        da.setValues(arr1,6,1)
        da2=da.invertArrayO2N2N2O(6)
        expected1=[1,3,0,5,2,4]
        for i in range(6):
            self.assertEqual(expected1[i],da2.getIJ(i,0))
            pass
#! [PySnippet_DataArrayInt_invertArrayO2N2N2O_1]
        return

    def testExample_DataArrayInt_invertArrayN2O2O2N(self):
#! [PySnippet_DataArrayInt_invertArrayN2O2O2N_1]
        arr1=[2,0,4,1,5,3]
        da=DataArrayInt()
        da.setValues(arr1,6,1)
        da2=da.invertArrayN2O2O2N(7)
        expected1=[1,3,0,5,2,4,-1]
        for i in range(6):
            self.assertEqual(expected1[i],da2.getIJ(i,0))
            pass
#! [PySnippet_DataArrayInt_invertArrayN2O2O2N_1]
        return


    def testExample_DataArrayDouble_getIdsInRange(self):
#! [PySnippet_DataArrayDouble_getIdsInRange_1]
        da=DataArrayDouble()
        da.alloc( 10, 1 )
        da[ :, :] = list(range(10))
        da2 = da.findIdsInRange( 2.5, 6 )
#! [PySnippet_DataArrayDouble_getIdsInRange_1]
        return

    def testExample_DataArrayDouble_setPartOfValues2(self):
#! [Snippet_DataArrayDouble_setPartOfValues2_1]
        da=DataArrayDouble()
        da.alloc( 4, 7 )
        #
        dv=DataArrayDouble()
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
#! [Snippet_DataArrayDouble_setPartOfValues2_1]
#! [Snippet_DataArrayDouble_setPartOfValues2_2]
        da.fillWithZero()
        da[ [0,1,2], [1,3] ] = dv
#! [Snippet_DataArrayDouble_setPartOfValues2_2]
#! [Snippet_DataArrayDouble_setPartOfValues2_3]
        da.fillWithZero()
        dv.rearrange( 6 )
        da[ [0,2,3], [0,2,3,4,5,6]] = dv
#! [Snippet_DataArrayDouble_setPartOfValues2_3]
        return

    def testExample_DataArrayInt_setPartOfValues2(self):
#! [Snippet_DataArrayInt_setPartOfValues2_1]
        da=DataArrayInt()
        da.alloc( 4, 7 )
        #
        dv=DataArrayInt()
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
#! [Snippet_DataArrayInt_setPartOfValues2_1]
#! [Snippet_DataArrayInt_setPartOfValues2_2]
        da.fillWithZero()
        da[ [0,1,2], [1,3] ] = dv
#! [Snippet_DataArrayInt_setPartOfValues2_2]
#! [Snippet_DataArrayInt_setPartOfValues2_3]
        da.fillWithZero()
        dv.rearrange( 6 )
        da[ [0,2,3], [0,2,3,4,5,6]] = dv
#! [Snippet_DataArrayInt_setPartOfValues2_3]
        return

    def testExample_DataArrayDouble_setPartOfValues3(self):
#! [Snippet_DataArrayDouble_setPartOfValues3_1]
        da=DataArrayDouble()
        da.alloc( 4, 7 )
        #
        dv=DataArrayDouble()
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
#! [Snippet_DataArrayDouble_setPartOfValues3_1]
#! [Snippet_DataArrayDouble_setPartOfValues3_2]
        da.fillWithZero()
        da[ 0:3, [1,3] ] = dv
#! [Snippet_DataArrayDouble_setPartOfValues3_2]
#! [Snippet_DataArrayDouble_setPartOfValues3_3]
        da.fillWithZero()
        dv.rearrange( 6 )
        da[ 0:4:2, [0,2,3,4,5,6]] = dv
#! [Snippet_DataArrayDouble_setPartOfValues3_3]
        return

    def testExample_DataArrayInt_setPartOfValues3(self):
#! [Snippet_DataArrayInt_setPartOfValues3_1]
        da=DataArrayInt()
        da.alloc( 4, 7 )
        #
        dv=DataArrayInt()
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
#! [Snippet_DataArrayInt_setPartOfValues3_1]
#! [Snippet_DataArrayInt_setPartOfValues3_2]
        da.fillWithZero()
        da[ 0:3, [1,3] ] = dv
#! [Snippet_DataArrayInt_setPartOfValues3_2]
#! [Snippet_DataArrayInt_setPartOfValues3_3]
        da.fillWithZero()
        dv.rearrange( 6 )
        da[ 0:4:2, [0,2,3,4,5,6]] = dv
#! [Snippet_DataArrayInt_setPartOfValues3_3]
        return

    def testExample_DataArrayDouble_setPartOfValues1(self):
#! [Snippet_DataArrayDouble_setPartOfValues1_1]
        da=DataArrayDouble()
        da.alloc( 4, 4 )
        da.setInfoOnComponents( ["v1","v2","v3","v4"])
        #
        dv=DataArrayDouble()
        dv.alloc( 4, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
        dv.setInfoOnComponents( ["a1","a2"])
#! [Snippet_DataArrayDouble_setPartOfValues1_1]
#! [Snippet_DataArrayDouble_setPartOfValues1_2]
        da.fillWithZero()
        da.setPartOfValues1( dv, 1,3,1, 1,3,1, True )
#! [Snippet_DataArrayDouble_setPartOfValues1_2]
#! [Snippet_DataArrayDouble_setPartOfValues1_3]
        da.fillWithZero()
        da.setPartOfValues1( dv, 0,4,1, 1,2,1, False )
#! [Snippet_DataArrayDouble_setPartOfValues1_3]
#! [Snippet_DataArrayDouble_setPartOfValues1_4]
        da.fillWithZero()
        da.setPartOfValues1( dv, 1,2,1, 0,4,1, False )
#! [Snippet_DataArrayDouble_setPartOfValues1_4]
#! [Snippet_DataArrayDouble_setPartOfValues1_5]
        da.fillWithZero()
        da.setPartOfValues1( dv, 0,3,2, 1,4,2, True )
#! [Snippet_DataArrayDouble_setPartOfValues1_5]
#! [Snippet_DataArrayDouble_setPartOfValues1_6]
        da2 = da.deepCopy()
        da2.fillWithZero()
        da2[ 0:3:2, 1:4:2 ] = dv
        self.assertTrue( da.isEqual( da2, 1e-20 ))
#! [Snippet_DataArrayDouble_setPartOfValues1_6]
        return

    def testExample_DataArrayInt_setPartOfValues1(self):
#! [Snippet_DataArrayInt_setPartOfValues1_1]
        da=DataArrayInt()
        da.alloc( 4, 4 )
        da.setInfoOnComponents( ["v1","v2","v3","v4"])
        #
        dv=DataArrayInt()
        dv.alloc( 4, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
        dv.setInfoOnComponents( ["a1","a2"])
#! [Snippet_DataArrayInt_setPartOfValues1_1]
#! [Snippet_DataArrayInt_setPartOfValues1_2]
        da.fillWithZero()
        da.setPartOfValues1( dv, 1,3,1, 1,3,1, True )
#! [Snippet_DataArrayInt_setPartOfValues1_2]
#! [Snippet_DataArrayInt_setPartOfValues1_3]
        da.fillWithZero()
        da.setPartOfValues1( dv, 0,4,1, 1,2,1, False )
#! [Snippet_DataArrayInt_setPartOfValues1_3]
#! [Snippet_DataArrayInt_setPartOfValues1_4]
        da.fillWithZero()
        da.setPartOfValues1( dv, 1,2,1, 0,4,1, False )
#! [Snippet_DataArrayInt_setPartOfValues1_4]
#! [Snippet_DataArrayInt_setPartOfValues1_5]
        da.fillWithZero()
        da.setPartOfValues1( dv, 0,3,2, 1,4,2, True )
#! [Snippet_DataArrayInt_setPartOfValues1_5]
#! [Snippet_DataArrayInt_setPartOfValues1_6]
        da2 = da.deepCopy()
        da2.fillWithZero()
        da2[ 0:3:2, 1:4:2 ] = dv
        self.assertTrue( da.isEqual( da2 ))
#! [Snippet_DataArrayInt_setPartOfValues1_6]
        return

    def testExample_DataArrayDouble_setPartOfValuesSimple1(self):
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_1]
        da=DataArrayDouble()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_1]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_2]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 1,3,1, 1,3,1 )
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_2]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_3]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 0,4,1, 1,2,1 )
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_3]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_4]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 1,2,1, 0,4,1 )
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_4]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_5]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 0,3,2, 1,4,2 )
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_5]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_6]
        da2 = da.deepCopy()
        da2.fillWithZero()
        da2[ 0:3:2, 1:4:2 ] = dv
        self.assertTrue( da.isEqual( da2, 1e-20 ))
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_6]
        return

    def testExample_DataArrayInt_setPartOfValuesSimple1(self):
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_1]
        da=DataArrayInt()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_1]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_2]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 1,3,1, 1,3,1 )
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_2]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_3]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 0,4,1, 1,2,1 )
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_3]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_4]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 1,2,1, 0,4,1 )
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_4]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_5]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 0,3,2, 1,4,2 )
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_5]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_6]
        da2 = da.deepCopy()
        da2.fillWithZero()
        da2[ 0:3:2, 1:4:2 ] = dv
        self.assertTrue( da.isEqual( da2 ))
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_6]
        return

    def testExample_DataArrayDouble_setPartOfValuesSimple2(self):
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_1]
        da=DataArrayDouble()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_1]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_2]
        da.fillWithZero()
        da[[1,2], [1,2]] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_2]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_3]
        da.fillWithZero()
        da[[0,1,2,3], [1]] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_3]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_4]
        da.fillWithZero()
        da[[1], [0,1,2,3]] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_4]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_5]
        da.fillWithZero()
        da[[0,2], [1,3]] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_5]
        return

    def testExample_DataArrayInt_setPartOfValuesSimple2(self):
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_1]
        da=DataArrayInt()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_1]
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_2]
        da.fillWithZero()
        da[[1,2], [1,2]] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_2]
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_3]
        da.fillWithZero()
        da[[0,1,2,3], [1]] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_3]
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_4]
        da.fillWithZero()
        da[[1], [0,1,2,3]] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_4]
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_5]
        da.fillWithZero()
        da[[0,2], [1,3]] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_5]
        return

    def testExample_DataArrayDouble_setPartOfValuesSimple3(self):
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_1]
        da=DataArrayDouble()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_1]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_2]
        da.fillWithZero()
        da[[1,2], 1:3] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_2]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_3]
        da.fillWithZero()
        da[[0,1,2,3], 1:2] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_3]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_4]
        da.fillWithZero()
        da[[1], 0:4] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_4]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_5]
        da.fillWithZero()
        da[[0,2], 1:4:2] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_5]
        return

    def testExample_DataArrayInt_setPartOfValuesSimple3(self):
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_1]
        da=DataArrayInt()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_1]
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_2]
        da.fillWithZero()
        da[[1,2], 1:3] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_2]
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_3]
        da.fillWithZero()
        da[[0,1,2,3], 1:2] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_3]
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_4]
        da.fillWithZero()
        da[[1], 0:4] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_4]
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_5]
        da.fillWithZero()
        da[[0,2], 1:4:2] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_5]
        return

    def testExample_DataArrayDouble_setSelectedComponents(self):
#! [Snippet_DataArrayDouble_setSelectedComponents1]
        array1=[1.,2., 3.,4., 5.,6.]
        da=DataArrayDouble(array1,3,2)
        da.setInfoOnComponents( ["a1","a2"])
#! [Snippet_DataArrayDouble_setSelectedComponents1]
#! [Snippet_DataArrayDouble_setSelectedComponents2]
        dv=DataArrayDouble()
        dv.alloc( 4, 4 )
        dv.fillWithZero()
        dv.setInfoOnComponents( ["v1","v2","v3","v4"])
        dv2 = dv.deepCopy()
        dv.setSelectedComponents( da, [1,0] )
#! [Snippet_DataArrayDouble_setSelectedComponents2]
#! [Snippet_DataArrayDouble_setSelectedComponents3]
        dv2[:3,[1,0]] = da
        self.assertTrue( dv.isEqualWithoutConsideringStr( dv2, 1e-20 ))
#! [Snippet_DataArrayDouble_setSelectedComponents3]
        return

    def testExample_DataArrayInt_setSelectedComponents(self):
#! [Snippet_DataArrayInt_setSelectedComponents1]
        da=DataArrayInt()
        array1=[1,2, 3,4, 5,6]
        da.setValues(array1,3,2)
        da.setInfoOnComponents( ["a1","a2"])
#! [Snippet_DataArrayInt_setSelectedComponents1]
#! [Snippet_DataArrayInt_setSelectedComponents2]
        dv=DataArrayInt()
        dv.alloc( 4, 4 )
        dv.fillWithZero()
        dv.setInfoOnComponents( ["v1","v2","v3","v4"])
        dv2 = dv.deepCopy()
        dv.setSelectedComponents( da, [1,0] )
#! [Snippet_DataArrayInt_setSelectedComponents2]
#! [Snippet_DataArrayInt_setSelectedComponents3]
        dv2[:3,[1,0]] = da
        self.assertTrue( dv.isEqualWithoutConsideringStr( dv2 ))
#! [Snippet_DataArrayInt_setSelectedComponents3]
        return

    def testExample_DataArrayDouble_getDifferentValues(self):
#! [Snippet_DataArrayDouble_getDifferentValues1]
        array1=[2.3,1.2,1.3,2.3,2.301,0.8]
        da=DataArrayDouble(array1,6,1)
        #
        dv=da.getDifferentValues(2e-1)
        expected2=[2.301,1.3,0.8]
        self.assertEqual(3,dv.getNbOfElems())
        for i in range(3):
            self.assertAlmostEqual(expected2[i],dv.getIJ(i,0),14)
            pass
#! [Snippet_DataArrayDouble_getDifferentValues1]
        return

    def testExample_DataArrayDouble_findCommonTuples1(self):
#! [PySnippet_DataArrayDouble_findCommonTuples1]
        array2=[2.3,2.3, 1.2,1.2, 1.3,1.3, 2.3,2.3, 2.301,2.301, 0.8,0.8]
        da=DataArrayDouble(array2,6,2)        
#! [PySnippet_DataArrayDouble_findCommonTuples1]
#! [PySnippet_DataArrayDouble_findCommonTuples2]
        c,cI=da.findCommonTuples(1.01e-1)
        expected3=[0,3,4,1,2]
        expected4=[0,3,5]
        self.assertEqual(expected3,c.getValues())
        self.assertEqual(expected4,cI.getValues())
#! [PySnippet_DataArrayDouble_findCommonTuples2]
        return

    def testExampleDataArrayDoubleMeldWith(self):
#! [PySnippet_DataArrayDouble_Meld1_1]
        da1=DataArrayDouble()
        da1.alloc(7,2)
        da2=DataArrayDouble()
        da2.alloc(7,1)
        #
        da1.fillWithValue(7.)
        da2.iota(0.)
        da3=da2.applyFunc(3,"10*x*IVec+100*x*JVec+1000*x*KVec")
        #
        da1.setInfoOnComponent(0,"c0da1")
        da1.setInfoOnComponent(1,"c1da1")
        da3.setInfoOnComponent(0,"c0da3")
        da3.setInfoOnComponent(1,"c1da3")
        da3.setInfoOnComponent(2,"c2da3")
        #
        da1C=da1.deepCopy()
        da1.meldWith(da3)
#! [PySnippet_DataArrayDouble_Meld1_1]

    def testExampleDataArrayIntMeldWith(self):
#! [PySnippet_DataArrayInt_Meld1_1]
        da1=DataArrayInt()
        da1.alloc(7,2)
        da2=DataArrayInt()
        da2.alloc(7,1)
        #
        da1.fillWithValue(7)
        da2.iota(0)
        #
        da1.setInfoOnComponent(0,"c0da1")
        da1.setInfoOnComponent(1,"c1da1")
        da2.setInfoOnComponent(0,"c0da2")
        #
        da1.meldWith(da2)
#! [PySnippet_DataArrayInt_Meld1_1]

    def testExampleDataArrayDoubleKeepSelectedComponents1(self):
#! [SnippeDataArrayDoubleKeepSelectedComponents1_1]
        arr1=[1.,2.,3.,4.,     # tuple 0
              11.,12.,13.,14., # tuple 1
              21.,22.,23.,24., # ...
              31.,32.,33.,34.,
              41.,42.,43.,44.]
        a1=DataArrayDouble(arr1,5,4)
        a1.setInfoOnComponent(0,"a")
        a1.setInfoOnComponent(1,"b")
        a1.setInfoOnComponent(2,"c")
        a1.setInfoOnComponent(3,"d")
#! [SnippeDataArrayDoubleKeepSelectedComponents1_1]
#! [SnippeDataArrayDoubleKeepSelectedComponents1_2]
        arr2V=[1,2,1,2,0,0]
        a2=a1.keepSelectedComponents(arr2V)
#! [SnippeDataArrayDoubleKeepSelectedComponents1_2]
        return

    def testExampleDataArrayIntKeepSelectedComponents1(self):
#! [SnippeDataArrayIntKeepSelectedComponents1_1]
        arr1=[1,2,3,4,     # tuple 0
              11,12,13,14, # tuple 1
              21,22,23,24, # 
              31,32,33,34,
              41,42,43,44]
        a1=DataArrayInt()
        a1.setValues(arr1,5,4)
        a1.setInfoOnComponent(0,"a")
        a1.setInfoOnComponent(1,"b")
        a1.setInfoOnComponent(2,"c")
        a1.setInfoOnComponent(3,"d")
#! [SnippeDataArrayIntKeepSelectedComponents1_1]
#! [SnippeDataArrayIntKeepSelectedComponents1_2]
        arr2V=[1,2,1,2,0,0]
        a2=a1.keepSelectedComponents(arr2V)
#! [SnippeDataArrayIntKeepSelectedComponents1_2]
#! [SnippeDataArrayIntKeepSelectedComponents1_3]
        a3=a1[:,arr2V ]
#! [SnippeDataArrayIntKeepSelectedComponents1_3]
        return

    def testExampleFieldDoubleBuildSubPart1(self):
        from MEDCouplingDataForTest import MEDCouplingDataForTest
#! [PySnippetFieldDoubleBuildSubPart1_1]
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1()
        f1=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
        f1.setTime(2.3,5,6)
        f1.setMesh(mesh1)
        arr1=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.]
        array=DataArrayDouble(arr1,mesh1.getNumberOfCells(),2)
        f1.setArray(array)
# ! [PySnippetFieldDoubleBuildSubPart1_1]
# ! [PySnippetFieldDoubleBuildSubPart1_2]
        part1=[2,1,4]
        f2=f1.buildSubPart(part1)
# ! [PySnippetFieldDoubleBuildSubPart1_2]
        f2.zipCoords()
        self.assertEqual(3,f2.getNumberOfTuples())
        self.assertEqual(2,f2.getNumberOfComponents())
        expected1=[5.,105.,4.,104.,7.,107.]
        for i in range(6):
            self.assertAlmostEqual(f2.getIJ(0,i),expected1[i],12)
            pass
        self.assertEqual(3,f2.getMesh().getNumberOfCells())
        self.assertEqual(6,f2.getMesh().getNumberOfNodes())
        self.assertEqual(2,f2.getMesh().getSpaceDimension())
        self.assertEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.assertEqual(13,m2C.getNodalConnectivityArrayLen())
        expected2=[0.2, -0.3, 0.7, -0.3, 0.2, 0.2, 0.7, 0.2, 0.2, 0.7, 0.7, 0.7]
        for i in range(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        expected3=[3,2,3,1,3,0,2,1,4,4,5,3,2]
        self.assertEqual(expected3,list(m2C.getNodalConnectivity().getValues()))
        expected4=[0,4,8,13]
        self.assertEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()))
        # Test with field on nodes.
# ! [PySnippetFieldDoubleBuildSubPart1_3]
        f1=MEDCouplingFieldDouble(ON_NODES,ONE_TIME)
        f1.setTime(2.3,5,6)
        f1.setMesh(mesh1)
        arr2=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.]
        array=DataArrayDouble(arr2,mesh1.getNumberOfNodes(),2)
        f1.setArray(array)
# ! [PySnippetFieldDoubleBuildSubPart1_3]
# ! [PySnippetFieldDoubleBuildSubPart1_4]
        part2=[1,2]
        f2=f1.buildSubPart(part2)
# ! [PySnippetFieldDoubleBuildSubPart1_4]
        self.assertEqual(4,f2.getNumberOfTuples())
        self.assertEqual(2,f2.getNumberOfComponents())
        expected5=[4.,104.,5.,105.,7.,107.,8.,108.]
        for i in range(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12)
            pass
        self.assertEqual(2,f2.getMesh().getNumberOfCells())
        self.assertEqual(4,f2.getMesh().getNumberOfNodes())
        self.assertEqual(2,f2.getMesh().getSpaceDimension())
        self.assertEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.assertEqual(8,m2C.getNodalConnectivityArrayLen())
        for i in range(8):  # 8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.assertEqual(expected3[:4],[int(i) for i in m2C.getNodalConnectivity()][4:])
        self.assertEqual(expected3[4:8],[int(i) for i in m2C.getNodalConnectivity()][:4])
        self.assertEqual(expected4[:3],[int(i) for i in m2C.getNodalConnectivityIndex()])
        #idem previous because nodes of cell#4 are not fully present in part3
        part3=[1,2]
        arrr=DataArrayInt()
        arrr.setValues(part3,2,1)
        f2=f1.buildSubPart(arrr)
        self.assertEqual(4,f2.getNumberOfTuples())
        self.assertEqual(2,f2.getNumberOfComponents())
        for i in range(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12)
            pass
        self.assertEqual(2,f2.getMesh().getNumberOfCells())
        self.assertEqual(4,f2.getMesh().getNumberOfNodes())
        self.assertEqual(2,f2.getMesh().getSpaceDimension())
        self.assertEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.assertEqual(8,m2C.getNodalConnectivityArrayLen())
        for i in range(8):  # 8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.assertEqual(expected3[:4],[int(i) for i in m2C.getNodalConnectivity()][4:8])
        self.assertEqual(expected3[4:8],[int(i) for i in m2C.getNodalConnectivity()][:4])
        self.assertEqual(expected4[:3],m2C.getNodalConnectivityIndex().getValues())
        part4=[1,2,4]
        f2=f1.buildSubPart(part4)
        self.assertEqual(6,f2.getNumberOfTuples())
        self.assertEqual(2,f2.getNumberOfComponents())
        expected6=[4.,104.,5.,105.,7.,107.,8.,108.,10.,110.,11.,111.]
        for i in range(12):
            self.assertAlmostEqual(f2.getIJ(0,i),expected6[i],12)
            pass
        self.assertEqual(3,f2.getMesh().getNumberOfCells())
        self.assertEqual(6,f2.getMesh().getNumberOfNodes())
        self.assertEqual(2,f2.getMesh().getSpaceDimension())
        self.assertEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.assertEqual(13,m2C.getNodalConnectivityArrayLen())
        for i in range(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.assertEqual(expected3[0:4],m2C.getNodalConnectivity().getValues()[4:8])
        self.assertEqual(expected3[4:8],m2C.getNodalConnectivity().getValues()[0:4])
        self.assertEqual(expected3[8:13],m2C.getNodalConnectivity().getValues()[8:13])
        self.assertEqual(expected4,m2C.getNodalConnectivityIndex().getValues())
        # previous line equivalent to
        self.assertEqual(expected4,[int(i) for i in m2C.getNodalConnectivityIndex()])
        return

    def testExampleUMeshStdBuild1(self):
# ! [PySnippetUMeshStdBuild1_1]
        coords=[-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                 0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. ]
        nodalConnPerCell=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
# ! [PySnippetUMeshStdBuild1_1]
# ! [PySnippetUMeshStdBuild1_2]
        mesh=MEDCouplingUMesh("My2DMesh",2)
# ! [PySnippetUMeshStdBuild1_2]
# ! [PySnippetUMeshStdBuild1_3]
        mesh.allocateCells(5)#You can put more than 5 if you want but not less.
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[:4])
        mesh.insertNextCell(NORM_TRI3,nodalConnPerCell[4:7])
        mesh.insertNextCell(NORM_TRI3,nodalConnPerCell[7:10])
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[10:14])
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[14:])
        mesh.finishInsertingCells()
# ! [PySnippetUMeshStdBuild1_3]
# ! [PySnippetUMeshStdBuild1_4]
        coordsArr=DataArrayDouble(coords,9,3)#here coordsArr are declared to have 3 components, mesh will deduce that its spaceDim==3. 
        mesh.setCoords(coordsArr)#coordsArr contains 9 tuples, that is to say mesh contains 9 nodes.
# ! [PySnippetUMeshStdBuild1_4]
# ! [PySnippetUMeshStdBuild1_5]
# ! [PySnippetUMeshStdBuild1_5]
        mesh.checkConsistencyLight()
        return

    def testExampleCMeshStdBuild1(self):
# ! [PySnippetCMeshStdBuild1_1]
        XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22] # 9 values along X
        YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007] # 7 values along Y
        arrX=DataArrayDouble(XCoords)
        arrX.setInfoOnComponent(0,"X [m]")
        arrY=DataArrayDouble(YCoords)
        arrY.setInfoOnComponent(0,"Y [m]")
# ! [PySnippetCMeshStdBuild1_1]
# ! [PySnippetCMeshStdBuild1_2]
        mesh=MEDCouplingCMesh("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetCMeshStdBuild1_2]
# ! [PySnippetCMeshStdBuild1_3]
        self.assertEqual(8*6,mesh.getNumberOfCells())
        self.assertEqual(9*7,mesh.getNumberOfNodes())
        self.assertEqual(2,mesh.getSpaceDimension())
        self.assertEqual(2,mesh.getMeshDimension())
# ! [PySnippetCMeshStdBuild1_3]
        mesh=MEDCouplingCMesh("My2D_CMesh")
# ! [PySnippetCMeshStdBuild1_2bis]
        mesh.setCoordsAt(0,arrX)
        mesh.setCoordsAt(1,arrY)
# ! [PySnippetCMeshStdBuild1_2bis]
        self.assertEqual(8*6,mesh.getNumberOfCells())
        self.assertEqual(9*7,mesh.getNumberOfNodes())
        self.assertEqual(2,mesh.getSpaceDimension())
        self.assertEqual(2,mesh.getMeshDimension())
# ! [PySnippetCMeshStdBuild1_4]
# ! [PySnippetCMeshStdBuild1_4]
        return

    def testExampleUMeshAdvBuild1(self):
# ! [PySnippetUMeshAdvBuild1_1]
        coords=[-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                 0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. ]
        nodalConnPerCell=[4,0,3,4,1, 3,1,4,2, 3,4,5,2, 4,6,7,4,3, 4,7,8,5,4]
        nodalConnPerCellIndex=[0,5,9,13,18,23]
# ! [PySnippetUMeshAdvBuild1_1]
# ! [PySnippetUMeshAdvBuild1_2]
        mesh=MEDCouplingUMesh("My2DMesh",2)
# ! [PySnippetUMeshAdvBuild1_2]
# ! [PySnippetUMeshAdvBuild1_3]
        nodalConn=DataArrayInt(nodalConnPerCell,23,1)
        nodalConnI=DataArrayInt(nodalConnPerCellIndex,6,1)
        mesh.setConnectivity(nodalConn,nodalConnI,True)
# ! [PySnippetUMeshAdvBuild1_3]
# ! [PySnippetUMeshAdvBuild1_4]
        coordsArr=DataArrayDouble(coords,9,3)#here coordsArr are declared to have 3 components, mesh will deduce that its spaceDim==3.
        mesh.setCoords(coordsArr)#coordsArr contains 9 tuples, that is to say mesh contains 9 nodes.
# ! [PySnippetUMeshAdvBuild1_4]
# ! [PySnippetUMeshAdvBuild1_5]
# ! [PySnippetUMeshAdvBuild1_5]
        mesh.checkConsistencyLight()
        return

    def testExampleDataArrayBuild1(self):
# ! [PySnippetDataArrayBuild1_0]
        dataDouble=[0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.]
# ! [PySnippetDataArrayBuild1_0]
# ! [PySnippetDataArrayBuild1_1]
        arrayDouble=DataArrayDouble()
        arrayDouble.setValues(dataDouble,5,3)# 5 tuples containing each 3 components
# ! [PySnippetDataArrayBuild1_1]
# ! [PySnippetDataArrayBuild1_1bis]
        arrayDouble=DataArrayDouble(dataDouble,5,3)
# ! [PySnippetDataArrayBuild1_1bis]
# ! [PySnippetDataArrayBuild1_2]
        dataInt=[0, 10, 20, 1, 11, 21, 2, 12, 22, 3, 13, 23, 4, 14, 24]
# ! [PySnippetDataArrayBuild1_2]
# ! [PySnippetDataArrayBuild1_3]
        arrayInt=DataArrayInt()
        arrayInt.setValues(dataInt,5,3)# 5 tuples containing each 3 components
# ! [PySnippetDataArrayBuild1_3]
# ! [PySnippetDataArrayBuild1_3bis]
        arrayInt=DataArrayInt(dataInt,5,3)
# ! [PySnippetDataArrayBuild1_3bis]
        return

    def testExampleFieldDoubleBuild1(self):
        XCoords=[-0.3,0.07,0.1,0.3,0.45,0.47,0.49,1.,1.22];  arrX=DataArrayDouble(XCoords)
        YCoords=[0.07,0.1,0.37,0.45,0.47,0.49,1.007]; arrY=DataArrayDouble(YCoords)
        mesh=MEDCouplingCMesh("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetFieldDoubleBuild1_1]
        fieldOnCells=MEDCouplingFieldDouble(ON_CELLS,NO_TIME)
        fieldOnCells.setName("MyTensorFieldOnCellNoTime")
        fieldOnCells.setMesh(mesh)
        array=DataArrayDouble()
        array.alloc(fieldOnCells.getMesh().getNumberOfCells(),9) # Implicitely fieldOnCells will be a 9 components field.
        array.fillWithValue(7.)
        fieldOnCells.setArray(array)
        # fieldOnCells is now usable
        # ...
# ! [PySnippetFieldDoubleBuild1_1]
# ! [PySnippetFieldDoubleBuild1_2]
        f1=mesh.fillFromAnalytic(ON_CELLS,1,"x*x+y*y*3+2.*x") # f1 is scalar
        f2=mesh.fillFromAnalytic(ON_CELLS,1,"cos(x+y/x)") # f2 is scalar too
        f2bis=mesh.fillFromAnalytic(ON_CELLS,2,"x*x*IVec+3*y*JVec") # f2bis is a vectors field
        f3=f1+f2 # f3 scalar
        f4=f3/f2 # f4 scalar
        f2bis.applyFunc(1,"sqrt(x*x+y*y)") # f2bis becomes scalar
        f5=f2bis*f4 # f5 scalar
        pos1=[0.48,0.38]
        res=f4.getValueOn(pos1) # f4 is scalar so the returned value is of size 1.
        # ...
# ! [PySnippetFieldDoubleBuild1_2]
# ! [PySnippetFieldDoubleBuild1_3]
# ! [PySnippetFieldDoubleBuild1_3]
        return

    def testExampleFieldDoubleBuild2(self):
        XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22];  arrX=DataArrayDouble(XCoords)
        YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007];  arrY=DataArrayDouble(YCoords)
        mesh=MEDCouplingCMesh("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetFieldDoubleBuild2_1]
        fieldOnNodes=MEDCouplingFieldDouble(ON_NODES,NO_TIME)
        fieldOnNodes.setName("MyScalarFieldOnNodeNoTime")
        fieldOnNodes.setMesh(mesh)
        array=DataArrayDouble()
        array.alloc(fieldOnNodes.getMesh().getNumberOfNodes(),1) # Implicitely fieldOnNodes will be a 1 component field.
        array.fillWithValue(7.)
        fieldOnNodes.setArray(array)
        # fieldOnNodes is now usable
        # ...
# ! [PySnippetFieldDoubleBuild2_1]
        return

    def testExampleFieldDoubleBuild3(self):
        XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22];  arrX=DataArrayDouble(XCoords)
        YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007];  arrY=DataArrayDouble(YCoords)
        mesh=MEDCouplingCMesh("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetFieldDoubleBuild3_1]
        fieldOnCells=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
        fieldOnCells.setName("MyTensorFieldOnCellNoTime")
        fieldOnCells.setTimeUnit("ms") # Time unit is ms.
        fieldOnCells.setTime(4.22,2,-1) # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
        fieldOnCells.setMesh(mesh)
        array=DataArrayDouble()
        array.alloc(fieldOnCells.getMesh().getNumberOfCells(),2) # Implicitely fieldOnCells will be a 2 components field.
        array.fillWithValue(7.)
        fieldOnCells.setArray(array)
        # fieldOnCells is now usable
        # ...
# ! [PySnippetFieldDoubleBuild3_1]
        return

    def testExampleFieldDoubleBuild4(self):
        XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22];  arrX=DataArrayDouble(XCoords)
        YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007];  arrY=DataArrayDouble(YCoords)
        mesh=MEDCouplingCMesh("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetFieldDoubleBuild4_1]
        fieldOnNodes=MEDCouplingFieldDouble(ON_NODES,CONST_ON_TIME_INTERVAL)
        fieldOnNodes.setName("MyVecFieldOnNodeWithConstTime")
        fieldOnNodes.setTimeUnit("ms") # Time unit is ms.
        fieldOnNodes.setStartTime(4.22,2,-1)
        fieldOnNodes.setEndTime(6.44,4,-1)# fieldOnNodes is defined in interval [4.22 ms,6.44 ms]
        fieldOnNodes.setMesh(mesh)
        array=DataArrayDouble()
        array.alloc(fieldOnNodes.getMesh().getNumberOfNodes(),3) # Implicitely fieldOnNodes will be a 3 components field.
        array.fillWithValue(7.)
        fieldOnNodes.setArray(array)
        # fieldOnNodes is now usable
        # ...
# ! [PySnippetFieldDoubleBuild4_1]
        return

    def testExampleDataArrayApplyFunc1(self):
# ! [PySnippetDataArrayApplyFunc1_1]
        d=DataArrayDouble([1.,2.,11.,12.,21.,22.,31.,41.],4,2)
        self.assertRaises(InterpKernelException,d.applyFunc,"x*y")
# ! [PySnippetDataArrayApplyFunc1_1]
# ! [PySnippetDataArrayApplyFunc1_2]
        d=DataArrayDouble([1.,2.,11.,12.,21.,22.,31.,41.],4,2)
        d1=d.applyFunc("smth*smth")
        self.assertTrue(d1.isEqual(DataArrayDouble([1.,4.,121.,144.,441.,484.,961.,1681.],4,2),1e-12))
# ! [PySnippetDataArrayApplyFunc1_2]
# ! [PySnippetDataArrayApplyFunc1_3]
        d2=d.applyFunc(2,"smth1*IVec+2*smth2*JVec")
        self.assertTrue(d2.isEqual(DataArrayDouble([1.,4.,11.,24.,21.,44.,31.,82.],4,2),1e-12))
# ! [PySnippetDataArrayApplyFunc1_3]
# ! [PySnippetDataArrayApplyFunc1_4]
        dd=DataArrayDouble([1.,4.,3.,11.,144.,13.,21.,484.,23.,31.,1024.,33.],4,3)
# ! [PySnippetDataArrayApplyFunc1_4]
# ! [PySnippetDataArrayApplyFunc1_5]
        dd1=dd.applyFunc(1,"f+sqrt(g)+h")
        self.assertTrue(dd1.isEqual(DataArrayDouble([6.,36.,66.,96.],4,1),1e-12))
# ! [PySnippetDataArrayApplyFunc1_5]
# ! [PySnippetDataArrayApplyFunc1_6]
        dd2=dd.applyFunc(1,"a+0.*b+c")
        self.assertTrue(dd2.isEqual(DataArrayDouble([4.,24.,44.,64.],4,1),1e-12))
# ! [PySnippetDataArrayApplyFunc1_6]
# ! [PySnippetDataArrayApplyFunc1_7]
        ddd=DataArrayDouble([1.,4.,3.,11.,144.,13.,21.,484.,23.,31.,1024.,33.],4,3)
        ddd.setInfoOnComponents(["Y [m]","AA [m/s]","GG [MW]"])
# ! [PySnippetDataArrayApplyFunc1_7]
# ! [PySnippetDataArrayApplyFunc1_8]
        ddd1=ddd.applyFuncCompo(1,"Y+GG")
        self.assertTrue(ddd1.isEqual(DataArrayDouble([4.,24.,44.,64.],4,1),1e-12))
# ! [PySnippetDataArrayApplyFunc1_8]
# ! [PySnippetDataArrayApplyFunc1_9]
        ddd1=ddd.applyFuncNamedCompo(1,["X","Y","Z"],"X+Z")
        self.assertTrue(ddd1.isEqual(DataArrayDouble([4.,24.,44.,64.],4,1),1e-12))
# ! [PySnippetDataArrayApplyFunc1_9]
        return

    pass

unittest.main()
