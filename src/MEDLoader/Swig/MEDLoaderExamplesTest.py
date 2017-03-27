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
# Author : Anthony Geay (CEA/DEN)

from MEDLoader import *
import unittest
import os

class MEDLoaderBasicsTest(unittest.TestCase):
    def testExampleReadFieldOnAllEntity1(self):
        from MEDLoaderDataForTest import MEDLoaderDataForTest
#! [PySnippetReadFieldOnAllEntity1_1]
        fname="PyExamples1.med"
        meshName="mesh"
        fieldName="FieldOnAll"
        iteration=3
        order=4
#! [PySnippetReadFieldOnAllEntity1_1]
#! [PySnippetWriteFieldOnAllEntity1_2]
        m=MEDLoaderDataForTest.build2DMesh_3()
        m=m[:10]
        m.setName(meshName)
        f=m.getMeasureField(False)
        f=f.buildNewTimeReprFromThis(ONE_TIME,False)
        f.setTime(5.5,iteration,order)
        f.setName(fieldName)
        # MEDCoupling finished, MEDLoader advanced API specific part starting from here
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,m)
        ff=MEDFileField1TS.New()
        ff.setFieldNoProfileSBT(f)
        mm.write(fname,2)
        ff.write(fname,0)
#! [PySnippetWriteFieldOnAllEntity1_2]
#! [PySnippetReadFieldOnAllEntity1_3]
        medfileField1TS=MEDFileField1TS.New(fname,fieldName,iteration,order)
        mm=MEDFileMesh.New(fname)
        fread=medfileField1TS.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
        fread2=medfileField1TS.getFieldAtLevel(ON_CELLS,0)
        self.assertTrue(fread.isEqual(f,1e-12,1e-12))
        self.assertTrue(fread2.isEqual(f,1e-12,1e-12))
#! [PySnippetReadFieldOnAllEntity1_3]
#! [PySnippetReadFieldOnAllEntity1_4]
        medfileFieldMTS=MEDFileFieldMultiTS.New(fname,fieldName)
        mm=MEDFileMesh.New(fname)
        fread=medfileFieldMTS.getFieldOnMeshAtLevel(ON_CELLS,iteration,order,0,mm)
        fread2=medfileFieldMTS.getFieldAtLevel(ON_CELLS,iteration,order,0)
        self.assertTrue(fread.isEqual(f,1e-12,1e-12))
        self.assertTrue(fread2.isEqual(f,1e-12,1e-12))
#! [PySnippetReadFieldOnAllEntity1_4]
#! [PySnippetReadFieldOnAllEntity1_5]
        medfileFieldMTS=MEDFileFieldMultiTS.New(fname,fieldName)
        for medfileField1TS in medfileFieldMTS:
            if medfileField1TS.getTime()[:2]==[iteration,order]:
                fread=medfileField1TS.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
                fread2=medfileField1TS.getFieldAtLevel(ON_CELLS,0)
                self.assertTrue(fread.isEqual(f,1e-12,1e-12))
                self.assertTrue(fread2.isEqual(f,1e-12,1e-12))
                pass
            pass
#! [PySnippetReadFieldOnAllEntity1_5]
        pass

    def testExampleReadFieldPartial1(self):
        from MEDLoaderDataForTest import MEDLoaderDataForTest
#! [PySnippetReadFieldPartial1_1]
        fname="PyExamples2.med"
        meshName="mesh"
        fieldName="FieldPartial"
        iteration=3
        order=4
#! [PySnippetReadFieldPartial1_1]
#! [PySnippetWriteFieldPartial1_2]
        m=MEDLoaderDataForTest.build2DMesh_1()
        m.sortCellsInMEDFileFrmt()
        m.setName(meshName)
        # end of generation of a mesh -> let's create a field on that mesh
        f=m.getMeasureField(False)
        f=f.buildNewTimeReprFromThis(ONE_TIME,False)
        f.setTime(5.5,iteration,order)
        f.setName(fieldName)
        # The MEDCoupling part is finished -> let's perform advanced API
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,m) # the MED file data structure is ready for writing. Of course mm could have been more complicated with groups, families, multilevel
        # Let's building a sub field
        pfl=DataArrayInt.New([1,3,4,5])
        pfl.setName("myPfl") # here it is necessary to give a name to be compliant with MED file
        f=f[pfl] ; f.getMesh().setName(m.getName()) # of course f should be in coherence with pfl -> f[pfl]
        #
        ff=MEDFileField1TS.New()
        tmp=f.getMesh() # useless line, only to show that mesh into f is not considered by MEDFileField1TS.setFieldProfile
        f.setMesh(None) # useless line, only to show that mesh into f is not considered by MEDFileField1TS.setFieldProfile
        ff.setFieldProfile(f,mm,0,pfl)
        f.setMesh(tmp) # useless line, only to show that mesh into f is not considered by MEDFileField1TS.setFieldProfile
        mm.write(fname,2)
        ff.write(fname,0)
#! [PySnippetWriteFieldPartial1_2]
#! [PySnippetReadFieldPartial1_3]
        mm=MEDFileMesh.New(fname)
        medfileField1TS=MEDFileField1TS.New(fname,fieldName,iteration,order)
        fread=medfileField1TS.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
        fread2=medfileField1TS.getFieldAtLevel(ON_CELLS,0)
        self.assertTrue(fread.isEqual(f,1e-12,1e-12))
        self.assertTrue(fread2.isEqual(f,1e-12,1e-12))
#! [PySnippetReadFieldPartial1_3]
#! [PySnippetReadFieldPartial1_4]
        medfileField1TS=MEDFileField1TS.New(fname,fieldName,iteration,order)
        mm=MEDFileMesh.New(fname)
        valsRead,pflRead=medfileField1TS.getFieldWithProfile(ON_CELLS,0,mm)
        self.assertEqual(valsRead.getName(),f.getName())
        valsRead.setName("")
        self.assertTrue(valsRead.isEqual(f.getArray(),1e-12))
        pflRead.setName(pfl.getName())
        self.assertTrue(pflRead.isEqual(pfl))
#! [PySnippetReadFieldPartial1_4]
#! [PySnippetReadFieldPartial1_5]
        firstApproachMesh=fread.getMesh()
        mm=MEDFileMesh.New(fname)
        wholeMesh=mm.getMeshAtLevel(0)
        wholeMesh.tryToShareSameCoords(firstApproachMesh,1e-12)
        isIncluded,pflComputed=wholeMesh.areCellsIncludedIn(firstApproachMesh,2)
        self.assertTrue(isIncluded)
        self.assertEqual(pflComputed.getName(),mm.getName())
        pflComputed.setName(pflRead.getName())
        self.assertTrue(pflComputed.isEqual(pflRead))
#! [PySnippetReadFieldPartial1_5]
#! [PySnippetReadFieldPartial1_6]
        mm=MEDFileMesh.New(fname)
        wholeMesh=mm.getMeshAtLevel(0)
        computedMesh=wholeMesh[pflRead] ; computedMesh.setName(mm.getName())
        self.assertTrue(computedMesh.isEqual(fread.getMesh(),1e-12))
        fieldFromSecondApproach=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
        fieldFromSecondApproach.setName(medfileField1TS.getName())
        fieldFromSecondApproach.setMesh(computedMesh)
        fieldFromSecondApproach.setArray(valsRead)
        fieldFromSecondApproach.setTime(medfileField1TS.getTime()[2],medfileField1TS.getTime()[0],medfileField1TS.getTime()[1])
        self.assertTrue(fieldFromSecondApproach.isEqual(fread,1e-12,1e-12))
#! [PySnippetReadFieldPartial1_6]
        pass

    def testExampleMeshAdvAPI1(self):
        da=DataArrayDouble.New([0.,1.1,2.3,3.6])
        meshName="Example2"
        cmesh=MEDCouplingCMesh.New() ; cmesh.setCoords(da,da,da)
        myMesh=cmesh.buildUnstructured()
#! [PySnippetMeshAdvAPI1_1]
        self.assertTrue(isinstance(myMesh,MEDCouplingUMesh))
        myMesh.setName(meshName)
        WriteUMesh("wFile1.med",myMesh,True)
#! [PySnippetMeshAdvAPI1_1]
        os.remove("wFile1.med")
#! [PySnippetMeshAdvAPI1_2]
        self.assertTrue(isinstance(myMesh,MEDCouplingUMesh))
        myMesh.setName(meshName)
        WriteUMesh("wFile1.med",myMesh,False)
#! [PySnippetMeshAdvAPI1_2]
        f=myMesh.getMeasureField(False)
        f=f.buildNewTimeReprFromThis(ONE_TIME,False)
        f.setName("myField")
#! [PySnippetMeshAdvAPI1_3]
        WriteUMesh("file3.med",f.getMesh(),True)
        f.setTime(1.2,1,0)
        fileNameMultiTimeStep="file3.med"
        WriteFieldUsingAlreadyWrittenMesh(fileNameMultiTimeStep,f)
        f.setTime(1.3,2,0)
        f.applyFunc("sqrt(x)");
        #Writing second time step with iteration==2 and order==0
        WriteFieldUsingAlreadyWrittenMesh("file3.med",f);
#! [PySnippetMeshAdvAPI1_3]
#! [PySnippetMeshAdvAPI1_11]
        timeStepsIds=GetCellFieldIterations("file3.med","Example2","myField")
        self.assertEqual([(1, 0),(2, 0)],timeStepsIds)
        fs=ReadFieldsOnSameMesh(ON_CELLS,"file3.med","Example2",0,"myField",timeStepsIds);
#! [PySnippetMeshAdvAPI1_11]
        ###
        myMesh0=myMesh[:] ; myMesh0.setName("Example2")
        myMesh1=myMesh0.buildDescendingConnectivity()[0] ; myMesh1.setName(myMesh0.getName())
        myMesh2=myMesh1.buildDescendingConnectivity()[0] ; myMesh2.setName(myMesh0.getName())
        myMesh3=myMesh2.buildDescendingConnectivity()[0] ; myMesh3.setName(myMesh0.getName())
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,myMesh0)
        mm.setMeshAtLevel(-1,myMesh1)
        mm.setMeshAtLevel(-2,myMesh2)
        mm.setMeshAtLevel(-3,myMesh3)
        mm.write("file2.med",2)
        F1Cell=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
        F1Cell.setMesh(myMesh0)
        F1Cell.setArray(myMesh0.getCoords()[:myMesh0.getNumberOfCells()])
        F1Cell.setTime(1000.,2,3)
        F1Cell.setName("F1Cell")
        WriteFieldUsingAlreadyWrittenMesh("file2.med",F1Cell)
        F1Cell1=F1Cell.deepCopy()
        F1Cell1.setMesh(myMesh1)
        F1Cell1.setArray(myMesh1.computeCellCenterOfMass())
        WriteFieldUsingAlreadyWrittenMesh("file2.med",F1Cell1)
#! [PySnippetMeshAdvAPI1_12]
        f1Cell_3D=ReadFieldCell("file2.med","Example2",0,"F1Cell",2,3)
#! [PySnippetMeshAdvAPI1_12]
#! [PySnippetMeshAdvAPI1_13]
        f1Cell_2D=ReadFieldCell("file2.med","Example2",-1,"F1Cell",2,3)
#! [PySnippetMeshAdvAPI1_13]
        self.assertTrue(F1Cell.isEqual(f1Cell_3D,1e-12,1e-12))
#! [PySnippetMeshAdvAPI1_8]
        self.assertEqual(3,ReadUMeshDimFromFile("file2.med","Example2"))
#! [PySnippetMeshAdvAPI1_8]
#! [PySnippetMeshAdvAPI1_7]
        m2D=ReadUMeshFromFile("file2.med","Example2",0)
#! [PySnippetMeshAdvAPI1_7]
#! [PySnippetMeshAdvAPI1_4]
        m2D=ReadUMeshFromFile("file2.med","Example2",-1)
#! [PySnippetMeshAdvAPI1_4]
#! [PySnippetMeshAdvAPI1_5]
        m1D=ReadUMeshFromFile("file2.med","Example2",-2)
#! [PySnippetMeshAdvAPI1_5]
#! [PySnippetMeshAdvAPI1_6]
        m0D=ReadUMeshFromFile("file2.med","Example2",-3)
#! [PySnippetMeshAdvAPI1_6]
        for i in range(4):
            mm.removeMeshAtLevel(-i)
            pass
        mm.setMeshAtLevel(0,myMesh1)
        mm.setMeshAtLevel(-1,myMesh2)
        mm.setName("MyMesh")
        mm.write("file1.med",2)
#! [PySnippetMeshAdvAPI1_9]
        m2D=ReadUMeshFromFile("file1.med","MyMesh",0)
#! [PySnippetMeshAdvAPI1_9]
#! [PySnippetMeshAdvAPI1_10]
        m1D=ReadUMeshFromFile("file1.med","MyMesh",-1)
#! [PySnippetMeshAdvAPI1_10]
        pass

    pass

unittest.main()
