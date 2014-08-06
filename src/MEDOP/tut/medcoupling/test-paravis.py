# Copyright (C) 2012-2014  CEA/DEN, EDF R&D
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

import CORBA
import PARAVIS_Gen_idl
import SALOME_ContainerManager_idl

from SALOME_NamingServicePy import SALOME_NamingServicePy_i

from MEDCouplingCorba import *

def createALocalMesh():
    targetCoords=[ 0., 0., 0., 50., 0., 0. , 200., 0., 0.  , 0., 50., 0., 50., 50., 0. , 200., 50., 0.,   0., 200., 0., 50., 200., 0. , 200., 200., 0. ,
                       0., 0., 50., 50., 0., 50. , 200., 0., 50.  , 0., 50., 50., 50., 50., 50. , 200., 50., 50.,   0., 200., 50., 50., 200., 50. , 200., 200., 50. ,
                       0., 0., 200., 50., 0., 200. , 200., 0., 200.  , 0., 50., 200., 50., 50., 200. , 200., 50., 200.,   0., 200., 200., 50., 200., 200. , 200., 200., 200. ];
    targetConn=[0,1,4,3,9,10,13,12, 1,2,5,4,10,11,14,13, 3,4,7,6,12,13,16,15, 4,5,8,7,13,14,17,16,
                9,10,13,12,18,19,22,21, 10,11,14,13,19,20,23,22, 12,13,16,15,21,22,25,24, 13,14,17,16,22,23,26,25];
    targetMesh=MEDCouplingUMesh.New();
    targetMesh.setMeshDimension(3);
    targetMesh.setName("MyMesh3D");
    targetMesh.setDescription("build3DMesh");
    targetMesh.allocateCells(12);
    for i in xrange(8):
        targetMesh.insertNextCell(NORM_HEXA8,8,targetConn[8*i:8*(i+1)]);
        pass
    targetMesh.finishInsertingCells();
    myCoords=DataArrayDouble.New();
    myCoords.setValues(targetCoords,27,3);
    targetMesh.setCoords(myCoords)
    myCoords.setName("check in case")
    return targetMesh;

def createALocalField1():
    m=createALocalMesh()
    field=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
    field.setMesh(m)
    da=DataArrayDouble.New()
    da.setValues([1.,11.,101.,1001., 2.,12.,102.,1002., 3.,13.,103.,1003., 4.,14.,104.,1004., 5.,15.,105.,1005., 6.,16.,106.,1006., 7.,17.,107.,1007., 8.,18.,108.,1008.,],8,4)
    field.setArray(da)
    field.setName("vitoo")
    field.setTime(4.5,3,4)
    return field

def createALocalField2():
    m=createALocalMesh()
    field=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME)
    field.setMesh(m)
    da=DataArrayDouble.New()
    da.setValues([float(3*i) for i in xrange(27)],27,1)
    field.setArray(da)
    field.setName("vitooNode")
    field.setTime(4.7,9,14)
    return field

def createALocalMultiField3():
    fName="FieldOverTimeCorba"
    m=createALocalMesh()
    nbOfFields=100
    fs=nbOfFields*[None]
    for i in xrange(nbOfFields):
        fs[i]=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
        fs[i].setMesh(m)
        da=DataArrayDouble.New()
        da.setValues([0.,1.,2.+i,3.,4.,5.,7.],8,1)
        fs[i].setArray(da)
        fs[i].setName(fName)
        fs[i].setTime(1.2+i,9,14)
        pass
    ret=MEDCouplingFieldOverTime.New(fs);
    return ret

def createALocalCMesh4():
    mesh=MEDCouplingCMesh.New();
    coordsX=DataArrayDouble.New();
    arrX=[ -1., 1., 2., 4. ]
    coordsX.setValues(arrX,4,1);
    coordsY=DataArrayDouble.New();
    arrY=[ -2., 4., 8. ]
    coordsY.setValues(arrY,3,1);
    coordsZ=DataArrayDouble.New();
    arrZ=[ -3., 3., 6., 12., 17. ]
    coordsZ.setValues(arrZ,5,1);
    mesh.setCoords(coordsX,coordsY,coordsZ);
    mesh.setName("CMeshSample")
    return mesh

def createALocalField5():
    m=createALocalCMesh4()
    field=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
    field.setMesh(m)
    da=DataArrayDouble.New()
    field.setTime(14.5,0,0)
    da.setValues([float(7*i) for i in xrange(24)],24,1)
    field.setName("MeshOnCMesh");
    field.setArray(da)
    return field;
    

######

orb = CORBA.ORB_init([], CORBA.ORB_ID)
poa=orb.resolve_initial_references("RootPOA");
mgr=poa._get_the_POAManager();
mgr.activate();

###### Searching for 

naming_service = SALOME_NamingServicePy_i(orb)
rp=SALOME_ContainerManager_idl._0_Engines.ResourceParameters("","","",["PARAVIS"],1,10,10,1,1,"first",[])
cp=SALOME_ContainerManager_idl._0_Engines.ContainerParameters("","get","",1,False,"",rp)
sm=naming_service.Resolve("/ContainerManager")
cont=sm.GiveContainer(cp)
paraviz=naming_service.Resolve("/Containers/%s/FactoryServer/PARAVIS_inst_1"%(cont.getHostName()))

######

meshCorba=MEDCouplingUMeshServant._this(createALocalMesh())
ior=orb.object_to_string(meshCorba)
print "mesh : ",ior

f1=MEDCouplingFieldDoubleServant._this(createALocalField1())
ior2=orb.object_to_string(f1)
print "Field on cell ",ior2

f2=MEDCouplingFieldDoubleServant._this(createALocalField2())
ior3=orb.object_to_string(f2)
print "Field on node ",ior3

fs3=MEDCouplingFieldOverTimeServant._this(createALocalMultiField3())
fs3.Register()
ior4=orb.object_to_string(fs3)
print "Fields over time ",ior4

m2=MEDCouplingCMeshServant._this(createALocalCMesh4())
ior5=orb.object_to_string(m2)
print "CMesh 2 : ",ior5

f5=MEDCouplingFieldDoubleServant._this(createALocalField5())
ior6=orb.object_to_string(f5)
print "Field on cell CMesh ",ior6

script="""
src1 = ParaMEDCorbaPluginSource()
src1.IORCorba = '%s'
asc=GetAnimationScene()
rw=GetRenderView()
dr=Show()\ndr.Visibility = 1
"""

content=script%(ior4)
paraviz.ExecuteScript(content)
