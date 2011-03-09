#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
#  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

from libMEDMEM_Swig import *
import os
import sys

##                                         ***************
##                                         *** TEST 2D ***
##                                         ***************
srcdir   = os.getenv("srcdir")
med_root = os.getenv("MED_ROOT_DIR")
if srcdir:
    # make check is being performed
    dir_renumber="./renumber"
    dir_mesh = os.path.join( srcdir, "../../resources")
elif med_root:
    # hope renumber has been already installed
    dir_renumber=os.path.join( med_root, "bin/salome/renumber")
    dir_mesh = os.path.join( med_root, "share/salome/resources/med")
else:
    # initial version
    dir_renumber="../../../MED_INSTALL/bin/salome/renumber"
    dir_mesh="../../resources"

filename="Test2D.med"
meshname="Mesh_1"

print "TEST 2D Boost"
method="BOOST"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)

field_ini=[2,3,12,13,14,15,4,5,6,7,8,9,16,17,0,1,10,11]
s = m.getSupportOnAll(MED_CELL)
f = FIELDDOUBLE(s,2)
id=f.addDriver(MED_DRIVER,dir_mesh+"/out_"+filename,"Test field")
f.read(id);
field=True
for i in range(9):
    field=field&(f.getValueIJ(i+1,1)==field_ini[i*2])
    field=field&(f.getValueIJ(i+1,2)==field_ini[i*2+1])
f.rmDriver(id)

nbcell2dboost=m.getNumberOfElementsWithPoly(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[2,6,13,11,11,13,14,12,6,5,15,13,12,14,10,4,13,15,16,14,5,1,7,15,14,16,9,10,15,7,8,16,16,8,3,9]
connectivite_index=[1,5,9,13,17,21,25,29,33,37]
conn=m.getConnectivity(MED_FULL_INTERLACE,MED_NODAL,MED_CELL,MED_QUAD4)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn2dboost=(len(conn)==len(connectivite))
if conn2dboost:
    for i in range(0,len(connectivite)):
        conn2dboost=conn2dboost&(conn[i]==connectivite[i])
conn_index2dboost=(len(conn_index)==len(connectivite_index))
if conn_index2dboost:
    for i in range(0,len(connectivite_index)):
        conn_index2dboost=conn_index2dboost&(conn_index[i]==connectivite_index[i])
Boost2D=conn2dboost&conn_index2dboost&(nbcell2dboost==9)&field
os.remove(dir_mesh+"/out_"+filename)


print "TEST 2D Metis"
method="METIS"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell2dmetis=m.getNumberOfElementsWithPoly(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[12,14,10,4,2,6,13,11,11,13,14,12,16,8,3,9,5,1,7,15,15,7,8,16,14,16,9,10,6,5,15,13,13,15,16,14]
connectivite_index=[1,5,9,13,17,21,25,29,33,37]
conn=m.getConnectivity(MED_FULL_INTERLACE,MED_NODAL,MED_CELL,MED_QUAD4)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn2dmetis=(len(conn)==len(connectivite))
if conn2dmetis:
    for i in range(0,len(connectivite)):
        conn2dmetis=conn2dmetis&(conn[i]==connectivite[i])
conn_index2dmetis=(len(conn_index)==len(connectivite_index))
if conn_index2dmetis:
    for i in range(0,len(connectivite_index)):
        conn_index2dmetis=conn_index2dmetis&(conn_index[i]==connectivite_index[i])
Metis2D=conn2dmetis&conn_index2dmetis&(nbcell2dmetis==9)
os.remove(dir_mesh+"/out_"+filename)

## *** Avec polygone ***

filename="Test2Dpoly.med"
meshname="Mesh_1"


print "TEST 2D Boost with polygons"
method="BOOST"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell2dpolyboost=m.getNumberOfElementsWithPoly(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[2,5,9,10,11,10,9,12,5,6,8,9,4,11,12,16,12,9,8,13,6,1,7,8,16,12,13,15,13,8,7,14,15,13,14,3]
connectivite_index=[1,5,9,13,17,21,25,29,33,37]
conn=m.getPolygonsConnectivity(MED_FULL_INTERLACE,MED_CELL)
conn_index=m.getPolygonsConnectivityIndex(MED_FULL_INTERLACE,MED_CELL);
conn2dpolyboost=(len(conn)==len(connectivite))
if conn2dpolyboost:
    for i in range(0,len(connectivite)):
        conn2dpolyboost=conn2dpolyboost&(conn[i]==connectivite[i])
conn_index2dpolyboost=(len(conn_index)==len(connectivite_index))
if conn_index2dpolyboost:
    for i in range(0,len(connectivite_index)):
        conn_index2dpolyboost=conn_index2dpolyboost&(conn_index[i]==connectivite_index[i])
PolyBoost2D=conn2dpolyboost&conn_index2dpolyboost&(nbcell2dpolyboost==9)
os.remove(dir_mesh+"/out_"+filename)

print "TEST 2D Metis with polygons"
method="METIS"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell2dpolymetis=m.getNumberOfElementsWithPoly(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[6,1,7,8,2,5,9,10,5,6,8,9,15,13,14,3,4,11,12,16,16,12,13,15,11,10,9,12,12,9,8,13,13,8,7,14]
connectivite_index=[1,5,9,13,17,21,25,29,33,37]
conn=m.getPolygonsConnectivity(MED_FULL_INTERLACE,MED_CELL)
conn_index=m.getPolygonsConnectivityIndex(MED_FULL_INTERLACE,MED_CELL);
conn2dpolymetis=(len(conn)==len(connectivite))
if conn2dpolymetis:
    for i in range(0,len(connectivite)):
        conn2dpolymetis=conn2dpolymetis&(conn[i]==connectivite[i])
conn_index2dpolymetis=(len(conn_index)==len(connectivite_index))
if conn_index2dpolymetis:
    for i in range(0,len(connectivite_index)):
        conn_index2dpolymetis=conn_index2dpolymetis&(conn_index[i]==connectivite_index[i])
PolyMetis2D=conn2dpolymetis&conn_index2dpolymetis&(nbcell2dpolymetis==9)
os.remove(dir_mesh+"/out_"+filename)


##                                         ***************
##                                         *** TEST 3D ***
##                                         ***************


filename="Test3D.med"
meshname="Mesh_1"


print "TEST 3D Boost"
method="BOOST"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell3dboost=m.getNumberOfElementsWithPoly(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[23,13,5,18,27,22,14,26,17,6,13,23,25,16,22,27,27,22,14,26,24,15,7,20,9,23,18,1,21,27,26,10,25,16,22,27,19,8,15,24,2,17,23,9,12,25,27,21,21,27,26,10,11,24,20,3,12,25,27,21,4,19,24,11]
connectivite_index=[1,9,17,25,33,41,49,57,65]
conn=m.getConnectivity(MED_FULL_INTERLACE,MED_NODAL,MED_CELL,MED_HEXA8)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn3dboost=(len(conn)==len(connectivite))
if conn3dboost:
    for i in range(0,len(connectivite)):
        conn3dboost=conn3dboost&(conn[i]==connectivite[i])
conn_index3dboost=(len(conn_index)==len(connectivite_index))
if conn_index3dboost:
    for i in range(0,len(connectivite_index)):
        conn_index3dboost=conn_index3dboost&(conn_index[i]==connectivite_index[i])
Boost3D=conn3dboost&conn_index3dboost&(nbcell3dboost==8)
os.remove(dir_mesh+"/out_"+filename)


print "TEST 3D Metis"
method="METIS"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell3dmetis=m.getNumberOfElementsWithPoly(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[12,25,27,21,4,19,24,11,27,22,14,26,24,15,7,20,17,6,13,23,25,16,22,27,9,23,18,1,21,27,26,10,23,13,5,18,27,22,14,26,25,16,22,27,19,8,15,24,2,17,23,9,12,25,27,21,21,27,26,10,11,24,20,3]
connectivite_index=[1,9,17,25,33,41,49,57,65]
conn=m.getConnectivity(MED_FULL_INTERLACE,MED_NODAL,MED_CELL,MED_HEXA8)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn3dmetis=(len(conn)==len(connectivite))
if conn3dmetis:
    for i in range(0,len(connectivite)):
        conn3dmetis=conn3dmetis&(conn[i]==connectivite[i])
conn_index3dmetis=(len(conn_index)==len(connectivite_index))
if conn_index3dmetis:
    for i in range(0,len(connectivite_index)):
        conn_index3dmetis=conn_index3dmetis&(conn_index[i]==connectivite_index[i])
Metis3D=conn3dmetis&conn_index3dmetis&(nbcell3dmetis==8)
os.remove(dir_mesh+"/out_"+filename)


## *** Avec polyedres ***

## 23,13,5,18,27,26,14,22,23,27,22,13,13,22,14,5,5,14,26,18,18,26,27,23,
## 21,27,26,10,11,3,20,24,21,11,24,27,27,24,20,26,26,20,3,10,10,3,11,21,
## 12,25,27,21,4,11,24,19,12,4,19,25,25,19,24,27,27,24,11,21,21,11,4,12,
## 9,23,18,1,21,10,26,27,9,21,27,23,23,27,26,18,18,26,10,1,1,10,21,9,
## 2,17,23,9,12,21,27,25,2,12,25,17,17,25,27,23,23,27,21,9,9,21,12,2,
## 25,16,22,27,19,24,15,8,25,19,8,16,16,8,15,22,22,15,24,27,27,24,19,25,
## 17,6,13,23,25,27,22,16,17,25,16,6,6,16,22,13,13,22,27,23,23,27,25,17,
## 27,22,14,26,24,20,7,15,27,24,15,22,22,15,7,14,14,7,20,26,26,20,24,27,


filename="Test3Dpoly.med"
meshname="Mesh_1"


print "TEST 3D Boost with polyhedra"
method="BOOST"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell3dpolyboost=m.getNumberOfElementsWithPoly(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[23,13,5,18,27,26,14,22,23,27,22,13,13,22,14,5,5,14,26,18,18,26,27,23,17,6,13,23,25,27,22,16,17,25,16,6,6,16,22,13,13,22,27,23,23,27,25,17,27,22,14,26,24,20,7,15,27,24,15,22,22,15,7,14,14,7,20,26,26,20,24,27,9,23,18,1,21,10,26,27,9,21,27,23,23,27,26,18,18,26,10,1,1,10,21,9,25,16,22,27,19,24,15,8,25,19,8,16,16,8,15,22,22,15,24,27,27,24,19,25,2,17,23,9,12,21,27,25,2,12,25,17,17,25,27,23,23,27,21,9,9,21,12,2,21,27,26,10,11,3,20,24,21,11,24,27,27,24,20,26,26,20,3,10,10,3,11,21,12,25,27,21,4,11,24,19,12,4,19,25,25,19,24,27,27,24,11,21,21,11,4,12]
connectivite_face_index=[1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93,97,101,105,109,113,117,121,125,129,133,137,141,145,149,153,157,161,165,169,173,177,181,185,189,193]
connectivite_index=[1,7,13,19,25,31,37,43,49]
conn=m.getPolyhedronConnectivity(MED_FULL_INTERLACE)
conn_index=m.getPolyhedronIndex(MED_FULL_INTERLACE)
conn_face_index=m.getPolyhedronFacesIndex()
conn3dpolyboost=(len(conn)==len(connectivite))
if conn3dpolyboost:
    for i in range(0,len(connectivite)):
        conn3dpolyboost=conn3dpolyboost&(conn[i]==connectivite[i])
conn_index3dpolyboost=(len(conn_index)==len(connectivite_index))
if conn3dpolyboost:
    for i in range(0,len(connectivite_index)):
        conn_index3dpolyboost=conn_index3dpolyboost&(conn_index[i]==connectivite_index[i])
conn_face_index3dpolyboost=(len(conn_face_index)==len(connectivite_face_index))
if conn_face_index3dpolyboost:
    for i in range(0,len(connectivite_face_index)):
        conn_face_index3dpolyboost=conn_face_index3dpolyboost&(conn_face_index[i]==connectivite_face_index[i])
PolyBoost3D=conn3dpolyboost&conn_index3dpolyboost&conn_face_index3dpolyboost&(nbcell3dpolyboost==8)
os.remove(dir_mesh+"/out_"+filename)


print "TEST 3D Metis with polyhedra"
method="METIS"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell3dpolymetis=m.getNumberOfElementsWithPoly(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[12,25,27,21,4,11,24,19,12,4,19,25,25,19,24,27,27,24,11,21,21,11,4,12,27,22,14,26,24,20,7,15,27,24,15,22,22,15,7,14,14,7,20,26,26,20,24,27,17,6,13,23,25,27,22,16,17,25,16,6,6,16,22,13,13,22,27,23,23,27,25,17,9,23,18,1,21,10,26,27,9,21,27,23,23,27,26,18,18,26,10,1,1,10,21,9,23,13,5,18,27,26,14,22,23,27,22,13,13,22,14,5,5,14,26,18,18,26,27,23,25,16,22,27,19,24,15,8,25,19,8,16,16,8,15,22,22,15,24,27,27,24,19,25,2,17,23,9,12,21,27,25,2,12,25,17,17,25,27,23,23,27,21,9,9,21,12,2,21,27,26,10,11,3,20,24,21,11,24,27,27,24,20,26,26,20,3,10,10,3,11,21]
connectivite_face_index=[1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93,97,101,105,109,113,117,121,125,129,133,137,141,145,149,153,157,161,165,169,173,177,181,185,189,193]
connectivite_index=[1,7,13,19,25,31,37,43,49]
conn=m.getPolyhedronConnectivity(MED_FULL_INTERLACE)
conn_index=m.getPolyhedronIndex(MED_FULL_INTERLACE)
conn_face_index=m.getPolyhedronFacesIndex()
conn3dpolymetis=(len(conn)==len(connectivite))
conn_index3dpolymetis=(len(conn_index)==len(connectivite_index))
conn_face_index3dpolymetis=(len(conn_face_index)==len(connectivite_face_index))
if conn3dpolymetis:
    for i in range(0,len(connectivite)):
        conn3dpolymetis=conn3dpolymetis&(conn[i]==connectivite[i])
if conn_index3dpolymetis:
    for i in range(0,len(connectivite_index)):
        conn_index3dpolymetis=conn_index3dpolymetis&(conn_index[i]==connectivite_index[i])
if conn_face_index3dpolymetis:
    for i in range(0,len(connectivite_face_index)):
        conn_face_index3dpolymetis=conn_face_index3dpolymetis&(conn_face_index[i]==connectivite_face_index[i])
PolyMetis3D=conn3dpolymetis&conn_index3dpolymetis&conn_face_index3dpolymetis&(nbcell3dpolymetis==8)
os.remove(dir_mesh+"/out_"+filename)




print ""
if Boost2D:
    print "Boost 2D ok"
else:
    print "ERROR Boost 2D"
if Metis2D:
    print "Metis 2D ok"
else:
    print "ERROR Metis 2D"
if PolyBoost2D:
    print "Poly Boost 2D ok"
else:
    print "ERROR Poly Boost 2D"
if PolyMetis2D:
    print "Poly Metis 2D ok"
else:
    print "ERROR Poly Metis 2D"
if Boost3D:
    print "Boost 3D ok"
else:
    print "ERROR Boost 3D"
if Metis3D:
    print "Metis 3D ok"
else:
    print "ERROR Metis 3D"
if PolyBoost3D:
    print "Poly Boost 3D ok"
else:
    print "ERROR Poly Boost 3D"
if PolyMetis3D:
    print "Poly Metis 3D ok"
else:
    print "ERROR Poly Metis 3D"


print ""
if Boost2D&Metis2D&PolyBoost2D&PolyMetis2D&Boost3D&Metis3D&PolyBoost3D&PolyMetis3D:
    print "Every mesh correctly renumbered"
    sys.exit()
else:
    print "Error"
    sys.exit("Error in the renumbering test")


