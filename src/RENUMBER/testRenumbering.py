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

nbcell2dboost=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[2,6,13,11,11,13,14,12,6,5,15,13,12,14,10,4,13,15,16,14,5,1,7,15,14,16,9,10,15,7,8,16,16,8,3,9]
connectivite_index=[1,5,9,13,17,21,25,29,33,37]
conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_QUAD4)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn2dboost=(list(conn)==connectivite) # convert numpy.ndarray to list
conn_index2dboost=(list(conn_index)==connectivite_index)
Boost2D=conn2dboost and conn_index2dboost and (nbcell2dboost==9) and field
os.remove(dir_mesh+"/out_"+filename)


print "TEST 2D Metis"
method="METIS"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell2dmetis=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[12,14,10,4,2,6,13,11,11,13,14,12,16,8,3,9,5,1,7,15,15,7,8,16,14,16,9,10,6,5,15,13,13,15,16,14]
connectivite_index=[1,5,9,13,17,21,25,29,33,37]
conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_QUAD4)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn2dmetis=(list(conn)==connectivite)
conn_index2dmetis=(list(conn_index)==connectivite_index)
Metis2D=conn2dmetis and conn_index2dmetis and (nbcell2dmetis==9)
os.remove(dir_mesh+"/out_"+filename)

## *** Avec polygone ***

filename="Test2Dpoly.med"
meshname="Mesh_1"


print "TEST 2D Boost with polygons"
method="BOOST"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell2dpolyboost=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[2,5,9,10,11,10,9,12,5,6,8,9,4,11,12,16,12,9,8,13,6,1,7,8,16,12,13,15,13,8,7,14,15,13,14,3]
connectivite_index=[1,5,9,13,17,21,25,29,33,37]
conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_POLYGON)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL)
conn2dpolyboost=(list(conn)==connectivite)
conn_index2dpolyboost=(list(conn_index)==connectivite_index)
PolyBoost2D=conn2dpolyboost and conn_index2dpolyboost and (nbcell2dpolyboost==9)
os.remove(dir_mesh+"/out_"+filename)

print "TEST 2D Metis with polygons"
method="METIS"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell2dpolymetis=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[6,1,7,8,2,5,9,10,5,6,8,9,15,13,14,3,4,11,12,16,16,12,13,15,11,10,9,12,12,9,8,13,13,8,7,14]
connectivite_index=[1,5,9,13,17,21,25,29,33,37]
conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_POLYGON)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL)
conn2dpolymetis=(list(conn)==connectivite)
conn_index2dpolymetis=(list(conn_index)==connectivite_index)
PolyMetis2D=conn2dpolymetis and conn_index2dpolymetis and (nbcell2dpolymetis==9)
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
nbcell3dboost=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[23,13,5,18,27,22,14,26,17,6,13,23,25,16,22,27,27,22,14,26,24,15,7,20,9,23,18,1,21,27,26,10,25,16,22,27,19,8,15,24,2,17,23,9,12,25,27,21,21,27,26,10,11,24,20,3,12,25,27,21,4,19,24,11]
connectivite_index=[1,9,17,25,33,41,49,57,65]
conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_HEXA8)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn3dboost=(list(conn)==connectivite)
conn_index3dboost=(list(conn_index)==connectivite_index)
Boost3D=conn3dboost and conn_index3dboost and (nbcell3dboost==8)
os.remove(dir_mesh+"/out_"+filename)


print "TEST 3D Metis"
method="METIS"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
nbcell3dmetis=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[12,25,27,21,4,19,24,11,27,22,14,26,24,15,7,20,17,6,13,23,25,16,22,27,9,23,18,1,21,27,26,10,23,13,5,18,27,22,14,26,25,16,22,27,19,8,15,24,2,17,23,9,12,25,27,21,21,27,26,10,11,24,20,3]
connectivite_index=[1,9,17,25,33,41,49,57,65]
conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_HEXA8)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn3dmetis=(list(conn)==connectivite)
conn_index3dmetis=(list(conn_index)==connectivite_index)
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
nbcell3dpolyboost=m.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
connectivite=[23,13,5,18,-1,27,26,14,22,-1,23,27,22,13,-1,13,22,14,5,-1,5,14,26,18,-1,18,26,27,23,
              17,6,13,23,-1,25,27,22,16,-1,17,25,16,6,-1,6,16,22,13,-1,13,22,27,23,-1,23,27,25,17,
              27,22,14,26,-1,24,20,7,15,-1,27,24,15,22,-1,22,15,7,14,-1,14,7,20,26,-1,26,20,24,27,
              9,23,18,1,-1,21,10,26,27,-1,9,21,27,23,-1,23,27,26,18,-1,18,26,10,1,-1,1,10,21,9,
              25,16,22,27,-1,19,24,15,8,-1,25,19,8,16,-1,16,8,15,22,-1,22,15,24,27,-1,27,24,19,25,
              2,17,23,9,-1,12,21,27,25,-1,2,12,25,17,-1,17,25,27,23,-1,23,27,21,9,-1,9,21,12,2,
              21,27,26,10,-1,11,3,20,24,-1,21,11,24,27,-1,27,24,20,26,-1,26,20,3,10,-1,10,3,11,21,
              12,25,27,21,-1,4,11,24,19,-1,12,4,19,25,-1,25,19,24,27,-1,27,24,11,21,-1,21,11,4,12]
connectivite_index=[1, 30, 59, 88, 117, 146, 175, 204, 233]
conn=m.getConnectivity(MED_NODAL,MED_CELL,MED_POLYHEDRA)
conn_index=m.getConnectivityIndex(MED_NODAL,MED_CELL);
conn3dpolyboost=(connectivite==list(conn))
conn_index3dpolyboost=(connectivite_index==list(conn_index))
PolyBoost3D=(conn3dpolyboost and conn_index3dpolyboost and (nbcell3dpolyboost==8))
os.remove(dir_mesh+"/out_"+filename)


print "TEST 3D Metis with polyhedra"
method="METIS"
string_to_execute="'"+dir_renumber+" "+dir_mesh+"/"+filename+" "+meshname+" "+method+" "+dir_mesh+"/out_"+filename+"'"
eval("os.system("+string_to_execute+")")
m = MESH(MED_DRIVER,dir_mesh+"/out_"+filename,meshname)
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


