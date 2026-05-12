#!/usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2025-2026  CEA, EDF
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

"""
EDF34966 : CFEMDEC test. 3 sources procs and 2 target procs. Test focusing on emission of part of meshes across processes. Also test of robustness in case of empty mesh.
           Source proc 1 emit noting to Target proc 1.
"""

# fmt: off
import medcoupling as mc
from mpi4py import MPI

def MyAssert( v ):
    if not v:
        raise RuntimeError( "Assertion failed !" )

def initParallel():
    globalComm = MPI.COMM_WORLD
    size = globalComm.size
    rank = globalComm.rank
    return rank, size

def GenerateSrcGlobalNodeIds( iproc : int, pts : mc.DataArrayDouble ) -> mc.DataArrayInt:
    ret = mc.DataArrayInt( len(pts) )
    ret.iota()
    ret += (iproc+1) * 10**7
    return ret

def FromDaToMeshOnly( nbCells : int, startPos : float ) -> mc.MEDCouplingUMesh:
    arr = mc.DataArrayDouble(nbCells+1)
    arr.iota()
    arr /= nbCells
    m = mc.MEDCouplingCMesh()
    m.setCoords( arr )
    m = m.buildUnstructured()
    m.translate([startPos])
    #m = m.buildDescendingConnectivity()[0]
    return m

def FromDaToMesh( iproc : int, nbCells : int, startPos : float ):
    mesh = FromDaToMeshOnly( nbCells, startPos )
    return mesh, GenerateSrcGlobalNodeIds( iproc, mesh.getCoords() )

def GetSrcProc0() -> mc.MEDCouplingUMesh:
    return FromDaToMesh( 0, nbCells = 20, startPos = 0.0 )

def GetSrcProc1() -> mc.MEDCouplingUMesh:
    return FromDaToMesh( 1, nbCells = 20, startPos = 2.0 )

def GetSrcProc2() -> mc.MEDCouplingUMesh:
    return FromDaToMesh( 2, nbCells = 20, startPos = 4.0 )

def GetTrgProcHelper( nbCells : int, start : float, end : float ):
    arr = mc.DataArrayDouble(nbCells+1)
    arr.iota()
    arr /= nbCells
    arr *= ( end - start )
    arr += start
    m = mc.MEDCouplingCMesh()
    m.setCoords( arr )
    m = m.buildUnstructured()
    return m

def GetTrgProc0() -> mc.MEDCouplingUMesh:
    mesh = mc.MEDCouplingUMesh.MergeUMeshes( [
        GetTrgProcHelper( 30, 0.8, 2.2 ) ,
        GetTrgProcHelper( 8, 4.6, 5.0 )
        ] )
    return mesh, GenerateSrcGlobalNodeIds( 3, mesh.getCoords() )

def GetTrgProc1() -> mc.MEDCouplingUMesh:
    mesh = mc.MEDCouplingUMesh.MergeUMeshes( [
        GetTrgProcHelper( 30, 0.4, 1.2 ) ,
        GetTrgProcHelper( 10, 3.8, 4.2 )
        ] )
    return mesh, GenerateSrcGlobalNodeIds( 4, mesh.getCoords() )

def BuildFieldFromMesh( i : int, mesh : mc.MEDCouplingUMesh ) -> mc.MEDCouplingFieldDouble:
    ret = mc.MEDCouplingFieldDouble(mc.ON_NODES_FE)
    arr = mc.DataArrayDouble( mesh.getNumberOfNodes() ) ; arr.iota()
    arr += i * 1000
    ret.setArray( arr )
    ret.setMesh( mesh )
    return ret

def WriteZone():
    for i in range(3):
        eval( f"GetSrcProc{i}" )()[0].writeVTK( f"src_{i}.vtu" )

    for i in range(2):
        eval( f"GetTrgProc{i}" )()[0].writeVTK( f"trg_{i}.vtu" )

def sequentialCase():
    srcMesh = mc.MEDCouplingUMesh.MergeUMeshes( [eval(f"GetSrcProc{i}")()[0] for i in range(3)] )
    srcFt = mc.MEDCouplingFieldTemplate(mc.ON_NODES_FE)
    srcFt.setMesh(srcMesh)

    trgMesh = mc.MEDCouplingUMesh.MergeUMeshes( [eval(f"GetTrgProc{i}")()[0] for i in range(2)] )
    trgFt = mc.MEDCouplingFieldTemplate(mc.ON_NODES_FE)
    trgFt.setMesh(trgMesh)

    rem = mc.MEDCouplingRemapper()
    rem.setIntersectionType(mc.PointLocator)
    #rem.setProjectionOnSurfStatus( True )
    print( rem.prepareEx(srcFt, trgFt) )

    fs = []
    for i in range(3):
        mesh,_ = eval(f"GetSrcProc{i}")()
        fs.append( BuildFieldFromMesh( i, mesh ) )
    f = mc.MEDCouplingFieldDouble.MergeFields( fs )
    f.setNature( mc.IntensiveMaximum )
    fRet = rem.transferField(f, 1e300)
    return fRet


procs_source = [0, 1, 2]
procs_target = [3, 4]
#fRet = sequentialCase()

rank, size = initParallel()

if size != 5:
    raise RuntimeError("Expected to be lanched with 5 procs !")

idec = mc.CFEMDEC(procs_source, procs_target)

if rank in procs_source:

    mesh, globalNodeIds  = eval( f"GetSrcProc{rank}" )()
    idec.attachLocalMesh(mesh, globalNodeIds)
    src_field_on_local = BuildFieldFromMesh( rank, mesh )
    idec.sendToTarget(src_field_on_local)

if rank in procs_target:
    ref_values = {
        3 : mc.DataArrayDouble( [16.0, 16.93333333333333, 17.866666666666667, 18.8, 19.733333333333334, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1000.2666666666668, 1001.2, 1002.1333333333334, 1003.0666666666666, 1004.0, 2012.0, 2013.0, 2014.0, 2015.0, 2016.0, 2017.0, 2018.0, 2019.0, 2020.0] ),
        4 : mc.DataArrayDouble( [8.0, 8.533333333333333, 9.066666666666668, 9.6, 10.133333333333333, 10.666666666666666, 11.2, 11.733333333333333, 12.266666666666666, 12.8, 13.333333333333332, 13.866666666666665, 14.399999999999999, 14.933333333333332, 15.466666666666665, 16.0, 16.53333333333333, 17.066666666666666, 17.599999999999998, 18.133333333333333, 18.666666666666664, 19.2, 19.73333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2000.0, 2000.8, 2001.6, 2002.4, 2003.2, 2004.0] )
        }
    #
    mesh, globalNodeIds = eval( f"GetTrgProc{rank-len(procs_source)}" )()
    idec.attachLocalMesh(mesh, globalNodeIds)
    zeResu = idec.receiveFromSource()
    MyAssert( zeResu.getArray().isEqual( ref_values[rank], 1e-12 ) )

# fmt: on
