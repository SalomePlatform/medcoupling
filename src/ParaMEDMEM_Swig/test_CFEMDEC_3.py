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
EDF34966 : CFEMDEC test. 3 sources procs and 2 target procs. Test checks BBTreeClosest ability to catch cells even when their bbox does not intercept.
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

def FromDaToMeshOnly( pts : mc.DataArrayDouble ) -> mc.MEDCouplingUMesh:
    m = mc.MEDCouplingUMesh.Build1DMeshFromCoords(pts)
    extrudePts = mc.DataArrayDouble([(0,0),(0,1)])
    m = m.buildExtrudedMesh( mc.MEDCouplingUMesh.Build1DMeshFromCoords(extrudePts),0 )
    m.getCoords()[:len(pts)]=pts
    m.getCoords()[len(pts):]=pts
    m.changeSpaceDimension(3,0.)
    m.getCoords()[len(pts):,2] = 1.0
    return m

def FromDaToMesh( iproc : int, pts : mc.DataArrayDouble ):
    mesh = FromDaToMeshOnly( pts )
    return mesh, GenerateSrcGlobalNodeIds( iproc, mesh.getCoords() )


def GetSrcProc0() -> mc.MEDCouplingUMesh:
    return FromDaToMesh( 0, mc.DataArrayDouble( [(2,4),(2,3),(2,2),(2,1),(3,1),(4,1),(5,1),(6,1),(6,2),(6,3),(6,4)] ) )

def GetSrcProc1() -> mc.MEDCouplingUMesh:
    return FromDaToMesh( 1, mc.DataArrayDouble( [(2,9),(2,8),(2,7),(2,6),(2,5),(3,5),(4,5),(5,5),(6,5),(6,6),(6,7),(6,8),(6,9)] ) )

def GetSrcProc2() -> mc.MEDCouplingUMesh:
    return FromDaToMesh( 2, mc.DataArrayDouble( [(2,15),(2,14),(2,13),(2,12),(2,11),(2,10),(3,10),(4,10),(5,10),(6,10),(6,11),(6,12),(6,13),(6,14),(6,15)] ) )

def GetTrgProc0() -> mc.MEDCouplingUMesh:
    mesh = mc.MEDCouplingUMesh.MergeUMeshes( [
        mc.MEDCouplingUMesh.Build1DMeshFromCoords( mc.DataArrayDouble([(0,15),(0,14)]) ) ,
        mc.MEDCouplingUMesh.Build1DMeshFromCoords( mc.DataArrayDouble([(4,0),(6,0)]) )
        ] )
    mesh.changeSpaceDimension(3,0.)
    mesh.getCoords()[:,2] = -1.0
    return mesh, GenerateSrcGlobalNodeIds( 3, mesh.getCoords() )

def GetTrgProc1() -> mc.MEDCouplingUMesh:
    mesh = mc.MEDCouplingUMesh.MergeUMeshes( [
        mc.MEDCouplingUMesh.Build1DMeshFromCoords( mc.DataArrayDouble([(1,10),(1,8)]) ) ,
        mc.MEDCouplingUMesh.Build1DMeshFromCoords( mc.DataArrayDouble([(7,1),(7,4)]) ) ,
        mc.MEDCouplingUMesh.Build1DMeshFromCoords( mc.DataArrayDouble([(7,12),(7,14)]) )
        ] )
    mesh.changeSpaceDimension(3,0.)
    mesh.getCoords()[:,2] = -1.0
    return mesh, GenerateSrcGlobalNodeIds( 4, mesh.getCoords() )

def BuildFieldFromMesh( i : int, mesh : mc.MEDCouplingUMesh ) -> mc.MEDCouplingFieldDouble:
    ret = mc.MEDCouplingFieldDouble(mc.ON_NODES_FE)
    arr = mc.DataArrayDouble( mesh.getNumberOfNodes() // 2 ) ; arr.iota()
    # deal with Z's
    arr2 = arr.deepCopy() ; arr2 += 200.0
    arr = mc.DataArrayDouble.Aggregate( [arr, arr2] )
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
    rem.prepareEx(srcFt, trgFt)

procs_source = [0, 1, 2]
procs_target = [3, 4]

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
        3 : mc.DataArrayDouble( [2000.0, 2001.0, 5.0, 7.0] ),
        4 : mc.DataArrayDouble( [2005.0, 1001.0, 7.0, 10.0, 2011., 2013.] )
        }
    #
    mesh, globalNodeIds = eval( f"GetTrgProc{rank-len(procs_source)}" )()
    idec.attachLocalMesh(mesh, globalNodeIds)
    zeResu = idec.receiveFromSource()
    if rank == 4:
        MyAssert( zeResu.getArray().isEqual( ref_values[rank], 1e-12 ) )

# fmt: on
