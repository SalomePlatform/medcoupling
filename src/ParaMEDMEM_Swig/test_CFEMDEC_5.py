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
EDF36048 : CFEMDEC test. 3 sources procs and 1 target proc.
           source proc#1 is voluntarary empty. Check that BBTreeClosest is resistant to this case
"""

# fmt: off
import medcoupling as mc
from mpi4py import MPI

nbNodes = 6

def MyAssert( v ):
    if not v:
        raise RuntimeError( "Assertion failed !" )

def initParallel():
    globalComm = MPI.COMM_WORLD
    size = globalComm.size
    rank = globalComm.rank
    return rank, size

def GetSrcProc_0_2_common( rank : int):
    s = mc.DataArray.GetSlice( slice(0, nbNodes-1, 1) , rank, 2 )
    arrX = mc.DataArrayDouble(nbNodes) ; arrX.iota()
    arrY = mc.DataArrayDouble([0,1])
    arrZ = mc.DataArrayDouble([0,1])
    m = mc.MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ)
    m = m.buildUnstructured()
    return m, s

def GetSrcProc_0() -> mc.MEDCouplingUMesh:
    m, s = GetSrcProc_0_2_common(0)
    ret = m[ mc.DataArrayInt.Range(s.start,s.stop,s.step) ]
    globalNodeIds = ret.computeFetchedNodeIds()
    ret.zipCoords()
    return ret, globalNodeIds

def GetSrcProc_1() -> mc.MEDCouplingUMesh:
    m = mc.MEDCouplingUMesh("",3)
    m.setCoords( mc.DataArrayDouble(0,3) )
    m.allocateCells(0)
    m.finishInsertingCells()
    gni = mc.DataArrayInt( m.getNumberOfNodes() ) ; gni.iota()
    m2 = GetSrcProc_0_2_common(0)[0]
    gni += m2.getNumberOfNodes()
    return m, gni

def GetSrcProc_2() -> mc.MEDCouplingUMesh:
    m, s = GetSrcProc_0_2_common(1)
    ret = m[ mc.DataArrayInt.Range(s.start,s.stop,s.step) ]
    globalNodeIds = ret.computeFetchedNodeIds()
    ret.zipCoords()
    return ret, globalNodeIds

def GetTrgProc_0() -> mc.MEDCouplingUMesh:
    arrX = mc.DataArrayDouble( [ 0.5, float(nbNodes-1) - 0.5 ] )
    arrY = mc.DataArrayDouble([0,1])
    m = mc.MEDCouplingCMesh() ; m.setCoords(arrX,arrY)
    m = m.buildUnstructured()
    globalNodeIds = mc.DataArrayInt( m.getNumberOfNodes() ) ; globalNodeIds.iota()
    m.changeSpaceDimension(3,0.)
    return m, globalNodeIds

def BuildFieldFromMesh( i : int, mesh : mc.MEDCouplingUMesh ) -> mc.MEDCouplingFieldDouble:
    ret = mc.MEDCouplingFieldDouble(mc.ON_NODES_FE)
    if mesh.getNumberOfNodes() >= 1:
        arr = mesh.getCoords().magnitude()
    else:
        arr = mc.DataArrayDouble([])
    ret.setArray( arr )
    ret.setMesh( mesh )
    return ret

procs_source = [0, 1, 2]
procs_target = [3]

rank, size = initParallel()

if size != 4:
    raise RuntimeError("Expected to be lanched with 4 procs !")

# EDF36048
# emulate Code_Saturne env to trap all overflows
mc.TrapHwOverflow()

idec = mc.CFEMDEC(procs_source, procs_target)

if rank in procs_source:
    mesh, globalNodeIds  = eval( f"GetSrcProc_{rank}" )()
    idec.attachLocalMesh(mesh, globalNodeIds)
    src_field_on_local = BuildFieldFromMesh( rank, mesh )
    idec.sendToTarget(src_field_on_local)
    resu = idec.receiveFromTarget()
    expect = {
        0:  [0.5, 1.0, 2.0, 1.2071067811865475, 1.6326012547388817, 2.4835902018435503, 0.5, 1.0, 2.0, 1.2071067811865475, 1.6326012547388817, 2.4835902018435503],
        1 : [],
        2 : [2.0, 3.0, 4.0, 4.5, 2.4835902018435503, 3.3345791489482193, 4.185568096052887, 4.611062569605222, 2.0, 3.0, 4.0, 4.5, 2.4835902018435503, 3.3345791489482193, 4.185568096052887, 4.611062569605222]
        }
    MyAssert( resu.getArray().isEqual( mc.DataArrayDouble(expect[rank]), 1e-12 ) )

if rank in procs_target:
    #
    mesh, globalNodeIds = eval( f"GetTrgProc_{rank-len(procs_source)}" )()
    idec.attachLocalMesh(mesh, globalNodeIds)
    zeResu = idec.receiveFromSource()
    expect = mc.DataArrayDouble([0.5, 4.5, 1.2071067811865475, 4.611062569605222])
    MyAssert( zeResu.getArray().isEqual( expect, 1e-2 ) )
    idec.sendToSource(zeResu)

# fmt: on
