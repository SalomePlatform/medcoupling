#! /usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2020-2021  CEA/DEN, EDF R&D
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
# Author : Anthony GEAY (EDF R&D)

import vtk
from vtk.util import numpy_support
import medcoupling as mc
import numpy as np

def mesh_convertor(fileName):
    #vtk.vtkDataSetReader()
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    ug = reader.GetOutput()
    return mesh_convertor_mem(ug),ug
    
def mesh_convertor_mem(ug):
    from distutils.version import LooseVersion
    #
    pts = numpy_support.vtk_to_numpy(ug.GetPoints().GetData())
    #
    cla = numpy_support.vtk_to_numpy(ug.GetCellLocationsArray())
    ctvtk = numpy_support.vtk_to_numpy(ug.GetCellTypesArray())
    conn = numpy_support.vtk_to_numpy(ug.GetCells().GetData())
    #
    ct=mc.DataArrayInt(np.array(ctvtk,dtype="int{}".format(mc.MEDCouplingSizeOfIDs())))[:]
    c=mc.DataArrayInt(conn)[:]
    ci=mc.DataArrayInt(cla)[:]
    # for pv580
    if LooseVersion(vtk.VTK_VERSION) >= LooseVersion("8.90.0"):
        ci = ci.deltaShiftIndex()+1
        ci.computeOffsetsFull()
    #
    vtk2med = mc.DataArrayInt(mc.vtk2med_cell_types())
    #
    ct.transformWithIndArr(vtk2med)
    c[ci]=ct
    ci = mc.DataArrayInt.Aggregate([ci,mc.DataArrayInt([len(c)])])
    #
    gtv = ct.getDifferentValues().getValues()
    dimArray = mc.DataArrayInt(len(gtv)) ; dimArray[:] = -1
    for i,gt in enumerate(gtv):
        dimArray[i] = mc.MEDCouplingUMesh.GetDimensionOfGeometricType(gt)
    dim = dimArray.getMaxValueInArray()
    if not dimArray.isUniform(dim):
        raise RuntimeError("The input vtkUnstructuredGrid instance is not a single level mesh ! need to split it !")
    #
    m=mc.MEDCouplingUMesh("mesh",dim)
    m.setCoords(mc.DataArrayDouble(np.array(pts,dtype=np.float64)))
    m.setConnectivity(c,ci,True)
    m.checkConsistencyLight()
    #
    return m
