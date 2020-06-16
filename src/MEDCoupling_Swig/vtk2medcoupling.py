import vtk
from vtk.util import numpy_support
import medcoupling as mc
import numpy as np

def mesh_convertor(fileName):
    reader = vtk.vtkDataSetReader()
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
    ct=mc.DataArrayInt(np.array(ctvtk,dtype=np.int32))[:]
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

    m=mc.MEDCouplingUMesh("mesh",2)
    m.setCoords(mc.DataArrayDouble(np.array(pts,dtype=np.float64)))
    m.setConnectivity(c,ci,True)
    m.checkConsistencyLight()
    #
    return m
