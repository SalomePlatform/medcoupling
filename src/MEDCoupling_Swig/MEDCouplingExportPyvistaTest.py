from importlib.util import find_spec

import medcoupling as mc
import unittest

def has_pyvista_and_numpy():
    n = find_spec("numpy")
    p = find_spec("pyvista")
    if n is not None and p is not None:
        return True
    return False

class MEDCouplingPyvistaTest(unittest.TestCase):
    @unittest.skipIf(not has_pyvista_and_numpy(), "Missing numpy or pyvista")
    def test1DCMesh(self):
        import numpy as np
        arr = mc.DataArrayDouble(np.logspace(1.0, 10.0))
        y = mc.DataArrayDouble([0.0])
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arr, y)
        mesh = mesh.buildUnstructured()
        mesh.toPyvista()

    @unittest.skipIf(not has_pyvista_and_numpy(), "Missing numpy or pyvista")
    def test2DCMesh(self):
        import numpy as np
        arr = mc.DataArrayDouble(np.logspace(1.0, 10.0))
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arr, arr)
        mesh = mesh.buildUnstructured()
        mesh.toPyvista()

    @unittest.skipIf(not has_pyvista_and_numpy(), "Missing numpy or pyvista")
    def test2D3DCMesh(self):
        import numpy as np
        arr = mc.DataArrayDouble(np.logspace(1.0, 10.0))
        z = mc.DataArrayDouble([0.0])
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arr, arr, z)
        mesh = mesh.buildUnstructured()
        mesh.toPyvista()

    @unittest.skipIf(not has_pyvista_and_numpy(), "Missing numpy or pyvista")
    def test3DCMesh(self):
        import numpy as np
        arr = mc.DataArrayDouble(np.logspace(1.0, 10.0))
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arr, arr, arr)
        mesh = mesh.buildUnstructured()
        mesh.toPyvista()

    @unittest.skipIf(not has_pyvista_and_numpy(), "Missing numpy or pyvista")
    def test1DField(self):
        import numpy as np
        arr = mc.DataArrayDouble(np.logspace(1.0, 10.0))
        y = mc.DataArrayDouble([0.0])
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arr, y)
        mesh = mesh.buildUnstructured()
        field = mesh.getMeasureField(True)
        field.toPyvista()

    @unittest.skipIf(not has_pyvista_and_numpy(), "Missing numpy or pyvista")
    def test2DField(self):
        import numpy as np
        arr = mc.DataArrayDouble(np.logspace(1.0, 10.0))
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arr, arr)
        mesh = mesh.buildUnstructured()
        field = mesh.getMeasureField(True)
        field.toPyvista()

    @unittest.skipIf(not has_pyvista_and_numpy(), "Missing numpy or pyvista")
    def test2D3DField(self):
        import numpy as np
        arr = mc.DataArrayDouble(np.logspace(1.0, 10.0))
        z = mc.DataArrayDouble([0.0])
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arr, arr, z)
        mesh = mesh.buildUnstructured()
        field = mesh.getMeasureField(True)
        field.toPyvista()

    @unittest.skipIf(not has_pyvista_and_numpy(), "Missing numpy or pyvista")
    def test3DField(self):
        import numpy as np
        arr = mc.DataArrayDouble(np.logspace(1.0, 10.0))
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arr, arr, arr)
        mesh = mesh.buildUnstructured()
        field = mesh.getMeasureField(True)
        field.toPyvista()

if __name__ == "__main__":
    unittest.main()
