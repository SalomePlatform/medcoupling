#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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
import sys
if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from medcoupling import *

from math import pi, sqrt

import numpy
if sys.platform == "win32":
    import MEDCouplingCompat as MEDCoupling
else:
    import MEDCoupling

#! [UG_DataArrayNumpy_0]
# NumPy is an optional pre-requisite!
assert(MEDCoupling.MEDCouplingHasNumPyBindings())
a=numpy.arange(20,dtype=numpy.int32)
d=DataArrayInt(a) # d owns data of a
e=DataArrayInt(a) # a not owned -> e only an access to chunk of a
a1=d.toNumPyArray()
#! [UG_DataArrayNumpy_0]
