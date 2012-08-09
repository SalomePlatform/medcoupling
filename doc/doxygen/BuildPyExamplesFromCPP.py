#!/usr/bin/env python
# Copyright (C) 2007-2012  CEA/DEN, EDF R&D, OPEN CASCADE
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
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
import re
import os

def Cpp2Python(st):
    st=st.replace("C++","Python")
    st=st.replace("Cxx","Py")
    st=st.replace("Cpp","Py")
    st=st.replace("cxx","py")
    st=st.replace("cpp","py")
    return st

fCpp=file(sys.argv[1],"r")
cppCont=fCpp.readlines() ; del fCpp
pyCont=cppCont[:]
pyCont=[elt.replace("medcouplingcppexamples","medcouplingpyexamples") for elt in pyCont]
pyCont=[Cpp2Python(elt) for elt in pyCont]

outFileName=os.path.join(sys.argv[2],os.path.basename(sys.argv[1]))

f=file(os.path.splitext(outFileName)[0]+".dox","w")
f.writelines(cppCont+pyCont) ; del f
