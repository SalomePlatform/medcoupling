#!/usr/bin/env python
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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
import re
import os

# :TRICKY:
# In input file, contents delimited by two lines containing BEGIN_CPP_ONLY and
# END_CPP_ONLY strings respectively, will only appear in C++ documentation.
# The same rule applies for Python-specific contents.

def Cpp2Python(contents):
    cpp_only = False
    output = []
    for st in contents:
        if "BEGIN_CPP_ONLY" in st:
            cpp_only = True
        elif "END_CPP_ONLY" in st:
            cpp_only = False
        elif not cpp_only:
            st=st.replace("C++","Python")
            st=st.replace("Cxx","Py")
            st=st.replace("Cpp","Py")
            st=st.replace("cxx","py")
            st=st.replace("cpp","py")
            output.append(st)
            pass
        pass

    return output
#
def discardPythonFrom(contents):
    python_only = False
    output = []
    for st in contents:
        if "BEGIN_PYTHON_ONLY" in st:
            python_only = True
        elif "END_PYTHON_ONLY" in st:
            python_only = False
        elif not python_only:
            output.append(st)
            pass
        pass

    return output
    #
#

# Usage: BuildPyExamplesFromCPP.py <examples.in> <output directory>

with open(sys.argv[1], "r") as fCpp:
    cppCont = fCpp.readlines()
pyCont=cppCont[:]
pyCont=[elt.replace("medcouplingcppexamples","medcouplingpyexamples") for elt in pyCont]
pyCont=Cpp2Python(pyCont)

s= "Be sure to take a look at the page \\ref python-api before proceeding."
pyCont.insert(2, s)
cppCont=discardPythonFrom(cppCont) # remove Python-only contents from Cpp

# Save CPP and PY examples in two separate dox files
outFileName=os.path.join(sys.argv[2],os.path.basename(sys.argv[1]))

with open(os.path.splitext(outFileName)[0] + "CPP.dox", "w") as f:
    f.writelines(cppCont)

with open(os.path.splitext(outFileName)[0] + "PY.dox", "w") as f:
    f.writelines(pyCont)
