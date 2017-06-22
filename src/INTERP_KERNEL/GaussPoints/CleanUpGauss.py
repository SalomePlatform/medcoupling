# Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (EDF R&D)

import re

s1=2709
s2=2848
f=file("InterpKernelGaussCoords.cxx","r")
lines=[elt[:-1] for elt in f.readlines()[s1:s2]]
pat0=re.compile("void[\s]+GaussInfo\:\:([^\(]+)\([\s]*\)[\s]*$")
pat1=re.compile("[\s]*\{[\s]*$")
pat2=re.compile("[\s]+LOCAL_COORD_MACRO_BEGIN[\s]*\;[\s]*$")
m0=pat0.match(lines[0])
m1=pat1.match(lines[1])
m2=pat2.match(lines[2])
if (not m0) or (not m1) or (not m2):
    raise Exception("Invalid first lines")
offsetLines=3
patEnd=re.compile("[\s]+LOCAL_COORD_MACRO_END[\s]*\;[\s]*$")
mEnd=patEnd.match(lines[-1])
if not mEnd:
    raise Exception("Invalid end lines")
#
nbLines=len(lines)-4
casePat=re.compile("[\s]+case[\s]+([\d]+)\:[\s]*$")
entries=[i_x for i_x in enumerate(lines[offsetLines:-1]) if casePat.match(i_x[1])]
#
nbPts=len(entries)
if nbLines%nbPts!=0:
    raise Exception("Invalid lines nb !")
dim=nbLines/nbPts-2
if dim<1 or dim>3:
    raise Exception("Ooops invalid dim !")
entries=[(i,int(casePat.match(elt).group(1))) for i,elt in entries]
assert({elt[1] for elt in entries} == set(range(nbPts)))
#
partEndEntries=re.compile("[\s]*break[\s]*\;[\s]*$")
zePat=re.compile("[\s]+coords\[([\d]+)\][\s]*=[\s]*([\-]?[\d]+[\.]?[\d]*)[\s]*\;[\s]*$")
zeTab=(nbPts*dim)*[None]
for lineId,ptId in entries:
    endLine=lines[offsetLines+lineId+1+dim]
    assert(partEndEntries.match(endLine))
    for j in range(dim):
        curLine=lines[offsetLines+lineId+1+j]
        m=zePat.match(curLine)
        assert(m)
        assert(int(m.group(1))==j)
        zeTab[ptId*dim+j]=m.group(2)
        pass
    pass
assert(None not in zeTab)
patInit="Init"
assert(m0.group(1)[-len(patInit):]==patInit)
varName="%s_REF"%((m0.group(1)[:-len(patInit)]).upper())
print(("const double %s[%d]={%s};"%(varName,len(zeTab),", ".join(zeTab))))
for i in range(nbPts):
    print(("  case %d:"%(i)))
    for j in range(dim):
        print(("    coords[%d] = %s[%d];"%(j,varName,i*dim+j)))
        pass
    print("    break;")
        
