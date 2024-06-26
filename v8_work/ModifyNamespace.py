# Copyright (C) 2016-2024  CEA, EDF
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

import os
rep=("namespace ParaMEDMEM","namespace MEDCoupling")
rep=("ParaMEDMEM::","MEDCoupling::")
#rep=("ParaMEDMEMImpl::","MEDCouplingImpl::")

#rep=("_ParaMEDMEM__","_MEDCoupling__")
#rep=("ParaMEDMEM_","MEDCoupling_")
#rep=("ParaMEDMEMData","MEDCouplingData")

#rep=("ParaMEDMEM_1","MEDCoupling_1")

def rep0(fi,rep):
    f=file(fi) ; lines=f.readlines() ; del f
    lines2=[line.replace(*rep) for line in lines]
    if lines2!=lines:
        f=file(fi,"w") ; f.writelines(lines2) ; f.flush()
        return 1
    else:
        return 0

def rep1(dirname,rep):
    i=0
    for fi in os.listdir(dirname):
        fi2=os.path.join(dirname,fi)
        if not os.path.isfile(fi2):
            continue
        i+=rep0(fi2,rep)
    return i

dirs=["MEDCoupling","MEDCoupling/Test","MEDLoader","MEDLoader/Swig","MEDLoader/Test","MEDPartitioner","MEDPartitioner/Test","MEDPartitioner_Swig","RENUMBER","RENUMBER_Swig","INTERP_KERNELTest","ParaMEDMEM","ParaMEDLoader","ParaMEDMEMTest","ParaMEDMEM_Swig","doc/user/doxygen/fakesources","doc/user/doxygen/doxy2swig","doc/user/doxygen/doxfiles","/home/H87074/salome/DEV/modules/src/MED/src/MEDCouplingCorba","/home/H87074/salome/DEV/modules/src/MED/src/MEDCouplingCorba/Client","/home/H87074/salome/DEV/modules/src/MED/src/MEDCouplingCorba/Test","/home/H87074/salome/DEV/modules/src/MED/src/MEDCalc/cmp","/home/H87074/salome/DEV/modules/src/MED/src/MEDCalculator","/home/H87074/salome/DEV/modules/src/MED/src/MEDCalculator/Swig","/home/H87074/salome/DEV/modules/src/MED/src/MEDCalculator/Test","/home/H87074/salome/DEV/modules/src/PARAVIS/src/Plugins/MEDReader/IO"]
dirname=dirs[-1]
i=0
print(rep1(dirname,rep))
"""for r,dirs,fis in os.walk(dirname):
    for fi in fis:
        if os.path.splitext(fi)[1] not in [".dox",".doxy"]:
            continue
        i+=rep0(os.path.join(r,fi),rep)

print i"""
