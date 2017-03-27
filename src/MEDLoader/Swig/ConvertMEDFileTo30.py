#! /usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

import MEDLoader as ml
import os

def ConvertTo30(nameOfMEDFile):
    fn,ext=os.path.splitext(nameOfMEDFile)
    assert(ext in [".med",".rmed"])
    realFnOut=fn+"_30.med"
    #
    initalVersion=ml.MEDFileVersionOfFileStr(nameOfMEDFile)
    #
    mfd=ml.MEDFileData(nameOfMEDFile)
    mfd.write30(realFnOut,2)
    #
    finalVersion=ml.MEDFileVersionOfFileStr(realFnOut)
    #
    print(("File \"%s\" has been converted to 3.0 successfuly ( %s -> %s ) !\nOutput file is here : \"%s\" !"%(fn,initalVersion,finalVersion,realFnOut)))
    pass

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Convert a MED file into a MED file with 3.0 version (3.0.8)')
    parser.add_argument('nameOfMEDFile', type=str, nargs=1,help='File name of the MED file to be converted into 3.0.')
    args=parser.parse_args()
    nameOfMEDFile=args.nameOfMEDFile[0]
    ConvertTo30(nameOfMEDFile)
    pass
