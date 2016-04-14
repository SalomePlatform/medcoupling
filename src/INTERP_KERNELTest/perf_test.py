#! /bin/env python
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

import re
import sys
import time
import subprocess

RES_FILE="perf_OPTIMIZE"

def test_pair(par1,par2):
    c=re.compile("[\D]*(?P<num>[\d]+)")
    val1=c.match(par1).group("num")
    val2=c.match(par2).group("num")
    st="%s %s"%(val1,val2)
    ret=subprocess.call(["./PerfTest",par1,par2])
    if ret!=0:
        f=file(RES_FILE,"w") ; f.write("Error on case %s %s !!!!\n"%(par1,par2)) ; del f
        sys.exit(ret)
        pass
    return st

def test_box_box():
    st="PerfBox PerfBox\n"
    st+=test_pair("PerfBox1495","PerfBox1495")
    st+=test_pair("PerfBox2506","PerfBox2506")
    st+=test_pair("PerfBox5708","PerfBox5708")
    st+=test_pair("PerfBox13461","PerfBox13461")
    st+=test_pair("PerfBox30808","PerfBox30808")
    st+=test_pair("PerfBox47176","PerfBox47176")
    st+=test_pair("PerfBox1495","PerfBox2506")
    st+=test_pair("PerfBox1495","PerfBox5708")
    st+=test_pair("PerfBox1495","PerfBox13461")
    st+=test_pair("PerfBox1495","PerfBox30808")
    st+=test_pair("PerfBox1495","PerfBox47176")
    st+=test_pair("PerfBox2506","PerfBox5708")
    st+=test_pair("PerfBox2506","PerfBox13461")
    st+=test_pair("PerfBox2506","PerfBox30808")
    st+=test_pair("PerfBox2506","PerfBox47176")
    st+=test_pair("PerfBox5708","PerfBox13461")
    st+=test_pair("PerfBox5708","PerfBox30808")
    st+=test_pair("PerfBox5708","PerfBox47176")
    st+=test_pair("PerfBox13461","PerfBox30808")
    st+=test_pair("PerfBox13461","PerfBox47176")
    st+=test_pair("PerfBox30808","PerfBox47176")
    pass

def test_cyl_cyl():
    st="PerfCyl PerfCyl\n"
    st+=test_pair("PerfCyl1047","PerfCyl1047")
    st+=test_pair("PerfCyl3020","PerfCyl3020")
    st+=test_pair("PerfCyl6556","PerfCyl6556")
    st+=test_pair("PerfCyl9766","PerfCyl9766")
    st+=test_pair("PerfCyl25745","PerfCyl25745")
    st+=test_pair("PerfCyl47601","PerfCyl47601")
    st+=test_pair("PerfCyl1047","PerfCyl3020")
    st+=test_pair("PerfCyl1047","PerfCyl6556")
    st+=test_pair("PerfCyl1047","PerfCyl9766")
    st+=test_pair("PerfCyl1047","PerfCyl25745")
    st+=test_pair("PerfCyl1047","PerfCyl47601")
    st+=test_pair("PerfCyl3020","PerfCyl6556")
    st+=test_pair("PerfCyl3020","PerfCyl9766")
    st+=test_pair("PerfCyl3020","PerfCyl25745")
    st+=test_pair("PerfCyl3020","PerfCyl47601")
    st+=test_pair("PerfCyl6556","PerfCyl9766")
    st+=test_pair("PerfCyl6556","PerfCyl25745")
    st+=test_pair("PerfCyl6556","PerfCyl47601")
    st+=test_pair("PerfCyl9766","PerfCyl25745")
    st+=test_pair("PerfCyl9766","PerfCyl47601")
    st+=test_pair("PerfCyl25745","PerfCyl47601")
    return st

def test_box_cyl():
    st="PerfBox PerfCyl\n"
    st+=test_pair("PerfBox1495","PerfCyl1047")
    st+=test_pair("PerfBox1495","PerfCyl3020")
    st+=test_pair("PerfBox1495","PerfCyl6556")
    st+=test_pair("PerfBox1495","PerfCyl9766")
    st+=test_pair("PerfBox1495","PerfCyl25745")
    st+=test_pair("PerfBox1495","PerfCyl47601")
    st+=test_pair("PerfBox2506","PerfCyl1047")
    st+=test_pair("PerfBox2506","PerfCyl3020")
    st+=test_pair("PerfBox2506","PerfCyl6556")
    st+=test_pair("PerfBox2506","PerfCyl9766")
    st+=test_pair("PerfBox2506","PerfCyl25745")
    st+=test_pair("PerfBox2506","PerfCyl47601")
    st+=test_pair("PerfBox5708","PerfCyl1047")
    st+=test_pair("PerfBox5708","PerfCyl3020")
    st+=test_pair("PerfBox5708","PerfCyl6556")
    st+=test_pair("PerfBox5708","PerfCyl9766")
    st+=test_pair("PerfBox5708","PerfCyl25745")
    st+=test_pair("PerfBox5708","PerfCyl47601")
    st+=test_pair("PerfBox13461","PerfCyl1047")
    st+=test_pair("PerfBox13461","PerfCyl3020")
    st+=test_pair("PerfBox13461","PerfCyl6556")
    st+=test_pair("PerfBox13461","PerfCyl9766")
    st+=test_pair("PerfBox13461","PerfCyl25745")
    st+=test_pair("PerfBox13461","PerfCyl47601")
    st+=test_pair("PerfBox30808","PerfCyl1047")
    st+=test_pair("PerfBox30808","PerfCyl3020")
    st+=test_pair("PerfBox30808","PerfCyl6556")
    st+=test_pair("PerfBox30808","PerfCyl9766")
    st+=test_pair("PerfBox30808","PerfCyl25745")
    st+=test_pair("PerfBox30808","PerfCyl47601")
    st+=test_pair("PerfBox47176","PerfCyl1047")
    st+=test_pair("PerfBox47176","PerfCyl3020")
    st+=test_pair("PerfBox47176","PerfCyl6556")
    st+=test_pair("PerfBox47176","PerfCyl9766")
    st+=test_pair("PerfBox47176","PerfCyl25745")
    st+=test_pair("PerfBox47176","PerfCyl47601")
    return st

def test_box_transbox():
    st=" PerfBox PerfBoxT\n"
    st+=test_pair("PerfBox1495","PerfBoxT1493")
    st+=test_pair("PerfBox2506","PerfBoxT2676")
    st+=test_pair("PerfBox5708","PerfBoxT5717")
    st+=test_pair("PerfBox13461","PerfBoxT12469")
    st+=test_pair("PerfBox30808","PerfBoxT29019")
    st+=test_pair("PerfBox47176","PerfBoxT47278")
    return st

gm=time.strftime("%a. %B %H:%M:%S",gm).lower()+time.strftime(" CET %Y",gm)
st="PerfTest execution on %s\n"%(time.strftime("%a. %B %H:%M:%S",gm).lower()+time.strftime(" CET %Y",gm))
st+=test_box_cyl()
st+=test_box_box()
st+=test_cyl_cyl()
st+=test_box_transbox()
f=file(RES_FILE,"w")
