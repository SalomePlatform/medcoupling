#!/bin/bash
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

# should be run from the build directory, so that ./PerfTest is available 
# output file
#
RES_FILE=perf_OPTIMIZE

#outputs lines of form :
#"no. source elems      no. target elems    user time"
function test_pair {
    echo -n $1 | sed 's/\(PerfCyl\)\([0-9]*\)/\2/' | sed 's/\(PerfBoxT\)\([0-9]*\)/\2/' | sed 's/\(PerfBox\)\([0-9]*\)/\2/' >> $RES_FILE
    echo -n " " >> $RES_FILE
    echo -n $2 | sed 's/\(PerfCyl\)\([0-9]*\)/\2/' | sed 's/\(PerfBoxT\)\([0-9]*\)/\2/' | sed 's/\(PerfBox\)\([0-9]*\)/\2/' >> $RES_FILE
    echo -n " " >> $RES_FILE
    time -o $RES_FILE --append -f"%U" ./PerfTest $1 $2 
    echo
}

function test_box_box {
echo PerfBox PerfBox >> $RES_FILE

test_pair PerfBox1495 PerfBox1495
test_pair PerfBox2506 PerfBox2506
test_pair PerfBox5708 PerfBox5708
test_pair PerfBox13461 PerfBox13461
test_pair PerfBox30808 PerfBox30808
test_pair PerfBox47176 PerfBox47176

test_pair PerfBox1495 PerfBox2506
test_pair PerfBox1495 PerfBox5708
test_pair PerfBox1495 PerfBox13461
test_pair PerfBox1495 PerfBox30808
test_pair PerfBox1495 PerfBox47176

test_pair PerfBox2506 PerfBox5708
test_pair PerfBox2506 PerfBox13461
test_pair PerfBox2506 PerfBox30808
test_pair PerfBox2506 PerfBox47176

test_pair PerfBox5708 PerfBox13461
test_pair PerfBox5708 PerfBox30808
test_pair PerfBox5708 PerfBox47176

test_pair PerfBox13461 PerfBox30808
test_pair PerfBox13461 PerfBox47176

test_pair PerfBox30808 PerfBox47176

}

function test_cyl_cyl {
echo PerfCyl PerfCyl >> $RES_FILE

test_pair PerfCyl1047 PerfCyl1047
test_pair PerfCyl3020 PerfCyl3020
test_pair PerfCyl6556 PerfCyl6556
test_pair PerfCyl9766 PerfCyl9766
test_pair PerfCyl25745 PerfCyl25745
test_pair PerfCyl47601 PerfCyl47601

test_pair PerfCyl1047 PerfCyl3020
test_pair PerfCyl1047 PerfCyl6556
test_pair PerfCyl1047 PerfCyl9766
test_pair PerfCyl1047 PerfCyl25745
test_pair PerfCyl1047 PerfCyl47601

test_pair PerfCyl3020 PerfCyl6556
test_pair PerfCyl3020 PerfCyl9766
test_pair PerfCyl3020 PerfCyl25745
test_pair PerfCyl3020 PerfCyl47601

test_pair PerfCyl6556 PerfCyl9766
test_pair PerfCyl6556 PerfCyl25745
test_pair PerfCyl6556 PerfCyl47601

test_pair PerfCyl9766 PerfCyl25745
test_pair PerfCyl9766 PerfCyl47601

test_pair PerfCyl25745 PerfCyl47601

}

function test_box_cyl {
    echo PerfBox PerfCyl >> $RES_FILE
    test_pair PerfBox1495 PerfCyl1047
    test_pair PerfBox1495 PerfCyl3020
    test_pair PerfBox1495 PerfCyl6556
    test_pair PerfBox1495 PerfCyl9766
    test_pair PerfBox1495 PerfCyl25745
    test_pair PerfBox1495 PerfCyl47601
    
    test_pair PerfBox2506 PerfCyl1047
    test_pair PerfBox2506 PerfCyl3020
    test_pair PerfBox2506 PerfCyl6556
    test_pair PerfBox2506 PerfCyl9766
    test_pair PerfBox2506 PerfCyl25745
    test_pair PerfBox2506 PerfCyl47601

    test_pair PerfBox5708 PerfCyl1047
    test_pair PerfBox5708 PerfCyl3020
    test_pair PerfBox5708 PerfCyl6556
    test_pair PerfBox5708 PerfCyl9766
    test_pair PerfBox5708 PerfCyl25745
    test_pair PerfBox5708 PerfCyl47601

    test_pair PerfBox13461 PerfCyl1047
    test_pair PerfBox13461 PerfCyl3020
    test_pair PerfBox13461 PerfCyl6556
    test_pair PerfBox13461 PerfCyl9766
    test_pair PerfBox13461 PerfCyl25745
    test_pair PerfBox13461 PerfCyl47601

    test_pair PerfBox30808 PerfCyl1047
    test_pair PerfBox30808 PerfCyl3020
    test_pair PerfBox30808 PerfCyl6556
    test_pair PerfBox30808 PerfCyl9766
    test_pair PerfBox30808 PerfCyl25745
    test_pair PerfBox30808 PerfCyl47601
    
    test_pair PerfBox47176 PerfCyl1047
    test_pair PerfBox47176 PerfCyl3020
    test_pair PerfBox47176 PerfCyl6556
    test_pair PerfBox47176 PerfCyl9766
    test_pair PerfBox47176 PerfCyl25745
    test_pair PerfBox47176 PerfCyl47601
}

function test_box_transbox {
    echo PerfBox PerfBoxT >> $RES_FILE
    test_pair PerfBox1495 PerfBoxT1493
    test_pair PerfBox2506 PerfBoxT2676
    test_pair PerfBox5708 PerfBoxT5717
    test_pair PerfBox13461 PerfBoxT12469
    test_pair PerfBox30808 PerfBoxT29019
    test_pair PerfBox47176 PerfBoxT47278
}
    

    
#functions to execute :     

echo PerfTest execution on `date` > $RES_FILE
test_box_cyl
test_box_box
test_cyl_cyl
test_box_transbox

cat $RES_FILE