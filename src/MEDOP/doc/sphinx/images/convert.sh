#!/bin/sh
# Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

factor="50%"
listfiles="\
    medop-gui-aliasfield.png \
    medop-gui-result.png \
    medop-gui-selectfield.png \
    medop-gui-visufield.png"

for file in $listfiles; do
    echo "Processing file $file ..."
    bn=$(basename $file .png)
    outfile=$bn"_scale.png"
    convert -scale $factor $file $outfile
done


