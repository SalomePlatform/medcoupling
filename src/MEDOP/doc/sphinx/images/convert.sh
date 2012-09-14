#!/bin/sh
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


