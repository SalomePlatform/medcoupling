#!/bin/bash

#
# Check that the ICoCo headers used in MEDCoupling are well synchronized with the official ICoCo version
# hosted at:
#   https://github.com/cea-trust-platform/icoco-coupling
#

rm -rf icoco-coupling
rm -rf tmp_compare

git clone https://github.com/cea-trust-platform/icoco-coupling.git

lst="ICoCo_DeclSpec.hxx ICoCoField.h ICoCoField.hxx ICoCoMEDDoubleField.h ICoCoMEDDoubleField.hxx ICoCoMEDIntField.h ICoCoMEDIntField.hxx"

mkdir tmp_compare
cd tmp_compare
for f in $lst; do
    tail -n+4 ../icoco-coupling/include/$f > "${f}_github"
    tail -n+20 ../../$f > "${f}_mc"
    diff "${f}_github" "${f}_mc"
    if [ ! $? -eq 0 ]; then
        echo "File $f is not the same in MEDCoupling repository and in official ICoCo GitHub repository!!"
        exit 1
    fi
done

cd ..
rm -rf icoco-coupling
rm -rf tmp_compare

