#!/bin/bash

# This script is a convenience for making BioExcel exascale co-design
# benchmark releases. Run it from the top level directory of the git
# repository.
#
# The resulting tarball contains code and inputs relevant for all
# BioExcel core applications. Details are found in the respective
# README files.
#
# The sha1sum and sha256sum for the tarball is computed and reported
# so we can put that on the BioExcel webpage easily.

version=0.1.0
tarball_prefix=BioExcel-co-design-benchmarks
tarball_name=$tarball_prefix-$version.tgz
dir_name=bioexcel-exascale-co-design-benchmarks
(
    rm -f $tarball_prefix*.tgz
    cd ..
    tar cfz /tmp/$tarball_name --exclude=.git* $dir_name
    mv /tmp/$tarball_name $dir_name
)

echo "sha1sum values"
for tarballs in *.tgz
do
    sha1sum $tarballs
done

echo
echo "sha256sum values"
for tarballs in *.tgz
do
    sha256sum $tarballs
done
