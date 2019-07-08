#!/bin/bash

# This script is a convenience for making BioExcel exascale co-design
# benchmark releases
#
# The resulting tarball contains code and inputs relevant for all
# BioExcel core applications. Details are found in the respective
# README files.
#
# The sha1sum and sha256sum for the tarball is computed and reported
# so we can put that on the BioExcel webpage easily.

version=0.1.0

tar cfz BioExcel-co-design-benchmarks-$version.tar.gz GROMACS HADDOCK CP2K

echo "sha1sum values"
for tarballs in *.tar.gz
do
    sha1sum $tarballs
done

echo
echo "sha256sum values"
for tarballs in *.tar.gz
do
    sha256sum $tarballs
done
