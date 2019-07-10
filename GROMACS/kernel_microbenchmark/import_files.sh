#!/bin/bash

srcdir=$HOME/git/features
subdir=src/gromacs/nbnxm/tests/kernelbenchmarks
builddir=$srcdir/build-cmake-clang-debug

for file in \
    gromacs/math/functions.h \
    gromacs/math/vectypes.h \
    gromacs/pbcutil/ishift.h \
    gromacs/utility/arrayref.h \
    gromacs/utility/basedefinitions.h \
    gromacs/utility/bitmask.h \
    gromacs/utility/current_function.h \
    gromacs/utility/gmxassert.h \
    gromacs/utility/real.h \
    ;
do
    rsync -av $srcdir/src/$file $file
done

rsync -av $srcdir/$subdir/*.cpp $srcdir/$subdir/*.h .
rsync -av $srcdir/$subdir/gromacs/nbnxm/ gromacs/nbnxm/
rsync -av $builddir/$subdir/*.h .
