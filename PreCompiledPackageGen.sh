#!/bin/bash
# OpenSMOKE++ Interface precompiled package generator for WSL 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script generates the package of precompiled OpenSMOKE++ libraries,
# including binaries and correct dependencies for WSL 2 environments.

VERSION=0.22.0  # your version number

pkgname=opensmoke-interface-$VERSION-WSL
ldpath=$PWD/$pkgname/lib
pathpath=$PWD/$pkgname/bin
aliases="alias set-$pkgname='export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$ldpath; export PATH=\$PATH:$pathpath'"

mkdir -p $pkgname/bin
mkdir -p $pkgname/lib
touch $pkgname/aliases

for exe in OpenSMOKE_* libopensmoke.so; do 
  cp $exe $pkgname/bin 
  list=$(ldd $exe | grep "=> " | awk '{print $3}')
  cp $list $pkgname/lib 
  for library in $list; do 
    cp $(ldd $library | grep "=> " | awk '{print $3}') $pkgname/lib 
  done 
done

echo $aliases > $pkgname/aliases

# Package it all up
zip -r $pkgname.zip $pkgname
echo "--> Prebuilt package generated: " $pkgname.zip
