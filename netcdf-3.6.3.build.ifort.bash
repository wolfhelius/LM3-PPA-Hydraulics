#!/bin/bash
set -e # exit on errors
#set -u # exit on undefines


#for arch in i386 x86_64 ; do
for arch in x86_64 ; do
   [[ -f Makefile ]] && make distclean
   export  FC=ifort
   export F77=ifort
   export F90=ifort
   export FFLAGS="-O2"
   export CFLAGS="-O2 -arch $arch"
   export CXXFLAGS="-O2 -arch $arch"
   case "$arch" in
      i386   ) source /opt/intel/composerxe-2011.3.167/bin/compilervars.sh ia32 ;;
      x86_64 ) source /opt/intel/composerxe-2011.3.167/bin/compilervars.sh intel64 ;;
   esac
   ./configure --prefix="/usr/local/ifort-2011.3.167/$arch" \
            --disable-fortran-type-check --disable-shared

#            --disable-cxx \
#            --disable-f90 \
#            --enable-docs-install \
#            --enable-large-file-tests

   make
   make check
   make install
done
