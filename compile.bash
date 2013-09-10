#!/bin/bash
set -e # quit on errors
set -u # quit on undefined vars

  code_dir="$PWD/src"
  exec_dir="$PWD/exec"
 pathnames="$exec_dir/path_names"
executable="LM3INDY.x"
   cppDefs="-Duse_netCDF -Duse_netCDF3"
      mkmf="$PWD/bin/mkmf"

# --- use path_names file if available ---
[[ -d $exec_dir ]] || mkdir -p $exec_dir
[[ -s $pathnames ]] ||
   find $code_dir -name '*90' -o -name '*.c' -o -name '*.inc' -o -name '*.h' | sort | sed "s:^$code_dir/::" > $pathnames

# --- execute mkmf ---
cd $exec_dir
cat > mkmf_template <<EOF
#from /home/slm/bin/mkmf.template.i686.noMPI
# template for the Intel fortran compiler
# typical use with mkmf
# mkmf -t template.ifc -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include
CPPFLAGS =
FFLAGS_BASE  = -fpp -Wp,-w -fno-alias -stack_temps -safe_cray_ptr -ftz -i_dynamic -assume byterecl -i4 -r8 -g -heap-arrays
FFLAGS_REPRO = -fltconsistency -pc64 -check bounds -check pointer -check uninit -traceback -fpe0
#FFLAGS_REPRO = -fltconsistency -check bounds -check pointer -check uninit -traceback -fpe0
FFLAGS = \$(FFLAGS_BASE) \$(FFLAGS_REPRO) -O2 -nowarn
FC = ifort
CC = cc
CFLAGS_BASE = -g -D__IFC
CFLAGS = \$(CFLAGS_BASE) -O2
LD = ifort
LDFLAGS = -L/usr/local/ifort/x86_64/lib -lnetcdf -lnetcdff
EOF

$mkmf -a "$code_dir" -t mkmf_template -p "$executable" -c "$cppDefs" \
   "$pathnames" \
   "$code_dir/shared/mpp/include" \
   "$code_dir/shared/include" \
   "/usr/local/ifort/x86_64/include" \
   || { echo "mkmf failed" ; exit 1 ; }

# --- execute make ---
make $executable || { echo "make failed" ; exit 1 ; }

echo "NOTE: make successful"
