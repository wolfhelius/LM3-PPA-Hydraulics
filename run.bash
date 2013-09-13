#!/bin/bash -x
#Minimal runscript for solo land single point model
set -e
set -u

forcing_dir="$PWD/../LM3-forcing/WCR"
 executable="$PWD/exec/LM3INDY.x"
    workdir="$PWD/work"
   namelist="$PWD/input.nml"

#[[ -d $workdir ]] ||  mkdir -p $workdir/{INPUT,RESTART,DIAGNOSTICS}
[[ -d $workdir ]] ||  mkdir -p $workdir/{INPUT,RESTART,DIAGNOSTICS,WOLF,WOLF/tmp}

# --- get land forcing data ---
ls -1 $forcing_dir/WCR_*.nc > "$workdir/file.table"

cd $workdir
cat > input.nml <<EOF
 &coupler_nml
     months = 120,
     days   = 0,
     seconds = 0,
     current_date = 1949,01,01,0,0,0,
     calendar = 'JULIAN'
     dt_atmos = 3600,
     dt_cpld  = 43200 
/
 &prescr_forcing_nml
     use_diurnal_solar = .true.,
     filelist = 'file.table'
     timeline = 'loop',
     start_loop = 1949,1,1,0,0,0,
     end_loop   = 2008,1,1,0,0,0,
     limit_swdn_by_TOA_flux = .true. 
/
EOF
cat $namelist >> input.nml

#------------------------------------------------------------------------------
# run the model
ulimit -n 1024 # this is necessary for Mac OS X because default max 
               # number of files is pretty small (256). Model opens a 
               # lot of output files
#for (( i=1; i<=30; i++ )) ; do 
for (( i=1; i<=2; i++ )) ; do 
   $executable 2>&1 | tee fms.out
   [[ ${PIPESTATUS[0]} == 0 ]] || die "Model failed"
   pushd RESTART
   	  if (i==1) then
	   	  cp $workdir/DIAGNOSTICS/vegn* $workdir/WOLF/.
	  else
	   	  cp $workdir/DIAGNOSTICS/vegn* $workdir/WOLF/tmp/.
	  fi
   	  #python append_files('vegn_cc*')
      tar cvf ../$((1949+i*10))0101.tar *
   popd
   mv -f RESTART/* INPUT/
done
