#!/bin/bash

# Load your modules here
module load compilers/intel/parallel_studio_xe_2015/15.0.1
module load tools/openmpi/intel/1.8.4
#export PERL5LIB=/usr/lib64/perl5/:/usr/lib64/perl5/:
#export PERL5LIB=/usr/lib64/perl5/:/usr/lib64/perl5/:${PERL5LIB}


export CASEROOT=/home/201803026n-2/cesm/$CASE
#export CCSMROOT=/home/201803026n-2/MODELOS/cesm1_2_2_1
export CASE=AMIP
export EXEROOT=/scratch/201803026n-2/cesm
export RUNDIR=/scratch/201803026n-2/cesm/$CASE/run
export DIN_LOC_ROOT_CSMDATA=/scratch/201803026n-2/inputdata
export CASEDIR=/home/201803026n-2/cesm/$CASE

task=240
nstop=nyears
years=27
#export LD_LIBRARY_PATH=:/home/201803026n-2/miniconda3/lib:/home/201803026n-2/lib

rm -fR $CASEDIR

echo "la ruta de CASEDIR es: "$CASEDIR

./create_newcase -case $CASEDIR -res 0.9x1.25_0.9x1.25 -compset F_AMIP_CAM5 -mach LNS

cd $CASEDIR

echo "Estamos en: "
pwd
echo " "
echo "cambiando numero de task a: "$task

#mod_pes.LNS
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $task
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $task
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $task
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $task
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $task
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $task
./xmlchange -file env_mach_pes.xml -id TOTALPES -val $task
./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE -val 24

echo "nstop in years: "$nstop
./xmlchange -file env_run.xml -id STOP_OPTION -val $nstop

echo "years: "$years
./xmlchange -file env_run.xml -id STOP_N -val $years



./cesm_setup
echo "./\${CASE}.build"
echo "sbatch \${CASE}.run"



#./xmlchange -file env_build.xml -id OS -val LINUX
#./xmlchange -file env_build.xml -id CESMSCRATCHROOT -val ' '
#./xmlchange -file env_build.xml -id MPILIB -val ' '
#./xmlchange -file env_run.xml -id RUNDIR -val /scratch/201803026n-2/cesm/$CASE/build
#./xmlchange -file env_run.xml -id DIN_LOC_ROOT -val /scratch/201803026n-2/inputdata
#./xmlchange -file env_build.xml -id COMPILER -val intel
#./xmlchange -file env_run.xml -id EXEROOT -val /scratch/201803026n-2/cesm/$CASE
#./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE -val 138


# NETCDF_PATH:=/home/201803026n-2/util/netcdf
# SLIBS+= -L/home/201803026n-2/util/netcdf/lib -lnetcdff  -lnetcdf -lm
# ./configure --prefix=/home/201803026n-2/util/netcdf --disable-dap --disable-netcdf-4 --disable-share
# make
# make install
# export PATH=/home/201803026n-2/util/netcdf/bin:$PATH
#
# source /opt/intel/bin/compilervars.sh intel64 

'
export CC=icc
export CXX=icpc
export CFLAGS='-O3 -xHost -ip -no-prec-div -static-intel'
export CXXFLAGS='-O3 -xHost -ip -no-prec-div -static-intel'
export F77=ifort
export FC=ifort
export F90=ifort
export FFLAGS='-O3 -xHost -ip -no-prec-div -static-intel'
export CPP='icc -E'
export CXXCPP='icpc -E'
export CPPFLAGS='-I/usr/include -I/usr/include/hdf'
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'/usr/lib64/usr/local/intel/composer_xe_2013_sp1.1.106/compiler/lib/intel64/usr/lib64/hdf'
'

nohup ./caso.build > log.build &

# /home/201803026n-2/cesm/input_data/atm/waccm/phot/temp_prs_GT200nm_jpl06_c080930.nc


#SBATCH -t 240:20:00
#SBATCH -t 240:20:00
#SBATCH -t 00:20:00








