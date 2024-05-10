#!/bin/sh

RES=ELMMOS_USRDAT
COMPSET=IELM
MACH=pm-cpu
COMPILER=gnu
PROJECT=m3780

FORCING=gfdl-esm4
SCENARIO=historical

SRC_DIR=/global/homes/d/donghui/e3sm_surface_water
CASE_DIR=${SRC_DIR}/cime/scripts

cd ${SRC_DIR}

GIT_HASH=`git log -n 1 --format=%h`
CASE_NAME=GLOBE_Surface_Water_Projection_test_${GIT_HASH}.`date "+%Y-%m-%d-%H%M%S"`

cd ${SRC_DIR}/cime/scripts

./create_newcase -case ${CASE_DIR}/${CASE_NAME} \
-res ${RES} -mach ${MACH} -compiler ${COMPILER} -compset ${COMPSET} --project ${PROJECT}


cd ${CASE_DIR}/${CASE_NAME}

./xmlchange -file env_run.xml -id DOUT_S             -val FALSE
./xmlchange -file env_run.xml -id INFO_DBUG          -val 2

./xmlchange DATM_MODE=CLMGSWP3v1
./xmlchange CLM_USRDAT_NAME=test_r05_r05
./xmlchange LND_DOMAIN_FILE=domain.lnd.r05_oEC60to30v3.190418.nc
./xmlchange ATM_DOMAIN_FILE=domain.lnd.r05_oEC60to30v3.190418.nc
./xmlchange LND_DOMAIN_PATH=/global/cfs/cdirs/e3sm/inputdata/share/domains
./xmlchange ATM_DOMAIN_PATH=/global/cfs/cdirs/e3sm/inputdata/share/domains
./xmlchange CIME_OUTPUT_ROOT=/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/outputs

./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_END -val 2000
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_START -val 1971
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_ALIGN -val 1

./xmlchange PIO_BUFFER_SIZE_LIMIT=67108864
./xmlchange STOP_N=5,STOP_OPTION=ndays
./xmlchange JOB_QUEUE=regular
./xmlchange NTASKS=384
./xmlchange JOB_WALLCLOCK_TIME=00:30:00

./preview_namelists

cat >> user_nl_mosart << EOF
frivinp_rtm = '/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/inputdata/MOSART_global_half_20180721b.nc'
rtmhist_fincl1='RIVER_DISCHARGE_OVER_LAND_LIQ','RIVER_DISCHARGE_TO_OCEAN_LIQ','FLOODED_FRACTION','FLOODPLAIN_FRACTION'
inundflag = .true.
opt_elevprof = 1
EOF

cat >> user_nl_elm << EOF
fsurdat = '/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/inputdata/surfdata_0.5x0.5_simyr2010_c200624_surfacewater.nc'
use_modified_infil = .true.
hist_empty_htapes = .true.
hist_fincl1 = 'QOVER', 'QDRAI', 'QH2OSFC', 'QRUNOFF', 'QINFL', 'QSNOMELT',    \
'FH2OSFC', 'EFLX_LH_TOT', 'RAIN', 'ZWT', 'ZWT_PERCH','FROST_TABLE','TSA',     \
'FSNO','FSAT','TWS','FSDS','FSA','FSR','FIRE','QSNOMELT','FPSN', 'FCTR'
EOF

cat >> user_nl_datm << EOF
tintalgo = 'coszen', 'nearest', 'linear', 'linear', 'lower'
dtlimit=2.0e0,2.0e0,2.0e0,2.0e0,2.0e0
EOF

./case.setup

./case.build

./case.submit
