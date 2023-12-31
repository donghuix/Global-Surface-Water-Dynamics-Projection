#!/bin/sh

RES=ELMMOS_USRDAT
COMPSET=IELM
MACH=compy
COMPILER=intel
PROJECT=esmd

SRC_DIR=/qfs/people/xudo627/e3sm_surface_water
CASE_DIR=${SRC_DIR}/cime/scripts

cd ${SRC_DIR}/cime/scripts

GIT_HASH=`git log -n 1 --format=%h`
CASE_NAME=ELMMOS_GLOBE_Surface_Water_${GIT_HASH}_cal04_compy.`date "+%Y-%m-%d-%H%M%S"`

./create_newcase \
-case ${CASE_NAME} \
-res ${RES} \
-mach ${MACH} \
-compiler ${COMPILER} \
-compset ${COMPSET} --project ${PROJECT}


cd ${CASE_DIR}/${CASE_NAME}

./xmlchange -file env_run.xml -id DOUT_S             -val FALSE
./xmlchange -file env_run.xml -id INFO_DBUG          -val 2

./xmlchange DATM_MODE=CLMGSWP3v1 #CLMNLDAS2
./xmlchange LND_DOMAIN_FILE=domain_lnd_GLOBE_cal.nc
./xmlchange ATM_DOMAIN_FILE=domain_lnd_GLOBE_cal.nc
./xmlchange LND_DOMAIN_PATH=/compyfs/xudo627/Global-Surface-Water-Dynamics-Projection/inputdata
./xmlchange ATM_DOMAIN_PATH=/compyfs/xudo627/Global-Surface-Water-Dynamics-Projection/inputdata
./xmlchange CIME_OUTPUT_ROOT=/compyfs/xudo627/Global-Surface-Water-Dynamics-Projection/outputs

./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_END -val 2014
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_START -val 1973
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_ALIGN -val 1

./xmlchange PIO_BUFFER_SIZE_LIMIT=67108864
./xmlchange STOP_N=42,STOP_OPTION=nyears
./xmlchange NTASKS=4000
./xmlchange JOB_WALLCLOCK_TIME=25:00:00

cat >> user_nl_elm << EOF
fsurdat = '/compyfs/xudo627/Global-Surface-Water-Dynamics-Projection/inputdata/surfdata_GLOBE_cal_04.nc'
use_modified_infil = .true.
hist_empty_htapes = .true.
hist_fincl1 = 'QOVER', 'QDRAI', 'QH2OSFC', 'QRUNOFF', 'QINFL', 'FH2OSFC', 'EFLX_LH_TOT', 'RAIN', 'ZWT', 'ZWT_PERCH','FROST_TABLE','TSA','FSNO','FSAT','TWS'
EOF

cat >> user_nl_mosart << EOF
frivinp_rtm = '/compyfs/xudo627/Global-Surface-Water-Dynamics-Projection/inputdata/MOSART_GLOBE_cal.nc'
rtmhist_fincl1='RIVER_DISCHARGE_OVER_LAND_LIQ','RIVER_DISCHARGE_TO_OCEAN_LIQ','FLOODED_FRACTION','FLOODPLAIN_FRACTION'
inundflag = .true.
routingmethod = 2
opt_elevprof = 1
EOF

./case.setup

./case.build

./case.submit

