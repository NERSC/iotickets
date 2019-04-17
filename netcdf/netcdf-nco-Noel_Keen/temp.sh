#!/bin/csh

#SBATCH --account=acme
#SBATCH --job-name=nct-knl-bb-stage
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -C knl,quad,cache

#SBATCH -t 29

#SBATCH --output=slurmoutput.%j 

#DW jobdw capacity=200GB access_mode=striped type=scratch

pwd
date

# There are different path to the NCO tools, but I found this doesn't matter
source /opt/modules/default/init/csh
module unload cray-hdf5-parallel
#module load nco
#setenv PATH ${PATH}:/global/u1/n/ndk/local/nco/4.7.4/intel/bin 
setenv PATH ${PATH}:/global/u1/n/ndk/local/ncoknl/4.7.4/intel/bin 
#/global/u1/n/ndk/local/nco/4.7.4/intel/bin
#setenv PATH ${PATH}:/global/u1/z/zender/bin_edison
which ncclimo

env >& envout-bb-stage.$SLURM_JOB_ID.txt

cd $DW_JOB_STRIPED
mkdir -p $DW_JOB_STRIPED/inputdir
mkdir -p $DW_JOB_STRIPED/outputdir

#set workdir='/global/cscratch1/sd/ndk/ncclimo_test-batch-cori-haswell'
set workdir=`pwd`

# scratch spaces:
#set drc_in='/global/cscratch1/sd/minxu/archive/F_acmev03_enso_ne30_knl_cesmmach_co2cyc_pmpet_25yr_climpac/lnd/hist/' # Input directory
set drc_in='/global/cscratch1/sd/ndk/Fstripe01'
#set drc_in='/global/cscratch1/sd/ndk/Fstripe08'
#set drc_in='/global/cscratch1/sd/ndk/Fstripe32'
#lfs getstripe -c $drc_in 

# staged in BB:
#set drc_in="$DW_JOB_STRIPED/inputdir"
#cp /global/cscratch1/sd/ndk/Fstripe01/* $drc_in
#/bin/ls -l $drc_in

# stage out BB
set drc_out="$DW_JOB_STRIPED/outputdir"
#set DATA=$workdir
#set drc_out="${DATA}/ne30/clm" # Native grid output directory
#set drc_rgr="${DATA}/ne30/rgr" # Regridded output directory
#set drc_tmp='/global/cscratch1/sd/zender/tmp' # Temporary/intermediate-file directory
set drc_tmp=$workdir/tmp
#lfs getstripe -c $drc_out

#set cmip6_opt='-7 --dfl_lvl=1 --no_cll_msr --no_frm_trm --no_stg_grd' # CMIP6-specific options
set cmip6_opt='-6 --no_cll_msr --no_frm_trm --no_stg_grd' # CMIP6-specific options

set spl_opt='--yr_srt=1 --yr_end=27 --ypf=27' # 2D+3D Splitter options
set vars='ALT,AR,BTRAN,CH4PROD,DENIT,EFLX_LH_TOT,ELAI,ER,ESAI,FAREA_BURNED,FCEV,FCH4,FCH4TOCO2,FCOV,FCTR,FGEV,FGR,FGR12,FH2OSFC,FINUNDATED,FIRA,FIRE,FLDS,FPG,FPI,FPSN,FROST_TABLE,FSA,FSAT,FSDS,FSH,FSM,FSNO,FSR,F_DENIT,F_NIT,GPP,GROSS_NMIN,H2OSFC,H2OSNO,HR,HTOP,LAND_USE_FLUX,LEAFC,FROOTC,NDEP_TO_SMINN,NBP,NEE,NEP,NET_NMIN,NFIX_TO_SMINN,NPP,Q2M,QCHARGE,QDRAI,QOVER,QRUNOFF,QRGWL,QSNOMELT,QSOIL,QVEGE,QVEGT,RAIN,RH2M,SMIN_NO3,SMIN_NH4,SNOW,SNOWDP,SNOWICE,SNOWLIQ,SNOW_DEPTH,SNOW_SINKS,SNOW_SOURCES,SOMHR,TG,TSA,TSAI,TLAI,TV,QBOT,TBOT,AGNPP,FROOTC_ALLOC,LEAFC_ALLOC,WOODC_ALLOC,WOOD_HARVESTC,CH4_SURF_AERE_SAT,CH4_SURF_AERE_UNSAT,CH4_SURF_DIFF_SAT,CH4_SURF_DIFF_UNSAT,CH4_SURF_EBUL_SAT,CONC_CH4_SAT,CONC_CH4_UNSAT,FCH4_DFSAT,MR,TOTCOLCH4,ZWT_CH4_UNSAT,FSDSND,FSDSNI,FSDSVD,FSDSVI,TWS,VOLR,WA,ZWT_PERCH,ZWT,WIND,COL_FIRE_CLOSS,F_DENIT_vr,F_NIT_vr,H2OSOI,O_SCALAR,SOILICE,SOILLIQ,SOILPSI,TLAKE,TSOI,T_SCALAR,W_SCALAR,SOIL1N,SOIL2N,SOIL3N,SOIL1C,SOIL2C,SOIL3C,TOTVEGC,TOTVEGN,TOTECOSYSC,TOTLITC,TOTLITC_1m,TOTLITN_1m,TOTSOMC,TOTSOMC_1m,TOTSOMN_1m,CWDC,PBOT'
#mkdir -p ${drc_out} ${drc_rgr} ${drc_tmp}
setenv TMPDIR ${drc_tmp}
#setenv HDF5_USE_FILE_LOCKING FALSE

set outprefix="batchout-stageout"
set exp='G'

set stripe="BB"

set np='136'
set output="$outprefix.$stripe.$exp.$np.txt"
/bin/ls ${drc_in}/*.clm2.h0.25[78][0-9]-*.nc ${drc_in}/*.clm2.h0.259[0-6]-*.nc | ncclimo --var=${vars} --job_nbr=${np} -d 0 ${cmip6_opt} ${spl_opt} --drc_out=${drc_out} >& $output

#set np='072'
#set output="$outprefix.$stripe.$exp.$np.txt"
#/bin/ls ${drc_in}/*.clm2.h0.25[78][0-9]-*.nc ${drc_in}/*.clm2.h0.259[0-6]-*.nc | ncclimo --var=${vars} --job_nbr=${np} -d 0 ${cmip6_opt} ${spl_opt} --drc_out=${drc_out} >& $output

#set np='064'
#set output="$outprefix.$stripe.$exp.$np.txt"
#/bin/ls ${drc_in}/*.clm2.h0.25[78][0-9]-*.nc ${drc_in}/*.clm2.h0.259[0-6]-*.nc | ncclimo --var=${vars} --job_nbr=${np} -d 0 ${cmip6_opt} ${spl_opt} --drc_out=${drc_out} >& $output

#set np='032'
#set output="$outprefix.$stripe.$exp.$np.txt"
#/bin/ls ${drc_in}/*.clm2.h0.25[78][0-9]-*.nc ${drc_in}/*.clm2.h0.259[0-6]-*.nc | ncclimo --var=${vars} --job_nbr=${np} -d 0 ${cmip6_opt} ${spl_opt} --drc_out=${drc_out} >& $output

#set np='136'
#set output="$outprefix.$stripe.$exp.$np.b.txt"
#/bin/ls ${drc_in}/*.clm2.h0.25[78][0-9]-*.nc ${drc_in}/*.clm2.h0.259[0-6]-*.nc | ncclimo --var=${vars} --job_nbr=${np} -d 0 ${cmip6_opt} ${spl_opt} --drc_out=${drc_out} >& $output

grep Elapsed batchout.*

date
cp -r $DW_JOB_STRIPED/outputdir $workdir/ne30/clm/
date
