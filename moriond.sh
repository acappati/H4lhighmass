#!/bin/bash

cd $YourCombine_CMSSW_src_Dir 
eval `scramv1 runtime -sh`
cd $YourCurrentWorkingDir

for dir in $3 
do
#for m in 120 124 126 130 
#do
mkdir -p "$dir/fig" 
root -q -b prepare_skim_mass.c"(\"htxs_stage1_reco_cat\",\"htxs_stage1_reco_catName\",\"$dir\",\"_moriond\",1,0, \"htxs_stage1_red_prod_cat\",\"$1\",$2, \"1000\")"
done

