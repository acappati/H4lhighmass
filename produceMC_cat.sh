## The code used fake rate root files from CJSLT ZZAnalysis framework, if you have the ZZAnalysis code somewhere, just setup cmsenv there. Otherwise, copy the fakerate root file from ZZAnalysis code
## and remove "$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/FakeRates/" in line 69 and 71 of moriond.c file
cd $YourZZanalysisCode_Dir 
cmsenv
cd -

#for i in AllData 
#do
##Zjet 
#	root -q -b moriond_stxs.c"(\"$i\",1,0,2016)"
#	root -q -b moriond_stxs.c"(\"$i\",1,0,2017)"
#	root -q -b moriond_stxs.c"(\"$i\",1,0,2018)"
## Data 
#	root -q -b moriond_stxs.c"(\"$i\",0,0,2018)"
#	root -q -b moriond_stxs.c"(\"$i\",0,0,2017)"
#	root -q -b moriond_stxs.c"(\"$i\",0,0,2016)"
#done

#for i in AllData 
#do
##Zjet 
#	root -q -b moriond_stxs.c"(\"$i\",1,0,2016)"
#	root -q -b moriond_stxs.c"(\"$i\",1,0,2017)"
#	root -q -b moriond.c"(\"$i\",1,0,2018)"
## Data 
#	root -q -b moriond_stxs.c"(\"$i\",0,0,2016)"
#	root -q -b moriond_stxs.c"(\"$i\",0,0,2017)"
#	root -q -b moriond.c"(\"$i\",0,0,2018)"
#done

## Rare backgrounds
#for j in TTZToLL_M1to10_MLM TTZToLLNuNu_M10 TTZJets_M10_MLM TTZZ TTWW ZZZ WWZ WZZ VBFToContinToZZ4l
#do
#        i=$j
#        root -q -b moriond_stxs.c"(\"$i\",0,1,2016)"
#        root -q -b moriond_stxs.c"(\"$i\",0,1,2017)"
#        root -q -b moriond_stxs.c"(\"$i\",0,1,2018)"
#done

## qqZZ and ggZZ backgrounds, with sys
#for j in VBFToContinToZZ4l #ggTo4mu_0MH125Contin_MCFM701 ggTo4e_0MH125Contin_MCFM701 ggTo2e2mu_0MH125Contin_MCFM701 VBFToHiggs0MToZZTo4l_M125_GaSM ZZTo4lext ggTo4mu_Contin_MCFM701 ggTo2e2mu_Contin_MCFM701 ggTo4e_Contin_MCFM701 ggTo4tau_Contin_MCFM701 ggTo2mu2tau_Contin_MCFM701 ggTo2e2tau_Contin_MCFM701 
#do
#	i=$j
#	root -q -b moriond_stxs.c"(\"$i\",0,1,2016)"
#	root -q -b moriond_stxs.c"(\"$i\",0,1,2017)"
#	root -q -b moriond.c"(\"$i\",0,1,2018)"
#        root -q -b skimmer_channel.c"(\"$i\",0,0,2018,0,0)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,0,1)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,0,2)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,1,0)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,1,1)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,1,2)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,2,0)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,2,1)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,2,2)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,3,0)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,3,1)"
#	root -q -b skimmer_channel.c"(\"$i\",0,0,2018,3,2)"
#done

## Mass variated points: used for fityield.c (Normalization parameter)
#for k in 115 120 125 130 135 140 145 150 155 160 165 170 175 190 200 210 230 250 270 350 400 450 500 550 600 700 750 800 900 1000 1500 2000 3000
for k in 250 270 350 400 450 500 550 600 700 750 800 900 1000 1500 2000 2500 3000
#for k in 2000 2500 3000
#for k in 125
do
	for m in ggH VBFH
	do
		i=$m$k
		#root -q -b moriond_stxs.c"(\"$i\",0,1,2016)"
		#root -q -b moriond_stxs.c"(\"$i\",0,1,2017)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,0,0)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,0,1)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,0,2)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,1,0)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,1,1)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,1,2)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,2,0)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,2,1)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,2,2)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,3,0)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,3,1)"
		root -q -b skimmer_channel.c"(\"$i\",0,0,2018,3,2)"
	done
done

#for k in 115 120 125 130 135 140 145 150 155 160 165 170 175 190 200 210 230 250 270 300 350 400 450 500 550 600 700 750 800 900 1000 1500 2000 2500 3000
#for k in 2000 2500 3000
#for k in 125 150 175
#do
#	for m in VBFH 
#	do
#		i=$m$k
		#root -q -b moriond_stxs.c"(\"$i\",0,1,2016)"
		#root -q -b moriond_stxs.c"(\"$i\",0,1,2017)"
#		root -q -b moriond.c"(\"$i\",0,0,2018)"
#	done
#done

## tune up/down points: used for tune.c (Systematics in up/dwn scales)
#for k in tunedown tuneup
#do
#	for m in ggZH ggH ttH WminusH WplusH VBFH ZH bbH tqH 
#	do
#		i=$m"125_"$k
#		root -q -b moriond_stxs.c"(\"$i\",0,1,2016)"
#		root -q -b moriond_stxs.c"(\"$i\",0,1,2017)"
#		root -q -b moriond_stxs.c"(\"$i\",0,1,2018)"
#	done
#done
