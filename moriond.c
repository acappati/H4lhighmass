// It would be nice to use the classes from this module
// TO DO: Find a way to link methods without having to use CJLST code.
// Similarly to stage1.cc
#include "ZZAnalysis/AnalysisStep/src/Discriminants.cc"
#include "ZZAnalysis/AnalysisStep/src/Category.cc"
#include "ZZAnalysis/AnalysisStep/src/cConstants.cc"
#include "src/stage1.cc"
#include "ZZAnalysis/AnalysisStep/interface/FinalStates.h"
#include "ZZAnalysis/AnalysisStep/src/bitops.cc"
#include "ZZAnalysis/AnalysisStep/test/Plotter_v2/src/FakeRates.cpp"
#include "ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Plotter.h"
#include "ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h"
#include "ZZAnalysis/AnalysisStep/src/LeptonSFHelper.cc"

// Could find a way to move this STXS ad interim uncertainties
// out from this code. In principle we also have this info in the
// CJLST code (hence in the TTrees) but atm is wrong (i.e. refers to STXS Stage 1.0)
float g_sig0 = 30.117;
float g_sig1 = 12.928;
float g_sig_ge2 = 5.475;
float g_sig_ge1 = g_sig1+g_sig_ge2;
float g_sig_tot = g_sig0+g_sig_ge1;
float g_sig_vbfTopo = 0.630;
float g_sig_ge2noVBF = g_sig_ge2-g_sig_vbfTopo;
float g_sig_ge1noVBF = g_sig_ge1-g_sig_vbfTopo;
//----
//---- Jet bin uncertainties 
std::vector<float> blptw(int Njets30) {
  
  static std::vector<float> sig({g_sig0,g_sig1,g_sig_ge2noVBF}); // NNLOPS subtracting VBF
  
  // BLPTW absolute uncertainties in pb
  static std::vector<float> yieldUnc({ 1.12, 0.66, 0.42});
  static std::vector<float> resUnc  ({ 0.03, 0.57, 0.42});
  static std::vector<float> cut01Unc({-1.22, 1.00, 0.21});
  static std::vector<float> cut12Unc({    0,-0.86, 0.86});
  
  // account for missing EW+quark mass effects by scaling BLPTW total cross section to sigma(N3LO)
  float sf = 48.52/47.4;
  int jetBin = (Njets30 > 1 ? 2 : Njets30);
  float normFact = sf/sig[jetBin];
  
  return { yieldUnc[jetBin]*normFact, 
    resUnc[jetBin]*normFact,
    cut01Unc[jetBin]*normFact,
    cut12Unc[jetBin]*normFact };
}

float vbf_2j(int STXS) {
  if (STXS==113 || STXS == 115 || (STXS >= 100 && STXS <105)) return 0.200; // 20.0%
  return 0.0; // Events with no VBF topology have no VBF uncertainty
}

float vbf_3j(int STXS) {
  if (STXS==116) return -0.320;//0.38 // GG2H_VBFTOPO_JET3VETO, tot unc 38%
  if (STXS==114) return  0.235;//0.304; // GG2H_VBFTOPO_JET3, tot unc 30.4%
  if ((STXS >= 100 && STXS <105)) return 0.235; 
  return 0.0; // Events with no VBF topology have no VBF uncertainty
}

float interpol(float x, float x1, float y1, float x2, float y2) {
  if (x<x1) return y1;
  if (x>x2) return y2;
  return y1+(y2-y1)*(x-x1)/(x2-x1);
}

// Difference between finite top mass dependence @NLO vs @LO evaluated using Powheg NNLOPS
// taken as uncertainty on the treamtment of top mass in ggF loop
float qm_t(float pT) {
  return interpol(pT,160,0.0,500,0.37);
}

// migration uncertaitny around the 120 GeV boundary
float pT120(float pT, int Njets30) {
  if (Njets30==0) return 0;
  return interpol(pT,90,-0.016,160,0.14);
}
// float pT120(float pT, int Njets30) {
//   if (Njets30==0) return 0;
//   else if (Njets30==1) return interpol(pT, 119, -0.0045, 121, 0.056);
//   else if (Njets30==2) return interpol(pT, 119, -0.012, 121, 0.034);
//   else return interpol(pT, 119, -0.006, 121, 0.02);

//aMC@NLO
/*
  if (Njets30==0) return 0;
  else if (Njets30==1) return interpol(pT, 119, -0.0032, 121, 0.041);
  else if (Njets30==2) return interpol(pT, 119, -0.0078, 121, 0.021);
  else return interpol(pT, 119, -0.032, 121, 0.05);
*/
// }

// migration uncertaitny around the 60 GeV boundary
float pT60(float pT, int Njets30) {
  if (Njets30==0) return 0;
  if (Njets30==1) return interpol(pT,20,-0.1,100,0.1);
  return interpol(pT,0,-0.1,180,0.10); // >=2 jets
}
// float pT60(float pT, int Njets30) {
//   if (Njets30==0) return 0;
//   else if (Njets30==1) return interpol(pT, 59, -0.02, 61, 0.021);//interpol(pT,20,-0.1,100,0.1);
//   else if (Njets30==2) return interpol(pT, 59, -0.034, 61, 0.014);
//   else return interpol(pT, 59, -0.03, 61, 0.01);//interpol(pT,0,-0.1,180,0.10); // >=2 jets

//aMC@NLO
/*
  if (Njets30==0) return 0;
  else if (Njets30==1) return interpol(pT, 59, -0.011, 61, 0.014);
  else if (Njets30==2) return interpol(pT, 59, -0.018, 61, 0.0068);
  else return interpol(pT, 59, -0.01, 61, 0.03);
*/
// }

// migration uncertaitny around the 10 GeV boundary (only for 0 jet bin)
float pT10(float pT, int Njets30) {
  if (Njets30!=0) return 0;
  return interpol(pT, 9, -0.052, 11, 0.167);
  // aMC@NLO
  //return interpol(pT, 9, -0.063, 11, 0.014);
}

std::vector<float> jetBinUnc(int Njets30, int STXS) {
  std::vector<float> result = blptw(Njets30);
  result.push_back(vbf_2j(STXS));
  result.push_back(vbf_3j(STXS));
  // set jet bin uncertainties to zero if we are in the VBF phase-space
  if (result.back()!=0.0) result[0]=result[1]=result[2]=result[3]=0.0;
  return result;
}

std::vector<float> qcd_ggF_uncert_2017_New (int Njets30, float pTH, int STXS) {
  std::vector<float> result = jetBinUnc(Njets30,STXS); // 6 nuisances
  result.push_back(pT10(pTH,Njets30));
  result.push_back(pT60(pTH,Njets30));
  result.push_back(pT120(pTH,Njets30));
  result.push_back(qm_t(pTH));
  return result; // tot 10 nuisances: m_top, mu, r, 01jet, 12jet, 10, 60, 120, VBF2j, VBF3j
}

std::vector<float> unc2sf(const std::vector<float> &unc, float Nsigma) {
  std::vector<float> sfs; 
  for (auto u:unc) sfs.push_back(1.0+Nsigma*u);
  return sfs;
}

std::vector<float> qcd_ggF_uncertSF_2017_New (int Njets30, float pTH, int STXS_Stage1, float Nsigma=1.0) {
	return unc2sf(qcd_ggF_uncert_2017_New(Njets30,pTH,STXS_Stage1),Nsigma);
}

static std::map<int, std::vector<float> > stxs_acc =
  {//STXS   TOT   ,  PTH200,  Mjj60 , Mjj120 , Mjj350 , Mjj700 ,Mjj1000 ,Mjj1500  ,  25       , JET01
   { 200 , {0.07  ,  0     , 0      , 0      , 0      , 0      , 0      , 0       ,        0  ,    0   }},
   { 201 , {0.0744,  0     , 0      , 0      , 0      , 0      , 0      , 0       ,        0  ,-0.1649 }}, // Jet0
   { 202 , {0.3367,  0     , 0      , 0      , 0      , 0      , 0      , 0       ,        0  ,-0.7464 }}, // Jet1
   { 203 , {0.014,   0     , -1     , 0      , 0      , 0      , 0      , 0       ,        0  , 0.0271 }}, // Mjj 0-60
   { 204 , {0.024,   0     , 0.0474 , -1     , 0      , 0      , 0      , 0       ,        0  , 0.0462 }}, // Mjj 60-120
   { 205 , {0.1201,  0     , 0.2379 , 0.1101 , -1     , 0      , 0      , 0       ,        0  , 0.2314 }}, // Mjj 120-350
   { 206 , {0.0391, -1     , 0.0775 , 0.0814 , 0.1083 , 0.1824 , 0.2804 , 0.5573  ,        0  , 0.0754 }}, // Mjj GT350
   { 207 , {0.0375, 0.1166 , 0.0743 , 0.078  , 0.1039 ,-0.2757 , 0      , 0       ,   -0.2306 ,  0.0723 }}, // Mjj 350-700,   PTHjj 0-25    , pTH 0-200
   { 208 , {0.0985, 0.3062 , 0.1951 , 0.2048 , 0.273  ,-0.7243 , 0      , 0       ,   +0.2306 ,  0.1898 }}, // Mjj 350-700,   PTHjj 25-inf  , pTH 0-200
   { 209 , {0.0408, 0.1268 , 0.0807 , 0.0849 , 0.1129 , 0.1903 , -0.0737, -0.0776 ,   -0.2508 , 0.0786 }}, // Mjj GT 700,   PTHjj 0-25  , pTH 0-200
   { 210 , {0.1449, 0.4505 , 0.2871 , 0.3014 , 0.4017 , 0.6763 , -0.0752, -0.089  ,    0.2508 , 0.2794 }} // Mjj GT 700,   PTHjj 25-inf  , pTH 0-200
 };

std::vector<float> uncert_deltas({14.867, 0.394, 9.762, 6.788, 7.276, 3.645, 2.638, 1.005, 20.073, 18.094});

std::map<int, float> powheg_xsec
  {{200,  273.952 },
   {201,  291.030 },
   {202, 1317.635 },
   {203,   54.934 },
   {204,   93.729 },
   {205,  470.017 },
   {206,   153.09 },
   {207,  146.782 },
   {208,  385.566 },
   {209,  159.624 },
   {210,  567.344 }};

float vbf_uncert_stage_1_1(int source, int event_STXS, double Nsigma=1.0){
  if(source < 10){
    float delta_var = stxs_acc[event_STXS][source] * uncert_deltas[source];
    return  1.0 + Nsigma * (delta_var/powheg_xsec[event_STXS]);
  }else{
    return 0.0;
  }
};


void moriond(TString inputname="VBFH1000",int isZX=0, int doSys=1, int year = 2018, bool splitVH=false, bool fail=false){

	LeptonSFHelper *lepSFHelper = new LeptonSFHelper();

	vector <float> _fs_ROS_SS;
	double _lumi;
   if(year == 2016) {
   	   _lumi = 35.9;
	   	_fs_ROS_SS.push_back(1.00229); // 4e
      _fs_ROS_SS.push_back(1.00008); // 4mu
	 	  _fs_ROS_SS.push_back(1.03601); // 2e2mu
	 	  _fs_ROS_SS.push_back(0.99873); // 2mu2e
   }
   else if (year == 2017) {
   	   _lumi = 41.5;
		  _fs_ROS_SS.push_back(1.01168); // 4e
      _fs_ROS_SS.push_back(1.03961); // 4mu
		  _fs_ROS_SS.push_back(1.01307); // 2e2mu
		  _fs_ROS_SS.push_back(1.00266); // 2mu2e
   }
   else {
   	   _lumi = 59.7;
    	_fs_ROS_SS.push_back(1.00635); // 4e
      _fs_ROS_SS.push_back(1.02763); // 4mu
		  _fs_ROS_SS.push_back(1.03170); // 2e2mu
		  _fs_ROS_SS.push_back(1.00492); // 2mu2e
   }

   // 4mu, 4e, 2e2mu
   // Matrix with comb/SS ratio from the 118-130 GeV mass range.
   // We select events in the 105-140 GeV mass range: the comb/SS
   // ratio scales differently, thus we use the reduced one.
  float cb_SS[3][3]={ // Matrix of Comb/SS yield for the three years
          0.9508, 1.1923, 0.968, 
          0.9267, 1.1568, 1.0307, 
          0.9296, 1.1538, 0.9898 
  };

	int rw_year = year;

	TFile *input_file;
	TString fullpath = Form("");

    TString input_dir = "/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All_WithoutGhost";
    if (inputname.Contains("AllData")) {
    		// path of the new Data ntuples, after 20/04 production (cleaning jets bug fixed)
    		input_dir = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200430_LegacyRun2";
            fullpath = Form("%s/Data_%i/%s/ZZ4lAnalysis.root", input_dir.Data(), year, inputname.Data());
    }
    if (!inputname.Contains("AllData"))  {
            if (year==2016) {
                    fullpath = Form("%s/MC_%i_CorrectBTag/%s/ZZ4lAnalysis.root", input_dir.Data(), year, inputname.Data());
            } else {
                    fullpath = Form("%s/MC_%i/%s/ZZ4lAnalysis.root", input_dir.Data(), year, inputname.Data());
            }
            if ((year == 2017) && (inputname.Contains("ZH125_tuneup") || inputname.Contains("ZH125_tunedown"))){
                    rw_year = 2018;
                    cout << " Taking ZH_tune* from 2018" << endl;
                    fullpath = Form("%s/MC_%i/%s/ZZ4lAnalysis.root", input_dir.Data(), 2018, inputname.Data());
            }
            if (year == 2017 && ((inputname.Contains("ggH125_tune") || (inputname.Contains("ggH120"))))) {
                    rw_year = 2018;
                    fullpath = Form("%s/MC_%i/%s/ZZ4lAnalysis.root", input_dir.Data(), 2018, inputname.Data());
            }
            if ((year == 2016) && (inputname.Contains("bbH") && !inputname.Contains("125"))){
                    rw_year = 2018;
                    cout << " Taking 2016 BBH from 2018" << endl;
                    fullpath = Form("%s/MC_%i/%s/ZZ4lAnalysis.root", input_dir.Data(), 2018, inputname.Data());
            }
            if (year==2018 && inputname.Contains("ZZTo4lext"))
                                    fullpath = Form("%s/MC_%i/%s1/ZZ4lAnalysis.root", input_dir.Data(), year, inputname.Data());
            if (year==2018 && (inputname.Contains("TTZJets_M10_MLM") || inputname.Contains("TTZToLLNuNu_M10")))
                                    fullpath = Form("%s/MC_%i/%sext1/ZZ4lAnalysis.root", input_dir.Data(), year, inputname.Data());

	    if (year==2018 && inputname.Contains("VBFToHiggs0MToZZTo4l_M125_GaSM"))
	                            fullpath = Form("%s/MC_%i/%s/ZZ4lAnalysis.root", input_dir.Data(), year, inputname.Data());
	    if (year==2018 && inputname.Contains("0MH125Contin_MCFM701"))
	                            fullpath = Form("%s/MC_%i/%s/ZZ4lAnalysis.root", input_dir.Data(), year, inputname.Data());
    }


	cout<< fullpath<<endl;
	input_file = TFile::Open(fullpath.Data());

	TH1F *hCounters;
	float gen_sum_weights ;
	if (!inputname.Contains("AllData")){
		hCounters= (TH1F*)input_file->Get("ZZTree/Counters");
		gen_sum_weights = hCounters->GetBinContent(40);
	}
	input_file->Close();

	FakeRates *FR;
	TChain *tqqzz;
	if(isZX){
		tqqzz=new TChain ("CRZLLTree/candTree");
		if (year==2016)
			// FR = new FakeRates( "$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/FakeRates/FakeRates_SS_2016_Legacy.root");
			// FR = new FakeRates( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/FRfiles/FakeRates_SS_2016.root");
			FR = new FakeRates( "$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/FakeRates/newData_FakeRates_SS_2016.root");
		else if (year==2017)
			// FR = new FakeRates( "$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/FakeRates/FakeRates_SS_2017_Legacy.root");
			// FR = new FakeRates( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/FRfiles/FakeRates_SS_2017.root");
			FR = new FakeRates( "$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/FakeRates/newData_FakeRates_SS_2017.root");
		else 
			// FR = new FakeRates( "$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/FakeRates/FakeRates_SS_2018_Legacy.root");
			// FR = new FakeRates( "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/FRfiles/FakeRates_SS_2018.root");
			FR = new FakeRates( "$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/FakeRates/newData_FakeRates_SS_2018.root");

	}
	else
		tqqzz= new TChain("ZZTree/candTree");
//		tqqzz= new TChain("ZZTree/candTree_failed");

	tqqzz->Add(fullpath.Data());


	TH1F *temp_zz_4e = new TH1F("vartemp_zz_4e","",30,0.,1.); 
	TH1F *temp_zz_4mu = new TH1F("vartemp_zz_4mu","",30,0.,1.); 
	TH1F *temp_zz_2e2mu = new TH1F("vartemp_zz_2e2mu","",30,0.,1.); 

	float ZZPt,ZZMass,ZZEta,ZZPhi;

	float xsec,KFactorEWKqqZZ,overallEventWeight,KFactorQCDqqZZ_M;
	float ggH_NNLOPS_weight;
	vector<float> *LepPt=new vector<float>;
	vector<float> *LepEta=new vector<float>;
	vector<float> *LepPhi=new vector<float>;
	vector<short> *LepLepId=new vector<short>;

	vector<float> *JetPt=new vector<float>;
	vector<float> *JetPt_JESUp=new vector<float>;
	vector<float> *JetPt_JESDown=new vector<float>;
	vector<float> *JetPt_JERUp=new vector<float>;
	vector<float> *JetPt_JERDown=new vector<float>;
	vector<float> *JetEta=new vector<float>;
	vector<float> *JetPhi=new vector<float>;
	vector<float> *JetMass=new vector<float>;
	vector<float> *JetIsBtagged=new vector<float>;

	vector<float> *ExtraLepPt=new vector<float>;
	vector<float> *ExtraLepEta=new vector<float>;
	vector<float> *ExtraLepPhi=new vector<float>;

	vector<float> *LHEMotherPz=new vector<float>;
	vector<float> *LHEMotherPx=new vector<float>;
	vector<float> *LHEMotherPy=new vector<float>;
	vector<short> *LHEMotherId=new vector<short>;

	vector<float> * LHEDaughterPt = new vector<float>;
	vector<float> * LHEDaughterEta = new vector<float>;
	vector<float> * LHEDaughterPhi = new vector<float>;
	vector<float> * LHEDaughterMass = new vector<float>;

	vector<float> *LHEAssociatedParticlePt=new vector<float>;
	vector<float> *LHEAssociatedParticleEta=new vector<float>;
	vector<float> *LHEAssociatedParticlePhi=new vector<float>;
	vector<short> *LHEAssociatedParticleId=new vector<short>;

	vector<float> *JetQGLikelihood=new vector<float>;

	short Z1Flav,Z2Flav;
	short nCleanedJetsPt30;
	float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
	float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
	float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_qq;
	float p_QQB_BKG_MCFM,p_GG_SIG_ghg2_1_ghz1_1_JHUGen;
	short ZZsel;
	short nExtraLep;
	short nCleanedJetsPt30BTagged_bTagSF;
	float p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal;
	float p_JJQCD_SIG_ghg4_1_JHUGen_JESUp;
	float p_JJQCD_SIG_ghg4_1_JHUGen_JESDn;
	float p_JJQCD_SIG_ghg4_1_JHUGen_JERUp;
	float p_JJQCD_SIG_ghg4_1_JHUGen_JERDn;
	float p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal_qq;
	float phjj_VAJHU_0m_Gen;
	float phjj_VAJHU_Gen;
	float pvbf_VAJHU_Gen;
	float p_m4l_SIG;
	float p_m4l_BKG;
	float GenHPt;

	float phjj_VAJHU_Gen_qq, phjj_VAJHU_0m_Gen_qq;
	int htxs_stage1p2_cat;
	int htxs_stage1_red_cat;
	int htxs_stage1_red_prod_cat;
	int htxs_stage1_reco_cat;
	int htxs_prodMode;
	int genFinalState;

	int CRflag;

	float p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal;
	float p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal; 
	float p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal;
	float p_JJVBF_BKG_MCFM_JECNominal;
	float p_HadZH_BKG_MCFM_JECNominal;
	float p_HadWH_BKG_MCFM_JECNominal;
	float p_JJQCD_BKG_MCFM_JECNominal;

	float PUWeight_Up;
	float PUWeight_Dn;
	float PUWeight;

	float LHEweight_QCDscale_muR2_muF1;
	float LHEweight_QCDscale_muR1_muF1;
	float LHEweight_QCDscale_muR0p5_muF1;
	float LHEweight_QCDscale_muR1_muF2;
	float LHEweight_QCDscale_muR1_muF0p5;
	float LHEweight_AsMZ_Up;
	float LHEweight_AsMZ_Dn;
	float LHEweight_PDFVariation_Up;
	float LHEweight_PDFVariation_Dn;
	float KFactor_EW_qqZZ_unc;
	// vector <float> *qcd_ggF_uncertSF=new vector<float>;
	float GenHEta, GenHPhi, GenHMass;

	float ZZMassCorr = ZZMass - GenHMass;

	float PythiaWeight_isr_muR4; 
	float PythiaWeight_isr_muR0p25;
	float PythiaWeight_fsr_muR4; 
	float PythiaWeight_fsr_muR0p25;

	short nCleanedJetsPt30_jesDn;
	short nCleanedJetsPt30_jerDn;

	short nCleanedJetsPt30BTagged_bTagSF_jesDn;
	short nCleanedJetsPt30BTagged_bTagSF_jerDn;
	short nCleanedJetsPt30BTagged_bTagSFDn;

	float p_JJQCD_SIG_ghg2_1_JHUGen_JESDn;
	float p_JQCD_SIG_ghg2_1_JHUGen_JESDn;
	float p_JJVBF_SIG_ghv1_1_JHUGen_JESDn;
	float p_JVBF_SIG_ghv1_1_JHUGen_JESDn;
	float pAux_JVBF_SIG_ghv1_1_JHUGen_JESDn;
	float p_HadWH_SIG_ghw1_1_JHUGen_JESDn;
	float p_HadZH_SIG_ghz1_1_JHUGen_JESDn;
	float p_HadWH_mavjj_JESDn;
	float p_HadWH_mavjj_true_JESDn;
	float p_HadZH_mavjj_JESDn;
	float p_HadZH_mavjj_true_JESDn;
	float p_JJQCD_SIG_ghg2_1_JHUGen_JERDn;
	float p_JQCD_SIG_ghg2_1_JHUGen_JERDn;
	float p_JJVBF_SIG_ghv1_1_JHUGen_JERDn;
	float p_JVBF_SIG_ghv1_1_JHUGen_JERDn;
	float pAux_JVBF_SIG_ghv1_1_JHUGen_JERDn;
	float p_HadWH_SIG_ghw1_1_JHUGen_JERDn;
	float p_HadZH_SIG_ghz1_1_JHUGen_JERDn;
	float p_HadWH_mavjj_JERDn;
	float p_HadWH_mavjj_true_JERDn;
	float p_HadZH_mavjj_JERDn;
	float p_HadZH_mavjj_true_JERDn;

	float p_Gen_JJEW_BSI_ghv1_1_MCFM;
	float p_Gen_JJEW_SIG_ghv1_1_MCFM;
	float p_Gen_JJEW_BKG_MCFM;
	float p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM;
	float p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM;
	float p_Gen_GG_BKG_MCFM;

	float PFMET_corrected_jesDn;
	float PFMET_corrected_jerDn;

	short nCleanedJetsPt30_jesUp;
	short nCleanedJetsPt30_jerUp;
	short nCleanedJetsPt30BTagged_bTagSF_jesUp;
	short nCleanedJetsPt30BTagged_bTagSF_jerUp;
	short nCleanedJetsPt30BTagged_bTagSFUp;

	float p_JJQCD_SIG_ghg2_1_JHUGen_JESUp;
	float p_JQCD_SIG_ghg2_1_JHUGen_JESUp;
	float p_JJVBF_SIG_ghv1_1_JHUGen_JESUp;
	float p_JVBF_SIG_ghv1_1_JHUGen_JESUp;
	float pAux_JVBF_SIG_ghv1_1_JHUGen_JESUp;
	float p_HadWH_SIG_ghw1_1_JHUGen_JESUp;
	float p_HadZH_SIG_ghz1_1_JHUGen_JESUp;
	float p_HadWH_mavjj_JESUp;
	float p_HadWH_mavjj_true_JESUp;
	float p_HadZH_mavjj_JESUp;
	float p_HadZH_mavjj_true_JESUp;
	float p_JJQCD_SIG_ghg2_1_JHUGen_JERUp;
	float p_JQCD_SIG_ghg2_1_JHUGen_JERUp;
	float p_JJVBF_SIG_ghv1_1_JHUGen_JERUp;
	float p_JVBF_SIG_ghv1_1_JHUGen_JERUp;
	float pAux_JVBF_SIG_ghv1_1_JHUGen_JERUp;
	float p_HadWH_SIG_ghw1_1_JHUGen_JERUp;
	float p_HadZH_SIG_ghz1_1_JHUGen_JERUp;
	float p_HadWH_mavjj_JERUp;
	float p_HadWH_mavjj_true_JERUp;
	float p_HadZH_mavjj_JERUp;
	float p_HadZH_mavjj_true_JERUp;

	float PFMET_corrected_jesUp;
	float PFMET_corrected_jerUp;

	float KFactor_QCD_ggZZ_Nominal;

	vector <float> *JetSigma=new vector <float>;
	int category_baysien;
	float genHEPMCweight;
	Float_t L1prefiringWeight = 0;
	Float_t dataMCWeight = 0;

	tqqzz->SetBranchStatus("*",0);
	tqqzz->SetBranchAddress("L1prefiringWeight", &L1prefiringWeight);

	if (!inputname.Contains("AllData")){
		tqqzz->SetBranchAddress("dataMCWeight", &dataMCWeight);
		// tqqzz->SetBranchAddress("category_baysien",				&category_baysien);
		tqqzz->SetBranchAddress("PUWeight_Up",				&PUWeight_Up);
		tqqzz->SetBranchAddress("PUWeight_Dn",                           & PUWeight_Dn);
		tqqzz->SetBranchAddress("PUWeight",                              & PUWeight);

		tqqzz->SetBranchAddress("LHEweight_QCDscale_muR2_muF1",          & LHEweight_QCDscale_muR2_muF1);
		tqqzz->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF1",        & LHEweight_QCDscale_muR0p5_muF1);
		tqqzz->SetBranchAddress("LHEweight_QCDscale_muR1_muF2",          & LHEweight_QCDscale_muR1_muF2);
		tqqzz->SetBranchAddress("LHEweight_QCDscale_muR1_muF1",          & LHEweight_QCDscale_muR1_muF1);
		tqqzz->SetBranchAddress("LHEweight_QCDscale_muR1_muF0p5",        & LHEweight_QCDscale_muR1_muF0p5);
		tqqzz->SetBranchAddress("LHEweight_AsMZ_Up",                     & LHEweight_AsMZ_Up);
		tqqzz->SetBranchAddress("LHEweight_AsMZ_Dn",                     & LHEweight_AsMZ_Dn);
		tqqzz->SetBranchAddress("LHEweight_PDFVariation_Up",             & LHEweight_PDFVariation_Up);
		tqqzz->SetBranchAddress("LHEweight_PDFVariation_Dn",             & LHEweight_PDFVariation_Dn);
		// tqqzz->SetBranchAddress("qcd_ggF_uncertSF",                      & qcd_ggF_uncertSF);
		tqqzz->SetBranchAddress("PythiaWeight_isr_muR4",                 & PythiaWeight_isr_muR4);
		tqqzz->SetBranchAddress("PythiaWeight_isr_muR0p25",              & PythiaWeight_isr_muR0p25);
		tqqzz->SetBranchAddress("PythiaWeight_fsr_muR4",                 & PythiaWeight_fsr_muR4);
		tqqzz->SetBranchAddress("PythiaWeight_fsr_muR0p25",              & PythiaWeight_fsr_muR0p25);
		tqqzz->SetBranchAddress("GenHPt",&GenHPt);
		tqqzz->SetBranchAddress("LHEMotherPz",&LHEMotherPz);
		// tqqzz->SetBranchAddress("LHEMotherPy",&LHEMotherPy);
		// tqqzz->SetBranchAddress("LHEMotherPx",&LHEMotherPx);


		tqqzz->SetBranchAddress("LHEMotherId",&LHEMotherId);

		tqqzz->SetBranchAddress("LHEDaughterPt",&LHEDaughterPt);
		tqqzz->SetBranchAddress("LHEDaughterEta",&LHEDaughterEta);
		tqqzz->SetBranchAddress("LHEDaughterPhi",&LHEDaughterPhi);
		tqqzz->SetBranchAddress("LHEDaughterMass",&LHEDaughterMass);
		tqqzz->SetBranchAddress("LHEAssociatedParticleEta",&LHEAssociatedParticleEta);
		tqqzz->SetBranchAddress("LHEAssociatedParticlePhi",&LHEAssociatedParticlePhi);
		tqqzz->SetBranchAddress("LHEAssociatedParticleId",&LHEAssociatedParticleId);
		tqqzz->SetBranchAddress("LHEAssociatedParticlePt",&LHEAssociatedParticlePt);
		tqqzz->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
	}
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSFDn",              & nCleanedJetsPt30BTagged_bTagSFDn);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF_jesDn",              & nCleanedJetsPt30BTagged_bTagSF_jesDn);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF_jerDn",              & nCleanedJetsPt30BTagged_bTagSF_jerDn);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSFDn",         &nCleanedJetsPt30BTagged_bTagSFDn);

	//pythia tune from seperate sample;
	tqqzz->SetBranchAddress("nCleanedJetsPt30_jesUp", 			&nCleanedJetsPt30_jesUp);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF_jesUp",         &nCleanedJetsPt30BTagged_bTagSF_jesUp);
	tqqzz->SetBranchAddress("nCleanedJetsPt30_jerUp", 			&nCleanedJetsPt30_jerUp);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF_jerUp",         &nCleanedJetsPt30BTagged_bTagSF_jerUp);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSFUp",         &nCleanedJetsPt30BTagged_bTagSFUp);

	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JESUp",              &p_JJQCD_SIG_ghg2_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JESUp",               &p_JQCD_SIG_ghg2_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JESUp",              &p_JJVBF_SIG_ghv1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JESUp",               &p_JVBF_SIG_ghv1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JESUp",            &pAux_JVBF_SIG_ghv1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JESUp",              &p_HadWH_SIG_ghw1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JESUp",              &p_HadZH_SIG_ghz1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_JESUp",                          &p_HadWH_mavjj_JESUp);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_true_JESUp",                     &p_HadWH_mavjj_true_JESUp);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_JESUp",                          &p_HadZH_mavjj_JESUp);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_true_JESUp",                     &p_HadZH_mavjj_true_JESUp);
	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JERUp",              &p_JJQCD_SIG_ghg2_1_JHUGen_JERUp);
	tqqzz->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JESUp",               &p_JQCD_SIG_ghg2_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JESUp",              &p_JJVBF_SIG_ghv1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JESUp",               &p_JVBF_SIG_ghv1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JESUp",            &pAux_JVBF_SIG_ghv1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JESUp",              &p_HadWH_SIG_ghw1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JESUp",              &p_HadZH_SIG_ghz1_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_JESUp",                          &p_HadWH_mavjj_JESUp);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_true_JESUp",                     &p_HadWH_mavjj_true_JESUp);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_JESUp",                          &p_HadZH_mavjj_JESUp);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_true_JESUp",                     &p_HadZH_mavjj_true_JESUp);

	tqqzz->SetBranchAddress("PFMET_corrected_jesUp",                                  &PFMET_corrected_jesUp);
	tqqzz->SetBranchAddress("PFMET_corrected_jerUp",                                  &PFMET_corrected_jerUp);

	tqqzz->SetBranchAddress("nCleanedJetsPt30_jesDn", 			&nCleanedJetsPt30_jesDn);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF_jesDn",         &nCleanedJetsPt30BTagged_bTagSF_jesDn);
	tqqzz->SetBranchAddress("nCleanedJetsPt30_jerDn", 			&nCleanedJetsPt30_jerDn);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF_jerDn",         &nCleanedJetsPt30BTagged_bTagSF_jerDn);
	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JESDn",              &p_JJQCD_SIG_ghg2_1_JHUGen_JESDn);
	tqqzz->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JESDn",               &p_JQCD_SIG_ghg2_1_JHUGen_JESDn);
	tqqzz->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JESDn",              &p_JJVBF_SIG_ghv1_1_JHUGen_JESDn);
	tqqzz->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JESDn",               &p_JVBF_SIG_ghv1_1_JHUGen_JESDn);
	tqqzz->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JESDn",            &pAux_JVBF_SIG_ghv1_1_JHUGen_JESDn);
	tqqzz->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JESDn",              &p_HadWH_SIG_ghw1_1_JHUGen_JESDn);
	tqqzz->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JESDn",              &p_HadZH_SIG_ghz1_1_JHUGen_JESDn);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_JESDn",                          &p_HadWH_mavjj_JESDn);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_true_JESDn",                     &p_HadWH_mavjj_true_JESDn);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_JESDn",                          &p_HadZH_mavjj_JESDn);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_true_JESDn",                     &p_HadZH_mavjj_true_JESDn);

	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JERDn",              &p_JJQCD_SIG_ghg2_1_JHUGen_JERDn);
	tqqzz->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JERDn",               &p_JQCD_SIG_ghg2_1_JHUGen_JERDn);
	tqqzz->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JERDn",              &p_JJVBF_SIG_ghv1_1_JHUGen_JERDn);
	tqqzz->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JERDn",               &p_JVBF_SIG_ghv1_1_JHUGen_JERDn);
	tqqzz->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JERDn",            &pAux_JVBF_SIG_ghv1_1_JHUGen_JERDn);
	tqqzz->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JERDn",              &p_HadWH_SIG_ghw1_1_JHUGen_JERDn);
	tqqzz->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JERDn",              &p_HadZH_SIG_ghz1_1_JHUGen_JERDn);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_JERDn",                          &p_HadWH_mavjj_JERDn);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_true_JERDn",                     &p_HadWH_mavjj_true_JERDn);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_JERDn",                          &p_HadZH_mavjj_JERDn);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_true_JERDn",                     &p_HadZH_mavjj_true_JERDn);

	tqqzz->SetBranchAddress("PFMET_corrected_jesDn",                                  &PFMET_corrected_jesDn);
	tqqzz->SetBranchAddress("PFMET_corrected_jerDn",                                  &PFMET_corrected_jerDn);

	tqqzz->SetBranchAddress("CRflag",&CRflag);
	tqqzz->SetBranchAddress("nExtraLep",&nExtraLep);
	tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF",&nCleanedJetsPt30BTagged_bTagSF);
	tqqzz->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",&p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",&p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal",&p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal);
	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg4_1_JHUGen_JESUp",&p_JJQCD_SIG_ghg4_1_JHUGen_JESUp);
	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg4_1_JHUGen_JESDn",&p_JJQCD_SIG_ghg4_1_JHUGen_JESDn);
	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg4_1_JHUGen_JERUp",&p_JJQCD_SIG_ghg4_1_JHUGen_JERUp);
	tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg4_1_JHUGen_JERDn",&p_JJQCD_SIG_ghg4_1_JHUGen_JERDn);

	tqqzz->SetBranchAddress("p_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal",&p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_qq);
	tqqzz->SetBranchAddress("p_JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECNominal",&p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal_qq);

	tqqzz->SetBranchAddress("p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal",&p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal);
	tqqzz->SetBranchAddress("p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal",&p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal);
	tqqzz->SetBranchAddress("p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal",&p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal);
	tqqzz->SetBranchAddress("p_JJVBF_BKG_MCFM_JECNominal",&p_JJVBF_BKG_MCFM_JECNominal);
	tqqzz->SetBranchAddress("p_HadZH_BKG_MCFM_JECNominal",&p_HadZH_BKG_MCFM_JECNominal);
	tqqzz->SetBranchAddress("p_HadWH_BKG_MCFM_JECNominal",&p_HadWH_BKG_MCFM_JECNominal);
	tqqzz->SetBranchAddress("p_JJQCD_BKG_MCFM_JECNominal",&p_JJQCD_BKG_MCFM_JECNominal);

	tqqzz->SetBranchAddress("p_Gen_HJJ_SIG_ghg2_1_JHUGen",&phjj_VAJHU_Gen);
	tqqzz->SetBranchAddress("p_Gen_JVBF_SIG_ghv1_1_JHUGen_JECNominal",&pvbf_VAJHU_Gen);
	tqqzz->SetBranchAddress("p_Gen_HJJ_SIG_ghg4_1_JHUGen",&phjj_VAJHU_0m_Gen);
	tqqzz->SetBranchAddress("p_Gen_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal",&phjj_VAJHU_Gen_qq);
	// tqqzz->SetBranchAddress("p_Gen_JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECNominal",&phjj_VAJHU_0m_Gen_qq);

	tqqzz->SetBranchAddress("p_m4l_SIG",&p_m4l_SIG);
	tqqzz->SetBranchAddress("p_m4l_BKG",&p_m4l_BKG);
	tqqzz->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",&p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
	tqqzz->SetBranchAddress("p_QQB_BKG_MCFM",&p_QQB_BKG_MCFM);
	tqqzz->SetBranchAddress("ZZPt",&ZZPt);
	tqqzz->SetBranchAddress("ZZEta",&ZZEta);
	tqqzz->SetBranchAddress("ZZPhi",&ZZPhi);
	tqqzz->SetBranchAddress("ZZMass",&ZZMass);
	tqqzz->SetBranchAddress("Z1Flav",&Z1Flav);
	tqqzz->SetBranchAddress("Z2Flav",&Z2Flav);
	tqqzz->SetBranchAddress("JetIsBtagged",&JetIsBtagged);
	tqqzz->SetBranchAddress("JetPt",&JetPt);
	tqqzz->SetBranchAddress("JetPt_JESUp",&JetPt_JESUp);
	tqqzz->SetBranchAddress("JetPt_JESDown",&JetPt_JESDown);
	tqqzz->SetBranchAddress("JetPt_JERUp",&JetPt_JERUp);
	tqqzz->SetBranchAddress("JetPt_JERDown",&JetPt_JERDown);
	tqqzz->SetBranchAddress("JetSigma",&JetSigma);
	tqqzz->SetBranchAddress("JetMass",&JetMass);
	tqqzz->SetBranchAddress("JetPhi",&JetPhi);
	tqqzz->SetBranchAddress("JetEta",&JetEta);
	tqqzz->SetBranchAddress("ExtraLepPt",&ExtraLepPt);
	tqqzz->SetBranchAddress("ExtraLepPhi",&ExtraLepPhi);
	tqqzz->SetBranchAddress("ExtraLepEta",&ExtraLepEta);
	
	tqqzz->SetBranchAddress("ZZMassCorr",&ZZMassCorr);
	tqqzz->SetBranchAddress("genFinalState",&genFinalState);


	if(inputname.Contains("0MH125Contin_MCFM701")){
	        tqqzz->SetBranchAddress("p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM",&p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM);
		tqqzz->SetBranchAddress("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM",&p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM);
		tqqzz->SetBranchAddress("p_Gen_GG_BKG_MCFM",&p_Gen_GG_BKG_MCFM);
	}
	if(inputname.Contains("VBFToHiggs0MToZZTo4l_M125_GaSM")){
	        tqqzz->SetBranchAddress("p_Gen_JJEW_BSI_ghv1_1_MCFM",&p_Gen_JJEW_BSI_ghv1_1_MCFM);
		tqqzz->SetBranchAddress("p_Gen_JJEW_SIG_ghv1_1_MCFM",&p_Gen_JJEW_SIG_ghv1_1_MCFM);
		tqqzz->SetBranchAddress("p_Gen_JJEW_BKG_MCFM",&p_Gen_JJEW_BKG_MCFM);;
	}

	if(inputname.Contains("ZZTo")){
		tqqzz->SetBranchAddress("KFactor_QCD_qqZZ_M",&KFactorQCDqqZZ_M);
		tqqzz->SetBranchAddress("KFactor_EW_qqZZ",&KFactorEWKqqZZ);	
		tqqzz->SetBranchAddress("KFactor_EW_qqZZ_unc",&KFactor_EW_qqZZ_unc);
	}

	if(inputname.Contains("VBFToContin")){
		tqqzz->SetBranchAddress("KFactor_QCD_qqZZ_Nominal",&KFactor_QCD_ggZZ_Nominal);
	}
	if (inputname.Contains("ggTo"))
		tqqzz->SetBranchAddress("KFactor_QCD_ggZZ_Nominal",&KFactor_QCD_ggZZ_Nominal);
	tqqzz->SetBranchAddress("xsec",&xsec);
	tqqzz->SetBranchAddress("ZZsel",&ZZsel);
	tqqzz->SetBranchAddress("LepLepId",&LepLepId);
	tqqzz->SetBranchAddress("LepPt",&LepPt);
	tqqzz->SetBranchAddress("LepEta",&LepEta);
	tqqzz->SetBranchAddress("LepPhi",&LepPhi);
	tqqzz->SetBranchAddress("overallEventWeight",&overallEventWeight);
	if (inputname.Contains("ggH") )
		tqqzz->SetBranchAddress("ggH_NNLOPS_weight",&ggH_NNLOPS_weight);
	tqqzz->SetBranchAddress("nCleanedJetsPt30",&nCleanedJetsPt30);
	tqqzz->SetBranchAddress("JetQGLikelihood",&JetQGLikelihood);
	tqqzz->SetBranchAddress("htxs_stage1p2_cat",&htxs_stage1p2_cat);
	tqqzz->SetBranchAddress("htxs_prodMode",&htxs_prodMode);
	short genExtInfo;
	tqqzz->SetBranchAddress("genExtInfo",&genExtInfo);
	short RunNumber;
	short EventNumber;
	tqqzz->SetBranchAddress("RunNumber",&RunNumber);
	tqqzz->SetBranchAddress("EventNumber",&EventNumber);

	// for reco category
	short nExtraZ;
	float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
	float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
	float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
	float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
	float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
	float p_HadWH_mavjj_JECNominal;
	float p_HadWH_mavjj_true_JECNominal;
	float p_HadZH_mavjj_JECNominal;
	float p_HadZH_mavjj_true_JECNominal;
	float PFMET_corrected;
	float pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal;
	Float_t                       pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal;
	Float_t                       pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal;
	Float_t                       pConst_JJVBF_BKG_MCFM_JECNominal;
	Float_t                       pConst_HadZH_BKG_MCFM_JECNominal;
	Float_t                       pConst_HadWH_BKG_MCFM_JECNominal;
	Float_t                       pConst_JJQCD_BKG_MCFM_JECNominal;
	tqqzz->SetBranchAddress("nExtraZ",&nExtraZ);
	tqqzz->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal",&p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
	tqqzz->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal",             &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
	tqqzz->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal",          &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
	tqqzz->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal",            &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
	tqqzz->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal",            &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_JECNominal",                        &p_HadWH_mavjj_JECNominal);
	tqqzz->SetBranchAddress("p_HadWH_mavjj_true_JECNominal",                   &p_HadWH_mavjj_true_JECNominal);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_JECNominal",                        &p_HadZH_mavjj_JECNominal);
	tqqzz->SetBranchAddress("p_HadZH_mavjj_true_JECNominal",                   &p_HadZH_mavjj_true_JECNominal);
	tqqzz->SetBranchAddress("PFMET_corrected",                                           &PFMET_corrected);

	tqqzz->SetBranchAddress("pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal",&pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal);
	tqqzz->SetBranchAddress("pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal",&pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal);
	tqqzz->SetBranchAddress("pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal",&pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal);
	tqqzz->SetBranchAddress("pConst_JJVBF_BKG_MCFM_JECNominal",         &pConst_JJVBF_BKG_MCFM_JECNominal);
	tqqzz->SetBranchAddress("pConst_HadZH_BKG_MCFM_JECNominal",         &pConst_HadZH_BKG_MCFM_JECNominal);
	tqqzz->SetBranchAddress("pConst_HadWH_BKG_MCFM_JECNominal",         &pConst_HadWH_BKG_MCFM_JECNominal);
	tqqzz->SetBranchAddress("pConst_JJQCD_BKG_MCFM_JECNominal",         &pConst_JJQCD_BKG_MCFM_JECNominal);
	tqqzz->SetBranchAddress("GenHMass",         &GenHMass);


	float weight;
	vector<float> *weight_up=new vector<float>;
	vector<float> *weight_dn=new vector<float>;
	float weight_vbf, weight_vbf_up, weight_vbf_dn;
	int chan;
	int vbfcate;
	float dbkg_kin;
	float dbkg;
	float d_2j;
	float d_2j_JESup, d_2j_JERup;
	float d_2j_JESdn, d_2j_JERdn;
	float d_2j_qq;
	float d_2j_other;
	float d_2jgen;
	float d_2j_qqgen;
	float dphi;
	float dphigen;
	float mjjgen;
	float mjj;
	// float mjj_up;
	float mjj_jes_up;
	float mjj_jes_dn;
	float mjj_jer_up;
	float mjj_jer_dn;
	float pt_hjj_jes_up;
	float pt_hjj_jes_dn;
	float pt_hjj_jer_up;
	float pt_hjj_jer_dn;
	// float pt_hjj_up;
	// float pt_hjj_dn;
	// float mjj_dn;
	float gene;
	float genpt;
	float vbfMela;
	float vbfMela_gen;
	float vbfMela_qg;
	float qgL;
	int njet;
	int category;
	// int category_stxs;
	vector <float> *dptj=new vector<float>;
	vector <float> *dphij=new vector<float>;
	vector <float> *drj=new vector<float>;
	vector <float> *genjpt=new vector<float>;
	vector <float> *reso=new vector<float>;
	vector <int> *reso_id=new vector<int>;
	vector <float> *opt=new vector<float>;
	float qqfrac;
	int nmatch;
	int ngenjet=0;
	int qg=0;
	if (isZX)
		inputname+="_ZX";
	TFile* fnew = new TFile(Form("/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All_WithoutGhost/ReducedTrees/%s_redTree_%d.root",inputname.Data(),year),"recreate");
	TTree *tnew = new TTree("SelectedTree","SelectedTree");
	tnew->Branch("dbkg_kin",&dbkg_kin,"dbkg_kin/F");
	tnew->Branch("dbkg",&dbkg,"dbkg/F");
	tnew->Branch("d_2j",&d_2j,"d_2j/F");
	tnew->Branch("d_2j_JESup",&d_2j_JESup,"d_2j_JESup/F");
	tnew->Branch("d_2j_JESdn",&d_2j_JESdn,"d_2j_JESdn/F");
	tnew->Branch("d_2j_JERup",&d_2j_JERup,"d_2j_JERup/F");
	tnew->Branch("d_2j_JERdn",&d_2j_JERdn,"d_2j_JERdn/F");
	tnew->Branch("d_2j_qq",&d_2j_qq,"d_2j_qq/F");
	tnew->Branch("d_2j_other",&d_2j_other,"d_2j_other/F");
	tnew->Branch("d_2j_qqgen",&d_2j_qqgen,"d_2j_qqgen/F");
	tnew->Branch("d_2jgen",&d_2jgen,"d_2jgen/F");
	tnew->Branch("dphi",&dphi,"dphi/F");
	tnew->Branch("dphigen",&dphigen,"dphigen/F");
	tnew->Branch("mjjgen",&mjjgen,"mjjgen/F");
	tnew->Branch("mjj",&mjj,"mjj/F");
	tnew->Branch("gene",&gene,"gene/F");
	tnew->Branch("genpt",&genpt,"genpt/F");
	tnew->Branch("ngenjet",&ngenjet,"ngenjet/I");
	tnew->Branch("vbfMela",&vbfMela,"vbfMela/F");
	tnew->Branch("vbfMela_gen",&vbfMela_gen,"vbfMela_gen/F");
	tnew->Branch("vbfMela_qg",&vbfMela_qg,"vbfMela_qg/F");
	tnew->Branch("weight",&weight,"weight/F");
	tnew->Branch("weight_up",&weight_up);
	tnew->Branch("weight_dn",&weight_dn);
	vector <TString> *weight_name=new vector<TString>;
	tnew->Branch("weight_name",&weight_name);
	tnew->Branch("nmatch",&nmatch,"nmatch/I");
	tnew->Branch("drj",&drj);
	tnew->Branch("dptj",&dptj);
	tnew->Branch("dphij",&dphij);
	tnew->Branch("GenHPt",&GenHPt);
	tnew->Branch("GenHEta",&GenHEta);
	tnew->Branch("GenHPhi",&GenHPhi);
	tnew->Branch("GenHMass",&GenHMass);
	tnew->Branch("ZZPt",&ZZPt);
	tnew->Branch("nCleanedJetsPt30BTagged_bTagSF",&nCleanedJetsPt30BTagged_bTagSF);
	tnew->Branch("ZZPhi",&ZZPhi);
	tnew->Branch("ZZEta",&ZZEta);
	tnew->Branch("ZZMass",&ZZMass);
	tnew->Branch("JetPt",&JetPt);
	tnew->Branch("JetPhi",&JetPhi);
	tnew->Branch("JetEta",&JetEta);
	tnew->Branch("JetMass",&JetMass);
	tnew->Branch("LepPt",&LepPt);
	tnew->Branch("LepPhi",&LepPhi);
	tnew->Branch("LepEta",&LepEta);
	tnew->Branch("ExtraLepPt",&ExtraLepPt);
	tnew->Branch("ExtraLepPhi",&ExtraLepPhi);
	tnew->Branch("ExtraLepEta",&ExtraLepEta);
	tnew->Branch("genjpt",&genjpt);
	tnew->Branch("reso",&reso);
	tnew->Branch("reso_id",&reso_id);
	tnew->Branch("opt",&opt);
	tnew->Branch("LHEDaughterMass",&LHEDaughterMass);
	tnew->Branch("LHEAssociatedParticleEta",&LHEAssociatedParticleEta);
	tnew->Branch("LHEAssociatedParticlePhi",&LHEAssociatedParticlePhi);
	tnew->Branch("LHEAssociatedParticlePt",&LHEAssociatedParticlePt);
	tnew->Branch("qgL",&qgL);
	tnew->Branch("chan",&chan,"chan/I");
	tnew->Branch("vbfcate",&vbfcate,"vbfcate/I");
	tnew->Branch("nCleanedJetsPt30",&njet,"nCleanedJetsPt30/I");
	tnew->Branch("nCleanedJetsPt30_jesUp",&nCleanedJetsPt30_jesUp);
	tnew->Branch("nCleanedJetsPt30_jesDn",&nCleanedJetsPt30_jesDn);
	tnew->Branch("nCleanedJetsPt30_jerUp",&nCleanedJetsPt30_jerUp);
	tnew->Branch("nCleanedJetsPt30_jerDn",&nCleanedJetsPt30_jerDn);
	tnew->Branch("nExtraLep",&nExtraLep,"nExtraLep/S");
	tnew->Branch("qg",&qg,"qg/I");
	tnew->Branch("qqfrac",&qqfrac,"qqfrac/F");
	tnew->Branch("category",&category,"category/I");
	char categoryName[200];
	char category_baysien7Name[200];
	tnew->Branch("categoryName",&categoryName,"categoryName/C");
	tnew->Branch("category_baysien7Name",&category_baysien7Name,"category_baysien7Name/C");

	int category_baysien7;
	tnew->Branch("category_baysien7",&category_baysien7,"category_baysien7/I");
	float mcweight;
	tnew->Branch("mcweight",&mcweight);
	tnew->Branch("xsec",&xsec);

	tnew->Branch("ZZMassCorr",&ZZMassCorr);
	tnew->Branch("genFinalState",&genFinalState,"genFinalState/I");

	float D2jet,D1jet,DWH,DZH,DVH,DVBFDEC,DVHDEC;
	float D2jet_JESup, D2jet_JESdn;
	float D2jet_JERup, D2jet_JERdn;
	TString htxs_stage1_catName;
	TString htxs_stage1_red_catName;
	TString htxs_stage1_red_prod_catName;
	TString htxs_stage1_reco_catName;
	TString htxs_stage1p2_catName;

	TString htxs_stage1_reco_catName_jes_up;
	TString htxs_stage1_reco_catName_jes_dn;
	TString htxs_stage1_reco_catName_jer_up;
	TString htxs_stage1_reco_catName_jer_dn;
	// TString htxs_stage1_reco_catName_jetPt_jes_up;
	// TString htxs_stage1_reco_catName_jetPt_jes_dn;
	// TString htxs_stage1_reco_catName_jetPt_jer_up;
	// TString htxs_stage1_reco_catName_jetPt_jer_dn;
	TString htxs_stage1_reco_catName_btag_up;
	TString htxs_stage1_reco_catName_btag_dn;

	int htxs_stage1_reco_cat_jes_up;
	int htxs_stage1_reco_cat_jes_dn;
	int htxs_stage1_reco_cat_jer_up;
	int htxs_stage1_reco_cat_jer_dn;
	// int htxs_stage1_reco_cat_jetPt_jes_up;
	// int htxs_stage1_reco_cat_jetPt_jes_dn;
	// int htxs_stage1_reco_cat_jetPt_jer_up;
	// int htxs_stage1_reco_cat_jetPt_jer_dn;
	int htxs_stage1_reco_cat_btag_up;
	int htxs_stage1_reco_cat_btag_dn;

	TString categoryName_jes_up;
	TString categoryName_jes_dn;
	TString categoryName_jer_up;
	TString categoryName_jer_dn;
	TString categoryName_btag_up;
	TString categoryName_btag_dn;

	int category_jes_up;
	int category_jes_dn;
	int category_jer_up;
	int category_jer_dn;
	// int category_jetPt_jes_up;
	// int category_jetPt_jes_dn;
	// int category_jetPt_jer_up;
	// int category_jetPt_jer_dn;
	int category_btag_up;
	int category_btag_dn;

	tnew->Branch("RunNumber",&RunNumber);
	tnew->Branch("EventNumber",&EventNumber,"EventNumber/I");
	tnew->Branch("D2jet",&D2jet,"D2jet/F");
	tnew->Branch("D2jet_JESdn",&D2jet_JESdn,"D2jet_JESdn/F");
	tnew->Branch("D2jet_JESup",&D2jet_JESup,"D2jet_JESup/F");
	tnew->Branch("D2jet_JERdn",&D2jet_JERdn,"D2jet_JERdn/F");
	tnew->Branch("D2jet_JERup",&D2jet_JERup,"D2jet_JERup/F");
	tnew->Branch("D1jet",&D1jet,"D1jet/F");
	tnew->Branch("DWH",&DWH,"DWH/F");
	tnew->Branch("DZH",&DZH,"DZH/F");
	tnew->Branch("DVH",&DVH,"DVH/F");
	tnew->Branch("DVBFDEC",&DVBFDEC,"DVBFDEC/F");
	tnew->Branch("DVHDEC",&DVHDEC,"DVHDEC/F");
	tnew->Branch("htxs_stage1p2_cat",&htxs_stage1p2_cat,"htxs_stage1p2_cat/I");
	tnew->Branch("htxs_stage1_red_cat",&htxs_stage1_red_cat,"htxs_stage1_red_cat/I");
	tnew->Branch("htxs_stage1_red_prod_cat",&htxs_stage1_red_prod_cat,"htxs_stage1_red_prod_cat/I");
	tnew->Branch("htxs_stage1_reco_cat",&htxs_stage1_reco_cat,"htxs_stage1_reco_cat/I");
	tnew->Branch("htxs_stage1p2_catName",&htxs_stage1p2_catName);

	tnew->Branch("htxs_stage1_reco_cat_jes_up",&htxs_stage1_reco_cat_jes_up,"htxs_stage1_reco_cat_jes_up/I");
	tnew->Branch("htxs_stage1_reco_cat_jes_dn",&htxs_stage1_reco_cat_jes_dn,"htxs_stage1_reco_cat_jes_dn/I");
	tnew->Branch("htxs_stage1_reco_cat_jer_up",&htxs_stage1_reco_cat_jer_up,"htxs_stage1_reco_cat_jer_up/I");
	tnew->Branch("htxs_stage1_reco_cat_jer_dn",&htxs_stage1_reco_cat_jer_dn,"htxs_stage1_reco_cat_jer_dn/I");

	// tnew->Branch("htxs_stage1_reco_cat_jetPt_jes_up",&htxs_stage1_reco_cat_jetPt_jes_up,"htxs_stage1_reco_cat_jetPt_jes_up/I");
	// tnew->Branch("htxs_stage1_reco_cat_jetPt_jes_dn",&htxs_stage1_reco_cat_jetPt_jes_dn,"htxs_stage1_reco_cat_jetPt_jes_dn/I");
	// tnew->Branch("htxs_stage1_reco_cat_jetPt_jer_up",&htxs_stage1_reco_cat_jetPt_jer_up,"htxs_stage1_reco_cat_jetPt_jer_up/I");
	// tnew->Branch("htxs_stage1_reco_cat_jetPt_jer_dn",&htxs_stage1_reco_cat_jetPt_jer_dn,"htxs_stage1_reco_cat_jetPt_jer_dn/I");

	tnew->Branch("htxs_stage1_reco_cat_btag_dn",&htxs_stage1_reco_cat_btag_dn,"htxs_stage1_reco_cat_btag_dn/I");
	tnew->Branch("htxs_stage1_reco_cat_btag_up",&htxs_stage1_reco_cat_btag_up,"htxs_stage1_reco_cat_btag_up/I");

	tnew->Branch("category_jes_up",&category_jes_up,"category_jes_up/I");
	tnew->Branch("category_jes_dn",&category_jes_dn,"category_jes_dn/I");
	tnew->Branch("category_jer_up",&category_jer_up,"category_jer_up/I");
	tnew->Branch("category_jer_dn",&category_jer_dn,"category_jer_dn/I");

	// tnew->Branch("category_jetPt_jes_up",&category,"category_jetPt_jes_up/I");
	// tnew->Branch("category_jetPt_jes_dn",&category,"category_jetPt_jes_dn/I");
	// tnew->Branch("category_jetPt_jer_up",&category,"category_jetPt_jer_up/I");
	// tnew->Branch("category_jetPt_jer_dn",&category,"category_jetPt_jer_dn/I");

	tnew->Branch("category_btag_dn",&category_btag_dn,"category_btag_dn/I");
	tnew->Branch("category_btag_up",&category_btag_up,"category_btag_up/I");

	tnew->Branch("htxs_stage1_catName",&htxs_stage1_catName);
	tnew->Branch("htxs_stage1_red_catName",&htxs_stage1_red_catName);
	tnew->Branch("htxs_stage1_red_prod_catName",&htxs_stage1_red_prod_catName);

	tnew->Branch("htxs_stage1_reco_catName",&htxs_stage1_reco_catName);
	tnew->Branch("htxs_stage1_reco_catName_jes_up",&htxs_stage1_reco_catName_jes_up);
	tnew->Branch("htxs_stage1_reco_catName_jes_dn",&htxs_stage1_reco_catName_jes_dn);
	tnew->Branch("htxs_stage1_reco_catName_jer_up",&htxs_stage1_reco_catName_jer_up);
	tnew->Branch("htxs_stage1_reco_catName_jer_dn",&htxs_stage1_reco_catName_jer_dn);

	//tnew->Branch("categoryHighMass",&categoryHighMass);

	// tnew->Branch("htxs_stage1_reco_catName_jetPt_jes_up",&htxs_stage1_reco_catName_jetPt_jes_up);
	// tnew->Branch("htxs_stage1_reco_catName_jetPt_jes_dn",&htxs_stage1_reco_catName_jetPt_jes_dn);
	// tnew->Branch("htxs_stage1_reco_catName_jetPt_jer_up",&htxs_stage1_reco_catName_jetPt_jer_up);
	// tnew->Branch("htxs_stage1_reco_catName_jetPt_jer_dn",&htxs_stage1_reco_catName_jetPt_jer_dn);

	tnew->Branch("htxs_stage1_reco_catName_btag_dn",&htxs_stage1_reco_catName_btag_dn);
	tnew->Branch("htxs_stage1_reco_catName_btag_up",&htxs_stage1_reco_catName_btag_up);
	tnew->Branch("htxs_prodMode",&htxs_prodMode,"htxs_prodMode/I");

	tnew->Branch("categoryName_jes_up",&categoryName_jes_up);
	tnew->Branch("categoryName_jes_dn",&categoryName_jes_dn);
	tnew->Branch("categoryName_jer_up",&categoryName_jer_up);
	tnew->Branch("categoryName_jer_dn",&categoryName_jer_dn);
	tnew->Branch("categoryName_btag_dn",&categoryName_btag_dn);
	tnew->Branch("categoryName_btag_up",&categoryName_btag_up);

	TString htxs_prodModeName;
	tnew->Branch("htxs_prodModeName",&htxs_prodModeName);
	float JetPt1,JetEta1,JetPhi1,JetMass1,JetIsBtagged1;
	float JetPt2,JetEta2,JetPhi2,JetMass2,JetIsBtagged2;
	float JetPt3,JetEta3,JetPhi3,JetMass3,JetIsBtagged3;
	float ExtraLepPt1,ExtraLepEta1,ExtraLepPhi1;
	float ExtraLepPt2,ExtraLepEta2,ExtraLepPhi2;

	tnew->Branch("JetPt1",&JetPt1,"JetPt1/F");
	tnew->Branch("JetEta1",&JetEta1,"JetEta1/F");
	tnew->Branch("JetPhi1",&JetPhi1,"JetPhi1/F");
	tnew->Branch("JetMass1",&JetMass1,"JetMass1/F");
	tnew->Branch("JetIsBtagged1",&JetIsBtagged1,"JetIsBtagged1/F");

	tnew->Branch("JetPt2",&JetPt2,"JetPt2/F");
	tnew->Branch("JetEta2",&JetEta2,"JetEta2/F");
	tnew->Branch("JetPhi2",&JetPhi2,"JetPhi2/F");
	tnew->Branch("JetMass2",&JetMass2,"JetMass2/F");
	tnew->Branch("JetIsBtagged2",&JetIsBtagged2,"JetIsBtagged2/F");

	tnew->Branch("JetPt3",&JetPt3,"JetPt3/F");
	tnew->Branch("JetEta3",&JetEta3,"JetEta3/F");
	tnew->Branch("JetPhi3",&JetPhi3,"JetPhi3/F");
	tnew->Branch("JetMass3",&JetMass3,"JetMass3/F");
	tnew->Branch("JetIsBtagged3",&JetIsBtagged3,"JetIsBtagged3/F");

	tnew->Branch("ExtraLepPt1",&ExtraLepPt1,"ExtraLepPt1/F");
	tnew->Branch("ExtraLepEta1",&ExtraLepEta1,"ExtraLepEta1/F");
	tnew->Branch("ExtraLepPhi1",&ExtraLepPhi1,"ExtraLepPhi1/F");

	tnew->Branch("ExtraLepPt2",&ExtraLepPt2,"ExtraLepPt2/F");
	tnew->Branch("ExtraLepEta2",&ExtraLepEta2,"ExtraLepEta2/F");
	tnew->Branch("ExtraLepPhi2",&ExtraLepPhi2,"ExtraLepPhi2/F");
	float gen_assoV_m; 
	float gen_assoV_pt; 
	float gen_assoV_eta; 
	float gen_assoV_phi; 
	float pt_hjj;
	tnew->Branch("gen_assoV_m",&gen_assoV_m);
	tnew->Branch("gen_assoV_pt",&gen_assoV_pt);
	tnew->Branch("gen_assoV_eta",&gen_assoV_eta);
	tnew->Branch("gen_assoV_phi",&gen_assoV_phi);
	tnew->Branch("pt_hjj",&pt_hjj);

	cout<< tqqzz->GetEntries()<<endl;

	std::vector <float> qcd_ggF_uncertSF;
	tnew->Branch("qcd_ggF_uncertSF",&qcd_ggF_uncertSF);

    std::vector <float> qcd_qqH_uncertSF;
    tnew->Branch("qcd_qqH_uncertSF", &qcd_qqH_uncertSF);

	for(int i=0;i<tqqzz->GetEntries();i++){
		//	for(int i=0;i<2000;i++){
		tqqzz->GetEntry(i);

		if (i%1000==0)
			cout<<i<<endl;
		if (ZZsel<0)
			continue;
		// if((ZZMass<105) || (ZZMass>140)) //(ZZMass<70);
		if(ZZMass<70)

			continue;
		njet = nCleanedJetsPt30;
		ZZMassCorr = ZZMass - GenHMass;

		// STXS theory uncertainties for ggH
		if(inputname.Contains("ggH")) {
			std::vector<float> qcd_ggF_uncertSF_tmp;
			qcd_ggF_uncertSF.clear();

			qcd_ggF_uncertSF_tmp = qcd_ggF_uncertSF_2017_New(njet, ZZPt, htxs_stage1p2_cat);
			qcd_ggF_uncertSF = std::vector<float>(qcd_ggF_uncertSF_tmp.begin(),qcd_ggF_uncertSF_tmp.end());
		}
		
		// STXS uncertainties for VBF
        if(inputname.Contains("VBF") && !inputname.Contains("VBFTo")) { 
                std::vector<float> qcd_qqH_uncertSF_tmp;
                qcd_qqH_uncertSF.clear();

                for(int s = 0; s < 10; s++) {
                        qcd_qqH_uncertSF_tmp.push_back(vbf_uncert_stage_1_1(s, htxs_stage1p2_cat));
                }

                qcd_qqH_uncertSF = std::vector<float>(qcd_qqH_uncertSF_tmp.begin(), qcd_qqH_uncertSF_tmp.end());
        }
		if(abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==121){
			chan=2;
		}
		else if (abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)!=121){
			chan=1;
		}
		else{
			chan=3;
		}
		if(xsec<=0)	
			xsec=1;
		if(isZX)
		{
			if ( !CRflag ) continue;
			if ( !test_bit(CRflag, CRZLLss) ) continue;

			int _current_final_state = -999; 
			if ( Z1Flav == -121 )
			{
				if ( Z2Flav == +121 )
					_current_final_state = 0; 
				else if ( Z2Flav == +169 )
					_current_final_state= 2;  
				else
					cerr << "[ERROR] in event "<<endl; 
			}
			else if ( Z1Flav == -169 )
			{
				if ( Z2Flav == +121 )
					_current_final_state= 3; 
				else if ( Z2Flav == +169 )
					_current_final_state= 1; 
				else
					cerr << "[ERROR] in event "<<endl; 
			}
			weight= cb_SS[year-2016][chan-1]*_fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2))*FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
		}
		else if (!inputname.Contains("AllData")){
			weight= L1prefiringWeight * xsec* overallEventWeight/gen_sum_weights*_lumi*1000;
			mcweight = genHEPMCweight;	
		}
		if (inputname.Contains("ggH") && (year !=2016 || (year==2016 && inputname=="ggH125"))){
			weight *= ggH_NNLOPS_weight;
			mcweight*= ggH_NNLOPS_weight;
		}

		if (inputname.Contains("VBFToContin")) {
			// VBF off-shell BKG has wrong units in xsec, hence the factor 1e3 is not needed
			weight= L1prefiringWeight * xsec* overallEventWeight/gen_sum_weights*_lumi;
		}

		if(inputname.Contains("ZZTo4l")){
			weight= L1prefiringWeight * xsec*KFactorEWKqqZZ * overallEventWeight*KFactorQCDqqZZ_M/gen_sum_weights*_lumi*1000;
			//			double rho = ZZPt/(LepPt->at(0)+LepPt->at(1)+LepPt->at(2)+LepPt->at(3)); 
			//			if(rho<0.3){
			//				weight_dn = weight*(1-abs((KFactorQCDqqZZ_M-1)*(KFactorEWKqqZZ-1))); 
			//				weight_up = weight*(1+abs((KFactorQCDqqZZ_M-1)*(KFactorEWKqqZZ-1))); 
			//			}
			//			else{
			//				weight_dn = weight*(1-abs((KFactorEWKqqZZ-1))); 
			//				weight_up = weight*(1+abs((KFactorEWKqqZZ-1))); 
			//			}
		}
		if (inputname.Contains("VBFToHiggs0MToZZTo4l_M125_GaSM")) {
			// VBF off-shell BKG has wrong units in xsec, hence the factor 1e3 is not needed
		  //weight= L1prefiringWeight * xsec * overallEventWeight * KFactorEWKqqZZ * KFactorQCDqqZZ_M * (p_Gen_JJEW_BSI_ghv1_1_MCFM - p_Gen_JJEW_SIG_ghv1_1_MCFM - p_Gen_JJEW_BKG_MCFM)/gen_sum_weights*_lumi;
		  //weight= L1prefiringWeight * xsec * overallEventWeight * (p_Gen_JJEW_BSI_ghv1_1_MCFM - p_Gen_JJEW_SIG_ghv1_1_MCFM - p_Gen_JJEW_BKG_MCFM)/gen_sum_weights*_lumi;
		  //weight= L1prefiringWeight * xsec * overallEventWeight * p_Gen_JJEW_SIG_ghv1_1_MCFM /gen_sum_weights*_lumi;
		  weight= L1prefiringWeight * xsec * overallEventWeight * p_Gen_JJEW_BKG_MCFM /gen_sum_weights*_lumi;

				  
		}
		

		if (inputname.Contains("_Contin_MCFM701")){
			weight = L1prefiringWeight * KFactor_QCD_ggZZ_Nominal*xsec*overallEventWeight/gen_sum_weights*_lumi*1000;
		}

		if (inputname.Contains("0MH125Contin_MCFM701")){
		  //weight = L1prefiringWeight * KFactor_QCD_ggZZ_Nominal*xsec*overallEventWeight*(p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM - p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM - p_Gen_GG_BKG_MCFM)/gen_sum_weights*_lumi*1000;
		  //weight = L1prefiringWeight * KFactor_QCD_ggZZ_Nominal*xsec*overallEventWeight* p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM /gen_sum_weights*_lumi*1000;
		  weight = L1prefiringWeight * KFactor_QCD_ggZZ_Nominal*xsec*overallEventWeight* p_Gen_GG_BKG_MCFM /gen_sum_weights*_lumi*1000;


		}

		if (!inputname.Contains("AllData")){
            float _updatedSF = ( lepSFHelper->getSF(year,LepLepId->at(0),LepPt->at(0),LepEta->at(0), LepEta->at(0), false) *
                   lepSFHelper->getSF(year,LepLepId->at(1),LepPt->at(1),LepEta->at(1), LepEta->at(1), false) *
                   lepSFHelper->getSF(year,LepLepId->at(2),LepPt->at(2),LepEta->at(2), LepEta->at(2), false) *
                   lepSFHelper->getSF(year,LepLepId->at(3),LepPt->at(3),LepEta->at(3), LepEta->at(3), false) );
            if(rw_year != year) {
                //cout << "Using SF from " << rw_year << endl;
	                _updatedSF = ( lepSFHelper->getSF(rw_year,LepLepId->at(0),LepPt->at(0),LepEta->at(0), LepEta->at(0), false) *
	                 lepSFHelper->getSF(rw_year,LepLepId->at(1),LepPt->at(1),LepEta->at(1), LepEta->at(1), false) *
	                lepSFHelper->getSF(rw_year,LepLepId->at(2),LepPt->at(2),LepEta->at(2), LepEta->at(2), false) *
	                lepSFHelper->getSF(rw_year,LepLepId->at(3),LepPt->at(3),LepEta->at(3), LepEta->at(3), false) );
            }
            weight *= _updatedSF/dataMCWeight;
		}
		drj->clear();
		genjpt->clear();
		dptj->clear();
		dphij->clear();
		mjj=-1;
		nmatch=0;
		qgL=-1;
		ngenjet=LHEAssociatedParticlePt->size();
		reso ->clear(); 
		opt->clear(); 
		reso_id->clear();

		//	TLorentzVector genH;
		//	GenHMass= genH.M();
		//	GenHPhi= genH.Phi();
		//	GenHEta= genH.Eta();
		qg=0;

		short ZZFlav = Z1Flav*Z2Flav;

		float jetPhi[nCleanedJetsPt30];
		float jetQGL[nCleanedJetsPt30];

		if (fail)
			goto linered;
		JetPt1=-999;JetEta1=-999;JetPhi1=-999;JetMass1=-999;JetIsBtagged1=-999;
		JetPt2=-999;JetEta2=-999;JetPhi2=-999;JetMass2=-999;JetIsBtagged2=-999;
		JetPt3=-999;JetEta3=-999;JetPhi3=-999;JetMass3=-999;JetIsBtagged3=-999;
		ExtraLepPt1=-999;ExtraLepEta1=-999;ExtraLepPhi1=-999;
		ExtraLepPt2=-999;ExtraLepEta2=-999;ExtraLepPhi2=-999;

		if(nCleanedJetsPt30>2){
			JetPt3= JetPt->at(2); 
			JetEta3= JetEta->at(2); 
			JetPhi3= JetPhi->at(2); 
			JetMass3= JetMass->at(2); 
			JetIsBtagged3= JetIsBtagged->at(2); 
		}
		if (nCleanedJetsPt30>1){
			JetPt2= JetPt->at(1); 
			JetEta2= JetEta->at(1); 
			JetPhi2= JetPhi->at(1); 
			JetMass2= JetMass->at(1); 
			JetIsBtagged2= JetIsBtagged->at(1); 
		}
		if (nCleanedJetsPt30>0){
			JetPt1= JetPt->at(0); 
			JetEta1= JetEta->at(0); 
			JetPhi1= JetPhi->at(0); 
			JetMass1= JetMass->at(0); 
			JetIsBtagged1= JetIsBtagged->at(0); 
		}
		if (nExtraLep>1){
			ExtraLepPt2 = ExtraLepPt->at(1); 
			ExtraLepEta2 = ExtraLepEta->at(1); 
			ExtraLepPhi2 = ExtraLepPhi->at(1); 
		}
		if (nExtraLep>0){
			ExtraLepPt1 = ExtraLepPt->at(0); 
			ExtraLepEta1 = ExtraLepEta->at(0); 
			ExtraLepPhi1 = ExtraLepPhi->at(0); 
		}


		for ( int j = 0; j < nCleanedJetsPt30; j++)
		{
			jetPhi[j] = JetPhi->at(j);
			if(JetQGLikelihood->size()!=0)
				jetQGL[j] = JetQGLikelihood->at(j);
			else
				jetQGL[j]=-1;
		}
		category= categoryMor18(nExtraLep,
				nExtraZ,
				nCleanedJetsPt30,
				nCleanedJetsPt30BTagged_bTagSF,
				jetQGL,
				p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
				p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
				p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
				p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
				pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
				p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
				p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
				p_HadWH_mavjj_JECNominal,
				p_HadWH_mavjj_true_JECNominal,
				p_HadZH_mavjj_JECNominal,
				p_HadZH_mavjj_true_JECNominal,
				jetPhi,
				ZZMass,
				PFMET_corrected,
				false,// Use VHMET category
				false);// Use QG tagging

		dbkg_kin = p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen + p_QQB_BKG_MCFM*getDbkgkinConstant(ZZFlav, ZZMass));
		dbkg = p_GG_SIG_ghg2_1_ghz1_1_JHUGen*p_m4l_SIG/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen*p_m4l_SIG + p_m4l_BKG*p_QQB_BKG_MCFM*getDbkgkinConstant(ZZFlav, ZZMass));
		if(chan==1)
			temp_zz_4mu->Fill(dbkg_kin,weight);
		else if(chan==2)
			temp_zz_4e->Fill(dbkg_kin,weight);
		else
			temp_zz_2e2mu->Fill(dbkg_kin,weight);



		DVBFDEC = (nCleanedJetsPt30>=2) ? D_bkg_VBFdec( p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
				p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
				p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
				p_JJVBF_BKG_MCFM_JECNominal,
				p_HadZH_BKG_MCFM_JECNominal,
				p_HadWH_BKG_MCFM_JECNominal,
				p_JJQCD_BKG_MCFM_JECNominal,
				p_HadZH_mavjj_JECNominal,
				p_HadZH_mavjj_true_JECNominal,
				p_HadWH_mavjj_JECNominal,
				p_HadWH_mavjj_true_JECNominal,
				pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
				pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
				pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
				pConst_JJVBF_BKG_MCFM_JECNominal,
				pConst_HadZH_BKG_MCFM_JECNominal,
				pConst_HadWH_BKG_MCFM_JECNominal,
				pConst_JJQCD_BKG_MCFM_JECNominal,
				Z1Flav*Z2Flav,
				ZZMass) : -2;

		DVHDEC = (nCleanedJetsPt30>=2) ?   D_bkg_VHdec( p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
				p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
				p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
				p_JJVBF_BKG_MCFM_JECNominal,
				p_HadZH_BKG_MCFM_JECNominal,
				p_HadWH_BKG_MCFM_JECNominal,
				p_JJQCD_BKG_MCFM_JECNominal,
				p_HadZH_mavjj_JECNominal,
				p_HadZH_mavjj_true_JECNominal,
				p_HadWH_mavjj_JECNominal,
				p_HadWH_mavjj_true_JECNominal,
				pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
				pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
				pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
				pConst_JJVBF_BKG_MCFM_JECNominal,
				pConst_HadZH_BKG_MCFM_JECNominal,
				pConst_HadWH_BKG_MCFM_JECNominal,
				pConst_JJQCD_BKG_MCFM_JECNominal,
				Z1Flav*Z2Flav,
				ZZMass) : -2;

		D2jet = ( nCleanedJetsPt30 >= 2)  ? DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass) : -2;
		D2jet_JESup = ( nCleanedJetsPt30_jesUp >= 2)  ? DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JESUp, p_JJQCD_SIG_ghg2_1_JHUGen_JESUp, ZZMass) : -2;
		D2jet_JESdn = ( nCleanedJetsPt30_jesDn >= 2)  ? DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JESDn, p_JJQCD_SIG_ghg2_1_JHUGen_JESDn, ZZMass) : -2;
		D2jet_JERup = ( nCleanedJetsPt30_jerUp >= 2)  ? DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JERUp, p_JJQCD_SIG_ghg2_1_JHUGen_JERUp, ZZMass) : -2;
		D2jet_JERdn = ( nCleanedJetsPt30_jerDn >= 2)  ? DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JERDn, p_JJQCD_SIG_ghg2_1_JHUGen_JERDn, ZZMass) : -2;
		D1jet = ( nCleanedJetsPt30 == 1 ) ? DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass) : -2;
		DWH =   ( nCleanedJetsPt30 >= 2 ) ? DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, p_HadWH_mavjj_true_JECNominal, ZZMass) : -2;
		DZH =   ( nCleanedJetsPt30 >= 2 ) ? DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, p_HadWH_mavjj_true_JECNominal, ZZMass) : -2;
linered:
		float oldCConstD2jet = getDVBF2jetsConstant(ZZMass);
		float oldCConstD1jet = getDVBF1jetConstant(ZZMass);
		float oldCConstDWH = getDWHhConstant(ZZMass);
		float oldCConstDZH = getDZHhConstant(ZZMass);
		float newCConstD2jet = getDVBF2jetsConstant_shiftWP(ZZMass,false,NEWWP2J);
		float newCConstD1jet = getDVBF1jetConstant_shiftWP(ZZMass,false,NEWWP1J);
		float newCConstDWH = getDWHhConstant_shiftWP(ZZMass,false,NEWWPWH);
		float newCConstDZH = getDZHhConstant_shiftWP(ZZMass,false,NEWWPZH);
		D2jet = 1/(newCConstD2jet/oldCConstD2jet*(1/D2jet-1)+1);
		D2jet_JESup = 1/(newCConstD2jet/oldCConstD2jet*(1/D2jet_JESup-1)+1);
		D2jet_JESdn = 1/(newCConstD2jet/oldCConstD2jet*(1/D2jet_JESdn-1)+1);
		D2jet_JERup = 1/(newCConstD2jet/oldCConstD2jet*(1/D2jet_JERup-1)+1);
		D2jet_JERdn = 1/(newCConstD2jet/oldCConstD2jet*(1/D2jet_JERdn-1)+1);
		D1jet = 1/(newCConstD1jet/oldCConstD1jet*(1/D1jet-1)+1);
		DWH = 1/(newCConstDWH/oldCConstDWH*(1/DWH-1)+1);
		DZH = 1/(newCConstDZH/oldCConstDZH*(1/DZH-1)+1);
		float DVH = max(DWH,DZH);
		//cout<< oldCConstD1jet<<"\t"<< newCConstD1jet<<"\t"<<newCConstD1jet/oldCConstD1jet<<endl;

		float WP_VBF2j = getDVBF2jetsWP(ZZMass, 0);
		d_2j= p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/(p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal+p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal); 
		d_2j_JESup= p_JJQCD_SIG_ghg2_1_JHUGen_JESUp/(p_JJQCD_SIG_ghg2_1_JHUGen_JESUp+p_JJQCD_SIG_ghg4_1_JHUGen_JESUp); 
		d_2j_JESdn= p_JJQCD_SIG_ghg2_1_JHUGen_JESDn/(p_JJQCD_SIG_ghg2_1_JHUGen_JESDn+p_JJQCD_SIG_ghg4_1_JHUGen_JESDn); 
		d_2j_JERup= p_JJQCD_SIG_ghg2_1_JHUGen_JERUp/(p_JJQCD_SIG_ghg2_1_JHUGen_JERUp+p_JJQCD_SIG_ghg4_1_JHUGen_JERUp); 
		d_2j_JERdn= p_JJQCD_SIG_ghg2_1_JHUGen_JERDn/(p_JJQCD_SIG_ghg2_1_JHUGen_JERDn+p_JJQCD_SIG_ghg4_1_JHUGen_JERDn); 

		d_2j_qq= p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_qq/(p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_qq+p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal_qq); 
		d_2j_other=  (p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal-p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_qq)/( p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal-p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_qq+p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal-p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal_qq); 
		qqfrac= p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_qq/(p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal); 
		d_2jgen= phjj_VAJHU_Gen/(phjj_VAJHU_Gen+phjj_VAJHU_0m_Gen); 
		d_2j_qqgen= phjj_VAJHU_Gen_qq/(phjj_VAJHU_0m_Gen_qq+phjj_VAJHU_Gen_qq); 

		float c_Mela2j = getDVBF2jetsConstant(ZZMass);
		vbfMela= 1./(1.+ c_Mela2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);

		vbfMela_gen= 1./(1.+ 0.05*c_Mela2j*phjj_VAJHU_Gen/pvbf_VAJHU_Gen);
		if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && vbfMela>WP_VBF2j )
			vbfcate=1;
		else
			vbfcate=0;
		TLorentzVector h_tlz ;
		if (fail)
			goto linest;
		h_tlz.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi,ZZMass);
		pt_hjj=-1;
		if(nCleanedJetsPt30>=2){
			//		if(JetPt->size()>=2){
			TLorentzVector j1,j2;
			j1.SetPtEtaPhiM(JetPt->at(0),JetEta->at(0),JetPhi->at(0),JetMass->at(0));
			j2.SetPtEtaPhiM(JetPt->at(1),JetEta->at(1),JetPhi->at(1),JetMass->at(1));
			if(j1.Pz()>j2.Pz())
				dphi = j1.Phi()-j2.Phi();
			else
				dphi = j2.Phi()-j1.Phi();
			if (dphi>TMath::Pi())
				dphi -= 2*TMath::Pi();
			if (dphi< -TMath::Pi())
				dphi += 2*TMath::Pi();
			mjj= (j1+j2).M();
			pt_hjj = (j1+j2+h_tlz).Pt();
			if(doSys){
				// TLorentzVector j1_up,j2_up;
				// TLorentzVector j1_dn,j2_dn;
				// j1_up.SetPtEtaPhiM(JetPt->at(0)*(1+JetSigma->at(0)),JetEta->at(0),JetPhi->at(0),JetMass->at(0));
				// j2_up.SetPtEtaPhiM(JetPt->at(1)*(1+JetSigma->at(1)),JetEta->at(1),JetPhi->at(1),JetMass->at(1));
				// j1_dn.SetPtEtaPhiM(JetPt->at(0)*(1-JetSigma->at(0)),JetEta->at(0),JetPhi->at(0),JetMass->at(0));
				// j2_dn.SetPtEtaPhiM(JetPt->at(1)*(1-JetSigma->at(1)),JetEta->at(1),JetPhi->at(1),JetMass->at(1));
				// mjj_up = (j1_up+j2_up).M();
				// mjj_dn = (j1_dn+j2_dn).M();
				// pt_hjj_up = (j1_up+j2_up+h_tlz).Pt();
				// pt_hjj_dn = (j1_dn+j2_dn+h_tlz).Pt();

				// Now in the Trees we already have JetPt->at(0)*(1+JetSigma->at(0)) stored in JetPt_JES(R)Up(Down)
				TLorentzVector j1pT_jes_up, j2pT_jes_up, j1pT_jer_up, j2pT_jer_up;
				TLorentzVector j1pT_jes_dn, j2pT_jes_dn, j1pT_jer_dn, j2pT_jer_dn;
				j1pT_jes_up.SetPtEtaPhiM(JetPt_JESUp->at(0),JetEta->at(0),JetPhi->at(0),JetMass->at(0));
				j2pT_jes_up.SetPtEtaPhiM(JetPt_JESUp->at(1),JetEta->at(1),JetPhi->at(1),JetMass->at(1));
				j1pT_jer_up.SetPtEtaPhiM(JetPt_JERUp->at(0),JetEta->at(0),JetPhi->at(0),JetMass->at(0));
				j2pT_jer_up.SetPtEtaPhiM(JetPt_JERUp->at(1),JetEta->at(1),JetPhi->at(1),JetMass->at(1));

				j1pT_jes_dn.SetPtEtaPhiM(JetPt_JESDown->at(0),JetEta->at(0),JetPhi->at(0),JetMass->at(0));
				j2pT_jes_dn.SetPtEtaPhiM(JetPt_JESDown->at(1),JetEta->at(1),JetPhi->at(1),JetMass->at(1));
				j1pT_jer_dn.SetPtEtaPhiM(JetPt_JERDown->at(0),JetEta->at(0),JetPhi->at(0),JetMass->at(0));
				j2pT_jer_dn.SetPtEtaPhiM(JetPt_JERDown->at(1),JetEta->at(1),JetPhi->at(1),JetMass->at(1));

				mjj_jes_up = (j1pT_jes_up+j2pT_jes_up).M();
				mjj_jes_dn = (j1pT_jes_dn+j1pT_jes_dn).M();
				mjj_jer_up = (j1pT_jer_up+j2pT_jer_up).M();
				mjj_jer_dn = (j1pT_jer_dn+j1pT_jer_dn).M();

				pt_hjj_jes_up = (j1pT_jes_up+j2pT_jes_up+h_tlz).Pt();
				pt_hjj_jes_dn = (j1pT_jes_dn+j2pT_jes_dn+h_tlz).Pt();
				pt_hjj_jer_up = (j1pT_jer_up+j2pT_jer_up+h_tlz).Pt();
				pt_hjj_jer_dn = (j1pT_jer_dn+j2pT_jer_dn+h_tlz).Pt();

			}

			//			qgL = JetQGLikelihood->at(0)*JetQGLikelihood->at(1);
			//			gOverq = (1./JetQGLikelihood->at(0)-1)*(1/JetQGLikelihood->at(1)-1);
			//			if(gOverq<0)
			//				d_2j_qg=-1;
			//			d_2j_qg= p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal*gOverq/(p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal*gOverq+p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal); 
		}
			else {
				mjj = 0;
				// mjj_up = 0;
				// mjj_dn = 0;
				// pt_hjj_up = 0;
				// pt_hjj_dn = 0;
				mjj_jes_up = 0;
				mjj_jes_dn = 0;
				mjj_jer_up = 0;
				mjj_jer_dn = 0;
				pt_hjj_jes_up = 0;
				pt_hjj_jes_dn = 0;
				pt_hjj_jer_up = 0;
				pt_hjj_jer_dn = 0;
			}

linest:
			switch(htxs_stage1p2_cat){
					  // ggH
				case 100: htxs_stage1_catName = "ggH_FWDH"; break; 
				case 101: htxs_stage1_catName = "ggH_200_300"; break;
				case 102: htxs_stage1_catName = "ggH_300_450"; break;
				case 103: htxs_stage1_catName = "ggH_450_650"; break;
				case 104: htxs_stage1_catName = "ggH_GT650"; break;				
				case 105: htxs_stage1_catName = "ggH_0j_0_10"; break;
				case 106: htxs_stage1_catName = "ggH_0j_10_200"; break;
				case 107: htxs_stage1_catName = "ggH_1j_0_60"; break;
				case 108: htxs_stage1_catName = "ggH_1j_60_120"; break;
				case 109: htxs_stage1_catName = "ggH_1j_120_200"; break;
				case 110: htxs_stage1_catName = "ggH_2j_0_60"; break;
				case 111: htxs_stage1_catName = "ggH_2j_60_120"; break;
				case 112: htxs_stage1_catName = "ggH_2j_120_200"; break;
				case 113: htxs_stage1_catName = "ggH_VBF_lt700_2j"; break;
				case 114: htxs_stage1_catName = "ggH_VBF_lt700_3j"; break;
				case 115: htxs_stage1_catName = "ggH_VBF_gt700_2j"; break;
				case 116: htxs_stage1_catName = "ggH_VBF_gt700_3j"; break;
					  // VBF
				case 200: htxs_stage1_catName = "VBF_FWDH"; break;
				case 201: htxs_stage1_catName = "VBF_0j"; break;
				case 202: htxs_stage1_catName = "VBF_1j"; break;
				case 203: htxs_stage1_catName = "VBF_2j_mjj_0_60"; break;
				case 204: htxs_stage1_catName = "VH_Had"; break;
 				case 205: htxs_stage1_catName = "VBF_2j_mjj_120_350"; break;
				case 206: htxs_stage1_catName = "VBF_GT200"; break; // QQ2HQQ_GE2J_MJJ_GT350_PTH_GT200
				case 207: htxs_stage1_catName = "VBF_2j_mjj_350_700_2j"; break;
				case 208: htxs_stage1_catName = "VBF_2j_mjj_350_700_3j"; break;
				case 209: htxs_stage1_catName = "VBF_2j_mjj_GT700_2j"; break;
				case 210: htxs_stage1_catName = "VBF_2j_mjj_GT700_3j"; break;
					  //VH
				case 300:htxs_stage1_catName = "WH_FWDH" ;break;
				case 301:htxs_stage1_catName = "WH_0_75" ;break;
				case 302:htxs_stage1_catName = "WH_75_150" ;break;
				case 303:htxs_stage1_catName = "WH_150_250_0J" ;break;
				case 304:htxs_stage1_catName = "WH_150_250_GE1J"; break;
				case 305:htxs_stage1_catName = "WH_GT250" ;break;
				case 400:htxs_stage1_catName = "ZH_FWDH" ;break;
				case 401:htxs_stage1_catName = "ZH_0_75" ;break;
				case 402:htxs_stage1_catName = "ZH_75_150" ;break;
				case 403:htxs_stage1_catName = "ZH_150_250_0J" ;break;
				case 404:htxs_stage1_catName = "ZH_150_250_GE1J" ;break;
				case 405:htxs_stage1_catName = "ZH_GT250" ;break;
				case 500:htxs_stage1_catName = "ggZH_FWDH" ;break;
				case 501:htxs_stage1_catName = "ggZH_0_75" ;break;
				case 502:htxs_stage1_catName = "ggZH_75_150" ;break;
				case 503:htxs_stage1_catName = "ggZH_GT150_0J" ;break;
				case 504:htxs_stage1_catName = "ggZH_GT150_GE1J"; break;
				case 505:htxs_stage1_catName = "ggZH_GT250"; break;
				case 600:htxs_stage1_catName = "TTH_FWDH" ;break;
				case 601:htxs_stage1_catName = "TTH_0_60" ;break;
				case 602:htxs_stage1_catName = "TTH_60_120" ;break;
				case 603:htxs_stage1_catName = "TTH_120_200" ;break;
				case 604:htxs_stage1_catName = "TTH_200_300" ;break;
				case 605:htxs_stage1_catName = "TTH_GT300" ;break;
				case 700:htxs_stage1_catName = "BBH_FWDH" ;break;
				case 701:htxs_stage1_catName = "BBH" ;break;
				case 800:htxs_stage1_catName = "TH_FWDH" ;break;
				case 801:htxs_stage1_catName = "TH" ;break;
			}
			htxs_stage1p2_catName= htxs_stage1_catName; 
			htxs_stage1_red_catName= htxs_stage1_catName; 
			htxs_stage1_red_cat = htxs_stage1p2_cat;
			if ((htxs_stage1p2_cat>112&&htxs_stage1p2_cat<117)){
				htxs_stage1_red_cat= 113;
				htxs_stage1_red_catName="ggH_VBF";
			}
			if ((htxs_stage1p2_cat>100&&htxs_stage1p2_cat<105)) {
				htxs_stage1_red_cat= 150;
				htxs_stage1_red_catName="ggH_GT200";
			}
			if (htxs_stage1p2_cat==208||htxs_stage1p2_cat==210){
				htxs_stage1_red_cat= 208;
				htxs_stage1_red_catName="VBF_2j_mjj_GT350_3j";
			}
			//Change 204 to 212 so that VH_Had_VBFori doesn't become 209 (same as VBF_2j_mjj_GT700_2j)
			if (splitVH) {
			if (htxs_stage1p2_cat==204) { 
			  if(inputname.Contains("ZH")) {	
					htxs_stage1_red_cat= 221;
					htxs_stage1_red_catName="ZH_Had";
				} 
				else if((inputname.Contains("WminusH") || inputname.Contains("WplusH"))) {
                                        htxs_stage1_red_cat= 222;
                                        htxs_stage1_red_catName="WH_Had";
				} else {	
					htxs_stage1_red_cat= 212;
                                        htxs_stage1_red_catName="VH_Had";
				}			
			}
			} else {

			if (htxs_stage1p2_cat==204) { 
				htxs_stage1_red_cat= 212;
				htxs_stage1_red_catName="VH_Had";				
			}

			}
			//Change 206 to 226 so that VBF_Rest_VHori doesn't become 206 (same as VBF_GT200)
			if (htxs_stage1p2_cat==206) { 
				htxs_stage1_red_cat= 226;
				htxs_stage1_red_catName="VBF_GT200";				
			}
				
			if ( (htxs_stage1p2_cat>200 && htxs_stage1p2_cat<204) || htxs_stage1p2_cat ==205){
				htxs_stage1_red_cat= 201;
				htxs_stage1_red_catName="VBF_Rest";
			}

			if (splitVH) {
			if (htxs_stage1p2_cat>300&&htxs_stage1p2_cat<400){
				if (htxs_stage1p2_cat%100==1 || htxs_stage1p2_cat%100==2){
					htxs_stage1_red_cat= 301;
					htxs_stage1_red_catName="WH_Lep_0_150";
				}
				else{
					htxs_stage1_red_cat= 303;
					htxs_stage1_red_catName="WH_Lep_GT150";
				}
			}
			if (htxs_stage1p2_cat>400&&htxs_stage1p2_cat<500){
				if (htxs_stage1p2_cat%100==1 || htxs_stage1p2_cat%100==2){
					htxs_stage1_red_cat= 401;
					htxs_stage1_red_catName="ZH_Lep_0_150";
				}
				else{
					htxs_stage1_red_cat= 403;
					htxs_stage1_red_catName="ZH_Lep_GT150";
				}
			}
			if (htxs_stage1p2_cat>500&&htxs_stage1p2_cat<600){
				if(htxs_stage1p2_cat<503){
					htxs_stage1_red_cat= 501;
					htxs_stage1_red_catName="ggZH_Lep_0_150";
				}
				else{
					htxs_stage1_red_cat= 503;
					htxs_stage1_red_catName="ggZH_Lep_GT150";
				}
			}

			} else {
			if (htxs_stage1p2_cat>300&&htxs_stage1p2_cat<500){
				if (htxs_stage1p2_cat%100==1 || htxs_stage1p2_cat%100==2){
					htxs_stage1_red_cat= 301;
					htxs_stage1_red_catName="VH_Lep_0_150";
				}
				else{
					htxs_stage1_red_cat= 303;
					htxs_stage1_red_catName="VH_Lep_GT150";
				}
			}
			if (htxs_stage1p2_cat>500&&htxs_stage1p2_cat<600){
				if(htxs_stage1p2_cat<503){
					htxs_stage1_red_cat= 301;
					htxs_stage1_red_catName="VH_Lep_0_150";
				}
				else{
					htxs_stage1_red_cat= 303;
					htxs_stage1_red_catName="VH_Lep_GT150";
				}
			}
			}

			if (htxs_stage1p2_cat>600&&htxs_stage1p2_cat<700) {
				htxs_stage1_red_cat= 601;
				htxs_stage1_red_catName="TTH";
			}

			if (inputname.Contains("ZZTo4l")){
				htxs_stage1_red_catName="qqZZ";
				htxs_stage1_red_cat = -1;
			}
			else if (isZX){
				htxs_stage1_red_catName="ZX";
				htxs_stage1_red_cat = -2;
			}
			else if (inputname.Contains("_Contin_MCFM701")){
				htxs_stage1_red_catName="ggZZ";
				htxs_stage1_red_cat = -3;
			}
			else if (inputname.Contains("TTZ") || inputname.Contains("ZZZ") || 
				inputname.Contains("WWZ") || inputname.Contains("WZZ") ||
				inputname.Contains("TTWW") || inputname.Contains("VBFToContinToZZ4l")){
				htxs_stage1_red_catName="EW_bkg";
				htxs_stage1_red_cat = -4;
			}

			htxs_stage1_red_prod_cat = htxs_stage1_red_cat;
			htxs_stage1_red_prod_catName = htxs_stage1_red_catName;

            // genExtInfo discriminates between Leptonic and Hadronic decay channels
			// we place this control here not to mistake a VH_Lep for VH_Had
      if(genExtInfo>10 && (htxs_prodMode==3 || htxs_prodMode==4))
              htxs_prodMode += 6;

			if (splitVH) {
			if ( (htxs_prodMode==3) && (htxs_stage1_red_cat == 226 || htxs_stage1_red_cat==201)){
				htxs_stage1_red_prod_cat +=5;
				htxs_stage1_red_prod_catName += "_WHori"; //"VBF_GT200_VHori","VBF_Rest_VHori"
			}
			else if ( (htxs_prodMode==4) && (htxs_stage1_red_cat == 226 || htxs_stage1_red_cat==201)){
                                htxs_stage1_red_prod_cat +=40;
                                htxs_stage1_red_prod_catName += "_ZHori"; //"VBF_GT200_VHori","VBF_Rest_VHori"
                        }
			else if (htxs_prodMode==2 && htxs_stage1_red_cat == 212){
				htxs_stage1_red_prod_cat +=5;
				htxs_stage1_red_prod_catName += "_VBFori"; //VH_Had_VBFori
			}
			} else {

			if ( (htxs_prodMode==3 || htxs_prodMode==4) && (htxs_stage1_red_cat == 226 || htxs_stage1_red_cat==201)){
				htxs_stage1_red_prod_cat +=5;
				htxs_stage1_red_prod_catName += "_VHori"; //"VBF_GT200_VHori","VBF_Rest_VHori"
			}
			else if (htxs_prodMode==2 && htxs_stage1_red_cat == 212){
				htxs_stage1_red_prod_cat +=5;
				htxs_stage1_red_prod_catName += "_VBFori"; //VH_Had_VBFori
			}

			}

			string reco_catName;
			string reco_1p1catName;

			htxs_stage1_reco_cat = stage1_reco_1p1_sync(njet,  mjj, ZZPt ,category,pt_hjj,reco_1p1catName);

			htxs_stage1_reco_catName = reco_1p1catName;
			// Not really needed
			if ( category > 7 ) htxs_stage1_reco_catName = "EW_bkg";

			TString categoryName_s;
			categoryName_s = convertStage0Name(category);
			strcpy(categoryName,categoryName_s.Data());
			switch(htxs_prodMode){
				case 0: htxs_prodModeName="TTV(V)";break; // Not really needed because htxs_prodModeName is used only for data
				case 1: htxs_prodModeName="ggH";break;
				case 2: htxs_prodModeName="VBF";break;
				case 3: htxs_prodModeName="WH_Had";break;
				case 4: htxs_prodModeName="ZH_Had";break;
				case 6: htxs_prodModeName="TTH";break;
				case 7: htxs_prodModeName="BBH";break;
				case 8: htxs_prodModeName="tqH";break;
				case 9: htxs_prodModeName="WH_Lep";break;
				case 10: htxs_prodModeName="ZH_Lep";break;
			}
			if(doSys){
				// sysname set using convertsys() method from stage1.cc
				TString sysname = convertsys(inputname);
				weight_up->clear();	
				weight_dn->clear();	
				weight_name->clear();	
				float weight_uptmp;
				float weight_dntmp;
				if(inputname.Contains("ggH")){
					for(int sfl=0;sfl<qcd_ggF_uncertSF.size();sfl++){
						weight_uptmp = weight* qcd_ggF_uncertSF.at(sfl);
						weight_dntmp = weight * (2-qcd_ggF_uncertSF.at(sfl)); 
						if(TMath::IsNaN(weight_uptmp))
							weight_uptmp=0;
						if(TMath::IsNaN(weight_dntmp))
							weight_dntmp=0;
						weight_up->push_back(weight_uptmp);
						weight_dn->push_back(weight_dntmp);
						weight_name->push_back(Form("qcd_ggH_%d",sfl));
					}
				}
        if(inputname.Contains("VBF")){
          for(int sfl=0;sfl<qcd_qqH_uncertSF.size();sfl++){
            weight_uptmp = weight * qcd_qqH_uncertSF.at(sfl);
            weight_dntmp = weight * (2-qcd_qqH_uncertSF.at(sfl));
            if(TMath::IsNaN(weight_uptmp))
                                        weight_uptmp = 0;
            if(TMath::IsNaN(weight_dntmp))
                                        weight_dntmp = 0;
            weight_up->push_back(weight_uptmp);
            weight_dn->push_back(weight_dntmp);
            weight_name->push_back(Form("qcd_qqH_%d", sfl));
          }
        }

				weight_uptmp = weight*PythiaWeight_isr_muR4    * PythiaWeight_fsr_muR4;
				weight_dntmp = weight*PythiaWeight_isr_muR0p25 * PythiaWeight_fsr_muR0p25;
				if(TMath::IsNaN(weight_uptmp))
					weight_uptmp=0;
				if(TMath::IsNaN(weight_dntmp))
					weight_dntmp=0;
				weight_up->push_back(weight_uptmp);
				weight_dn->push_back(weight_dntmp);
				weight_name->push_back("pythia");

				weight_uptmp = weight*LHEweight_QCDscale_muR2_muF1/LHEweight_QCDscale_muR1_muF1;
				weight_dntmp = weight*LHEweight_QCDscale_muR0p5_muF1/LHEweight_QCDscale_muR1_muF1;
				if(TMath::IsNaN(weight_uptmp))
					weight_uptmp=0;
				if(TMath::IsNaN(weight_dntmp))
					weight_dntmp=0;
				weight_up->push_back(weight_uptmp);
				weight_dn->push_back(weight_dntmp);
				weight_name->push_back("QCDscale_muR_"+sysname);

				weight_uptmp = weight*LHEweight_QCDscale_muR1_muF2/LHEweight_QCDscale_muR1_muF1;
				weight_dntmp = weight*LHEweight_QCDscale_muR1_muF0p5/LHEweight_QCDscale_muR1_muF1;
				if(TMath::IsNaN(weight_uptmp))
					weight_uptmp=0;
				if(TMath::IsNaN(weight_dntmp))
					weight_dntmp=0;
				weight_up->push_back(weight_uptmp);
				weight_dn->push_back(weight_dntmp);
				weight_name->push_back("QCDscale_muF_"+sysname);

				weight_uptmp = weight*LHEweight_PDFVariation_Up;
				weight_dntmp = weight*LHEweight_PDFVariation_Dn;
				if(TMath::IsNaN(weight_uptmp))
					weight_uptmp=0;
				if(TMath::IsNaN(weight_dntmp))
					weight_dntmp=0;
				weight_up->push_back(weight_uptmp);
				weight_dn->push_back(weight_dntmp);
				TString pdfsysname = "Higgs_qqbar";
				if (sysname.Contains("ggH"))
					pdfsysname = "Higgs_gg";
				else if (sysname.Contains("ttH")) 
				        pdfsysname = "Higgs_ttH";
				else if (sysname.Contains("qqZZ"))
					pdfsysname = "qqbar";
				else if (sysname.Contains("ggZZ"))
					pdfsysname = "gg";
				weight_name->push_back("pdf_"+pdfsysname);

				weight_uptmp = weight*LHEweight_AsMZ_Up;
				weight_dntmp = weight*LHEweight_AsMZ_Dn;
				if(TMath::IsNaN(weight_uptmp))
					weight_uptmp=0;
				if(TMath::IsNaN(weight_dntmp))
					weight_dntmp=0;
				weight_up->push_back(weight_uptmp);
				weight_dn->push_back(weight_dntmp);
				weight_name->push_back("pdf_As_"+pdfsysname);

				if(inputname.Contains("ZZTo")){
					weight_uptmp =  weight* ( 1. + KFactor_EW_qqZZ_unc ) ;
					weight_dntmp =  weight* ( 1. - KFactor_EW_qqZZ_unc ) ;
					if(TMath::IsNaN(weight_uptmp))
						weight_uptmp=0;
					if(TMath::IsNaN(weight_dntmp))
						weight_dntmp=0;
					weight_up->push_back(weight_uptmp);
					weight_dn->push_back(weight_dntmp);
					weight_name->push_back("EWCorr_"+sysname);
				}
				int _current_category_JES_UP = categoryMor18(nExtraLep,
						nExtraZ,
						nCleanedJetsPt30_jesUp,
						nCleanedJetsPt30BTagged_bTagSF_jesUp,
						jetQGL,
						p_JJQCD_SIG_ghg2_1_JHUGen_JESUp,
						p_JQCD_SIG_ghg2_1_JHUGen_JESUp,
						p_JJVBF_SIG_ghv1_1_JHUGen_JESUp,
						p_JVBF_SIG_ghv1_1_JHUGen_JESUp,
						pAux_JVBF_SIG_ghv1_1_JHUGen_JESUp,
						p_HadWH_SIG_ghw1_1_JHUGen_JESUp,
						p_HadZH_SIG_ghz1_1_JHUGen_JESUp,
						p_HadWH_mavjj_JESUp,
						p_HadWH_mavjj_true_JESUp,
						p_HadZH_mavjj_JESUp,
						p_HadZH_mavjj_true_JESUp,
						jetPhi,
						ZZMass,
						PFMET_corrected_jesUp,
						false,
						false);
				int _current_category_JER_UP = categoryMor18(nExtraLep,
						nExtraZ,
						nCleanedJetsPt30_jerUp,
						nCleanedJetsPt30BTagged_bTagSF_jerUp,
						jetQGL,
						p_JJQCD_SIG_ghg2_1_JHUGen_JERUp,
						p_JQCD_SIG_ghg2_1_JHUGen_JERUp,
						p_JJVBF_SIG_ghv1_1_JHUGen_JERUp,
						p_JVBF_SIG_ghv1_1_JHUGen_JERUp,
						pAux_JVBF_SIG_ghv1_1_JHUGen_JERUp,
						p_HadWH_SIG_ghw1_1_JHUGen_JERUp,
						p_HadZH_SIG_ghz1_1_JHUGen_JERUp,
						p_HadWH_mavjj_JERUp,
						p_HadWH_mavjj_true_JERUp,
						p_HadZH_mavjj_JERUp,
						p_HadZH_mavjj_true_JERUp,
						jetPhi,
						ZZMass,
						PFMET_corrected_jerUp,
						false,
						false);

				float ptmax_up=30;
				float ptmax_dn=30;
				vector <float> jetpt_up;
				vector <float> jetpt_dn;
				jetpt_up.clear();
				jetpt_dn.clear();
				for(int jecl =0;jecl<nCleanedJetsPt30;jecl++){
					float pt_up = JetPt->at(jecl)*(1+JetSigma->at(jecl));
					float pt_dn = JetPt->at(jecl)*(1-JetSigma->at(jecl));
					jetpt_up.push_back(pt_up);
					jetpt_dn.push_back(pt_dn);
					if(pt_up<ptmax_up)
						ptmax_up = pt_up;
					if(pt_dn<ptmax_dn)
						ptmax_dn = pt_dn;
					//cout <<"jet pt "<< pt_up<<endl;
				}
				std::sort(jetpt_up.begin(), jetpt_up.end(), std::greater<float>());
				std::sort(jetpt_dn.begin(), jetpt_dn.end(), std::greater<float>());

				//			if(jetpt_up.size()>1){
				//				cout <<"1st pt "<< jetpt_up .at(0)<<endl;
				//				cout <<"2nd pt "<< jetpt_up .at(1)<<endl;
				//			}

				int _current_category_JES_DN = categoryMor18(nExtraLep,
						nExtraZ,
						nCleanedJetsPt30_jesDn,
						nCleanedJetsPt30BTagged_bTagSF_jesDn,
						jetQGL,
						p_JJQCD_SIG_ghg2_1_JHUGen_JESDn,
						p_JQCD_SIG_ghg2_1_JHUGen_JESDn,
						p_JJVBF_SIG_ghv1_1_JHUGen_JESDn,
						p_JVBF_SIG_ghv1_1_JHUGen_JESDn,
						pAux_JVBF_SIG_ghv1_1_JHUGen_JESDn,
						p_HadWH_SIG_ghw1_1_JHUGen_JESDn,
						p_HadZH_SIG_ghz1_1_JHUGen_JESDn,
						p_HadWH_mavjj_JESDn,
						p_HadWH_mavjj_true_JESDn,
						p_HadZH_mavjj_JESDn,
						p_HadZH_mavjj_true_JESDn,
						jetPhi,
						ZZMass,
						PFMET_corrected_jesDn,
						false,
						false);

				int _current_category_JER_DN = categoryMor18(nExtraLep,
						nExtraZ,
						nCleanedJetsPt30_jerDn,
						nCleanedJetsPt30BTagged_bTagSF_jerDn,
						jetQGL,
						p_JJQCD_SIG_ghg2_1_JHUGen_JERDn,
						p_JQCD_SIG_ghg2_1_JHUGen_JERDn,
						p_JJVBF_SIG_ghv1_1_JHUGen_JERDn,
						p_JVBF_SIG_ghv1_1_JHUGen_JERDn,
						pAux_JVBF_SIG_ghv1_1_JHUGen_JERDn,
						p_HadWH_SIG_ghw1_1_JHUGen_JERDn,
						p_HadZH_SIG_ghz1_1_JHUGen_JERDn,
						p_HadWH_mavjj_JERDn,
						p_HadWH_mavjj_true_JERDn,
						p_HadZH_mavjj_JERDn,
						p_HadZH_mavjj_true_JERDn,
						jetPhi,
						ZZMass,
						PFMET_corrected_jerDn,
						false,
						false);

				int _current_category_BTag_UP = categoryMor18(nExtraLep,
						nExtraZ,
						nCleanedJetsPt30,
						nCleanedJetsPt30BTagged_bTagSFUp,
						jetQGL,
						p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
						p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
						p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
						p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
						pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
						p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
						p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
						p_HadWH_mavjj_JECNominal,
						p_HadWH_mavjj_true_JECNominal,
						p_HadZH_mavjj_JECNominal,
						p_HadZH_mavjj_true_JECNominal,
						jetPhi,
						ZZMass,
						PFMET_corrected,
						false,// Use VHMET category
						false);// Use QG tagging
				int _current_category_BTag_DN = categoryMor18(nExtraLep,
						nExtraZ,
						nCleanedJetsPt30,
						nCleanedJetsPt30BTagged_bTagSFDn,
						jetQGL,
						p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
						p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
						p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
						p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
						pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
						p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
						p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
						p_HadWH_mavjj_JECNominal,
						p_HadWH_mavjj_true_JECNominal,
						p_HadZH_mavjj_JECNominal,
						p_HadZH_mavjj_true_JECNominal,
						jetPhi,
						ZZMass,
						PFMET_corrected,
						false,// Use VHMET category
						false);// Use QG tagging

				//float D1jet_up = ( nCleanedJetsPt30_jecUp == 1 ) ? DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECUp, pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp, p_JQCD_SIG_ghg2_1_JHUGen_JECUp, ZZMass) : -2;
				//D1jet_up = 1/(newCConstD1jet/oldCConstD1jet*(1/D1jet_up-1)+1);
				//float D1jet_dn = ( nCleanedJetsPt30_jecDn == 1 ) ? DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECDn, pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn, p_JQCD_SIG_ghg2_1_JHUGen_JECDn, ZZMass) : -2;
				//D1jet_dn = 1/(newCConstD1jet/oldCConstD1jet*(1/D1jet_dn-1)+1);

				string reco_catName_jes_up; 
				string reco_catName_jes_dn; 
				string reco_catName_jer_up; 
				string reco_catName_jer_dn; 				
				// string reco_catName_jetpT_jes_up; 
				// string reco_catName_jetpT_jes_dn; 
				// string reco_catName_jetpT_jer_up; 
				// string reco_catName_jetpT_jer_dn; 
				string reco_catName_btag_dn; 
				string reco_catName_btag_up; 

				categoryName_jes_up=convertStage0Name(_current_category_JES_UP); 
				categoryName_jes_dn=convertStage0Name(_current_category_JES_DN); 
				categoryName_jer_up=convertStage0Name(_current_category_JER_UP); 
				categoryName_jer_dn=convertStage0Name(_current_category_JER_DN); 
				categoryName_btag_up=convertStage0Name(_current_category_BTag_UP); 
				categoryName_btag_dn=convertStage0Name(_current_category_BTag_DN); 

				category_jes_up = _current_category_JES_UP;
				category_jes_dn = _current_category_JES_DN;
				category_jer_up = _current_category_JER_UP;
				category_jer_dn = _current_category_JER_DN;
				category_btag_up = _current_category_BTag_UP;
				category_btag_dn = _current_category_BTag_DN;

				htxs_stage1_reco_cat_jes_up = stage1_reco_1p1_sync(nCleanedJetsPt30_jesUp, mjj_jes_up , ZZPt, _current_category_JES_UP, pt_hjj_jes_up, reco_catName_jes_up);
				htxs_stage1_reco_cat_jes_dn = stage1_reco_1p1_sync(nCleanedJetsPt30_jesDn, mjj_jes_dn , ZZPt, _current_category_JES_DN, pt_hjj_jes_dn, reco_catName_jes_dn);
				htxs_stage1_reco_cat_jer_up = stage1_reco_1p1_sync(nCleanedJetsPt30_jerUp, mjj_jer_up , ZZPt, _current_category_JER_UP, pt_hjj_jer_up, reco_catName_jer_up);
				htxs_stage1_reco_cat_jer_dn = stage1_reco_1p1_sync(nCleanedJetsPt30_jerDn, mjj_jer_dn , ZZPt, _current_category_JER_DN, pt_hjj_jer_dn, reco_catName_jer_dn);

				htxs_stage1_reco_cat_btag_up =stage1_reco_1p1_sync(njet, mjj, ZZPt, _current_category_BTag_UP, pt_hjj, reco_catName_btag_up);
				htxs_stage1_reco_cat_btag_dn =stage1_reco_1p1_sync(njet, mjj, ZZPt, _current_category_BTag_DN, pt_hjj, reco_catName_btag_dn);

				// int njet_jer_up=0;
				// int njet_jer_dn=0;
				// for (int ljer = 0;ljer<JetPt_JERUp->size();ljer++){
				// 	if (JetPt_JERUp->at(ljer)>30)
				// 		njet_jer_up++;
				// 	if (JetPt_JERDown->at(ljer)>30)
				// 		njet_jer_dn++;
				// }
				// int njet_jes_up=0;
				// int njet_jes_dn=0;
				// for (int ljes = 0;ljes<JetPt_JESUp->size();ljes++){
				// 	if (JetPt_JESUp->at(ljes)>30)
				// 		njet_jes_up++;
				// 	if (JetPt_JESDown->at(ljes)>30)
				// 		njet_jes_dn++;
				// }

				// htxs_stage1_reco_cat_jetPt_jes_up = stage1_reco_1p1_sync(njet_jes_up, mjj_jes_up , ZZPt, category,pt_hjj_jes_up,reco_catName_jetpT_jes_up);
				// htxs_stage1_reco_cat_jetPt_jes_dn = stage1_reco_1p1_sync(njet_jes_dn, mjj_jes_dn , ZZPt, category,pt_hjj_jes_dn,reco_catName_jetpT_jes_dn);
				// htxs_stage1_reco_cat_jetPt_jer_up = stage1_reco_1p1_sync(njet_jer_up, mjj_jer_up , ZZPt, category,pt_hjj_jer_up,reco_catName_jetpT_jer_up);
				// htxs_stage1_reco_cat_jetPt_jer_dn = stage1_reco_1p1_sync(njet_jer_dn, mjj_jer_dn , ZZPt, category,pt_hjj_jer_dn,reco_catName_jetpT_jer_dn);


				htxs_stage1_reco_catName_jes_up=reco_catName_jes_up;	
				htxs_stage1_reco_catName_jes_dn=reco_catName_jes_dn;	
				htxs_stage1_reco_catName_jer_up=reco_catName_jer_up;	
				htxs_stage1_reco_catName_jer_dn=reco_catName_jer_dn;	

				// htxs_stage1_reco_catName_jetPt_jes_up=reco_catName_jetpT_jes_up;	
				// htxs_stage1_reco_catName_jetPt_jes_dn=reco_catName_jetpT_jes_dn;	
				// htxs_stage1_reco_catName_jetPt_jer_up=reco_catName_jetpT_jer_up;	
				// htxs_stage1_reco_catName_jetPt_jer_dn=reco_catName_jetpT_jer_dn;	
				
				htxs_stage1_reco_catName_btag_up=reco_catName_btag_up;	
				htxs_stage1_reco_catName_btag_dn=reco_catName_btag_dn;	

			}
			tnew->Fill();
		}
		tqqzz->SetLineColor(2);
		tqqzz->SetMarkerColor(2);
		fnew->cd();
		fnew->Write();
		fnew->Close();
	}


