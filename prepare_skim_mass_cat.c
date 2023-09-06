using namespace std;
using namespace RooFit;
//      chan1: 4m
//      chan2: 4e
//      chan3: 2e2m
#include <fstream>      // std::ofstream
//#include "HiggsAnalysis/CombinedLimit/src/HZZ2L2QRooPdfs.cc" 

const int extraSysn_zjet = 3;
const int extraSysn_ggzz = 1;
const int extraSysn_all = 3;
const int extraCat= 4; //6;
// TString extraSysName[extraCat]={"jec","jer","btag","pythia_tune"};
// New nuisances for the new branches defined for JES and JER
// TString extraSysName[extraCat]={"jes","jer", "jetPt_jes", "jetPt_jer", "btag","pythia_tune"};
TString extraSysName[extraCat]={"jes", "jer", "btag","pythia_tune"};
//TString extraSys_shape[2]={"Res","Scale"};
TString extraSys_handName_ggzz [extraSysn_ggzz]={"ggzz_QCD_scale"};
TString extraSys_handName_zjet [extraSysn_zjet]={"zjet_4mu","zjet_4e","zjet_2e2mu"};
TString extraSys_handName_all[extraSysn_all]={"lumi_13TeV","CMS_eff_m","CMS_eff_e"};

// ZX yield uncertainties. 
// These values are the sum of three components:
// * Stat uncertainty
// * FR up/do variations
// * Bkg composition
// Updated values from: 
// https://elfontan.web.cern.ch/elfontan/CutBased_analysis/ZX_ESTIMATE/ZX_SystematicUncertainties.txt

float extraSys_zjet_up_2016 [extraSysn_zjet][3]= {
	1.30542,	0,	0, // 1.104,   0,      0,
	0,	1.42863,	0, // 0,     1.314,   0,
	0,	0,	1.35484 // 0,      0,      1.152
};
float extraSys_zjet_up [extraSysn_zjet][3]= {
	1.30685,	0,	0, // 1.32,	0,	0,  
	0,	1.37505,	0, // 0,	1.38,	0,
	0,	0,	1.33282 // 0,	0,	1.33
};

float extraSys_zjet_up_2018 [extraSysn_zjet][3]= {
	1.30459,	0,	0, // 1.30,	0,	0,  
	0,	1.36539,	0, // 0,	1.37,	0,
	0,	0,	1.32828 // 0,	0,	1.24
};
float extraSys_zjet_dn [extraSysn_zjet][3]= {
	0.69350,	0,	0, // 0.68, 	0,	0,
	0,	0.63816,	0, // 0., 	0.64,	0,
	0,	0,	0.67262 // 0., 	0,	0.67
};

float extraSys_zjet_dn_2016 [extraSysn_zjet][3]= {
	0.69481,	0,	0, // 0.899,   0,      0,
	0,	0.60745,	0, // 0.,     0.728,   0,
	0,	0,	0.65673 // 0.,     0,      0.868
};
float extraSys_zjet_dn_2018 [extraSysn_zjet][3]= {
	0.69559,	0,	0, // 0.70,   0,      0,
	0,	0.64540,	0, // 0.,     0.63,   0,
	0,	0,	0.67618 // 0.,     0,      0.76
};

// Columns: 4mu, 4e, 2e2mu
// Rows : LUMI, mu sys, ele sys
// Set to 1. the values corresponding to 
// ele sys in the 4mu channel and mu sys
// in the 4e channel
float extraSys_all_up_2016 [extraSysn_all][3]= { 
	1.026,  1.026,  1.026,
	1.016,		1.,	1.012, //1.046,  1.,     1.025,
	1.,		1.155,	1.105  //1.,     1.082,  1.039
};
float extraSys_all_up [extraSysn_all][3]= {
	1.023,	1.023,	1.023,
	1.011,	1.,		1.008, //1.056,	1.,	1.03,
	1.,		1.121,	1.082 //1.,	1.125,	1.058
};
float extraSys_all_up_2018 [extraSysn_all][3]= {
	1.025,	1.025,	1.025,
	1.01,	1.,		1.007, //1.016,	1.,	1.011,
	1.,		1.11,	1.074 //1.,	1.161,	1.074
};
float extraSys_all_dn_2016 [extraSysn_all][3]= {
	0.974,  0.974,  0.974,
	0.977,	1.,		0.983, //0.953,  1.,     0.975,
	1.,		0.835,  0.893  //1.,     0.914,  0.96
};
float extraSys_all_dn [extraSysn_all][3]= {
	0.977,	0.977,	0.977,
	0.98,	1.,		0.985, //0.937,	1.,	0.968,
	1.,		0.867,	0.915  //1.,	0.862,	0.939
};
float extraSys_all_dn_2018 [extraSysn_all][3]= {
	0.975,	0.975,	0.975,
	0.981,	1.,		0.986, //0.978,	1.,	0.992,
	1.,		0.877,	0.923 //1.,	0.850,	0.928
};
float extraSys_ggzz [extraSysn_ggzz][3]= {
	1.1, 1.1, 1.1
};

float res_unc_year[3][3] ={
	//4mu, 4e, 2e2mu
	0.2,0.2,0.2, // 2016
	0.2,0.2,0.2, // 2017
	0.12,0.17,0.17 // 2018
};	
float scale_unc_m_year[3][3] ={
	0.0004,0.0002,0,
	0.0004,0.0002,0,
	0.0004,0.0002,0
};	
float scale_unc_e_year[3][3] ={
	0.,0.003,0.0015,
	0.,0.003,0.0015,
	0.,0.003,0.0015,
};	
void Floor(TH2F* histo,int newb){
	TH2F *histo_new = (TH2F*)histo->Clone(Form("%s_new",histo->GetName())); 
	histo_new->RebinX(newb);
	histo_new->Smooth();
	for (int i=0;i<histo_new->GetNbinsX();i++){
		for (int j=0;j<histo_new->GetNbinsY();j++){
			if(histo_new->GetBinContent(i+1,j+1)==0)
				histo_new->SetBinContent(i+1,j+1,1.E-10);
		}
		float integral = histo_new->Integral(i+1,i+1);
		//histo->SetBinContent(i+1,j+1,histo->GetBinContent(i+1,j+1)/integral);
		for (int k =0;k<newb;k++){
			//cout<< i*newb+k+1<<endl;
			for (int j=0;j<histo_new->GetNbinsY();j++){
				histo->SetBinContent(i*newb+k+1,j+1,histo_new->GetBinContent(i+1,j+1)/integral);
			}
		}
	}
	//histo->Smooth();
}

void writeline(vector<TString> arr , ofstream &card){
	for (int i=0;i<arr.size();i++)
		card<< arr[i]<<"\t";
	card<<endl;
}

void writeline(vector<int> arr , ofstream &card){
  for (int i=0;i<arr.size();i++)
		card<< arr[i]<<"\t";
	card<<endl;
}
void writeline(vector<float> arr , ofstream &card){
	//	gStyle->SetPaintTextFormat("2.2f");
	card<< std::fixed;
	card<< std::setprecision(3);
	for (int i=0;i<arr.size();i++)
		card<< arr[i]<<"\t";
	card<<endl;
}

void prepare_skim_mass_cat(TString category="htxs_stage1_reco_cat", TString categoryName="htxs_stage1_reco_catName",TString outputDir="/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All",TString file_app = "_redTree",int doSys=1, int readPara=0, TString stageName = "htxs_stage1_red_prod_cat", TString year="2018", int flav=0, TString mass="125", int cate_vbfNum=1, int channelNum=1){

        int ch = channelNum-1;

	TString channelName[4]={"0","1","2","3"};
	
	TString channel=channelName[channelNum];

	TString cate_Name[2]={"0","1"};
	
	TString cate_vbf=cate_Name[cate_vbfNum];
	
	TString chanName[3] = {"4mu" , "4e", "2e2mu"}; 
	//TString stageName[2] = {"ggH", "VBF"};
	TString scale_e="CMS_scale_e";
	TString scale_m="CMS_scale_m";
	TString res_ch="CMS_res_"+chanName[flav];
	gStyle->SetOptStat(0);
	float lumi;
	if (year=="2016"){lumi = 35.9; 
		for (int i=0;i<extraSysn_zjet;i++){ 
			for (int k =0;k<3;k++) { 
				extraSys_zjet_up[i][k] = extraSys_zjet_up_2016[i][k];
				extraSys_zjet_dn[i][k] = extraSys_zjet_dn_2016[i][k];
			}
		}
		for (int i=0;i<extraSysn_all;i++){ 
			for (int k =0;k<3;k++) { 
				extraSys_all_up[i][k] = extraSys_all_up_2016[i][k];
				extraSys_all_dn[i][k] = extraSys_all_dn_2016[i][k];
			}
		}
	}
	else if(year=="2017") lumi = 41.5; 
	else if(year=="2018"){ 
		lumi = 59.7; 
		for (int i=0;i<extraSysn_zjet;i++){ 
			for (int k =0;k<3;k++) { 
				extraSys_zjet_up[i][k] = extraSys_zjet_up_2018[i][k];
				extraSys_zjet_dn[i][k] = extraSys_zjet_dn_2018[i][k];
			}
		}
		for (int i=0;i<extraSysn_all;i++){ 
			for (int k =0;k<3;k++) { 
				extraSys_all_up[i][k] = extraSys_all_up_2018[i][k];
				extraSys_all_dn[i][k] = extraSys_all_dn_2018[i][k];
			}
		}
	}
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
	//	chan1: 4m
	//	chan2: 4e
	//	chan3: 2e2m

	TChain *tsig = new TChain ("SelectedTree");
	TChain *tzx = new TChain ("SelectedTree");
	TChain *tqqzz = new TChain ("SelectedTree");
	TChain *tggzz = new TChain ("SelectedTree");
	TChain *tggh = new TChain ("SelectedTree");
	TChain *tvbf = new TChain ("SelectedTree");
	TChain *tggh125 = new TChain ("SelectedTree");
	TChain *tvbf125 = new TChain ("SelectedTree");
	TChain *tgghint = new TChain ("SelectedTree");
	TChain *tvbfint = new TChain ("SelectedTree");
	TChain *tqqzzew = new TChain ("SelectedTree");
	
	TString inputDir = "/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All/ReducedTrees";//"input_newProd";
	TChain *tdata = new TChain("SelectedTree");
	tdata->Add(inputDir+"/AllData"+file_app+"_"+year+".root");
	tzx->Add(inputDir+"/AllData_ZX"+file_app+"*"+year+".root");

	/* tzx->Add(inputDir+"/ZZZ"+file_app+"*"+year+".root"); */
	/* tzx->Add(inputDir+"/WWZ"+file_app+"*"+year+".root"); */
	/* tzx->Add(inputDir+"/WZZ"+file_app+"*"+year+".root"); */
	/* tzx->Add(inputDir+"/TTZZ"+file_app+"*"+year+".root"); */
	/* tzx->Add(inputDir+"/TTWW"+file_app+"*"+year+".root"); */
	/* tzx->Add(inputDir+"/TTZJets_M10_MLM"+file_app+"*"+year+".root"); */
	/* tzx->Add(inputDir+"/TTZToLL*"+file_app+"*"+year+".root"); */
	/* tzx->Add(inputDir+"/VBFTo*"+file_app+"*"+year+".root"); */

	tqqzz->Add(inputDir+"/ZZTo4lext"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tggzz->Add(inputDir+"/ggTo*_Contin_MCFM701"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tsig->Add(inputDir+"/VBFH*"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tsig->Add(inputDir+"/ggH*"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tvbf->Add(inputDir+"/VBFH*"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tggh->Add(inputDir+"/ggH*"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tvbf125->Add(inputDir+"/VBFH125"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tggh125->Add(inputDir+"/ggH125"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tvbfint->Add(inputDir+"/VBFToHiggs0MToZZTo4l_M125_GaSM"+file_app+"_"+year+".root");
	tgghint->Add(inputDir+"/ggTo*_0MH125Contin_MCFM701"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	tqqzzew->Add(inputDir+"/VBFToContinToZZ4l"+file_app+"_chan"+channel+"_cat"+cate_vbf+"_"+year+".root");
	/* tsig->Add(inputDir+"/WplusH125"+file_app+"*"+year+".root"); */
	/* tsig->Add(inputDir+"/WminusH125"+file_app+"*"+year+".root"); */
	/* tsig->Add(inputDir+"/ZH125"+file_app+"*"+year+".root"); */
	/* tsig->Add(inputDir+"/ggZH125"+file_app+"*"+year+".root"); */
	/* tsig->Add(inputDir+"/ttH125"+file_app+"*"+year+".root"); */
	/* tsig->Add(inputDir+"/tqH125"+file_app+"*"+year+".root"); */
	/* tsig->Add(inputDir+"/bbH125"+file_app+"*"+year+".root"); */
	/* tsig->Add(inputDir+"/tHW125"+file_app+"*"+year+".root"); */
	
	tsig->Draw(categoryName+">>h_reco");
	TH1F *h=(TH1F*)gROOT->FindObject("h_reco")->Clone();
	int reco_cat_all = h->GetNbinsX();
	vector<int> reco_cat_num;
	cout<< "reconstructed categories "<<reco_cat_all <<endl;
	for(int i=0;i<reco_cat_all;i++){
		cout<< h->GetXaxis()->GetBinLabel(i+1)<<endl;
		TString reco_name = h->GetXaxis()->GetBinLabel(i+1);
		tsig->Draw(Form("%s",category.Data()),Form("%s==\"%s\"",categoryName.Data(),reco_name.Data()));
		int reco_cat = *(tsig->GetV1());
		reco_cat_num.push_back(reco_cat);
	}
	cout<<stageName.Data()<<endl;
	tsig->Draw(Form("%sName>>h_stage1",stageName.Data()),Form("%s>0 && %s%100!=0",stageName.Data(),stageName.Data()));
	TH1F *h_stage1=(TH1F*)gROOT->FindObject("h_stage1")->Clone();
	const int stage1_cat = h_stage1->GetNbinsX();
	for(int i=0;i<stage1_cat;i++){
		cout<< h_stage1->GetXaxis()->GetBinLabel(i+1)<<endl;
	}
	cout<< "stage1 categories "<< stage1_cat <<endl;
	//tsig->Draw("weight_name>>h_sys","!weight_name.Contains(\"qcd_ggH\")");
	//tsig->Draw("weight_name>>h_sys","!weight_name.Contains(\"qcd_ggH\") && !weight_name.Contains(\"qcd_qqH\") && !weight_name.Contains(\"ggH_scale\")");
	
	tsig->Draw("weight_name>>h_sys","!weight_name.Contains(\"qcd_ggH\") && !weight_name.Contains(\"qcd_qqH\") && !weight_name.Contains(\"_scale\")");
	TH1F *h_sys=(TH1F*)gROOT->FindObject("h_sys")->Clone();
	vector <TString> sysnames ;
	for(int i=0;i<h_sys->GetNbinsX();i++){
		TString obl = h_sys->GetXaxis()->GetBinLabel(i+1);
		obl+="_cat";
		sysnames.push_back(obl);
	}
	tqqzz->Draw("weight_name>>h_sys");
	tsig->Draw("weight_name>>+h_sys");
	h_sys=(TH1F*)gROOT->FindObject("h_sys")->Clone();
	for(int i=0;i<h_sys->GetNbinsX();i++){
		sysnames.push_back( h_sys->GetXaxis()->GetBinLabel(i+1));
	}
	int nsys = sysnames.size();

	cout<< "Systematic uncertainties "<<nsys<<endl;

	RooRealVar* ZZMass=new RooRealVar("ZZMass","",100,3500);
	RooRealVar* GenHMass=new RooRealVar("GenHMass","",100,3500);
	RooRealVar* MH=new RooRealVar("MH","",100,100,3500);
	//	RooConstVar* MH=new RooConstVar("MH","",125);
	RooRealVar* dbkg=new RooRealVar("dbkg_kin","",0,1);
	RooRealVar* D1jet=new RooRealVar("D1jet","",0.7,1);
	RooRealVar* D1jet_norm=new RooRealVar("D1jet_norm","",0.5,1);
	RooRealVar* ZZPt=new RooRealVar("ZZPt","",0,4500);
	RooRealVar* DVHDEC=new RooRealVar("DVHDEC","",0,1);
	RooRealVar* DVBFDEC=new RooRealVar("DVBFDEC","",0,1);

	dbkg->setBins(20);
	D1jet->setBins(10);
	D1jet_norm->setBins(15);
	ZZPt->setBins(6);
	DVHDEC->setBins(20);
	DVBFDEC->setBins(20);

	RooRealVar *Ddis=dbkg;
	cout<<"NBins "<<Ddis->getBins()<<endl;

	RooRealVar* weight=new RooRealVar("weight","",0,10000);
	RooRealVar* chan=new RooRealVar("chan","",1,5);
	//RooRealVar* vbfcate=new RooRealVar("vbfcate","",0,1000);
	RooRealVar* cat=new RooRealVar(category,"",0,1000);
	RooRealVar* stage=new RooRealVar(stageName,"",0,1000);
	

	RooArgSet ntupleVarSet(*GenHMass,*weight,*chan,*cat,*stage);
	RooArgSet ntupleVarSet_bkg(*GenHMass,*weight,*chan,*cat);

	RooDataSet *dataset = new RooDataSet("mc","mc",tsig,ntupleVarSet,"","weight");
	RooDataSet *dataset_stage[3][stage1_cat];
	RooDataSet *dataset_qqzz = new RooDataSet("qqzz","qqzz",tqqzz,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_ggzz = new RooDataSet("ggzz","ggzz",tggzz,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_qqzzew = new RooDataSet("qqzzew","qqzzew",tqqzzew,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_ggh = new RooDataSet("ggh","ggh",tggh,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_vbf = new RooDataSet("vbf","vbf",tvbf,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_ggh125 = new RooDataSet("ggh125","ggh125",tggh125,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_vbf125 = new RooDataSet("vbf125","vbf125",tvbf125,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_gghint = new RooDataSet("gghint","gghint",tgghint,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_vbfint = new RooDataSet("vbfint","vbfint",tvbfint,ntupleVarSet_bkg,"","weight");
	RooDataSet *dataset_zx = new RooDataSet("zx","zx",tzx,ntupleVarSet_bkg,"","weight");
	RooAbsPdf *pdf_stage[3][stage1_cat];
	RooAbsPdf *pdf_stage_mean_e_up[3][stage1_cat];
	RooAbsPdf *pdf_stage_mean_e_dn[3][stage1_cat];
	RooAbsPdf *pdf_stage_mean_m_up[3][stage1_cat];
	RooAbsPdf *pdf_stage_mean_m_dn[3][stage1_cat];
	RooAbsPdf *pdf_stage_sigmaup[3][stage1_cat];
	RooAbsPdf *pdf_stage_sigmadn[3][stage1_cat];

	TH2F *template2d[3][stage1_cat];
	RooHistPdf *pdf2d[3][stage1_cat];
	RooDataHist *datahist[3][stage1_cat];
	int stage_cat_arr[stage1_cat];

	TH2F *template2d_zx [3];
	TH2F *template2d_qqzz [3];
	TH2F *template2d_qqzzew [3];
	TH2F *template2d_ggzz [3];
	TH2F *template2d_ggh [3];
	TH2F *template2d_vbf [3];
	TH2F *template2d_ggh125 [3];
	TH2F *template2d_vbf125 [3];
	TH2F *template2d_gghint [3];
	TH2F *template2d_vbfint [3];
	RooHistPdf *pdfzx2d[3];
	RooDataHist *datahistzx[3];
	RooHistPdf *pdfqqzz2d[3];
	RooHistPdf *pdfqqzzew2d[3];
	RooHistPdf *pdfggzz2d[3];
	RooHistPdf *pdfggh2d[3];
	RooHistPdf *pdfvbf2d[3];
	RooHistPdf *pdfggh125_2d[3];
	RooHistPdf *pdfvbf125_2d[3];
	RooHistPdf *pdfgghint_2d[3];
	RooHistPdf *pdfvbfint_2d[3];
	RooDataHist *datahistqqzz[3];
	RooDataHist *datahistqqzzew[3];
	RooDataHist *datahistggzz[3];
	RooDataHist *datahistggh[3];
	RooDataHist *datahistvbf[3];
	RooDataHist *datahistggh125[3];
	RooDataHist *datahistvbf125[3];
	RooDataHist *datahistgghint[3];
	RooDataHist *datahistvbfint[3];
	RooBernstein *pdf_qqzz_chan[3];
	RooBernstein *pdf_qqzzew_chan[3];
	RooBernstein *pdf_ggzz_chan[3];
	RooBernstein *pdf_ggh_chan[3];
	RooBernstein *pdf_vbf_chan[3];
	RooBernstein *pdf_ggh125_chan[3];
	RooBernstein *pdf_vbf125_chan[3];
	RooBernstein *pdf_gghint_chan[3];
	RooBernstein *pdf_vbfint_chan[3];
	RooGenericPdf *pdf_zx_chan[3];

	/* RooRealVar* ber1=new RooRealVar("ber1_","ber1_",1.21409,-5,5); */
	/* RooRealVar* ber2=new RooRealVar("ber2_","ber2_",1.54396,-5,5); */
	/* RooRealVar* ber3=new RooRealVar("ber3_","ber3_",1.71314,-5,5); */

	/* RooRealVar* ber4=new RooRealVar("ber4_","ber4_",1.21409,-5,5); */
	/* RooRealVar* ber5=new RooRealVar("ber5_","ber5_",1.54396,-5,5); */
	/* RooRealVar* ber6=new RooRealVar("ber6_","ber6_",1.71314,-5,5); */

	/* RooRealVar* land1=new RooRealVar("land1_","land1_",105, 140); */
	/* RooRealVar* land2=new RooRealVar("land2_","land2_",10,5,20); */
	/* RooRealVar* frac_land=new RooRealVar("frac_land","frac_land",0,1); */
	/* RooRealVar* land1_2=new RooRealVar("land1_2_","land1_2",105, 140); */
	/* RooRealVar* land2_2=new RooRealVar("land2_2","land2_2",10,5,20); */

	/* RooBernstein *pdf_qqzz = new RooBernstein("pdfm4l_qqzz_", "pdfm4l_qqzz_",*ZZMass,RooArgList(*ber1,*ber2,*ber3)); */
	/* RooBernstein *pdf_ggzz = new RooBernstein("pdfm4l_ggzz_", "pdfm4l_ggzz_",*ZZMass,RooArgList(*ber4,*ber5,*ber6));	 */
	/* RooGenericPdf* pdf_zx= new RooGenericPdf("pdfm4l_zjets_", "pdfm4l_zjets_","TMath::Landau(@0,@1,@2)",RooArgList(*ZZMass,*land1,*land2));	 */
	/* RooGenericPdf* pdf_zx_2land= new RooGenericPdf("pdfm4l_zjets_2lan", "pdfm4l_zjets_2lan","TMath::Landau(@0,@1,@2)+@3*TMath::Landau(@0,@4,@5)",RooArgList(*ZZMass,*land1,*land2,*frac_land,*land1_2,*land2_2));	 */

	ifstream massPara;
	ifstream yieldPara;
	ifstream sysPara;
	//	ofstream massPara_output;
	massPara.open(inputDir+"/sim_massParam_"+category+chanName[flav]+year+".txt");
	yieldPara.open(inputDir+"/Param_"+year+".txt");
	sysPara.open(inputDir+"/Sys_"+year+".txt");
	//	massPara_output.open(outputDir+"/fig/massParam_"+category+chanName[flav]+year+".txt");

	map <TString, TString> massParaMap;
	string line;
	while( getline( massPara, line ) ){
		istringstream iss( line );
		TString paraName;
		TString paraValue;
		iss>> paraName >> paraValue;
		massParaMap.insert(pair<TString, TString>(paraName, paraValue));
		//cout<< paraName<<"\t"<<paraValue<<endl;
	}
	map <TString, TString> yieldParaMap;
	while( getline( yieldPara, line ) ){
		istringstream iss( line );
		TString paraName;
		TString paraValue;
		iss>> paraName >> paraValue;
		yieldParaMap.insert(pair<TString, TString>(paraName, paraValue));
		//cout<< paraName<<"\t"<<paraValue<<endl;
	}
	map <TString, float> sysParaMap;
	while( getline( sysPara, line ) ){
		istringstream iss( line );
		TString paraName;
		float paraValue;
		iss>> paraName >> paraValue;
		sysParaMap.insert(pair<TString, float>(paraName, paraValue));
	}

	
	//for (int ch=0;ch<3;ch++){
	  //for (int ch=flav;ch<flav+1;ch++){
		//for (int ch=2;ch<3;ch++){
	        int nbins=29;
	        double xbin[30]={
		  100,110,120,130,140,150,160,170,180,190,200,
		  210,220,230,240,250,270,290,310,330,350,
		  370,390,420,460,500,700,1000,1500,3500
		};
	       	  
		template2d_zx[ch] = new TH2F("template_zx_"+chanName[ch],"",nbins,xbin,30,0,1.);
		template2d_qqzz[ch] = new TH2F("template_qqzz_"+chanName[ch],"",nbins,xbin,30,0,1.);
		template2d_ggzz[ch] = new TH2F("template_ggzz_"+chanName[ch],"",nbins,xbin,30,0,1.);
		template2d_qqzzew[ch] = new TH2F("template_qqzzew_"+chanName[ch],"",nbins,xbin,30,0,1.);

		tqqzz->Draw("dbkg_kin:ZZMass>>template_qqzz_"+chanName[ch],Form("weight*(chan==%d)",ch+1));
		tggzz->Draw("dbkg_kin:ZZMass>>template_ggzz_"+chanName[ch],Form("weight*(chan==%d)",ch+1));
		tqqzzew->Draw("dbkg_kin:ZZMass>>template_qqzzew_"+chanName[ch],Form("weight*(chan==%d)",ch+1));

		template2d_ggh[ch] = new TH2F("template_ggh_"+chanName[ch],"",nbins,xbin,30,0,1.);
		template2d_vbf[ch] = new TH2F("template_vbf_"+chanName[ch],"",nbins,xbin,30,0,1.);

		tggh->Draw("dbkg_kin:ZZMass>>template_ggh_"+chanName[ch],Form("weight*(chan==%d)",ch+1));
		tvbf->Draw("dbkg_kin:ZZMass>>template_vbf_"+chanName[ch],Form("weight*(chan==%d)",ch+1));
		tzx->Draw("dbkg_kin:ZZMass>>template_zx_"+chanName[ch],Form("weight*(chan==%d)",ch+1));

		template2d_ggh125[ch] = new TH2F("template_ggh125_"+chanName[ch],"",nbins,xbin,30,0,1.);
		template2d_vbf125[ch] = new TH2F("template_vbf125_"+chanName[ch],"",nbins,xbin,30,0,1.);

		tggh125->Draw("dbkg_kin:ZZMass>>template_ggh125_"+chanName[ch],Form("weight*(chan==%d)",ch+1));
		tvbf125->Draw("dbkg_kin:ZZMass>>template_vbf125_"+chanName[ch],Form("weight*(chan==%d)",ch+1));

		template2d_gghint[ch] = new TH2F("template_gghint_"+chanName[ch],"",nbins,xbin,30,0,1.);
		template2d_vbfint[ch] = new TH2F("template_vbfint_"+chanName[ch],"",nbins,xbin,30,0,1.);

		tgghint->Draw("dbkg_kin:ZZMass>>template_gghint_"+chanName[ch],Form("weight*(chan==%d)",ch+1));
		tvbfint->Draw("dbkg_kin:ZZMass>>template_vbfint_"+chanName[ch],Form("weight*(chan==%d)",ch+1));

		//		template2d_zx[ch]->Smooth();
		//		template2d_qqzz[ch]->Smooth();
		Floor(template2d_zx[ch],1);
		Floor(template2d_qqzz[ch],1);
		Floor(template2d_qqzzew[ch],1);
		Floor(template2d_ggzz[ch],1);
		Floor(template2d_ggh[ch],1);
		Floor(template2d_vbf[ch],1);
		Floor(template2d_ggh125[ch],1);
		Floor(template2d_vbf125[ch],1);
		Floor(template2d_gghint[ch],1);
		Floor(template2d_vbfint[ch],1);
		template2d_zx[ch]->Smooth();
		template2d_qqzz[ch]->Smooth();
		template2d_qqzzew[ch]->Smooth();
		template2d_ggzz[ch]->Smooth();
		template2d_ggh[ch]->Smooth();
		template2d_vbf[ch]->Smooth();
		template2d_ggh125[ch]->Smooth();
		template2d_vbf125[ch]->Smooth();
		template2d_gghint[ch]->Smooth();
		template2d_vbfint[ch]->Smooth();
		//template2d_zx[ch]->Draw("colz");
		template2d_zx[ch]->SetMaximum(0.15);
		template2d_zx[ch]->SetTitle("Z+X "+chanName[ch]);
		template2d_zx[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_zx[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_zx[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/zx_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/zx_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");
		//template2d_ggzz[ch]->Draw("colz");
		template2d_ggzz[ch]->SetMaximum(0.15);
		template2d_ggzz[ch]->SetTitle("ggZZ "+chanName[ch]);
		template2d_ggzz[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_ggzz[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_ggzz[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/ggzz_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/ggzz_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");
		//template2d_qqzz[ch]->Draw("colz");
		template2d_qqzz[ch]->SetMaximum(0.15);
		template2d_qqzz[ch]->SetTitle("qqZZ "+chanName[ch]);
		template2d_qqzz[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_qqzz[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_qqzz[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/qqzz_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/qqzz_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");

		template2d_qqzzew[ch]->SetMaximum(0.15);
		template2d_qqzzew[ch]->SetTitle("qqZZ EW"+chanName[ch]);
		template2d_qqzzew[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_qqzzew[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_qqzzew[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/qqzzew_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/qqzzew_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");

		datahistqqzz[ch]=new RooDataHist("datahistqqzz_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistqqzz_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_qqzz[ch]);
		pdfqqzz2d[ch]=new RooHistPdf("histpdfqqzz_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfqqzz_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistqqzz[ch]);

		datahistqqzzew[ch]=new RooDataHist("datahistqqzzew_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistqqzzew_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_qqzzew[ch]);
		pdfqqzzew2d[ch]=new RooHistPdf("histpdfqqzzew_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfqqzzew_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistqqzzew[ch]);

		datahistggzz[ch]=new RooDataHist("datahistggzz_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistggzz_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_ggzz[ch]);
		pdfggzz2d[ch]=new RooHistPdf("histpdfggzz_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfggzz_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistggzz[ch]);

		//template2d_ggh[ch]->Draw("colz");
		template2d_ggh[ch]->SetMaximum(0.15);
		template2d_ggh[ch]->SetTitle("ggH "+chanName[ch]);
		template2d_ggh[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_ggh[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_ggh[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/ggh_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/ggh_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");
		//template2d_vbf[ch]->Draw("colz");
		template2d_vbf[ch]->SetMaximum(0.15);
		template2d_vbf[ch]->SetTitle("VBFH "+chanName[ch]);
		template2d_vbf[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_vbf[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_vbf[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/vbf_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/vbf_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");

		datahistvbf[ch]=new RooDataHist("datahistvbf_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistvbf_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_vbf[ch]);
		pdfvbf2d[ch]=new RooHistPdf("histpdfvbf_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfvbf_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistvbf[ch]);

		datahistggh[ch]=new RooDataHist("datahistggh_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistggh_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_ggh[ch]);
		pdfggh2d[ch]=new RooHistPdf("histpdfggh_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfggh_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistggh[ch]);

		datahistzx[ch]=new RooDataHist("datahistzx_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistzx_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_zx[ch]);
		pdfzx2d[ch]=new RooHistPdf("histpdfzx_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfzx_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistzx[ch]);

		//template2d_ggh125[ch]->Draw("colz");
		template2d_ggh125[ch]->SetMaximum(0.15);
		template2d_ggh125[ch]->SetTitle("ggH 125 "+chanName[ch]);
		template2d_ggh125[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_ggh125[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_ggh125[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/ggh125_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/ggh125_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");
		//template2d_vbf125[ch]->Draw("colz");
		template2d_vbf125[ch]->SetMaximum(0.15);
		template2d_vbf125[ch]->SetTitle("VBFH 125 "+chanName[ch]);
		template2d_vbf125[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_vbf125[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_vbf125[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/vbf125_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/vbf125_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");

		datahistvbf125[ch]=new RooDataHist("datahistvbf125_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistvbf125_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_vbf125[ch]);
		pdfvbf125_2d[ch]=new RooHistPdf("histpdfvbf125_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfvbf125_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistvbf125[ch]);

		datahistggh125[ch]=new RooDataHist("datahistggh125_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistggh125_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_ggh125[ch]);
		pdfggh125_2d[ch]=new RooHistPdf("histpdfggh125_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfggh125_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistggh125[ch]);

		template2d_gghint[ch]->SetMaximum(0.15);
		template2d_gghint[ch]->SetTitle("ggH int "+chanName[ch]);
		template2d_gghint[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_gghint[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_gghint[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/gghint_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/gghint_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");
		//template2d_vbfint[ch]->Draw("colz");
		template2d_vbfint[ch]->SetMaximum(0.15);
		template2d_vbfint[ch]->SetTitle("VBFH int "+chanName[ch]);
		template2d_vbfint[ch]->GetXaxis()->SetTitle("ZZ Mass [GeV]");
		template2d_vbfint[ch]->GetYaxis()->SetTitle("Dkin");
		template2d_vbfint[ch]->Draw("colz");
		gPad->Print(outputDir+"/fig/vbfint_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.png");
		gPad->Print(outputDir+"/fig/vbfint_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+"template.pdf");

		datahistvbfint[ch]=new RooDataHist("datahistvbfint_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistvbfint_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_vbfint[ch]);
		pdfvbfint_2d[ch]=new RooHistPdf("histpdfvbfint_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfvbfint_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistvbfint[ch]);

		datahistgghint[ch]=new RooDataHist("datahistgghint_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"datahistgghint_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),template2d_gghint[ch]);
		pdfgghint_2d[ch]=new RooHistPdf("histpdfgghint_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year,"histpdfgghint_"+chanName[ch],RooArgSet(*GenHMass,*dbkg),*datahistgghint[ch]);

		TFile* templates_int_HM = new TFile(outputDir+"/template/templates_int_HM_"+chanName[ch]+"_chan"+channel+"_cat"+cate_vbf+year+".root","recreate");

		templates_int_HM->WriteObject(template2d_zx[ch],"template2d_zx");
		templates_int_HM->WriteObject(template2d_qqzz[ch],"template2d_qqzz");
		templates_int_HM->WriteObject(template2d_qqzzew[ch],"template2d_qqzzew");
		templates_int_HM->WriteObject(template2d_ggzz[ch],"template2d_ggzz");
		templates_int_HM->WriteObject(template2d_ggh[ch],"template2d_ggh");
		templates_int_HM->WriteObject(template2d_vbf[ch],"template2d_vbf");
	       	templates_int_HM->WriteObject(template2d_gghint[ch],"template2d_gghint");
		templates_int_HM->WriteObject(template2d_vbfint[ch],"template2d_vbfint");
		templates_int_HM->WriteObject(template2d_ggh125[ch],"template2d_ggh125");
		templates_int_HM->WriteObject(template2d_vbf125[ch],"template2d_vbf125");
		//templates_int_HM->Write();

		templates_int_HM->Close();		
	
	/* 	RooDataSet *red_qqzz = (RooDataSet*)dataset_qqzz->reduce(Form("chan==%d",ch+1)); */
	/* 	RooDataSet *red_ggzz = (RooDataSet*)dataset_ggzz->reduce(Form("chan==%d",ch+1)); */
	/* 	RooDataSet *red_zx = (RooDataSet*)dataset_zx->reduce(Form("chan==%d",ch+1)); */
	/* 	cout<<red_qqzz->sumEntries()<<endl; */
	/* 	cout<<red_ggzz->sumEntries()<<endl; */
	/* 	cout<<red_zx->sumEntries()<<endl; */

	/* 	int lineN = 0; */
	/* 	pdf_qqzz->fitTo(*red_qqzz); */
	/* 	pdf_ggzz->fitTo(*red_ggzz); */
	/* 	// if(ch==3) */
	/* 	// 	pdf_zx_2land->fitTo(*red_zx); */
	/* 	// else */
	/* 	pdf_zx->fitTo(*red_zx); */
	/* 	//massPara_output << "ber1_"+chanName[ch]<< "\t"<< ber1->getVal()<<endl; */
	/* 	//massPara_output << "ber2_"+chanName[ch]<< "\t"<< ber2->getVal()<<endl; */
	/* 	//massPara_output << "ber3_"+chanName[ch]<< "\t"<< ber3->getVal()<<endl; */

	/* 	//massPara_output << "ber4_"+chanName[ch]<< "\t"<< ber4->getVal()<<endl; */
	/* 	//massPara_output << "ber5_"+chanName[ch]<< "\t"<< ber5->getVal()<<endl; */
	/* 	//massPara_output << "ber6_"+chanName[ch]<< "\t"<< ber6->getVal()<<endl; */

	/* 	//massPara_output << "land1_"+chanName[ch]<< "\t"<< land1->getVal()<<endl; */
	/* 	//massPara_output << "land2_"+chanName[ch]<< "\t"<< land2->getVal()<<endl; */
	/* 	//massPara_output << "land1_2_"+chanName[ch]<< "\t"<< land1_2->getVal()<<endl; */
	/* 	//massPara_output << "land2_2_"+chanName[ch]<< "\t"<< land2_2->getVal()<<endl; */
	/* 	//massPara_output << "frac_land"+chanName[ch]<< "\t"<< frac_land->getVal()<<endl; */
	/* 	RooConstVar* cber1=new RooConstVar("ber1_"+chanName[ch]+year,"ber1_"+chanName[ch²]+year,ber1->getVal()); */
	/* 	RooConstVar* cber2=new RooConstVar("ber2_"+chanName[ch]+year,"ber2_"+chanName[ch]+year,ber2->getVal()); */
	/* 	RooConstVar* cber3=new RooConstVar("ber3_"+chanName[ch]+year,"ber3_"+chanName[ch]+year,ber3->getVal()); */

	/* 	RooConstVar* cber4=new RooConstVar("ber4_"+chanName[ch]+year,"ber4_"+chanName[ch]+year,ber4->getVal()); */
	/* 	RooConstVar* cber5=new RooConstVar("ber5_"+chanName[ch]+year,"ber5_"+chanName[ch]+year,ber5->getVal()); */
	/* 	RooConstVar* cber6=new RooConstVar("ber6_"+chanName[ch]+year,"ber6_"+chanName[ch]+year,ber6->getVal()); */

	/* 	RooConstVar* cland1=new RooConstVar("land1_"+chanName[ch]+year,"land1_"+chanName[ch]+year,land1->getVal()); */
	/* 	RooConstVar* cland2=new RooConstVar("land2_"+chanName[ch]+year,"land2_"+chanName[ch]+year,land2->getVal()); */
	/* 	RooConstVar* cland1_2=new RooConstVar("land1_2_"+chanName[ch]+year,"land1_2_"+chanName[ch]+year,land1_2->getVal()); */
	/* 	RooConstVar* cland2_2=new RooConstVar("land2_2_"+chanName[ch]+year,"land2_2_"+chanName[ch]+year,land2_2->getVal()); */
	/* 	RooConstVar* cfrac_land=new RooConstVar("frac_land_"+chanName[ch]+year,"frac_land_"+chanName[ch]+year,frac_land->getVal()); */

	/* 	pdf_qqzz_chan[ch] = new RooBernstein("pdfm4l_qqzz_"+chanName[ch]+year, "pdfm4l_qqzz_"+chanName[ch]+year,*ZZMass,RooArgList(*cber1,*cber2,*cber3)); */
	/* 	pdf_ggzz_chan[ch] = new RooBernstein("pdfm4l_ggzz_"+chanName[ch]+year, "pdfm4l_ggzz_"+chanName[ch]+year,*ZZMass,RooArgList(*cber4,*cber5,*cber6)); */
	/* 	// if(ch!=3) */
	/* 	pdf_zx_chan [ch]= new RooGenericPdf("pdfm4l_zjets_"+chanName[ch]+year, "pdfm4l_zjets_"+chanName[ch]+year,"TMath::Landau(@0,@1,@2)",RooArgList(*ZZMass,*cland1,*cland2)); */
	/* 	// else */
	/* 	// 	pdf_zx_chan [ch]= new RooGenericPdf("pdfm4l_zjets_"+chanName[ch]+year, "pdfm4l_zjets_"+chanName[ch]+year,"TMath::Landau(@0,@1,@2)+@3*TMath::Landau(@0,@4,@5)",RooArgList(*ZZMass,*cland1,*cland2,*cfrac_land,*cland1_2,*cland2_2)); */


	/* 	RooPlot *frame = ZZMass->frame(); */
	/* 	red_qqzz->plotOn(frame,RooFit::Binning(50)); */
	/* 	pdf_qqzz_chan[ch]->plotOn(frame); */
	/* 	pdf_qqzz_chan[ch]->paramOn(frame,Layout(0.1,0.4,0.9)); */
	/* 	red_ggzz->plotOn(frame,RooFit::Binning(50)); */
	/* 	pdf_ggzz_chan[ch]->plotOn(frame); */
	/* 	pdf_ggzz_chan[ch]->paramOn(frame,Layout(0.1,0.4,0.3)); */
	/* 	red_zx->plotOn(frame,RooFit::Binning(50)); */
	/* 	pdf_zx_chan[ch]->plotOn(frame); */
	/* 	pdf_zx_chan[ch]->paramOn(frame,Layout(0.1,0.4,0.6)); */

	/* 	frame->getAttText()->SetTextSize(0.03) ; */
	/* 	frame->Draw(); */
	/* 	gPad->Print(outputDir+"/fig/fit_bkg_"+chanName[ch]+".png"); */
	/* 	gPad->Print(outputDir+"/fig/fit_bkg_"+chanName[ch]+".pdf"); */


	/* 	for (int loop=0;loop<stage1_cat;loop++){ */
	/* 		TString stage1_name = h_stage1->GetXaxis()->GetBinLabel(loop+1); */
	/* 		TString cut = Form("%sName==\"%s\" && chan==%d",stageName.Data(),stage1_name.Data(),ch+1); */
	/* 		//TString cut_weight = Form("weight*(htxs_stage1_red_catName==\"%s\" && chan==%d)",stage1_name.Data(),ch+1); */
	/* 		TString cut_weight = Form("weight*(chan==%d)",ch+1); */
	/* 		//TString cut_weight = Form("weight*(%sName==\"%s\"&&chan==%d)",stageName.Data(),stage1_name.Data(),ch+1); */
	/* 		TString histname = "template_"+stage1_name+"_"+chanName[ch]+year; */
	/* 		template2d[ch][loop] = new TH2F (histname,histname, 200,100, 4000,20,0,1); */
	/* 		tsig->Draw(Form("dbkg_kin:ZZMass>>%s",histname.Data()),cut_weight); */
	/* 		//template2d[ch][loop]->Smooth(); */
	/* 		//tsigothers->Draw(Form("dbkg_kin:ZZMass>>+%s",histname.Data()),cut_weight); */
	/* 		//template2d[ch][loop]->Smooth(); */
	/* 		Floor(template2d[ch][loop],10); */
	/* 		template2d[ch][loop]->Draw("colz"); */
	/* 		template2d[ch][loop]->SetMaximum(0.15); */
	/* 		gPad->Print(outputDir+"/fig/"+histname+"template.png"); */
	/* 		gPad->Print(outputDir+"/fig/"+histname+"template.pdf"); */
	/* 		//gPad->Print(histname+".png"); */
	/* 		//return; */
	/* 		datahist[ch][loop]=new RooDataHist(histname.ReplaceAll("template","datahist"),histname,RooArgSet(*ZZMass,*dbkg),template2d[ch][loop]); */
	/* 		pdf2d[ch][loop]=new RooHistPdf(histname.ReplaceAll("datahist","histpdf"),histname,RooArgSet(*ZZMass,*dbkg),*datahist[ch][loop]); */

	/* 		tsig->Draw(stageName,cut); */
	/* 		int stage_cat = *(tsig->GetV1()); */
	/* 		stage_cat_arr[loop] = stage_cat; */
	/* 		dataset_stage[ch][loop]=(RooDataSet*)dataset->reduce(Form("%s==%d && chan==%d",stageName.Data(),stage_cat,ch+1)); */
	/* 		TString cat_name = stage1_name+"_"+chanName[ch]+year; */

	/* 		RooFormulaVar* ca1 = new RooFormulaVar("a1_"+cat_name,"",massParaMap["a1_"+cat_name],*MH); */
	/* 		RooFormulaVar* ca2 = new RooFormulaVar("a2_"+cat_name,"",massParaMap["a2_"+cat_name],*MH); */
	/* 		RooFormulaVar* cn1 = new RooFormulaVar("n1_"+cat_name,"",massParaMap["n1_"+cat_name],*MH); */
	/* 		RooFormulaVar* cn2 = new RooFormulaVar("n2_"+cat_name,"",massParaMap["n2_"+cat_name],*MH); */
	/* 		RooFormulaVar* cmean = new RooFormulaVar("mean_"+cat_name,"",massParaMap["mean_"+cat_name],*MH); */
	/* 		RooFormulaVar* csigma=new RooFormulaVar("sigma_"+cat_name,"",massParaMap["sigma_"+cat_name],*MH); */

	/* 		int yearN = year.Atoi(); */
	/* 		RooFormulaVar* cmean_e_up = new RooFormulaVar("meaneup_"+cat_name,"",Form("(%s)*%.4f",massParaMap["mean_"+cat_name].Data(),(scale_unc_e_year[yearN-2016][ch]+1)),*MH); */
	/* 		RooFormulaVar* cmean_e_dn = new RooFormulaVar("meanedn_"+cat_name,"",Form("(%s)*%.4f",massParaMap["mean_"+cat_name].Data(),(1-scale_unc_e_year[yearN-2016][ch])),*MH); */
	/* 		RooFormulaVar* cmean_m_up = new RooFormulaVar("meanmup_"+cat_name,"",Form("(%s)*%.4f",massParaMap["mean_"+cat_name].Data(),(scale_unc_m_year[yearN-2016][ch]+1)),*MH); */
	/* 		RooFormulaVar* cmean_m_dn = new RooFormulaVar("meanmdn_"+cat_name,"",Form("(%s)*%.4f",massParaMap["mean_"+cat_name].Data(),(1-scale_unc_m_year[yearN-2016][ch])),*MH); */

	/* 		RooFormulaVar* csigma_up=new RooFormulaVar("sigmaup_"+cat_name,"",Form("(%s)*%.3f",massParaMap["sigma_"+cat_name].Data(),(res_unc_year[yearN-2016][ch]+1)),*MH); */
	/* 		RooFormulaVar* csigma_dn=new RooFormulaVar("sigmadn_"+cat_name,"",Form("(%s)*%.3f",massParaMap["sigma_"+cat_name].Data(),(1-res_unc_year[yearN-2016][ch])),*MH); */
	/* 		RooFormulaVar* cmean_l = new RooFormulaVar("landau_mean_"+cat_name,"",massParaMap["mean_l_"+cat_name],*MH); */
	/* 		RooFormulaVar* csigma_l = new RooFormulaVar("landau_sigma_"+cat_name,"",massParaMap["sigma_l_"+cat_name],*MH); */
	/* 		RooConstVar* cfrac = new RooConstVar("frac_"+cat_name,"",massParaMap["frac_"+cat_name].Atof()); */

	/* 		RooDoubleCBFast *cpdf = new RooDoubleCBFast("DCBall_"+cat_name,"Double Crystal ball",*ZZMass,*cmean,*csigma,*ca1,*cn1,*ca2,*cn2); */
	/* 		RooLandau *cpdf_landau = new RooLandau("landau_"+cat_name,"",*ZZMass, *cmean_l, *csigma_l); */
	/* 		RooAddPdf *caddpdf = new RooAddPdf("cpdf"+cat_name,"",RooArgList(*cpdf,*cpdf_landau),*cfrac); */

	/* 		RooDoubleCBFast *cpdf_mean_e_up = new RooDoubleCBFast("DCBallmean_e_up_"+cat_name,"",*ZZMass,*cmean_e_up,*csigma,*ca1,*cn1,*ca2,*cn2); */
	/* 		RooDoubleCBFast *cpdf_mean_e_dn = new RooDoubleCBFast("DCBallmean_e_dn_"+cat_name,"",*ZZMass,*cmean_e_dn,*csigma,*ca1,*cn1,*ca2,*cn2); */
	/* 		RooDoubleCBFast *cpdf_mean_m_up = new RooDoubleCBFast("DCBallmean_m_up_"+cat_name,"",*ZZMass,*cmean_m_up,*csigma,*ca1,*cn1,*ca2,*cn2); */
	/* 		RooDoubleCBFast *cpdf_mean_m_dn = new RooDoubleCBFast("DCBallmean_m_dn_"+cat_name,"",*ZZMass,*cmean_m_dn,*csigma,*ca1,*cn1,*ca2,*cn2); */

	/* 		RooDoubleCBFast *cpdf_sigmaup = new RooDoubleCBFast("DCBallsigmaup_"+cat_name,"",*ZZMass,*cmean,*csigma_up,*ca1,*cn1,*ca2,*cn2); */
	/* 		RooDoubleCBFast *cpdf_sigmadn = new RooDoubleCBFast("DCBallsigmadn_"+cat_name,"",*ZZMass,*cmean,*csigma_dn,*ca1,*cn1,*ca2,*cn2); */
	/* 		RooAddPdf *caddpdf_mean_e_up = new RooAddPdf("cpdfmean_e_up"+cat_name,"",RooArgList(*cpdf_mean_e_up,*cpdf_landau),*cfrac); */
	/* 		RooAddPdf *caddpdf_mean_e_dn = new RooAddPdf("cpdfmean_e_dn"+cat_name,"",RooArgList(*cpdf_mean_e_dn,*cpdf_landau),*cfrac); */
	/* 		RooAddPdf *caddpdf_mean_m_up = new RooAddPdf("cpdfmean_m_up"+cat_name,"",RooArgList(*cpdf_mean_m_up,*cpdf_landau),*cfrac); */
	/* 		RooAddPdf *caddpdf_mean_m_dn = new RooAddPdf("cpdfmean_m_dn"+cat_name,"",RooArgList(*cpdf_mean_m_dn,*cpdf_landau),*cfrac); */
	/* 		RooAddPdf *caddpdf_sigmaup = new RooAddPdf("cpdfsigmaup"+cat_name,"",RooArgList(*cpdf_sigmaup,*cpdf_landau),*cfrac); */
	/* 		RooAddPdf *caddpdf_sigmadn = new RooAddPdf("cpdfsigmadn"+cat_name,"",RooArgList(*cpdf_sigmadn,*cpdf_landau),*cfrac); */


	/* 		if (stage1_name.Contains("VH_Lep")||stage1_name.Contains("ZH_Lep")||stage1_name.Contains("WH_Lep")|| stage1_name.Contains("TTH")){ */
	/* 			pdf_stage[ch][loop]=caddpdf; */
	/* 			pdf_stage_mean_m_up[ch][loop]=caddpdf_mean_m_up; */
	/* 			pdf_stage_mean_m_dn[ch][loop]=caddpdf_mean_m_dn; */
	/* 			pdf_stage_mean_e_up[ch][loop]=caddpdf_mean_e_up; */
	/* 			pdf_stage_mean_e_dn[ch][loop]=caddpdf_mean_e_dn; */
	/* 			pdf_stage_sigmaup[ch][loop]=caddpdf_sigmaup; */
	/* 			pdf_stage_sigmadn[ch][loop]=caddpdf_sigmadn; */
	/* 		} */
	/* 		else{ */
	/* 			pdf_stage[ch][loop]=cpdf; */
	/* 			pdf_stage_mean_m_up[ch][loop]=cpdf_mean_m_up; */
	/* 			pdf_stage_mean_m_dn[ch][loop]=cpdf_mean_m_dn; */
	/* 			pdf_stage_mean_e_up[ch][loop]=cpdf_mean_e_up; */
	/* 			pdf_stage_mean_e_dn[ch][loop]=cpdf_mean_e_dn; */
	/* 			pdf_stage_sigmaup[ch][loop]=cpdf_sigmaup; */
	/* 			pdf_stage_sigmadn[ch][loop]=cpdf_sigmadn; */
	/* 		} */

	/* 		pdf_stage[ch][loop]->SetName("pdfm4l_"+cat_name); */
	/* 		pdf_stage[ch][loop]->SetTitle("pdfm4l_"+cat_name); */

	/* 		RooPlot *frame = ZZMass->frame(); */
	/* 		cout<< dataset_stage[ch][loop]->sumEntries()<<endl; */
	/* 		dataset_stage[ch][loop]->plotOn(frame); */
	/* 		ca1->Print("v"); */
	/* 		ca2->Print("v"); */
	/* 		cn1->Print("v"); */
	/* 		cn2->Print("v"); */
	/* 		cmean->Print("v"); */
	/* 		csigma->Print("v"); */
	/* 		//			pdf_stage[ch][loop]->Print("v"); */

	/* 		MH->setVal(125); */
	/* 		//			MH->setConstant(); */
	/* 		pdf_stage[ch][loop]->plotOn(frame); */
	/* 		//pdfcat->paramOn(frame,Layout(0.1,0.4)); */
	/* 		//frame->getAttText()->SetTextSize(0.03) ; */
	/* 		frame->Draw(); */
	/* 		gPad->Print(outputDir+"/fig/fit_stage_"+cat_name+".png"); */
	/* 		gPad->Print(outputDir+"/fig/fit_stage_"+cat_name+".pdf"); */

	/* 	} */
        /* } */



	/* for (int i = 0; i < reco_cat_all; i++){ */
	/* 	//for (int i = 0; i < 2; i++){ */
	/* 	TString reco_name = h->GetXaxis()->GetBinLabel(i+1); */

	/* 	cout<<reco_name<<endl; */
	/* 	RooProdPdf *pdfprod[3][stage1_cat]; */
	/* 	RooProdPdf *pdfprod_mean_e_up[3][stage1_cat]; */
	/* 	RooProdPdf *pdfprod_mean_e_dn[3][stage1_cat]; */
	/* 	RooProdPdf *pdfprod_mean_m_up[3][stage1_cat]; */
	/* 	RooProdPdf *pdfprod_mean_m_dn[3][stage1_cat]; */
	/* 	RooProdPdf *pdfprod_sigmaup[3][stage1_cat]; */
	/* 	RooProdPdf *pdfprod_sigmadn[3][stage1_cat]; */
	/* 	for (int k=flav;k<flav+1;k++){ */

	/* 		RooRealVar* res_unc=new RooRealVar(res_ch,"",0,-7,7); */
	/* 		RooRealVar* scale_unc_e=new RooRealVar(scale_e,"",0,-7,7); */
	/* 		RooRealVar* scale_unc_m=new RooRealVar(scale_m,"",0,-7,7); */
	/* 		ofstream card; */
	/* 		TString cat_name = reco_name + "_"+chanName[k]; */
	/* 		card.open (outputDir+"/"+cat_name+year+".txt"); */
	/* 		card <<"Datacard for event category: "<< cat_name<<endl; */
	/* 		card<< "imax 1 number of bins"<<endl; */
	/* 		//				card<< "jmax " <<stage1_cat + 2<< " number of processes minus 1"<<endl; */
	/* 		card<< "jmax * number of processes minus 1"<<endl; */
	/* 		//			card<< "kmax " << nsys+3 <<" number of nuisance parameters"<<endl; */
	/* 		card<< "kmax * number of nuisance parameters"<<endl; */
	/* 		card<<"---------------------------------"<<endl; */
	/* 		card<<endl; */
	/* 		card<< "shapes * "<< cat_name << " "<< cat_name+year+".root w:$PROCESS w:$PROCESS_$SYSTEMATIC" <<endl; */
	/* 		card<<"---------------------------------"<<endl; */
	/* 		card<< "bin "<< cat_name <<endl; */
	/* 		card<< "observation -1 "<<endl; */
	/* 		card<<"---------------------------------"<<endl; */
	/* 		card<<endl; */

	/* 		vector<TString> bin_arr; */
	/* 		vector<TString> processName_arr; */
	/* 		vector<int> process_arr; */
	/* 		vector<float> yield_arr; */
	/* 		//			vector <vector<float>> sys_arr_up; */
	/* 		//			vector <vector<float>> sys_arr_dn; */

	/* 		TFile *f=new TFile(outputDir+"/"+cat_name+year+".root","recreate"); */
	/* 		RooWorkspace w("w"); */
	/* 		int reco_cat=-1; */

	/* 		TString varname ="dbkg_kin"; */
	/* 		int nbin = 20; */
	/* 		float low=0; */
	/* 		float high=1; */
	/* 		Ddis = dbkg; */

	/* 		if(reco_name.Contains("VBF_1j")){ */
	/* 			varname = "D1jet"; */
	/* 			Ddis = D1jet; */
	/* 		} */

	/* 		else if(reco_name.Contains("VBF_ptj_GT200_1J")){ */
	/* 			varname = "D1jet"; */
	/* 			Ddis = D1jet_norm; */
	/* 		} */
	/* 		else if (reco_name.Contains("VBF_2j")|| reco_name.Contains("VBF_VBF")|| reco_name.Contains("VBF_ptj_GT200_2J")){ */
	/* 			//if (reco_name.Contains("VBF_2j")|| reco_name.Contains("VBF_VBF")|| reco_name.Contains("VBF_ptj_GT200_2J")|| reco_name.Contains("ggH_2j_")){ */
	/* 			varname = "DVBFDEC"; */
	/* 			Ddis = DVBFDEC; */
	/* 		} */
	/* 		else if (reco_name.Contains("VH_Had")){ */
	/* 			varname = "DVHDEC"; */
	/* 			Ddis = DVHDEC; */
	/* 		} */

	/* 		low  = Ddis->getMin(); */
	/* 		high = Ddis->getMax(); */
	/* 		nbin = Ddis->getBins(); */
	/* 		cout<<"NBins "<<nbin<<reco_name<<endl; */
	/* 		//TString cutwtf = Form("%s==\"%s\"",categoryName.Data(),reco_name.Data()); */
	/* 		//cout<<cutwtf<<endl; */
	/* 		//tsig->Draw(Form("%s>>%s",category.Data(),cat_name.Data()),cutwtf); */
	/* 		//reco_cat = *(tsig->GetV1()); */
	/* 		reco_cat = reco_cat_num.at(i); */
	/* 		vector <map<TString,float>> sys_arr_up; */
	/* 		vector <map<TString,float>> sys_arr_dn; */
	/* 		for (int j = 0; j < stage1_cat; j++){ */
	/* 			TString stage1_name = h_stage1->GetXaxis()->GetBinLabel(j+1); */
	/* 			cout<<stage1_name<<endl; */
	/* 			TString histname = "template_"+cat_name+"_"+stage1_name+year; */
	/* 			TH2F *template2d_tmp = new TH2F (histname,histname, 140,105, 140,nbin,low,high); */
	/* 			TString cut = Form("weight*(%sName==\"%s\" && %s==\"%s\" && chan==%d)",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 			cout<<cut<<endl; */
	/* 			TString cut_stage = Form("weight*(%sName==\"%s\" && chan==%d)",stageName.Data(),stage1_name.Data(), k+1); */
	/* 			cout<<cut_stage<<endl; */
	/* 			TTree *subsig = tsig->CopyTree(cut); */
	/* 			int cutentry = subsig->Draw(Form("%s:ZZMass>>%s",varname.Data(),histname.Data()),cut); */
	/* 			//				tsig->Draw(Form("%s:%s",stageName.Data(),category.Data()),cut); */
	/* 			float yield = template2d_tmp->Integral(); */
	/* 			if(yield<0.005) */
	/* 				continue; */
	/* 			bin_arr.push_back(cat_name); */
	/* 			processName_arr.push_back(stage1_name); */
	/* 			process_arr.push_back ( 0-stage1_cat+j); */
	/* 			yield_arr.push_back(yield); */
	/* 			int stage_cat=stage_cat_arr[j]; */
	/* 			cout<< stage_cat<<"\t"<<reco_cat<<"\t"<<yield<<endl; */
	/* 			TString cut_num = Form("%s==%d && %s==%d && chan==%d",stageName.Data(),stage_cat,category.Data(),reco_cat, k+1); */
	/* 			//TString cut_num_red = Form("htxs_stage1_red_cat==%d && chan==%d",stage_cat, k+1); */

	/* 			TTree *subsig_stage = tsig->CopyTree(cut_stage); */
	/* 			//vector <float> sys_arr_stage_up ; */
	/* 			//vector <float> sys_arr_stage_dn ; */
	/* 			map<TString, float> sys_arr_stage_up; */
	/* 			map<TString, float> sys_arr_stage_dn; */

	/* 			if(doSys){ */
	/* 				if( cutentry !=0){ */
	/* 					// TString cut_jec_up = Form("weight*(%sName==\"%s\" && %s_jec_up==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					// TString cut_jec_dn = Form("weight*(%sName==\"%s\" && %s_jec_dn==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					// TString cut_jer_up = Form("weight*(%sName==\"%s\" && %s_jer_up==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					// TString cut_jer_dn = Form("weight*(%sName==\"%s\" && %s_jer_dn==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */

	/* 					TString cut_jes_up = Form("weight*(%sName==\"%s\" && %s_jes_up==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					TString cut_jes_dn = Form("weight*(%sName==\"%s\" && %s_jes_dn==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					TString cut_jer_up = Form("weight*(%sName==\"%s\" && %s_jer_up==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					TString cut_jer_dn = Form("weight*(%sName==\"%s\" && %s_jer_dn==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */

	/* 					// TString cut_jetPt_jes_up = Form("weight*(%sName==\"%s\" && %s_jetPt_jes_up==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					// TString cut_jetPt_jes_dn = Form("weight*(%sName==\"%s\" && %s_jetPt_jes_dn==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					// TString cut_jetPt_jer_up = Form("weight*(%sName==\"%s\" && %s_jetPt_jer_up==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					// TString cut_jetPt_jer_dn = Form("weight*(%sName==\"%s\" && %s_jetPt_jer_dn==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */

	/* 					TString cut_btag_up = Form("weight*(%sName==\"%s\" && %s_btag_up==\"%s\" && chan==%d )",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */
	/* 					TString cut_btag_dn = Form("weight*(%sName==\"%s\" && %s_btag_dn==\"%s\" && chan==%d )", stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1); */

	/* 					TH1F *jet_sys =new TH1F("jet_sys","",1,105, 140); */
	/* 					subsig_stage->Draw("ZZMass>>jet_sys",cut_jes_up); */
	/* 					float yield_jes_up = jet_sys->Integral(); */
	/* 					subsig_stage->Draw("ZZMass>>jet_sys",cut_jes_dn); */
	/* 					float yield_jes_dn = jet_sys->Integral(); */
	/* 					subsig_stage->Draw("ZZMass>>jet_sys",cut_jer_up); */
	/* 					float yield_jer_up = jet_sys->Integral(); */
	/* 					subsig_stage->Draw("ZZMass>>jet_sys",cut_jer_dn); */
	/* 					float yield_jer_dn = jet_sys->Integral(); */
	/* 					// subsig_stage->Draw("ZZMass>>jet_sys",cut_jetPt_jes_up); */
	/* 					// float yield_jetPt_jes_up = jet_sys->Integral(); */
	/* 					// subsig_stage->Draw("ZZMass>>jet_sys",cut_jetPt_jes_dn); */
	/* 					// float yield_jetPt_jes_dn = jet_sys->Integral(); */
	/* 					// subsig_stage->Draw("ZZMass>>jet_sys",cut_jetPt_jer_up); */
	/* 					// float yield_jetPt_jer_up = jet_sys->Integral(); */
	/* 					// subsig_stage->Draw("ZZMass>>jet_sys",cut_jetPt_jer_dn); */
	/* 					// float yield_jetPt_jer_dn = jet_sys->Integral(); */
	/* 					subsig_stage->Draw("ZZMass>>jet_sys",cut_btag_up); */
	/* 					float yield_btag_up = jet_sys->Integral(); */
	/* 					subsig_stage->Draw("ZZMass>>jet_sys",cut_btag_dn); */
	/* 					float yield_btag_dn = jet_sys->Integral(); */

	/* 					sys_arr_stage_up["jes"]=yield_jes_up/yield; */
	/* 					sys_arr_stage_dn["jes"]=yield_jes_dn/yield; */
	/* 					sys_arr_stage_up["jer"]=yield_jer_up/yield; */
	/* 					sys_arr_stage_dn["jer"]=yield_jer_dn/yield; */
	/* 					// sys_arr_stage_up["jetPt_jes"]=yield_jetPt_jes_up/yield; */
	/* 					// sys_arr_stage_dn["jetPt_jes"]=yield_jetPt_jes_dn/yield; */
	/* 					// sys_arr_stage_up["jetPt_jer"]=yield_jetPt_jer_up/yield; */
	/* 					// sys_arr_stage_dn["jetPt_jer"]=yield_jetPt_jer_dn/yield; */
	/* 					sys_arr_stage_up["hzz_br"]=1.02; */
	/* 					sys_arr_stage_dn["hzz_br"]=0.98; */
	/* 					sys_arr_stage_up["btag"]=yield_btag_up/yield; */
	/* 					sys_arr_stage_dn["btag"]=yield_btag_dn/yield; */
	/* 					//						TString res="CMS_res_"+chanName[k]; */
	/* 					sys_arr_stage_up[res_ch]=1; */
	/* 					sys_arr_stage_dn[res_ch]=1; */

	/* 					TString lookup_tune= reco_name+"_"+stage1_name+"_total_125_tune"; */
						
	/* 						cout<<lookup_tune+"up"<<"\t"<< sysParaMap[lookup_tune+"up"]<<endl; */
	/* 					if(sysParaMap.find(lookup_tune+"up") != sysParaMap.end()){ */
	/* 						sys_arr_stage_up["pythia_tune"]=sysParaMap[lookup_tune+"up"]; */
	/* 						sys_arr_stage_dn["pythia_tune"]=sysParaMap[lookup_tune+"down"]; */
	/* 					} */

	/* 					if(k==0 || k==2){ */
	/* 						sys_arr_stage_up[scale_m]=1; */
	/* 						sys_arr_stage_dn[scale_m]=1; */
	/* 					} */
	/* 					if (k==1||k==2){ */
	/* 						sys_arr_stage_up[scale_e]=1; */
	/* 						sys_arr_stage_dn[scale_e]=1; */
	/* 					} */
	/* 					cout<<"jes "<<sys_arr_stage_up["jes"]<<"\t"<< sys_arr_stage_dn["jes"]<<endl; */
	/* 					cout<<"jer "<<sys_arr_stage_up["jer"]<<"\t"<< sys_arr_stage_dn["jer"]<<endl; */
	/* 					// cout<<"jetPt jes "<<sys_arr_stage_up["jetPt_jes"]<<"\t"<< sys_arr_stage_dn["jetPt_jes"]<<endl; */
	/* 					// cout<<"jetPt jer "<<sys_arr_stage_up["jetPt_jer"]<<"\t"<< sys_arr_stage_dn["jetPt_jer"]<<endl; */
	/* 					cout<<"btag "<<sys_arr_stage_up["btag"]<<"\t"<< sys_arr_stage_dn["btag"]<<endl; */
	/* 					for (int nextra = 0; nextra < extraSysn_all; nextra++){ */
	/* 						sys_arr_stage_up[extraSys_handName_all[nextra]]=extraSys_all_up[nextra][k]; */
	/* 						sys_arr_stage_dn[extraSys_handName_all[nextra]]=extraSys_all_dn[nextra][k]; */
	/* 					} */


	/* 					TString syshist_name ="sys_hist_"+stage1_name; */
	/* 					subsig->Draw("weight_name>>"+syshist_name,cut); */

	/* 					TH1F *systemp = (TH1F*)gROOT->FindObject(syshist_name); */
	/* 					int nbins_sys= systemp->GetNbinsX(); */


			/* 			for (int sysl = 0; sysl< nbins_sys;sysl++){ */
			/* 				TString sysname = systemp->GetXaxis()->GetBinLabel(sysl+1); */
			/* 				TString sysname_cat = sysname +"_cat"; */
			/* 				cout<<sysname_cat<<endl; */
			/* 				bool qcd_ggh_sys= false; */
			/* 				if(sysname=="QCDscale_muR_ggH"||sysname=="QCDscale_muF_ggH") */
			/* 					qcd_ggh_sys=true; */
			/* 				if(yield==0){ */
			/* 					if(!qcd_ggh_sys){ */
			/* 						sys_arr_stage_up[sysname]=0; */
			/* 						sys_arr_stage_dn[sysname]=0; */
			/* 					} */
			/* 					sys_arr_stage_up[sysname_cat]=0; */
			/* 					sys_arr_stage_dn[sysname_cat]=0; */
			/* 				} */
			/* 				else{ */
			/* 					TString hsys_tmpName_up = "hsys_tmp_"+histname+"_"+sysname+"_up"; */
			/* 					TString hsys_tmpName_dn = "hsys_tmp_"+histname+"_"+sysname+"_dn"; */
			/* 					TString hsys_tmpName_norm = "hsys_tmp_"+histname+"_"+sysname+"_norm"; */
			/* 					TH1F *hsys_tmp_up = new TH1F (hsys_tmpName_up,hsys_tmpName_up,1,105, 140); */
			/* 					TH1F *hsys_tmp_dn = new TH1F (hsys_tmpName_dn,hsys_tmpName_dn,1,105, 140); */
			/* 					TH1F *hsys_tmp_norm = new TH1F (hsys_tmpName_norm,hsys_tmpName_norm,1,105, 140); */
			/* 					TString cut_norm = Form("weight*(%sName==\"%s\" && %s==\"%s\" && chan==%d && weight_name==\"%s\")",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1,sysname.Data()); */
			/* 					TString cut_up = Form("weight_up*(%sName==\"%s\" && %s==\"%s\" && chan==%d && weight_name==\"%s\")",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1,sysname.Data()); */
			/* 					TString cut_dn = Form("weight_dn*(%sName==\"%s\" && %s==\"%s\" && chan==%d && weight_name==\"%s\")",stageName.Data(),stage1_name.Data(),categoryName.Data(),reco_name.Data(), k+1,sysname.Data()); */
			/* 					subsig->Draw("ZZMass>>"+hsys_tmpName_up,cut_up); */
			/* 					subsig->Draw("ZZMass>>"+hsys_tmpName_dn,cut_dn); */
			/* 					subsig->Draw("ZZMass>>"+hsys_tmpName_norm,cut_norm); */

			/* 					cout<<"up tree "<<hsys_tmp_up->GetEntries()<<"\t"<< sysname<<endl; */
			/* 					float yield_sys_up = hsys_tmp_up->Integral(); */
			/* 					float yield_sys_dn = hsys_tmp_dn->Integral(); */
			/* 					float yield_sys_norm = hsys_tmp_norm->Integral(); */

			/* 					TString allcut_norm = Form("weight*(weight_name==\"%s\")",sysname.Data()); */
			/* 					TString allcut_up = Form("weight_up*(weight_name==\"%s\")",sysname.Data()); */
			/* 					TString allcut_dn = Form("weight_dn*(weight_name==\"%s\")",sysname.Data()); */
			/* 					subsig_stage->Draw("ZZMass>>"+hsys_tmpName_up,allcut_up); */
			/* 					subsig_stage->Draw("ZZMass>>"+hsys_tmpName_dn,allcut_dn); */
			/* 					subsig_stage->Draw("ZZMass>>"+hsys_tmpName_norm,allcut_norm); */
			/* 					float allyield_sys_up = hsys_tmp_up->Integral(); */
			/* 					float allyield_sys_dn = hsys_tmp_dn->Integral(); */
			/* 					float allyield_sys_norm = hsys_tmp_norm->Integral(); */

			/* 					// We don't have scale up/do anymore */
			/* 					// Always use yields scale variation */
			/* 					if(sysname.Contains("pythia")){ */
			/* 						//if(year!="2016"){ */
			/* 						sys_arr_stage_up[sysname]=(yield-yield_sys_norm+yield_sys_up)/yield; */
			/* 						sys_arr_stage_dn[sysname]=(yield-yield_sys_norm+yield_sys_dn)/yield; */
			/* 						//} */
			/* 						//else{ */
			/* 						//TString lookup= reco_name+"_"+stage1_name+"_total_125_scale"; */
			/* 						//cout<< lookup<<"\t"<<sysParaMap[lookup+"up"]<<endl; */
			/* 						//sys_arr_stage_up[sysname]=sysParaMap[lookup+"up"]; */
			/* 						//sys_arr_stage_dn[sysname]=sysParaMap[lookup+"down"]; */
			/* 						//} */
			/* 					} */
			/* 					else{ */
			/* 						if(!qcd_ggh_sys){ */
			/* 							sys_arr_stage_up[sysname]=(allyield_sys_up)/allyield_sys_norm; */
			/* 							sys_arr_stage_dn[sysname]=(allyield_sys_dn)/allyield_sys_norm; */
			/* 						} */
			/* 						sys_arr_stage_up[sysname_cat]=(yield-yield_sys_norm+yield_sys_up)/yield/allyield_sys_up*allyield_sys_norm; */
			/* 						sys_arr_stage_dn[sysname_cat]=(yield-yield_sys_norm+yield_sys_dn)/yield/allyield_sys_dn*allyield_sys_norm; */
			/* 					} */

			/* 					cout<< "up and down  "<<sys_arr_stage_up[sysname_cat]<<"\t"<<sys_arr_stage_dn[sysname_cat]<<"\t"<<yield<<endl; */
			/* 				} */
			/* 			} */
			/* 		} */
			/* 		sys_arr_up.push_back(sys_arr_stage_up); */
			/* 		sys_arr_dn.push_back(sys_arr_stage_dn); */
			/* 		cout<< "sys "<<sys_arr_dn.size()<<"\t"<< sys_arr_stage_dn.size()<<endl; */
			/* 	} */
			/* 	//				for (int test =0;test<sys_arr_dn.size();) */

			/* 	RooDataHist *datahist_tmp; */
			/* 	RooHistPdf *pdf2d_tmp; */
			/* 	if(reco_name.Contains("VH_Had")||reco_name.Contains("WH_Had")||reco_name.Contains("ZH_Had")|| reco_name.Contains("VBF_1j")|| reco_name.Contains("VBF_VBF")||reco_name.Contains("VBF_2j")||reco_name.Contains("VBF_ptj_GT200_1J")||reco_name.Contains("VBF_ptj_GT200_2J") ){ */
			/* 		//					template2d_tmp->Smooth(); */
			/* 		Floor(template2d_tmp,20); */
			/* 		template2d_tmp->Draw("colz"); */
			/* 		template2d_tmp->SetMaximum(0.2); */
			/* 		gPad->Print(outputDir+"/fig/"+histname+"template.png"); */
			/* 		gPad->Print(outputDir+"/fig/"+histname+"template.pdf"); */
			/* 		gPad->Print(); */
			/* 		datahist_tmp=new RooDataHist(histname.ReplaceAll("template","datahist"),histname,RooArgSet(*ZZMass,*Ddis),template2d_tmp); */
			/* 		pdf2d_tmp=new RooHistPdf(histname.ReplaceAll("datahist","histpdf"),histname,RooArgSet(*ZZMass,*Ddis),*datahist_tmp); */
			/* 	} */

			/* 	RooDataSet *red = (RooDataSet*)dataset->reduce(cut_num); */
			/* 	cout<<red->sumEntries()<<endl; */
			/* 	RooAbsPdf *pdfcat ; */
			/* 	RooAbsPdf *pdfcat_mean_e_up; */
			/* 	RooAbsPdf *pdfcat_mean_e_dn; */
			/* 	RooAbsPdf *pdfcat_mean_m_up; */
			/* 	RooAbsPdf *pdfcat_mean_m_dn; */
			/* 	RooAbsPdf *pdfcat_sigmaup; */
			/* 	RooAbsPdf *pdfcat_sigmadn; */

			/* 	pdfcat = pdf_stage[k][j]; */
			/* 	pdfcat_mean_e_up = pdf_stage_mean_e_up[k][j]; */
			/* 	pdfcat_mean_e_dn = pdf_stage_mean_e_dn[k][j]; */
			/* 	pdfcat_mean_m_up = pdf_stage_mean_m_up[k][j]; */
			/* 	pdfcat_mean_m_dn = pdf_stage_mean_m_dn[k][j]; */
			/* 	pdfcat_sigmaup = pdf_stage_sigmaup[k][j]; */
			/* 	pdfcat_sigmadn = pdf_stage_sigmadn[k][j]; */

			/* 	pdfcat->SetName("pdfm4l_"+cat_name+"_"+stage1_name+year); */
			/* 	pdfcat->SetTitle("pdfm4l_"+cat_name+"_"+stage1_name+year); */
			/* 	pdfcat_mean_e_up->SetName("pdfm4l_mean_e_up_"+cat_name+"_"+stage1_name+year); */
			/* 	pdfcat_mean_e_dn->SetName("pdfm4l_mean_e_dn_"+cat_name+"_"+stage1_name+year); */
			/* 	pdfcat_mean_m_up->SetName("pdfm4l_mean_m_up_"+cat_name+"_"+stage1_name+year); */
			/* 	pdfcat_mean_m_dn->SetName("pdfm4l_mean_m_dn_"+cat_name+"_"+stage1_name+year); */
			/* 	pdfcat_sigmaup->SetName("pdfm4l_sigmaup_"+cat_name+"_"+stage1_name+year); */
			/* 	pdfcat_sigmadn->SetName("pdfm4l_sigmadn_"+cat_name+"_"+stage1_name+year); */

			/* 	//VerticalInterpPdf *pdfMorph = new VerticalInterpPdf(stage1_name+"morph",stage1_name+"morph",RooArgList(*pdfcat,*pdfcat_meanup,*pdfcat_meandn,*pdfcat_sigmaup,*pdfcat_sigmadn),RooArgList(*scale_unc,*res_unc)); */

			/* 	if(cat_name.Contains ("VH_Had")||cat_name.Contains ("WH_Had")||cat_name.Contains ("ZH_Had")|| cat_name.Contains("VBF_1j")|| cat_name.Contains ("VBF_VBF")|| cat_name.Contains("VBF_2j")||cat_name.Contains("VBF_ptj_GT200_1J")||cat_name.Contains("VBF_ptj_GT200_2J")){ */

			/* 		//if(cat_name.Contains ("VH_Had")|| cat_name.Contains ("VBF_VBF")|| cat_name.Contains("VBF_2j")||cat_name.Contains("VBF_ptj_GT200_1J")||cat_name.Contains("VBF_ptj_GT200_2J")) */
			/* 		//					pdfprod[k][j]=new RooProdPdf (stage1_name,stage1_name,*pdfMorph,Conditional(*pdf2d_tmp,*Ddis)); */
			/* 		pdfprod[k][j]=new RooProdPdf (stage1_name,stage1_name,*pdfcat,Conditional(*pdf2d_tmp,*Ddis)); */
			/* 		pdfprod_mean_e_up[k][j]=new RooProdPdf (stage1_name+"_"+scale_e+"Up",stage1_name+"_"+scale_e+"Up",*pdfcat_mean_e_up,Conditional(*pdf2d_tmp,*Ddis)); */
			/* 		pdfprod_mean_e_dn[k][j]=new RooProdPdf (stage1_name+"_"+scale_e+"Down",stage1_name+"_"+scale_e+"Down",*pdfcat_mean_e_dn,Conditional(*pdf2d_tmp,*Ddis)); */
			/* 		pdfprod_mean_m_up[k][j]=new RooProdPdf (stage1_name+"_"+scale_m+"Up",stage1_name+"_"+scale_m+"Up",*pdfcat_mean_m_up,Conditional(*pdf2d_tmp,*Ddis)); */
			/* 		pdfprod_mean_m_dn[k][j]=new RooProdPdf (stage1_name+"_"+scale_m+"Down",stage1_name+"_"+scale_m+"Down",*pdfcat_mean_m_dn,Conditional(*pdf2d_tmp,*Ddis)); */
			/* 		pdfprod_sigmaup[k][j]=new RooProdPdf (stage1_name+"_"+res_ch+"Up",stage1_name+"_"+res_ch+"Up",*pdfcat_sigmaup,Conditional(*pdf2d_tmp,*Ddis)); */
			/* 		pdfprod_sigmadn[k][j]=new RooProdPdf (stage1_name+"_"+res_ch+"Down",stage1_name+"_"+res_ch+"Down",*pdfcat_sigmadn,Conditional(*pdf2d_tmp,*Ddis)); */

			/* 	} */
			/* 	else{ */
			/* 		//pdfprod[k][j]=new RooProdPdf (stage1_name,stage1_name,*pdfMorph,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		pdfprod[k][j]=new RooProdPdf (stage1_name,stage1_name,*pdfcat,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		pdfprod_mean_e_up[k][j]=new RooProdPdf (stage1_name+"_"+scale_e+"Up",stage1_name+"_"+scale_e+"Up",*pdfcat_mean_e_up,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		pdfprod_mean_e_dn[k][j]=new RooProdPdf (stage1_name+"_"+scale_e+"Down",stage1_name+"_"+scale_e+"Down",*pdfcat_mean_e_dn,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		pdfprod_mean_m_up[k][j]=new RooProdPdf (stage1_name+"_"+scale_m+"Up",stage1_name+"_"+scale_m+"Up",*pdfcat_mean_m_up,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		pdfprod_mean_m_dn[k][j]=new RooProdPdf (stage1_name+"_"+scale_m+"Down",stage1_name+"_"+scale_m+"Down",*pdfcat_mean_m_dn,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		pdfprod_sigmaup[k][j]=new RooProdPdf (stage1_name+"_"+res_ch+"Up",stage1_name+"_"+res_ch+"Up",*pdfcat_sigmaup,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		pdfprod_sigmadn[k][j]=new RooProdPdf (stage1_name+"_"+res_ch+"Down",stage1_name+"_"+res_ch+"Down",*pdfcat_sigmadn,Conditional(*pdf2d[k][j],*Ddis)); */

			/* 		//pdfprod_mean_e_up[k][j]=new RooProdPdf (stage1_name+"_CMS_scale_e_"+year+"Up",stage1_name+"_CMS_scale_e_"+year+"Up",*pdfcat_mean_e_up,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		//pdfprod_mean_e_dn[k][j]=new RooProdPdf (stage1_name+"_CMS_scale_e_"+year+"Down",stage1_name+"_CMS_scale_e_"+year+"Down",*pdfcat_mean_e_dn,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		//pdfprod_mean_m_up[k][j]=new RooProdPdf (stage1_name+"_CMS_scale_m_"+year+"Up",stage1_name+"_CMS_scale_m_"+year+"Up",*pdfcat_mean_m_up,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		//pdfprod_mean_m_dn[k][j]=new RooProdPdf (stage1_name+"_CMS_scale_m_"+year+"Down",stage1_name+"_CMS_scale_m_"+year+"Down",*pdfcat_mean_m_dn,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		//pdfprod_sigmaup[k][j]=new RooProdPdf (stage1_name+"_CMS_res_"+chanName[k]+"_"+year+"Up",stage1_name+"_Res"+chanName[k]+year+"Up",*pdfcat_sigmaup,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 		//pdfprod_sigmadn[k][j]=new RooProdPdf (stage1_name+"_CMS_res_"+chanName[k]+"_"+year+"Down",stage1_name+"_Res"+chanName[k]+year+"Down",*pdfcat_sigmadn,Conditional(*pdf2d[k][j],*Ddis)); */
			/* 	} */
			/* 	//VerticalInterpPdf *pdfMorph = new VerticalInterpPdf(stage1_name,stage1_name,*ZZMass,RooArgList(*pdfprod[k][j],*pdfprod_meanup[k][j],*pdfprod_meandn[k][j],*pdfprod_sigmaup[k][j],*pdfprod_sigmadn[k][j]),RooArgList(*scale_unc,*res_unc)); */
			/* 	//w.import(*pdfMorph,RecycleConflictNodes()); */
			/* 	w.import(*pdfprod[k][j],RecycleConflictNodes()); */
			/* 	if(k==0 || k==2){ */
			/* 		w.import(*pdfprod_mean_m_up[k][j],RecycleConflictNodes()); */
			/* 		w.import(*pdfprod_mean_m_dn[k][j],RecycleConflictNodes()); */
			/* 	} */
			/* 	if (k==1 || k==2){ */
			/* 		w.import(*pdfprod_mean_e_up[k][j],RecycleConflictNodes()); */
			/* 		w.import(*pdfprod_mean_e_dn[k][j],RecycleConflictNodes()); */
			/* 	} */
			/* 	w.import(*pdfprod_sigmaup[k][j],RecycleConflictNodes()); */
			/* 	w.import(*pdfprod_sigmadn[k][j],RecycleConflictNodes()); */
			/* 	TString yieldName = reco_name+"_"+stage1_name+"_"+chanName[k]; */
			/* 	if ( yieldParaMap.find(yieldName) != yieldParaMap.end() ){ */
			/* 		RooFormulaVar *stage_norm = new RooFormulaVar(stage1_name+"_norm","",yieldParaMap[yieldName],*MH); */
			/* 		w.import(*stage_norm,RecycleConflictNodes()); */
			/* 	} */
			/* } */

			/* TString cut_reco= Form("%s==%d && chan==%d",category.Data(),reco_cat, k+1); */

			/* RooDataSet *red_qqzz = (RooDataSet*)dataset_qqzz->reduce(cut_reco); */
			/* RooDataSet *red_ggzz = (RooDataSet*)dataset_ggzz->reduce(cut_reco); */
			/* RooDataSet *red_zx = (RooDataSet*)dataset_zx->reduce(cut_reco); */

			/* RooBernstein *pdf_qqzz_cat = new RooBernstein("pdfm4l_qqzz_"+cat_name, "pdfm4l_qqzz_"+cat_name,*ZZMass,RooArgList(*ber1,*ber2,*ber3)); */
			/* RooBernstein *pdf_ggzz_cat = new RooBernstein("pdfm4l_ggzz_"+cat_name, "pdfm4l_ggzz_"+cat_name,*ZZMass,RooArgList(*ber4,*ber5,*ber6)); */
			/* RooGenericPdf *pdf_zx_cat ; */
			/* // if (k!=2) */
			/* pdf_zx_cat = new RooGenericPdf("pdfm4l_zjets_"+cat_name, "pdfm4l_zjets_"+cat_name,"TMath::Landau(@0,@1,@2)",RooArgList(*ZZMass,*land1,*land2)); */
			/* // else */
			/* // 	pdf_zx_cat = new RooGenericPdf("pdfm4l_zjets_"+cat_name, "pdfm4l_zjets_"+cat_name,"TMath::Landau(@0,@1,@2)+@3*TMath::Landau(@0,@4,@5)",RooArgList(*ZZMass,*land1,*land2,*frac_land,*land1_2,*land2_2)); */

			/* TH2F *template2d_zx_tmp = new TH2F("template_zx"+cat_name,"", 140,105, 140,nbin,low,high); */
			/* TH2F *template2d_qqzz_tmp = new TH2F("template_qqzz"+cat_name,"", 140,105, 140,nbin,low,high); */
			/* TH2F *template2d_ggzz_tmp = new TH2F("template_ggzz"+cat_name,"", 140,105, 140,nbin,low,high); */
			/* TH2F *template2d_zx_tmp_all = new TH2F("template_zx"+cat_name+"_all","", 140,105, 140,nbin,low,high); */
			/* TH2F *template2d_qqzz_tmp_all = new TH2F("template_qqzz"+cat_name+"_all","", 140,105, 140,nbin,low,high); */
			/* TH2F *template2d_ggzz_tmp_all = new TH2F("template_ggzz"+cat_name+"_all","", 140,105, 140,nbin,low,high); */
			/* tzx ->Draw(varname+":ZZMass>>template_zx"+cat_name,Form("weight*(%s==\"%s\"&&chan==%d)",categoryName.Data(),reco_name.Data(),k+1)); */
			/* float yield_zx = template2d_zx_tmp->Integral(); */
			/* tqqzz ->Draw(varname+":ZZMass>>template_qqzz"+cat_name,Form("weight*(%s==\"%s\"&&chan==%d)",categoryName.Data(),reco_name.Data(),k+1)); */
			/* float yield_qqzz = template2d_qqzz_tmp->Integral(); */
			/* tggzz ->Draw(varname+":ZZMass>>template_ggzz"+cat_name,Form("weight*(%s==\"%s\"&&chan==%d)",categoryName.Data(),reco_name.Data(),k+1)); */
			/* float yield_ggzz = template2d_ggzz_tmp->Integral(); */

			/* tzx ->Draw(varname+":ZZMass>>template_zx"+cat_name+"_all",Form("weight*(%s==\"%s\")",categoryName.Data(),reco_name.Data())); */
			/* tqqzz ->Draw(varname+":ZZMass>>template_qqzz"+cat_name+"_all",Form("weight*(%s==\"%s\")",categoryName.Data(),reco_name.Data())); */
			/* tggzz ->Draw(varname+":ZZMass>>template_ggzz"+cat_name+"_all",Form("weight*(%s==\"%s\")",categoryName.Data(),reco_name.Data())); */


			/* map<TString, float> sys_arr_qqzz_up ; */
			/* map<TString, float> sys_arr_qqzz_dn ; */
			/* map<TString, float> sys_arr_ggzz_up ; */
			/* map<TString, float> sys_arr_ggzz_dn ; */
			/* map<TString, float> sys_arr_zjet_up ; */
			/* map<TString, float> sys_arr_zjet_dn ; */
			/* if(doSys){ */

			/* 	// TString cut_jec_up = Form("weight*(%s_jec_up==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	// TString cut_jec_dn = Form("weight*(%s_jec_dn==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	// TString cut_jer_up = Form("weight*(%s_jer_up==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	// TString cut_jer_dn = Form("weight*(%s_jer_dn==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */

			/* 	TString cut_jes_up = Form("weight*(%s_jes_up==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	TString cut_jes_dn = Form("weight*(%s_jes_dn==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	TString cut_jer_up = Form("weight*(%s_jer_up==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	TString cut_jer_dn = Form("weight*(%s_jer_dn==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */

			/* 	// TString cut_jetPt_jes_up = Form("weight*(%s_jetPt_jes_up==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	// TString cut_jetPt_jes_dn = Form("weight*(%s_jetPt_jes_dn==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	// TString cut_jetPt_jer_up = Form("weight*(%s_jetPt_jer_up==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	// TString cut_jetPt_jer_dn = Form("weight*(%s_jetPt_jer_dn==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */

			/* 	TString cut_btag_up = Form("weight*( %s_btag_up==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */
			/* 	TString cut_btag_dn = Form("weight*( %s_btag_dn==\"%s\" && chan==%d )",categoryName.Data(),reco_name.Data(), k+1); */

			/* 	TH1F *jet_sys =new TH1F("jet_sys","",1,105, 140); */
			/* 	tqqzz->Draw("ZZMass>>jet_sys",cut_jes_up); */
			/* 	float yield_jes_up = jet_sys->Integral(); */
			/* 	tqqzz->Draw("ZZMass>>jet_sys",cut_jes_dn); */
			/* 	float yield_jes_dn = jet_sys->Integral(); */
			/* 	tqqzz->Draw("ZZMass>>jet_sys",cut_jer_up); */
			/* 	float yield_jer_up = jet_sys->Integral(); */
			/* 	tqqzz->Draw("ZZMass>>jet_sys",cut_jer_dn); */
			/* 	float yield_jer_dn = jet_sys->Integral(); */
			/* 	// tqqzz->Draw("ZZMass>>jet_sys",cut_jetPt_jes_up); */
			/* 	// float yield_jetPt_jes_up = jet_sys->Integral(); */
			/* 	// tqqzz->Draw("ZZMass>>jet_sys",cut_jetPt_jes_dn); */
			/* 	// float yield_jetPt_jes_dn = jet_sys->Integral(); */
			/* 	// tqqzz->Draw("ZZMass>>jet_sys",cut_jetPt_jer_up); */
			/* 	// float yield_jetPt_jer_up = jet_sys->Integral(); */
			/* 	// tqqzz->Draw("ZZMass>>jet_sys",cut_jetPt_jer_dn); */
			/* 	// float yield_jetPt_jer_dn = jet_sys->Integral(); */
			/* 	tqqzz->Draw("ZZMass>>jet_sys",cut_btag_up); */
			/* 	float yield_btag_up = jet_sys->Integral(); */
			/* 	tqqzz->Draw("ZZMass>>jet_sys",cut_btag_dn); */
			/* 	float yield_btag_dn = jet_sys->Integral(); */

			/* 	if(yield_qqzz!=0){ */
			/* 		// sys_arr_qqzz_up["jec"]=yield_jec_up/yield_qqzz; */
			/* 		// sys_arr_qqzz_dn["jec"]=yield_jec_dn/yield_qqzz; */
			/* 		// sys_arr_qqzz_up["jer"]=yield_jer_up/yield_qqzz; */
			/* 		// sys_arr_qqzz_dn["jer"]=yield_jer_dn/yield_qqzz; */

			/* 		sys_arr_qqzz_up["jes"]=yield_jes_up/yield_qqzz; */
			/* 		sys_arr_qqzz_dn["jes"]=yield_jes_dn/yield_qqzz; */
			/* 		sys_arr_qqzz_up["jer"]=yield_jer_up/yield_qqzz; */
			/* 		sys_arr_qqzz_dn["jer"]=yield_jer_dn/yield_qqzz; */
			/* 		// sys_arr_qqzz_up["jetPt_jes"]=yield_jetPt_jes_up/yield_qqzz; */
			/* 		// sys_arr_qqzz_dn["jetPt_jes"]=yield_jetPt_jes_dn/yield_qqzz; */
			/* 		// sys_arr_qqzz_up["jetPt_jer"]=yield_jetPt_jer_up/yield_qqzz; */
			/* 		// sys_arr_qqzz_dn["jetPt_jer"]=yield_jetPt_jer_dn/yield_qqzz; */
			/* 		sys_arr_qqzz_up["btag"]=yield_btag_up/yield_qqzz; */
			/* 		sys_arr_qqzz_dn["btag"]=yield_btag_dn/yield_qqzz; */
			/* 	} */
			/* 	for (int nextra = 0; nextra < extraSysn_all; nextra++){ */
			/* 		sys_arr_qqzz_up[extraSys_handName_all[nextra]]=extraSys_all_up[nextra][k]; */
			/* 		sys_arr_qqzz_dn[extraSys_handName_all[nextra]]=extraSys_all_dn[nextra][k]; */
			/* 	} */
			/* 	// tggzz->Draw("ZZMass>>jec_up",cut_jec_up); */
			/* 	// float yield_jec_up_gg = jec_up->Integral(); */
			/* 	// tggzz->Draw("ZZMass>>jec_up",cut_jec_dn); */
			/* 	// float yield_jec_dn_gg = jec_up->Integral(); */
			/* 	// tggzz->Draw("ZZMass>>jec_up",cut_jer_up); */
			/* 	// float yield_jer_up_gg = jec_up->Integral(); */
			/* 	// tggzz->Draw("ZZMass>>jec_up",cut_jer_dn); */
			/* 	// float yield_jer_dn_gg = jec_up->Integral(); */

			/* 	tggzz->Draw("ZZMass>>jet_sys",cut_jes_up); */
			/* 	float yield_jes_up_gg = jet_sys->Integral(); */
			/* 	tggzz->Draw("ZZMass>>jet_sys",cut_jes_dn); */
			/* 	float yield_jes_dn_gg = jet_sys->Integral(); */
			/* 	tggzz->Draw("ZZMass>>jet_sys",cut_jer_up); */
			/* 	float yield_jer_up_gg = jet_sys->Integral(); */
			/* 	tggzz->Draw("ZZMass>>jet_sys",cut_jer_dn); */
			/* 	float yield_jer_dn_gg = jet_sys->Integral(); */
			/* 	// tggzz->Draw("ZZMass>>jet_sys",cut_jetPt_jes_up); */
			/* 	// float yield_jetPt_jes_up_gg = jet_sys->Integral(); */
			/* 	// tggzz->Draw("ZZMass>>jet_sys",cut_jetPt_jes_dn); */
			/* 	// float yield_jetPt_jes_dn_gg = jet_sys->Integral(); */
			/* 	// tggzz->Draw("ZZMass>>jet_sys",cut_jetPt_jer_up); */
			/* 	// float yield_jetPt_jer_up_gg = jet_sys->Integral(); */
			/* 	// tggzz->Draw("ZZMass>>jet_sys",cut_jetPt_jer_dn); */
			/* 	// float yield_jetPt_jer_dn_gg = jet_sys->Integral(); */
			/* 	tggzz->Draw("ZZMass>>jet_sys",cut_btag_up); */
			/* 	float yield_btag_up_gg = jet_sys->Integral(); */
			/* 	tggzz->Draw("ZZMass>>jet_sys",cut_btag_dn); */
			/* 	float yield_btag_dn_gg = jet_sys->Integral(); */
			/* 	if(yield_ggzz!=0){ */
			/* 		// sys_arr_ggzz_up["jec"]=yield_jec_up_gg/yield_ggzz; */
			/* 		// sys_arr_ggzz_dn["jec"]=yield_jec_dn_gg/yield_ggzz; */
			/* 		// sys_arr_ggzz_up["jer"]=yield_jer_up_gg/yield_ggzz; */
			/* 		// sys_arr_ggzz_dn["jer"]=yield_jer_dn_gg/yield_ggzz; */

			/* 		sys_arr_ggzz_up["jes"]=yield_jes_up_gg/yield_ggzz; */
			/* 		sys_arr_ggzz_dn["jes"]=yield_jes_dn_gg/yield_ggzz; */
			/* 		sys_arr_ggzz_up["jer"]=yield_jer_up_gg/yield_ggzz; */
			/* 		sys_arr_ggzz_dn["jer"]=yield_jer_dn_gg/yield_ggzz; */
			/* 		// sys_arr_ggzz_up["jetPt_jes"]=yield_jetPt_jes_up_gg/yield_ggzz; */
			/* 		// sys_arr_ggzz_dn["jetPt_jes"]=yield_jetPt_jes_dn_gg/yield_ggzz; */
			/* 		// sys_arr_ggzz_up["jetPt_jer"]=yield_jetPt_jer_up_gg/yield_ggzz; */
			/* 		// sys_arr_ggzz_dn["jetPt_jer"]=yield_jetPt_jer_dn_gg/yield_ggzz; */
			/* 		sys_arr_ggzz_up["btag"]=yield_btag_up_gg/yield_ggzz; */
			/* 		sys_arr_ggzz_dn["btag"]=yield_btag_dn_gg/yield_ggzz; */
			/* 	} */
			/* 	for (int nextra=0; nextra< extraSysn_ggzz;nextra++){ */
			/* 		sys_arr_ggzz_up[extraSys_handName_ggzz[nextra]]= extraSys_ggzz[nextra][k]; */
			/* 		sys_arr_ggzz_dn[extraSys_handName_ggzz[nextra]]= 2-extraSys_ggzz[nextra][k]; */
			/* 	} */
			/* 	for (int nextra = 0; nextra < extraSysn_all; nextra++){ */
			/* 		sys_arr_ggzz_up[extraSys_handName_all[nextra]]=extraSys_all_up[nextra][k]; */
			/* 		sys_arr_ggzz_dn[extraSys_handName_all[nextra]]=extraSys_all_dn[nextra][k]; */
			/* 	} */
			/* 	TString syshist_name ="sys_hist_qqzz"; */
			/* 	tqqzz->Draw("weight_name>>"+syshist_name,cut_reco); */
			/* 	TH1F *systemp = (TH1F*)gROOT->FindObject(syshist_name); */
			/* 	int nbins_sys= systemp->GetNbinsX(); */

			/* 	for (int sysl = 0; sysl< nbins_sys;sysl++){ */
			/* 		TString sysname = systemp->GetXaxis()->GetBinLabel(sysl+1); */
			/* 		TString hsys_tmpName_up = "hsys_tmp_qqzz_"+sysname+"_up"; */
			/* 		TString hsys_tmpName_dn = "hsys_tmp_qqzz_"+sysname+"_dn"; */
			/* 		TH1F *hsys_tmp_up = new TH1F (hsys_tmpName_up,hsys_tmpName_up, 140,105, 140); */
			/* 		TH1F *hsys_tmp_dn = new TH1F (hsys_tmpName_dn,hsys_tmpName_dn, 140,105, 140); */
			/* 		tqqzz->Draw("ZZMass>>"+hsys_tmpName_up,Form("weight_up*(%s==\"%s\"&&chan==%d && weight_name == \"%s\")",categoryName.Data(),reco_name.Data(),k+1,sysname.Data())); */
			/* 		tqqzz->Draw("ZZMass>>"+hsys_tmpName_dn,Form("weight_dn*(%s==\"%s\"&&chan==%d && weight_name == \"%s\")",categoryName.Data(),reco_name.Data(),k+1,sysname.Data())); */
			/* 		float yield_sys_up = hsys_tmp_up->Integral(); */
			/* 		float yield_sys_dn = hsys_tmp_dn->Integral(); */
			/* 		sys_arr_qqzz_up[sysname]=yield_sys_up/yield_qqzz; */
			/* 		sys_arr_qqzz_dn[sysname]=yield_sys_dn/yield_qqzz; */
			/* 	} */

			/* 	for (int nextra = 0; nextra < extraSysn_zjet; nextra++){ */
			/* 		sys_arr_zjet_up[extraSys_handName_zjet[nextra]]=extraSys_zjet_up[nextra][k]; */
			/* 		sys_arr_zjet_dn[extraSys_handName_zjet[nextra]]=extraSys_zjet_dn[nextra][k]; */
			/* 	} */

			/* 	sys_arr_up.push_back(sys_arr_zjet_up); */
			/* 	sys_arr_dn.push_back(sys_arr_zjet_dn); */
			/* 	sys_arr_up.push_back(sys_arr_qqzz_up); */
			/* 	sys_arr_dn.push_back(sys_arr_qqzz_dn); */
			/* 	sys_arr_up.push_back(sys_arr_ggzz_up); */
			/* 	sys_arr_dn.push_back(sys_arr_ggzz_dn); */
			/* } */

			/* // Setting at 50 the min entries for which ZX histo is fitted */
			/* // See: https://indico.cern.ch/event/889021/contributions/3748912/attachments/1987318/3311893/CutBasedAnalysis_ZX_jets.pdf */
			/* // if(red_zx->numEntries()>500){ */
			/* if(red_zx->numEntries()>50){ */
			/* 	RooPlot *frame = ZZMass->frame(); */
			/* 	pdf_zx_cat->fitTo(*red_zx); */
			/* 	red_zx->plotOn(frame,RooFit::Binning(35)); */
			/* 	pdf_zx_cat->plotOn(frame); */
			/* 	pdf_zx_cat->paramOn(frame,Layout(0.1,0.4,0.6)); */
			/* 	frame->getAttText()->SetTextSize(0.03) ; */
			/* 	frame->Draw(); */
			/* 	gPad->Print(outputDir+"/fig/fit_zx_"+cat_name+year+".png"); */
			/* 	gPad->Print(outputDir+"/fig/fit_zx_"+cat_name+year+".pdf"); */

			/* 	RooConstVar* cland1=new RooConstVar("land1_"+cat_name+year,"land1_"+cat_name+year,land1->getVal()); */
			/* 	RooConstVar* cland2=new RooConstVar("land2_"+cat_name+year,"land2_"+cat_name+year,land2->getVal()); */
			/* 	RooConstVar* cland1_2=new RooConstVar("land1_2_"+cat_name+year,"land1_2_"+cat_name+year,land1_2->getVal()); */
			/* 	RooConstVar* cland2_2=new RooConstVar("land2_2_"+cat_name+year,"land2_2_"+cat_name+year,land2_2->getVal()); */
			/* 	RooConstVar* cfrac_land=new RooConstVar("frac_land_"+cat_name+year,"frac_land_"+cat_name+year,frac_land->getVal()); */
			/* 	// if(k!=2) */
			/* 	pdf_zx_cat=new RooGenericPdf("pdfm4l_zjets_"+cat_name+year, "pdfm4l_zjets_"+cat_name+year,"TMath::Landau(@0,@1,@2)",RooArgList(*ZZMass,*cland1,*cland2)); */
			/* 	// else */
			/* 	// 	pdf_zx_cat= new RooGenericPdf("pdfm4l_zjets_"+cat_name+year, "pdfm4l_zjets_"+cat_name+year,"TMath::Landau(@0,@1,@2)+@3*TMath::Landau(@0,@4,@5)",RooArgList(*ZZMass,*cland1,*cland2,*cfrac_land,*cland1_2,*cland2_2)); */

			/* } */
			/* else{ */
			/* 	pdf_zx_cat = pdf_zx_chan[k]; */
			/* } */
			/* if(red_qqzz->numEntries()>1000){ */
			/* 	RooPlot *frame = ZZMass->frame(); */
			/* 	pdf_qqzz_cat->fitTo(*red_qqzz); */
			/* 	red_qqzz->plotOn(frame,RooFit::Binning(35)); */
			/* 	pdf_qqzz_cat->plotOn(frame); */
			/* 	pdf_qqzz_cat->paramOn(frame,Layout(0.1,0.4,0.9)); */
			/* 	frame->getAttText()->SetTextSize(0.03) ; */
			/* 	frame->Draw(); */
			/* 	gPad->Print(outputDir+"/fig/fit_qqzz_"+cat_name+".png"); */
			/* 	gPad->Print(outputDir+"/fig/fit_qqzz_"+cat_name+".pdf"); */
			/* 	RooConstVar* cber1=new RooConstVar("ber1_"+cat_name+year,"ber1_"+cat_name+year,ber1->getVal()); */
			/* 	RooConstVar* cber2=new RooConstVar("ber2_"+cat_name+year,"ber2_"+cat_name+year,ber2->getVal()); */
			/* 	RooConstVar* cber3=new RooConstVar("ber3_"+cat_name+year,"ber3_"+cat_name+year,ber3->getVal()); */
			/* 	pdf_qqzz_cat = new RooBernstein("pdfm4l_qqzz_"+cat_name+year, "pdfm4l_qqzz_"+cat_name+year,*ZZMass,RooArgList(*cber1,*cber2,*cber3)); */
			/* } */
			/* else */
			/* 	pdf_qqzz_cat = pdf_qqzz_chan[k]; */

			/* pdf_ggzz_cat = pdf_ggzz_chan[k]; */

			/* RooProdPdf *pdfprod_qqzz; */
			/* RooProdPdf *pdfprod_ggzz; */
			/* //RooProdPdf *pdfprod_zx; */



			/* if(reco_name.Contains("VH_Had")|| reco_name.Contains("VBF_1j")|| reco_name.Contains("VBF_2j")|| reco_name.Contains("VBF_VBF")|| reco_name.Contains("VBF_ptj_GT200_1J")||reco_name.Contains("VBF_ptj_GT200_2J") ){ */

			/* 	//			if(reco_name.Contains("VH_Had") || reco_name.Contains("VBF_2j")|| reco_name.Contains("VBF_VBF")|| reco_name.Contains("VBF_ptj_GT200_1J")||reco_name.Contains("VBF_ptj_GT200_2J")){ */


			/* 	//Floor(template2d_zx_tmp_all,20); */
			/* 	Floor(template2d_qqzz_tmp_all,20); */
			/* 	Floor(template2d_ggzz_tmp_all,20); */
			/* 	//template2d_zx_tmp_all->Draw("colz"); */
			/* 	//template2d_zx_tmp_all->SetMaximum(0.2); */
			/* 	//gPad->Print(outputDir+"/fig/"+cat_name+"_zx_template.png"); */
			/* 	//gPad->Print(outputDir+"/fig/"+cat_name+"_zx_template.pdf"); */
			/* 	template2d_qqzz_tmp_all->Draw("colz"); */
			/* 	template2d_qqzz_tmp_all->SetMaximum(0.2); */
			/* 	gPad->Print(outputDir+"/fig/"+cat_name+"_qqzz_template.png"); */
			/* 	gPad->Print(outputDir+"/fig/"+cat_name+"_qqzz_template.pdf"); */
			/* 	template2d_ggzz_tmp_all->Draw("colz"); */
			/* 	template2d_ggzz_tmp_all->SetMaximum(0.2); */
			/* 	gPad->Print(outputDir+"/fig/"+cat_name+"_ggzz_template.png"); */
			/* 	gPad->Print(outputDir+"/fig/"+cat_name+"_ggzz_template.pdf"); */

			/* 	RooDataHist *datahist_qqzz_tmp; */
			/* 	RooDataHist *datahist_ggzz_tmp; */
			/* 	//RooDataHist *datahist_zx_tmp; */

			/* 	RooHistPdf *pdf2d_qqzz_tmp; */
			/* 	RooHistPdf *pdf2d_ggzz_tmp; */
			/* 	//RooHistPdf *pdf2d_zx_tmp; */

			/* 	datahist_qqzz_tmp=new RooDataHist("datahist_qqzz"+cat_name+year,"datahist_qqzz"+cat_name+year,RooArgSet(*ZZMass,*Ddis),template2d_qqzz_tmp_all); */
			/* 	pdf2d_qqzz_tmp=new RooHistPdf("histpdf_aqqzz"+cat_name+year,"histpdf_qqzz"+cat_name+year,RooArgSet(*ZZMass,*Ddis),*datahist_qqzz_tmp); */

			/* 	datahist_ggzz_tmp=new RooDataHist("datahist_ggzz"+cat_name+year,"datahist_ggzz"+cat_name+year,RooArgSet(*ZZMass,*Ddis),template2d_ggzz_tmp_all); */
			/* 	pdf2d_ggzz_tmp=new RooHistPdf("histpdf_aggzz"+cat_name+year,"histpdf_ggzz"+cat_name+year,RooArgSet(*ZZMass,*Ddis),*datahist_ggzz_tmp); */

			/* 	/\* datahist_zx_tmp=new RooDataHist("datahist_zx"+cat_name+year,"datahist_zx"+cat_name+year,RooArgSet(*ZZMass,*Ddis),template2d_zx_tmp_all); *\/ */
			/* 	/\* pdf2d_zx_tmp=new RooHistPdf("histpdf_azx"+cat_name+year,"histpdf_zx"+cat_name+year,RooArgSet(*ZZMass,*Ddis),*datahist_zx_tmp); *\/ */

			/* 	pdfprod_qqzz=new RooProdPdf ("qqzz","qqzz",*pdf_qqzz_cat,Conditional(*pdf2d_qqzz_tmp,*Ddis)); */
			/* 	pdfprod_ggzz=new RooProdPdf ("ggzz","ggzz",*pdf_ggzz_cat,Conditional(*pdf2d_ggzz_tmp,*Ddis)); */
			/* 	//pdfprod_zx=new RooProdPdf ("zjets","zjets",*pdf_zx_cat,Conditional(*pdf2d_zx_tmp,*Ddis)); */
			/* } */
			/* else{ */

			/* 	pdfprod_qqzz=new RooProdPdf ("qqzz","qqzz",*pdf_qqzz_cat,Conditional(*pdfqqzz2d[k],*dbkg)); */
			/* 	pdfprod_ggzz=new RooProdPdf ("ggzz","ggzz",*pdf_ggzz_cat,Conditional(*pdfggzz2d[k],*dbkg)); */
			/* 	//pdfprod_zx=new RooProdPdf ("zjets","zjets",*pdf_zx_cat,Conditional(*pdfzx2d[k],*dbkg)); */
			/* } */

			/* w.import(*pdfprod_qqzz,RecycleConflictNodes()); */
			/* w.import(*pdfprod_ggzz,RecycleConflictNodes()); */
			/* //w.import(*pdfprod_zx,RecycleConflictNodes()); */

			/* /\* bin_arr.push_back(cat_name); *\/ */
			/* /\* processName_arr.push_back("zjets"); *\/ */
			/* /\* process_arr.push_back (1); *\/ */
			/* /\* yield_arr.push_back(yield_zx); *\/ */

			/* bin_arr.push_back(cat_name); */
			/* processName_arr.push_back("qqzz"); */
			/* process_arr.push_back (2); */
			/* yield_arr.push_back(yield_qqzz); */

			/* bin_arr.push_back(cat_name); */
			/* processName_arr.push_back("ggzz"); */
			/* process_arr.push_back (3); */
			/* yield_arr.push_back(yield_ggzz); */

			/* card<<"bin \t"; */
			/* writeline(bin_arr,card); */
			/* card<<"process \t"; */
			/* writeline(processName_arr,card); */
			/* card<<"process \t"; */
			/* writeline(process_arr,card); */
			/* card<<"rate \t"; */
			/* writeline(yield_arr,card); */
			/* int extrashape = 3; */
			/* if(doSys){ */
			/* 	TString extraSys_shape[2]={scale_e,scale_m}; */
			/* 	for (int sysl = 0; sysl< nsys+extraSysn_zjet+extraSysn_all+extraSysn_ggzz+extraCat+1+extrashape;sysl++){ */
			/* 		TString sysname; */
			/* 		if(sysl<nsys){ */
			/* 			sysname	= sysnames.at(sysl); */
			/* 		} */
			/* 		else if (sysl<nsys+extraCat) */
			/* 			sysname = extraSysName[sysl-nsys]; */
			/* 		else if(sysl<nsys+extraSysn_zjet+extraCat) */
			/* 			sysname = extraSys_handName_zjet[sysl-nsys-extraCat]; */
			/* 		else if (sysl <nsys+extraSysn_zjet+extraSysn_all+extraCat) */
			/* 			sysname = extraSys_handName_all[sysl-nsys-extraSysn_zjet-extraCat]; */
			/* 		else if (sysl <nsys+extraSysn_zjet+extraSysn_all+extraCat+extraSysn_ggzz) */
			/* 			sysname = extraSys_handName_ggzz[sysl-nsys-extraSysn_zjet-extraSysn_all-extraCat]; */
			/* 		else if (sysl <nsys+extraSysn_zjet+extraSysn_all+extraCat+extraSysn_ggzz+1) */
			/* 			sysname = "hzz_br"; */
			/* 		else if (sysl < nsys+extraSysn_zjet+extraSysn_all+extraCat+extraSysn_ggzz+1+1) */
			/* 			sysname = res_ch; */
			/* 		else */
			/* 			//sysname = extraSys_shape[sysl-nsys-extraSysn_zjet-extraSysn_all-extraSysn_ggzz-4]+"_"+year; */
			/* 			sysname = extraSys_shape[sysl-nsys-extraSysn_zjet-extraSysn_all-extraSysn_ggzz-extraCat-2]; */

			/* 		cout<< sysname<<endl; */
			/* 		if(sysname =="QCDscale_muR_ggH" || sysname =="QCDscale_muF_ggH"||sysname=="QCDscale_muF_VBF"||sysname=="QCDscale_muR_VBF"||sysname=="pythia_cat") */
			/* 			continue; */
			/* 		if (sysname.Contains("CMS_scale_e") && (k==0)) */
			/* 			continue; */
			/* 		if (sysname.Contains("CMS_scale_m") && (k==1)) */
			/* 			continue; */
			/* 		if(sysl <nsys+extraSysn_zjet+extraSysn_all+extraCat+extraSysn_ggzz+1) */
			/* 			card<<sysname<<"\t lnN \t"; */
			/* 		else */
			/* 			card<<sysname<<"\t shape \t"; */

			/* 		for (int nstage = 0; nstage < yield_arr.size(); nstage++){ */
			/* 			std::map<TString,float>::iterator it = sys_arr_up[nstage].find(sysname); */
			/* 			if(it!=sys_arr_up[nstage].end()){ */
			/* 				std::map<TString,float>::iterator it_dn = sys_arr_dn[nstage].find(sysname); */
			/* 				float up = it->second >0 ? it->second : 1.; */
			/* 				float dn = it_dn->second >0 ? it_dn->second :1; */
			/* 				if(sysl <nsys+extraSysn_zjet+extraSysn_all+extraCat) */
			/* 					card<< up<<"/"<<dn<<"\t"; */
			/* 				else */
			/* 					card<< up <<"\t"; */
			/* 			} */
			/* 			else */
			/* 				card<< "- \t"; */
			/* 		} */
			/* 		card<<endl; */
			/* 	} */
			/* } */
			/* //int yearN = year.Atoi(); */
			/* //TString res = Form("Res_%s%s param 0 %.2f [-5,5]",chanName[k].Data(), year.Data(),res_unc_year[yearN-2016][k]); */
			/* //TString scale = Form("Scale_%s%s param 0 %.4f [-0.015,0.015]",chanName[k].Data(), year.Data(),scale_unc_year[yearN-2016][k]); */
			/* //card<<res<<endl; */
			/* //card<<scale<<endl; */
			/* //card<<"lumiscale rateParam * * 1"<<endl; */
			/* //card<<"nuisance edit freeze lumiscale"<<endl; */
			/* card.close(); */
			/* f->cd(); */


			/* TTree* reducetree= tdata->CopyTree(Form("chan==%d&&%s==%d",k+1,category.Data(),reco_cat)); */
			/* RooDataSet *data =new RooDataSet("data_obs","",reducetree,RooArgList(*ZZMass,*Ddis)); */
			/* cout<< data->numEntries()<<endl; */
			/* w.import(*data); */
			/* w.Write(); */
			/* f->Close(); */
			/* } */

			//}

		}




		//Datacard for event category: hzz_13TeV_2e2mu_Untagged
		//
		//
		//------------------------------------------------------------
		//imax 1 number of bins
		//jmax 12 number of processes minus 1
		//kmax 41 number of nuisance parameters
		//------------------------------------------------------------
		//shapes *    cat_hzz_13TeV_2e2mu_Untagged  hzz4lcard_2e2mu_0.input.root w:$PROCESS
		//------------------------------------------------------------
		//bin          cat_hzz_13TeV_2e2mu_Untagged
		//observation  -1
		//------------------------------------------------------------
		//bin          cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged cat_hzz_13TeV_2e2mu_Untagged
		//process      WH_had_hzz WH_lep_hzz ZH_had_hzz ZH_lep_hzz bbH_hzz ggH_hzz qqH_hzz tqH_hzz ttH_had_hzz ttH_lep_hzz ggZZ_hzz qqZZ_hzz zjets_hzz
		//process      -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3
		//rate         1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.9454 27.3966 11.8287
		//------------------------------------------------------------
