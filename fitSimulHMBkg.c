void fitSimulHMBkg(int year=2018){

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	TString chanName[3]={"4mu","4e","2e2mu"};
	int mass[1]={200};

	const int nstage=2;
	TString stageName [nstage]={"qqZZ","ggZZ"};
	int stage [nstage]= {-1,-3};
	//TString stageName [nstage]={"ggZZ"};
	//int stage [nstage]= {-3};
	
        RooRealVar *ZZMass =new RooRealVar("ZZMass","",100,1100);
	RooRealVar *stagev=new RooRealVar("htxs_stage1_red_prod_cat","",-4,1000);
	RooRealVar *chan=new RooRealVar("chan","",0,1000);
	RooRealVar *weight=new RooRealVar("weight","",1.e-9,1.e9);

	TString inputTrees = "/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All/ReducedTrees/";
	TString path = "/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All/";
	
	for (int j =0;j<3;j++){
		ofstream outmass ;
		//outmass.open(Form(path + "shapeParam/sim_massParam_qqZZ_%s_%d_%d.txt",chanName[j].Data(),mass[0],year));
		RooDataSet *data[1];
		for (int i =0;i<1;i++){
    			TChain *t = new TChain ("SelectedTree");
			//outmass.open(Form(path + "shapeParam/sim_massParam_%s_%d_%d.txt",chanName[j].Data(),mass[0],year));
			t->Add(Form(inputTrees + "*1_redTree_*%d*.root",year));
			//t->Add(Form(inputTrees + "*_Contin_MCFM701_redTree_*%d*.root",year));
			data[i]=new RooDataSet(Form("data%d_%d",i,mass[i]),"",t,RooArgSet(*ZZMass,*weight,*chan,*stagev),Form("chan==%d",j+1),"weight");
		}

		for (int k =0;k<nstage;k++){

		        outmass.open(Form(path + "shapeParamNew/sim_massParam_%s_%s_%d_%d.txt",stageName[k].Data(),chanName[j].Data(),mass[0],year));

			RooCategory massrc("massrc","");
			RooSimultaneous simPdf("simPdf","simultaneous pdf",massrc) ;
			RooDataSet* data_stage[1]; 
			RooDataSet* data_stage_cat[1]; 
			RooDataSet* data_all;
			RooDoubleCBFast* cpdf[1]; 
			RooAddPdf* addpdf[1]; 
		        RooqqZZPdf_v2* qqzzpdf[1];
			RooggZZPdf_v2* ggzzpdf[1];
        

			for (int i=0;i<1;i++){
				TString cat_type =Form("mh%d",mass[i]);
				massrc.defineType(cat_type,mass[i]);
			}

			RooRealVar *landau_mean_1=new RooRealVar("landau_mean_1","",0,-5,5);
			RooRealVar *landau_mean_2=new RooRealVar("landau_mean_2","",0,-5,5);
			RooRealVar *landau_mean_0=new RooRealVar("landau_mean_0","",150,100,300);
			RooRealVar *landau_sigma_1=new RooRealVar("landau_sigma_1","",0,-5,5);
			RooRealVar *landau_sigma_2=new RooRealVar("landau_sigma_2","",0,-5,5);
			RooRealVar *landau_sigma_0=new RooRealVar("landau_sigma_0","",20,0,100);
			

			RooRealVar *a0_1=new RooRealVar("a0_1","",0.001,-1.,1.);
			RooRealVar *a0_2=new RooRealVar("a0_2","",0.001,-0.1,0.1);
			RooRealVar *a0_0=new RooRealVar("a0_0","",177.9,50.,250.);
			RooRealVar *a1_1=new RooRealVar("a1_1","",0.001,-1.,1.);
			RooRealVar *a1_2=new RooRealVar("a1_2","",0.0,-0.1,0.1);
			RooRealVar *a1_0=new RooRealVar("a1_0","",70.38,0.,150.);
			RooRealVar *a2_1=new RooRealVar("a2_1","",0.001,-1.,1.);
			RooRealVar *a2_2=new RooRealVar("a2_2","",0.001,-0.1,0.1);
			RooRealVar *a2_0=new RooRealVar("a2_0","",153.1,0.,250.);
			RooRealVar *a3_1=new RooRealVar("a3_1","",0.001,-2.,2.);
			RooRealVar *a3_2=new RooRealVar("a3_2","",0.001,-0.1,0.1);
			RooRealVar *a3_0=new RooRealVar("a3_0","",0.041,0.,1.);
			RooRealVar *a4_1=new RooRealVar("a4_1","",0.001,-1.,1.);
			RooRealVar *a4_2=new RooRealVar("a4_2","",0.001,-0.1,0.1);
			RooRealVar *a4_0=new RooRealVar("a4_0","",183.5,0.,400.);
			RooRealVar *a5_1=new RooRealVar("a5_1","",0.001,-1.,1.);
			RooRealVar *a5_2=new RooRealVar("a5_2","",0.001,-0.1,0.1);
			RooRealVar *a5_0=new RooRealVar("a5_0","",11.67,0.,100.);
			RooRealVar *a6_1=new RooRealVar("a6_1","",0.001,-1.,1.);
			RooRealVar *a6_2=new RooRealVar("a6_2","",0.001,-0.1,0.1);
			RooRealVar *a6_0=new RooRealVar("a6_0","",35.24,0.,70.);
			RooRealVar *a7_1=new RooRealVar("a7_1","",0.001,-1.,1.);
			RooRealVar *a7_2=new RooRealVar("a7_2","",0.0,-0.1,0.1);
			RooRealVar *a7_0=new RooRealVar("a7_0","",0.54,0.,1.);
			RooRealVar *a8_1=new RooRealVar("a8_1","",0.001,-1.,1.);
			RooRealVar *a8_2=new RooRealVar("a8_2","",0.001,-0.1,0.1);
			RooRealVar *a8_0=new RooRealVar("a8_0","",28.94,0.,300.);
			RooRealVar *a9_1=new RooRealVar("a9_1","",0.001,-1.,1.);
			RooRealVar *a9_2=new RooRealVar("a9_2","",0.001,-0.1,0.1);
			RooRealVar *a9_0=new RooRealVar("a9_0","",-0.322,-1.,1.);

			RooRealVar *a10_1=new RooRealVar("a10_1","",0.001,-1.,1.);
			RooRealVar *a10_2=new RooRealVar("a10_2","",0.00146,-0.1,0.1);
			RooRealVar *a10_0=new RooRealVar("a10_0","",100.,0.,2000.);
			RooRealVar *a11_1=new RooRealVar("a11_1","",0.001,-1.,1.);
			RooRealVar *a11_2=new RooRealVar("a11_2","",0.00146,-0.1,0.1);
			RooRealVar *a11_0=new RooRealVar("a11_0","",-5.,-20.,20.);
			RooRealVar *a12_1=new RooRealVar("a12_1","",0.001,-1.,1.);
			RooRealVar *a12_2=new RooRealVar("a12_2","",0.00146,-0.1,0.1);
			RooRealVar *a12_0=new RooRealVar("a12_0","",100,0.,2000.);
			RooRealVar *a13_1=new RooRealVar("a13_1","",0.001,-1.,1.);
			RooRealVar *a13_2=new RooRealVar("a13_2","",0.00146,-0.1,0.1);
			RooRealVar *a13_0=new RooRealVar("a13_0","",0.2,0.,1.);
			 
						
			for (int i=0;i<1;i++){


				RooConstVar *MH=new RooConstVar("MH","",mass[i]);
		        
				RooFormulaVar* a0=new RooFormulaVar(Form("a0_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a0_0,*a0_1,*MH));
				RooFormulaVar* a1=new RooFormulaVar(Form("a1_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a1_0,*a1_1,*MH));
				RooFormulaVar* a2=new RooFormulaVar(Form("a2_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a2_0,*a2_1,*MH));
				RooFormulaVar* a3=new RooFormulaVar(Form("a3_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a3_0,*a3_1,*MH));
				RooFormulaVar* a4=new RooFormulaVar(Form("a4_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a4_0,*a4_1,*MH));
				RooFormulaVar* a5=new RooFormulaVar(Form("a5_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a5_0,*a5_1,*MH));
				RooFormulaVar* a6=new RooFormulaVar(Form("a6_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a6_0,*a6_1,*MH));
				RooFormulaVar* a7=new RooFormulaVar(Form("a7_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a7_0,*a7_1,*MH));
				RooFormulaVar* a8=new RooFormulaVar(Form("a8_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a8_0,*a8_1,*MH));
				RooFormulaVar* a9=new RooFormulaVar(Form("a9_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a9_0,*a9_1,*MH));

				RooFormulaVar* a10=new RooFormulaVar(Form("a10_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a10_0,*a10_1,*MH));
				RooFormulaVar* a11=new RooFormulaVar(Form("a11_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a11_0,*a11_1,*MH));
				RooFormulaVar* a12=new RooFormulaVar(Form("a12_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a12_0,*a12_1,*MH));
				RooFormulaVar* a13=new RooFormulaVar(Form("a13_%d",mass[i]),"","@0+@1*(MH-MH)",RooArgList(*a13_0,*a13_1,*MH));

				

				TString cat_type =Form("mh%d",mass[i]);
				data_stage[i] = (RooDataSet*)data[i]->reduce(Form("htxs_stage1_red_prod_cat==%d",stage[k])); 
				data_stage_cat[i] = new RooDataSet(Form("data_%d",i),"",RooArgSet(*ZZMass,*weight),Index(massrc),Import(cat_type,*data_stage[i]),WeightVar("weight")); 

				//qqzzpdf[i] = new RooqqZZPdf_v2(Form("qqZZpdf_%d",i),"",*ZZMass,*a0,*a1,*a2,*a3,*a4,*a5,*a6,*a7,*a8,*a9,*a10,*a11,*a12,*a13);
				//ggzzpdf[i] = new RooggZZPdf_v2(Form("ggZZpdf_%d",i),"",*ZZMass,*a0,*a1,*a2,*a3,*a4,*a5,*a6,*a7,*a8,*a9);
				cout<< data_stage_cat[i]->sumEntries()<<endl;
				if (stageName[k]=="qqZZ"){
				  qqzzpdf[i] = new RooqqZZPdf_v2(Form("qqZZpdf_%d",i),"",*ZZMass,*a0,*a1,*a2,*a3,*a4,*a5,*a6,*a7,*a8,*a9,*a10,*a11,*a12,*a13);
				  simPdf.addPdf(*qqzzpdf[i],cat_type);
				}
				if (stageName[k]=="ggZZ"){
				  ggzzpdf[i] = new RooggZZPdf_v2(Form("ggZZpdf_%d",i),"",*ZZMass,*a0,*a1,*a2,*a3,*a4,*a5,*a6,*a7,*a8,*a9);
				  simPdf.addPdf(*ggzzpdf[i],cat_type);
				}
				if(i==0)
					data_all = data_stage_cat[i];
				else
					data_all->append(*data_stage_cat[i]);
				if(i==0){ 
				  //RooqqZZPdf_v2 *qqzzpdftmp = new RooqqZZPdf_v2("qqZZpdf","",*ZZMass,*a0_0,*a1_0,*a2_0,*a3_0,*a4_0,*a5_0,*a6_0,*a7_0,*a8_0,*a9_0,*a10_0,*a11_0,*a12_0,*a13_0);
				  //RooggZZPdf_v2 *ggzzpdftmp = new RooggZZPdf_v2(Form("ggZZpdf_%d",i),"",*ZZMass,*a0_0,*a1_0,*a2_0,*a3_0,*a4_0,*a5_0,*a6_0,*a7_0,*a8_0,*a9_0);
					

					if (stageName[k]=="qqZZ"){
					  RooqqZZPdf_v2 *qqzzpdftmp = new RooqqZZPdf_v2("qqZZpdf","",*ZZMass,*a0_0,*a1_0,*a2_0,*a3_0,*a4_0,*a5_0,*a6_0,*a7_0,*a8_0,*a9_0,*a10_0,*a11_0,*a12_0,*a13_0);
					  qqzzpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));
					}
					if (stageName[k]=="ggZZ"){
					  RooggZZPdf_v2 *ggzzpdftmp = new RooggZZPdf_v2(Form("ggZZpdf_%d",i),"",*ZZMass,*a0_0,*a1_0,*a2_0,*a3_0,*a4_0,*a5_0,*a6_0,*a7_0,*a8_0,*a9_0);
					  ggzzpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));
					}
				}
			}
	        
			simPdf.Print("v");
			data_all->Print("v");
			simPdf.fitTo(*data_all,InitialHesse(true),Strategy(2)) ;
			TString a0outs= Form ("a0_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a0_0->getVal(),a0_1->getVal());
			TString a1outs= Form ("a1_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a1_0->getVal(),a1_1->getVal());
			TString a2outs= Form ("a2_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a2_0->getVal(),a2_1->getVal());
			TString a3outs= Form ("a3_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a3_0->getVal(),a3_1->getVal());
			TString a4outs= Form ("a4_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a4_0->getVal(),a4_1->getVal());
			TString a5outs= Form ("a5_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a5_0->getVal(),a5_1->getVal());
			TString a6outs= Form ("a6_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a6_0->getVal(),a6_1->getVal());
			TString a7outs= Form ("a7_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a7_0->getVal(),a7_1->getVal());
			TString a8outs= Form ("a8_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a8_0->getVal(),a8_1->getVal());
			TString a9outs= Form ("a9_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a9_0->getVal(),a9_1->getVal());

			TString a10outs= Form ("a10_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a10_0->getVal(),a10_1->getVal());
			TString a11outs= Form ("a11_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a11_0->getVal(),a11_1->getVal());
			TString a12outs= Form ("a12_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a12_0->getVal(),a12_1->getVal());
			TString a13outs= Form ("a13_%s_%d_%s_%d \t %.4f+%.5f*(MH-MH)",stageName[k].Data(),mass[0],chanName[j].Data(),year,a13_0->getVal(),a13_1->getVal());
			

			outmass<< a0outs<<endl;
			outmass<< a1outs<<endl;
			outmass<< a2outs<<endl;
			outmass<< a3outs<<endl;
			outmass<< a4outs<<endl;
			outmass<< a5outs<<endl;
			outmass<< a6outs<<endl;
			outmass<< a7outs<<endl;
			outmass<< a8outs<<endl;
			outmass<< a9outs<<endl;
			
			if (stageName[k]=="qqZZ"){			
			  outmass<< a10outs<<endl;
			  outmass<< a11outs<<endl;
			  outmass<< a12outs<<endl;
			  outmass<< a13outs<<endl;
			}
			cout<<"where1"<<endl;
			for (int i =0;i<1;i++){
			        RooPlot *frame = ZZMass->frame();
				TString cat_type =Form("mh%d",mass[i]);
				data_all->plotOn(frame,Cut(Form("massrc==massrc::mh%d",mass[i]))) ;
				simPdf.plotOn(frame,Slice(massrc,cat_type),ProjWData(massrc,*data_all)) ;
				frame->chiSquare();
				frame->pullHist() ;
				frame->residHist() ;
				frame->SetTitle(Form("%s -> %s",stageName[k].Data(),chanName[j].Data()));
				frame->GetXaxis()->SetTitle("ZZ mass (GeV)");
			    	frame->Draw();
				//gPad->SetLogy();
				gPad->Print(Form(path + "shapeNew/fig_%d/simFit_%s_%s_%d.png", year, stageName[k].Data(),chanName[j].Data(),year));
				
			}
		}
		outmass.close();
	}
}
