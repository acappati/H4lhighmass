void fitSimulHM_cat(int year=2018){

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TString chanName[3]={"4mu","4e","2e2mu"};
       	int mass [21]={125,150,175,200,250,300,350,400,450,500,550,600,700,750,800,900,1000,1500,2000,2500,3000};
	//int mass [11]={200,250,300,350,400,450,500,550,600,700,750};
	//int mass [6]={800,900,1000,1500,2000,3000};
	int N = 21;
       
	const int nstage=2;
	TString stageName [nstage]={"ggH","VBF_2jets"};
	int stage [nstage]= {0,1};;

	float a1_DCB [6][21]={
	  1.2419,1.15,1.2,1.2194,1.3056,1.1834,1.3413,1.3367,1.2401,1.2318,1.3577,1.1751,1.1597,1.1928,1.2527,1.1341,1.2005,1.2003,1.2233,12122,1.2064, //4mu ggH
	  0.9713,1.02,1.04,1.0670,1.1198,1.1782,1.2477,1.1772,1.3321,1.2215,1.0133,1.3804,1.1538,1.2304,1.3627,1.3488,1.3540,1.2640,1.4640,1.5,1.5446, //4e ggH
	  0.9684,1.12,1.14,1.1755,1.1485,1.1330,1.1965,1.2451,1.2556,1.2068,1.3711,1.3008,1.3171,1.3458,1.2813,1.2628,1.2317,1.0120,1.0343,1.07,1.1026, //2mu2e ggH
	  1.2,1.2,1.2,1.2470,1.2669,1.3139,1.2892,1.2717,1.1509,1.1685,1.1076,1.2109,1.2654,1.3116,1.3301,1.0984,1.0898,1.0002,1.0367,1.2,1.4186, //4mu VBFH
	  1.0,1.,1.,1.0024,1.3677,1.0246,1.1189,1.2159,1.0011,1.1527,1.0551,1.3614,1.2521,1.1731,1.2292,1.1205,1.1301,1.2138,1.2611,1.4,1.5926, //4e VBFH
	  1.,1.,1.,1.0590,1.1080,1.1260,1.2288,1.1138,1.1882,1.1735,1.0918,1.2937,1.2089,1.1811,1.1615,1.2028,1.1740,1.1746,1.1886,1.16,1.1474  //2mu2e VBFH
	};
	float a2_DCB [6][21]={
	  1.737,1.8,1.7,1.7329,1.8515,1.6352,1.7223,1.6797,1.5217,1.5876,1.7580,1.1891,1.1888,1.2232,1.1869,1.2630,1.2961,1.3098,1.3309,1.4,1.4513,
	  1.4817,1.8,1.7,1.7445,1.5366,1.7465,1.5923,1.5026,1.6559,1.4645,1.1241,1.6973,1.0017,1.3534,1.6169,1.4211,1.4492,1.2599,1.3479,1.33,1.3210,
	  1.49,1.8,1.7,1.6979,1.7005,1.5763,1.5858,1.6247,1.6256,1.5809,1.6747,1.5011,1.4856,1.4784,1.3876,1.3233,1.3701,1.0423,1.0289,1.05,1.0788,
	  1.9,1.9,1.9,1.9183,1.8129,1.8694,1.6360,1.5863,1.0010,1.0129,1.0608,1.5351,1.1812,2.0776,1.4164,1.2640,1.2893,1.0035,1.0106,1.2,1.4869,
	  1.9,1.9,1.9,1.9585,1.8763,1.7642,1.4552,1.5432,1.3811,1.4278,1.4475,1.5333,1.4776,1.5569,1.3670,1.4236,1.3001,1.4630,1.4887,1.47,1.4653,
	  1.9,1.9,1.9,1.7349,1.5972,1.5839,1.5865,1.4639,1.4984,1.4808,1.4942,1.4835,1.3565,1.2604,1.2841,1.2649,1.3415,1.2723,1.2403,1.3,1.3707
	};
	float n1_DCB [6][21]={
	  2.037,1.4,1.6,1.7812,1.8153,2.2590,2.0530,2.1682,2.7169,2.8426,2.4775,2.9803,3.2299,3.4444,3.0335,3.7047,4.2011,3.0141,1.5893,1.6,1.6273,
	  5.4568,2.3,2.4,2.5282,2.3005,2.2804,2.0976,2.3546,1.8135,2.1955,2.9664,1.7682,1.9978,2.0053,1.7789,1.7352,1.7147,1.8324,1.2500,1.1,1.0010,
	  3.6264,1.7,1.8,1.8120,2.1907,2.3989,2.3008,2.2780,2.4135,2.7077,2.0415,2.3285,2.2460,2.1987,2.2848,2.3302,2.3305,2.6363,2.2268,1.9,1.5629,
	  1.4,1.5,1.6,1.6126,1.7532,1.8159,1.9957,2.2397,2.3804,2.4510,3.1545,2.9462,2.6360,2.6051,2.5435,2.6779,2.6625,2.5020,2.5142,2.52,2.5247,
	  2.2,2.4,2.6,2.7588,1.4158,2.9465,2.3930,2.1092,3.1218,2.4703,2.7062,1.8016,2.1265,2.3799,2.0897,2.6818,2.5226,2.4442,2.3063,2.0,1.4456,
	  1.7,1.8,1.9,1.9433,2.0555,2.1605,2.1170,2.7849,2.4604,2.2141,2.7054,2.0967,2.4989,2.3064,2.5088,2.9055,2.8528,2.5623,4.5470,3.2,2.5025
	};
	float n2_DCB [6][21]={
	  3.375,2.3,2.5,2.7300,3.1547,3.9043,3.7195,4.1315,5.6250,4.3445,5.7724,7.5556,9.8253,9.9259,9.6037,7.1887,6.9848,5.7123,3.9931,2.3,1.9656,
	  8.4376,2.5,2.8,3.0247,5.8120,4.7525,5.7845,6.1304,4.8972,7.8430,7.9764,6.7244,7.0848,6.7780,7.1157,6.9672,5.5329,5.2339,4.6659,4.7,4.7416,
	  5.7589,2.5,3.,3.7065,3.5213,4.7145,5.1935,5.0584,5.0224,5.7709,4.7919,7.3961,6.5437,6.7708,7.5857,7.2731,6.2044,6.5999,5.2371,3.5,2.7961,
	  1.5,1.7,1.8,1.7967,2.4162,2.4376,3.5257,3.5147,4.1033,5.2645,7.3337,7.5293,6.5710,6.0014,4.8617,3.9502,3.7235,3.7503,3.9706,4.0,4.0411,
	  1.5,1.5,1.5,1.4292,2.3610,2.6239,5.3727,4.1560,6.0071,5.5596,4.9926,5.6982,5.5362,4.5224,6.3557,5.8163,7.2659,5.6072,5.0098,5.2,5.3683,
	  1.5,1.8,2.1,2.1439,3.1216,3.4119,3.6354,4.3947,4.3994,5.6987,6.9944,7.6535,8.8642,8.9118,8.5872,7.2481,5.6786,5.0373,5.1350,5.15,5.2018
	};
	float mean_DCB [6][21]={
	  -0.1497,-0.15,-0.2,-0.2733,-0.4556,-0.5071,-0.8003,-0.8387,-0.8950,-0.9835,-1.0939,-1.1994,-1.4743,-1.6215,-1.7332,-1.6752,-1.3437,-1.8231,-2.1106,-1.78,-1.1403,
	  -1.4902,-2.,-2.,-2.0441,-2.3380,-2.4606,-2.6021,-2.6051,-2.7330,-2.7124,-2.6962,-2.7921,-2.9619,-2.8934,-2.9749,-3.0776,-3.2397,-3.2648,-3.8052,-3.7,-3.6949,
	  -0.7704,-1.2,-1.2,-1.2179,-1.3987,-1.4794,-1.6126,-1.6873,-1.7742,-1.8225,-1.9961,-2.0185,-2.1402,-2.4026,-2.2819,-2.5162,-2.4708,-2.4248,-2.6631,-2.5,-2.3918,
	  -0.2,-0.2,-0.3,-0.3190,-0.4010,-0.5921,-0.6866,-0.7632,-1.0320,-1.3180,-1.3115,-1.2150,-1.5638,-1.6059,-1.9173,-1.7971,-1.3286,-2.4428,-2.4847,-2.5,-4.5943,
	  -2.,-2.,-2.,-2.0672,-2.4014,-2.3934,-2.6989,-2.7380,-2.8086,-2.7428,-2.7260,-3.0907,-3.1382,-2.9275,-3.1296,-3.1150,-3.2262,-3.6326,-4.2784,-5.0,-5.8606,
	  -1.,-1.05,-1.1,-1.1232,-1.3046,-1.4514,-1.6463,-1.7011,-1.8138,-2.1489,-2.0556,-2.2590,-2.7241,-2.8792,-3.3554,-2.8010,-2.4750,-3.8210,-3.6401,-4.2,-4.8466
	};
	float sigma_DCB [6][21]={
	  1.1249,1.2,1.5,1.8439,2.6977,3.3445,4.3044,5.1506,5.8175,6.7200,7.9383,7.9957,9.6456,10.6918,11.4362,12.3759,12.7858,14.8817,16.6634,17.3,18.9136,
	  2.3557,2.5,2.8,3.5003,4.0150,4.7134,5.0622,5.2898,5.8487,6.0417,6.5309,7.0486,7.2767,7.4056,8.1837,8.4639,9.0455,9.9154,11.6088,12.0,12.3301,
	  1.7477,1.5,2.,2.7603,3.4052,3.9622,4.5838,5.2583,5.9739,6.5516,7.3605,7.8096,8.9956,9.6403,9.9994,10.8533,11.4922,11.9462,12.2392,12.5,12.9958,
	  1.,1.3,1.7,1.9129,2.6722,3.5879,4.3313,5.2615,5.5533,6.5403,7.6873,9.1676,11.1969,13.2381,13.7754,15.8128,18.2346,26.7205,34.5603,50.0,59.6160,
	  1.5,2.3,3.,3.4220,4.3056,4.4942,4.8157,5.4075,5.4011,5.9669,6.2477,7.3166,8.0090,8.3197,8.7633,9.4238,10.0955,14.0834,16.4086,19.0,22.0820,
	  1.,1.5,2.,2.6199,3.2837,3.9145,4.6630,5.2155,6.0262,6.0252,6.9880,8.4331,9.3758,9.6942,12.1845,13.0429,14.7118,22.7513,30.4955,40.0,48.6552
	};


        //RooRealVar *ZZMassCorr =new RooRealVar("ZZMassCorr","",-mass[i]/8,mass[i]/8);
        RooRealVar *ZZMassCorr=new RooRealVar("ZZMassCorr","",-100,100);
	RooRealVar *stagev=new RooRealVar("htxs_prodMode","",0,1000);
	RooRealVar *chan=new RooRealVar("chan","",0,1000);
	RooRealVar *vbfcate=new RooRealVar("vbfcate","",0,1000);
	RooRealVar *weight=new RooRealVar("weight","",1.e-9,1.e9);

	TString inputTrees = "/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All/ReducedTrees/";
	TString path = "/eos/user/a/axbuchot/SWAN_projects/HighMass/ZZ4lNutuples/All/";
	
	for (int j =0;j<3;j++){
	        ofstream outmass ;
	        //outmass.open(Form(path + "shapeParamNew/sim_massParam_prodMode%s_%d.txt",chanName[j].Data(),year));
		RooDataSet *data[N];
		//RooRealVar *ZZMassCorr[17];
		for (int i =0;i<N;i++){
		  //ofstream outmass ;
		  //outmass.open(Form(path + "shapeParamNew/sim_massParam_prodMode%s_%d.txt",chanName[j].Data(),year));

		  //RooRealVar *ZZMassCorr - *GenHMass = new RooRealVar("ZZMassCorr","",100,2*mass[i]);
    			TChain *t = new TChain ("SelectedTree");
			t->Add(Form(inputTrees + "*H%d_redTree_*%d*.root",mass[i],year));
			// Uncomment to use only TTrees for which we have all the samples (mass variated and tuneup/down)
		  // for(int l = 0; l < sigName.size(); l++) {
		  //   	TString name = sigName.at(l);
		  //	t->Add(Form(inputTrees+name+"%d_redTree_*%d*.root", mass[i],year));
		  //    cout << sigName.at(l) << " added to TChain" << endl;
    		  // }
		        //ZZMassCorr[i] =new RooRealVar("ZZMassCorr","",-mass[i]/8,mass[i]/8);
			data[i]=new RooDataSet(Form("data%d_%d",i,mass[i]),"",t,RooArgSet(*ZZMassCorr,*weight,*chan,*vbfcate),Form("chan==%d",j+1),"weight");
		}

		for (int k =0;k<nstage;k++){
			//	MH->setVal(mass[i]);
			//	MH->setConstant(kTRUE);
		        outmass.open(Form(path + "shapeParamCategory/sim_massParam_Cat%s_%s_%d.txt",stageName[k].Data(),chanName[j].Data(),year));

			RooCategory massrc("massrc","");
			RooSimultaneous simPdf("simPdf","simultaneous pdf",massrc) ;
			RooDataSet* data_stage[N]; 
			RooDataSet* data_stage_cat[N]; 
			RooDataSet* data_all;
			RooDoubleCBFast* cpdf[N]; 
			RooAddPdf* addpdf[N]; 
		              

			for (int i=0;i<N;i++){
				TString cat_type =Form("mh%d",mass[i]);
				massrc.defineType(cat_type,mass[i]);
			}
			/* RooRealVar *a1_1=new RooRealVar("a1_1","",0.1,-10.,10.); */
			/* RooRealVar *a1_2=new RooRealVar("a1_2","",0.01,-5.,5.); */
			/* RooRealVar *a1_0=new RooRealVar("a1_0","",2.,-50.,50.); */
			/* RooRealVar *a2_1=new RooRealVar("a2_1","",0.1,-10.,10.); */
			/* RooRealVar *a2_2=new RooRealVar("a2_2","",0.05,-5.,5.); */
			/* RooRealVar *a2_0=new RooRealVar("a2_0","",2.,-50.,50.); */
			/* RooRealVar *n1_1=new RooRealVar("n1_1","",0.01,-5.,5.); */
			/* RooRealVar *n1_2=new RooRealVar("n1_2","",0.01,-5.,5.); */
			/* RooRealVar *n1_0=new RooRealVar("n1_0","",2.,1.,50.); */
			/* RooRealVar *n2_1=new RooRealVar("n2_1","",0.01,-5.,5.); */
			/* RooRealVar *n2_2=new RooRealVar("n2_2","",-0.01,-5.,5.); */
			/* RooRealVar *n2_0=new RooRealVar("n2_0","",2.,1.,50.); */

			/* RooRealVar *mean_1=new RooRealVar("mean_1","",1.,-20.,20.); */
			/* RooRealVar *mean_2=new RooRealVar("mean_2","",0.01,-5.,5.); */
			/* RooRealVar *mean_0=new RooRealVar("mean_0","",3000,190,210); */

			/* RooRealVar *sigma_1=new RooRealVar("sigma_1","",1.,-20.,20.); */
			/* RooRealVar *sigma_2=new RooRealVar("sigma_2","",0.01,-10.,10.); */
			/* RooRealVar *sigma_0=new RooRealVar("sigma_0","",10.,0.,100.); */

			/* RooRealVar* mean_l_0 = new RooRealVar("landau_mean_0","",3000,180,220); */
			/* RooRealVar* mean_l_1 = new RooRealVar("landau_mean_1","",1.,-1.5,1.5); */
			/* RooRealVar* mean_l_2 = new RooRealVar("landau_mean_2","",0.,-1.5,1.5); */
			/* RooRealVar* sigma_l_0 = new RooRealVar("landau_sigma_0","",1.,0.,100.); */
			/* RooRealVar* sigma_l_1 = new RooRealVar("landau_sigma_1","",1.,-1.5,1.5); */
			/* RooRealVar* sigma_l_2 = new RooRealVar("landau_sigma_2","",0.,-1.5,1.5); */
			/* RooRealVar* frac_0 = new RooRealVar("frac_0","",0.5,0,1); */
			/* RooRealVar* frac_1 = new RooRealVar("frac_1","",-0.1,0.1); */
			/* RooRealVar* frac_2 = new RooRealVar("frac_2","",-0.1,0.1); */

			/* RooRealVar *a1_1=new RooRealVar("a1_1","",0.00111,-0.2,0.2); */
			/* RooRealVar *a1_2=new RooRealVar("a1_2","",0.00146,-0.1,0.1); */
			/* RooRealVar *a1_0=new RooRealVar("a1_0","",1.2917,1.,5.); */
			/* RooRealVar *a2_1=new RooRealVar("a2_1","",0.00130,-0.2,0.2); */
			/* RooRealVar *a2_2=new RooRealVar("a2_2","",0.00127,-0.1,0.1); */
			/* RooRealVar *a2_0=new RooRealVar("a2_0","",1.8125,1.,5.); */
			/* RooRealVar *n1_1=new RooRealVar("n1_1","",0.00138,-0.2,0.2); */
			/* RooRealVar *n1_2=new RooRealVar("n1_2","",0.00185,-0.1,0.1); */
			/* RooRealVar *n1_0=new RooRealVar("n1_0","",2.0312,1.,7.); */
			/* RooRealVar *n2_1=new RooRealVar("n2_1","",0.001,-0.2,0.2); */
			/* RooRealVar *n2_2=new RooRealVar("n2_2","",-0.001,-0.1,0.1); */
			/* RooRealVar *n2_0=new RooRealVar("n2_0","",3.0369,1.,5.); */

			/* RooRealVar *mean_1=new RooRealVar("mean_1","",1,0,1.5); */
			/* RooRealVar *mean_2=new RooRealVar("mean_2","",0.1,-1.,1.); */
			/* RooRealVar *mean_0=new RooRealVar("mean_0","",0.1,-100.,100); */

			/* RooRealVar *sigma_1=new RooRealVar("sigma_1","",0.01,-0.5,0.5); */
			/* RooRealVar *sigma_2=new RooRealVar("sigma_2","",0.001,-1.,1.); */
			/* RooRealVar *sigma_0=new RooRealVar("sigma_0","",5.,0.5,100.); */


			/* RooRealVar *a1_1=new RooRealVar("a1_1","",0.1,-1.,1.); */
			/* RooRealVar *a1_2=new RooRealVar("a1_2","",0.05,-0.1,0.1); */
			/* RooRealVar *a1_0=new RooRealVar("a1_0","",1.4,0.,10.); */
			/* RooRealVar *a2_1=new RooRealVar("a2_1","",0.1,-1.,1.); */
			/* RooRealVar *a2_2=new RooRealVar("a2_2","",0.05,-0.1,0.1); */
			/* RooRealVar *a2_0=new RooRealVar("a2_0","",1.2,0.,10.); */
			/* RooRealVar *n1_1=new RooRealVar("n1_1","",0.1,-1.,1.); */
			/* RooRealVar *n1_2=new RooRealVar("n1_2","",0.05,-0.1,0.1); */
			/* RooRealVar *n1_0=new RooRealVar("n1_0","",3.,0.,20.); */
			/* RooRealVar *n2_1=new RooRealVar("n2_1","",0.1,-1.,1.); */
			/* RooRealVar *n2_2=new RooRealVar("n2_2","",0.05,-0.1,0.1); */
			/* RooRealVar *n2_0=new RooRealVar("n2_0","",5.,0.,20.); */

			/* RooRealVar *mean_1=new RooRealVar("mean_1","",0.01,-10.,10.); */
			/* RooRealVar *mean_2=new RooRealVar("mean_2","",0.01,-1.,1.); */
			/* RooRealVar *mean_0=new RooRealVar("mean_0","",-1.,-20.,20.); */

			/* RooRealVar *sigma_1=new RooRealVar("sigma_1","",0.01,-5.,5.); */
			/* RooRealVar *sigma_2=new RooRealVar("sigma_2","",0.026,-1.,1.); */
			/* RooRealVar *sigma_0=new RooRealVar("sigma_0","",50.,0.,100); */

			/* RooRealVar* mean_l_0 = new RooRealVar("landau_mean_0","",130,110,140); */
			/* RooRealVar* mean_l_1 = new RooRealVar("landau_mean_1","",0,-1.5,1.5); */
			/* RooRealVar* mean_l_2 = new RooRealVar("landau_mean_2","",0,-1,1); */
			/* RooRealVar* sigma_l_0 = new RooRealVar("landau_sigma_0","",15,2,20); */
			/* RooRealVar* sigma_l_1 = new RooRealVar("landau_sigma_1","",0.,-1,1); */
			/* RooRealVar* sigma_l_2 = new RooRealVar("landau_sigma_2","",0,-1,1); */
			/* RooRealVar* frac_0 = new RooRealVar("frac_0","",0.65,0,1); */
			/* RooRealVar* frac_1 = new RooRealVar("frac_1","",-0.1,0.1); */
			/* RooRealVar* frac_2 = new RooRealVar("frac_2","",-0.1,0.1); */
			
			
			for (int i=0;i<N;i++){

			        int i2 = i;

				//RooRealVar *a1_1_mh=new RooRealVar("a1_1_mh","",0.1,-1.,1.);
				//RooRealVar *a1_2_mh=new RooRealVar("a1_2_mh","",0.05,-0.1,0.1);
				RooRealVar *a1_0_mh=new RooRealVar("a1_0_mh","",a1_DCB[k*3+j][i2],-3.,3.);
				//RooRealVar *a1_0_mh=new RooRealVar("a1_0_mh","",0,-3.,3.);
				//RooRealVar *a2_1_mh=new RooRealVar("a2_1_mh","",0.1,-1.,1.);
				//RooRealVar *a2_2_mh=new RooRealVar("a2_2_mh","",0.05,-0.1,0.1);
				RooRealVar *a2_0_mh=new RooRealVar("a2_0_mh","",a2_DCB[k*3+j][i2],-3.,3.);
				//RooRealVar *a2_0_mh=new RooRealVar("a2_0_mh","",0,-3.,3.);
				//RooRealVar *n1_1_mh=new RooRealVar("n1_1_mh","",0.1,-1.,1.);
				//RooRealVar *n1_2_mh=new RooRealVar("n1_2_mh","",0.05,-0.1,0.1);
				RooRealVar *n1_0_mh=new RooRealVar("n1_0_mh","",n1_DCB[k*3+j][i2],-5.,30.);
				//RooRealVar *n2_1_mh=new RooRealVar("n2_1_mh","",0.1,-1.,1.);
				//RooRealVar *n2_2_mh=new RooRealVar("n2_2_mh","",0.05,-0.1,0.1);
				RooRealVar *n2_0_mh=new RooRealVar("n2_0_mh","",n2_DCB[k*3+j][i2],-5.,30.);
				//RooRealVar *mean_1_mh=new RooRealVar("mean_1_mh","",0.01,-10.,10.);
				//RooRealVar *mean_2_mh=new RooRealVar("mean_2_mh","",0.01,-1.,1.);
				RooRealVar *mean_0_mh=new RooRealVar("mean_0_mh","",mean_DCB[k*3+j][i2],-mass[i]/20,mass[i]/20);
				//RooRealVar *sigma_1_mh=new RooRealVar("sigma_1_mh","",0.01,-5.,5.);
				//RooRealVar *sigma_2_mh=new RooRealVar("sigma_2_mh","",0.026,-1.,1.);
				RooRealVar *sigma_0_mh=new RooRealVar("sigma_0_mh","",sigma_DCB[k*3+j][i2],0.,100);


				RooConstVar *MH=new RooConstVar("MH","",mass[i]);
				// As of AN: the CB parameters are extracted by a fit p = a1+a11(mh-125)+a12(mh-125)**2
				// RooFormulaVar* a1=new RooFormulaVar(Form("a1_%d",mass[i]),"","@0+@1*(MH-125)+@2*(MH-125)*(MH-125)",RooArgList(*a1_0,*a1_1,*a1_2,*MH));
				// RooFormulaVar* a2=new RooFormulaVar(Form("a2_%d",mass[i]),"","@0+@1*(MH-125)+@2*(MH-125)*(MH-125)",RooArgList(*a2_0,*a2_1,*a2_2,*MH));
				// RooFormulaVar* n1=new RooFormulaVar(Form("n1_%d",mass[i]),"","@0+@1*(MH-125)+@2*(MH-125)*(MH-125)",RooArgList(*n1_0,*n1_1,*n1_2,*MH));
				// RooFormulaVar* n2=new RooFormulaVar(Form("n2_%d",mass[i]),"","@0+@1*(MH-125)+@2*(MH-125)*(MH-125)",RooArgList(*n2_0,*n2_1,*n2_2,*MH));
				RooFormulaVar* a1=new RooFormulaVar(Form("a1_%d",mass[i]),"","@0",RooArgList(*a1_0_mh));
				RooFormulaVar* a2=new RooFormulaVar(Form("a2_%d",mass[i]),"","@0",RooArgList(*a2_0_mh));
				RooFormulaVar* n1=new RooFormulaVar(Form("n1_%d",mass[i]),"","@0",RooArgList(*n1_0_mh));
				RooFormulaVar* n2=new RooFormulaVar(Form("n2_%d",mass[i]),"","@0",RooArgList(*n2_0_mh));
				RooFormulaVar* mean=new RooFormulaVar(Form("mean_%d",mass[i]),"","@0",RooArgList(*mean_0_mh));
				RooFormulaVar* sigma=new RooFormulaVar(Form("sigma_%d",mass[i]),"","@0",RooArgList(*sigma_0_mh));
				/* RooFormulaVar* sigma_l=new RooFormulaVar(Form("sigma_l_%d",mass[i]),"","@0+@1*(MH)",RooArgList(*sigma_l_0_mh,*sigma_l_1_mh,*MH)); */
				/* RooFormulaVar* mean_l=new RooFormulaVar(Form("mean_l_%d",mass[i]),"","@0+@1*(MH)",RooArgList(*mean_l_0_mh,*mean_l_1_mh,*MH)); */
				/* RooFormulaVar* frac=new RooFormulaVar(Form("frac_%d",mass[i]),"","@0",RooArgList(*frac_0_mh));*/
				/* if((stageName[k]=="ggH_VBF" && chanName[j] == "4e" && year == 2016) || stageName[k].Contains("BBH")|| (stageName[k].Contains("TH") && !stageName[k].Contains("TTH"))){ */
				/* 	a1=new RooFormulaVar(Form("a1_%d",mass[i]),"","@0",RooArgList(*a1_0_mh)); */
				/* 	a2=new RooFormulaVar(Form("a2_%d",mass[i]),"","@0",RooArgList(*a2_0_mh)); */
				/* 	n1=new RooFormulaVar(Form("n1_%d",mass[i]),"","@0",RooArgList(*n1_0_mh)); */
				/* 	n2=new RooFormulaVar(Form("n2_%d",mass[i]),"","@0",RooArgList(*n2_0_mh)); */
				/* 	// mean=new RooFormulaVar(Form("mean_%d",mass[i]),"","@0",RooArgList(*mean_0)); */
				/* 	// sigma=new RooFormulaVar(Form("sigma_%d",mass[i]),"","@0",RooArgList(*sigma_0)); */
				/* } */

				TString cat_type =Form("mh%d",mass[i]);
				data_stage[i] = (RooDataSet*)data[i]->reduce(Form("vbfcate==%d",stage[k])); 
				data_stage_cat[i] = new RooDataSet(Form("data_%d",i),"",RooArgSet(*ZZMassCorr,*weight),Index(massrc),Import(cat_type,*data_stage[i]),WeightVar("weight")); 
				cpdf[i] = new RooDoubleCBFast(Form("DCBall_%d",i),"",*ZZMassCorr,*mean,*sigma,*a1,*n1,*a2,*n2);
				/* RooLandau *pdf_landau = new RooLandau("landau","",*ZZMassCorr, *mean_l, *sigma_l); */
				/* addpdf[i] = new RooAddPdf(Form("apdf_%d",i),"",RooArgList(*cpdf[i],*pdf_landau),*frac); */
				cout<< data_stage_cat[i]->sumEntries()<<endl;
				/* if (stageName[k].Contains("TTH")|| stageName[k].Contains("VH_Lep") || stageName[k].Contains("WH_Lep") || stageName[k].Contains("ZH_Lep")) // WH, ZH commented */
				/*         simPdf.addPdf(*addpdf[i],cat_type); // Crystal Ball + Landau in TTH and VH to cope with the rising tail */
				//else
				simPdf.addPdf(*cpdf[i],cat_type); // Crystal Ball in all categories
				if(i==0)
					data_all = data_stage_cat[i];
				else
				        data_all->append(*data_stage_cat[i]);
				if(i==0){ // mH = 125 GeV
					// Quoting from analysis note: "The initial value for the parmsCB0 is obtained by fitting 125 GeV mass sample alone,
					// and refitted in the simultaneous fits. "
				  //RooDoubleCBFast *cpdftmp = new RooDoubleCBFast("DCBall","",*ZZMassCorr,*mean_0_mh,*sigma_0_mh,*a1_0_mh,*n1_0_mh,*a2_0_mh,*n2_0_mh);
				     RooDoubleCBFast *cpdftmp = new RooDoubleCBFast("DCBall","",*ZZMassCorr,*mean,*sigma,*a1,*n1,*a2,*n2);
					//RooLandau *pdf_landautmp = new RooLandau("landautmp","",*ZZMassCorr, *mean_l_0, *sigma_l_0);
					//RooAddPdf *addpdftmp = new RooAddPdf(Form("apdf_%d",i),"",RooArgList(*cpdftmp,*pdf_landautmp),*frac);
					//if (stageName[k].Contains("TTH")|| stageName[k].Contains("VH_Lep") || stageName[k].Contains("WH_Lep") || stageName[k].Contains("ZH_Lep"))
					//   cc     addpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));
				        //else
				     cpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));
					//addpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));
					

			

					//a1_0->setConstant(1);
					//a2_0->setConstant(1);
					//n1_0->setConstant(1);
					//n2_0->setConstant(1);
					// The CB and landau means, as well as CB sigma are set to constants for the simultaneous fit,
					//using the parameters obtained from 125GeV fit.
					//sigma_0->setConstant(1);
					//mean_0->setConstant(1);
					//mean_l_0->setConstant(1);
					//sigma_l_0->setConstant(1);
				     }

				/* TString a1outs= Form ("a1_TString a1outs= Form ("a1_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],a1_0_mh->getVal());
			     TString a2outs= Form ("a2_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],a2_0_mh->getVal());
			     TString n1outs= Form ("n1_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],n1_0_mh->getVal());
			     TString n2outs= Form ("n2_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],n2_0_mh->getVal());
			     TString meanouts= Form ("mean_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],mean_0_mh->getVal());
			     TString sigmaouts= Form ("sigma_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],sigma_0_mh->getVal());%s_%s_%d_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,mass[i],a1_0_mh->getVal(),a1_1_mh->getVal()); */
				/* TString a2outs= Form ("a2_%s_%s_%d_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,mass[i],a2_0_mh->getVal(),a2_1_mh->getVal()); */
				/* TString n1outs= Form ("n1_%s_%s_%d_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,mass[i],n1_0_mh->getVal(),n1_1_mh->getVal()); */
				/* TString n2outs= Form ("n2_%s_%s_%d_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,mass[i],n2_0_mh->getVal(),n2_1_mh->getVal()); */
				/* TString meanouts= Form ("mean_%s_%s_%d_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,mass[i],mean_0_mh->getVal(),mean_1_mh->getVal()); */
				/* TString sigmaouts= Form ("sigma_%s_%s_%d_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,mass[i],sigma_0_mh->getVal(),sigma_1_mh->getVal()); */
				
				/* outmass<< a1outs<<endl; */
				/* outmass<< a2outs<<endl; */
				/* outmass<< meanouts<<endl; */
				/* outmass<< sigmaouts<<endl; */
				/* outmass<< n1outs<<endl; */
				/* outmass<< n2outs<<endl; */
			
			     //for (int i=0;i<17;i++){

		 	     //RooRealVar *ZZMassCorr - *GenHMass = new RooRealVar("ZZMassCorr","",100,2*mass[i]);
			     simPdf.Print("v");
			     data_all->Print("v");
			     //			MH->setConstant(1);
			     simPdf.fitTo(*data_all,InitialHesse(true),Strategy(2)) ;
			     /* TString a1outs= Form ("a1_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],a1_0_mh->getVal()); */
			     /* TString a2outs= Form ("a2_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],a2_0_mh->getVal()); */
			     /* TString n1outs= Form ("n1_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],n1_0_mh->getVal()); */
			     /* TString n2outs= Form ("n2_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],n2_0_mh->getVal()); */
			     /* TString meanouts= Form ("mean_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],mean_0_mh->getVal()); */
			     /* TString sigmaouts= Form ("sigma_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],sigma_0_mh->getVal()); */

			     TString a1outs= Form ("a1_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],a1->getVal());
			     TString a2outs= Form ("a2_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],a2->getVal());
			     TString n1outs= Form ("n1_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],n1->getVal());
			     TString n2outs= Form ("n2_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],n2->getVal());
			     TString meanouts= Form ("mean_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],mean->getVal());
			     TString sigmaouts= Form ("sigma_%s_%s_%d_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mass[i],sigma->getVal());


			     /* TString a1outs= Form ("a1_%s_%s_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,a1_0_mh->getVal(),a1_1_mh->getVal()); */
			     /* TString a2outs= Form ("a2_%s_%s_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,a2_0_mh->getVal(),a2_1_mh->getVal()); */
			     /* TString n1outs= Form ("n1_%s_%s_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,n1_0_mh->getVal(),n1_1_mh->getVal()); */
			     /* TString n2outs= Form ("n2_%s_%s_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,n2_0_mh->getVal(),n2_1_mh->getVal()); */
			     /* TString meanouts= Form ("mean_%s_%s_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,mean_0_mh->getVal(),mean_1_mh->getVal()); */
			     /* TString sigmaouts= Form ("sigma_%s_%s_%d \t %.4f+%.5f",stageName[k].Data(),chanName[j].Data(),year,sigma_0_mh->getVal(),sigma_1_mh->getVal()); */
			     /* if((stageName[k]=="ggH_VBF" && chanName[j] == "4e" && year == 2016) || stageName[k].Contains("BBH")|| (stageName[k].Contains("TH") && !stageName[k].Contains("TTH"))){ */
			     /*      a1outs= Form ("a1_%s_%s_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,a1_0->getVal()); */
			     /*      a2outs= Form ("a2_%s_%s_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,a2_0->getVal()); */
			     /*      n1outs= Form ("n1_%s_%s_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,n1_0->getVal()); */
			     /*      n2outs= Form ("n2_%s_%s_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,n2_0->getVal()); */
			     /*      meanouts= Form ("mean_%s_%s_%d \t %.4f+%.5f*(MH)",stageName[k].Data(),chanName[j].Data(),year,mean_0->getVal(),mean_1->getVal()); */
			     /*      sigmaouts= Form ("sigma_%s_%s_%d \t %.4f+%.5f*(MH)",stageName[k].Data(),chanName[j].Data(),year,sigma_0->getVal(),sigma_1->getVal()); */
			     /*      //				meanouts= Form ("mean_%s_%s%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,mean_0->getVal()); */
			     /*      //	 			sigmaouts= Form ("sigma_%s_%s%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,sigma_0->getVal());*/
			     //}
			     
			     outmass<< a1outs<<endl;
			     outmass<< a2outs<<endl;
			     outmass<< meanouts<<endl;
			     outmass<< sigmaouts<<endl;
			     outmass<< n1outs<<endl;
			     outmass<< n2outs<<endl;
			     /* if (stageName[k].Contains("TTH")|| stageName[k].Contains("VH_Lep") || stageName[k].Contains("WH_Lep") || stageName[k].Contains("ZH_Lep")){ */
			     /*          outmass<<Form ("mean_l_%s_%s_%d \t %.4f+%.5f*(MH)",stageName[k].Data(),chanName[j].Data(),year,mean_l_0->getVal(),mean_l_1->getVal())<<endl; */
			     /*          outmass<<Form ("sigma_l_%s_%s_%d \t %.4f+%.5f*(MH)",stageName[k].Data(),chanName[j].Data(),year,sigma_l_0->getVal(),sigma_l_1->getVal())<<endl; */
			     /*          outmass<<Form ("frac_%s_%s_%d \t %.4f",stageName[k].Data(),chanName[j].Data(),year,frac_0->getVal())<<endl;*/
			     //}
			     }
			cout<<"where1"<<endl;
			for (int i =0;i<N;i++){
			        RooPlot *frame = ZZMassCorr->frame();
				TString cat_type =Form("mh%d",mass[i]);
				data_all->plotOn(frame,Cut(Form("massrc==massrc::mh%d",mass[i]))) ;
				simPdf.plotOn(frame,Slice(massrc,cat_type),ProjWData(massrc,*data_all)) ;
				frame->chiSquare();
				frame->residHist();
				frame->SetTitle(Form("%s -> %s (%d GeV)",stageName[k].Data(),chanName[j].Data(),mass[i]));
				frame->GetYaxis()->SetTitle("Events");
				frame->GetXaxis()->SetTitle("m_{4l}^{reco} - m_{4l}^{gen} (GeV)");
				frame->Draw();
				gPad->Print(Form(path + "shapeCategory/fig_%d/simFit_%s_%d_%s_%d.png", year, stageName[k].Data(),mass[i],chanName[j].Data(),year));
			}
		}
		outmass.close();
	}
}
