#include <vector>
//Kevin Eisenberg, April 2026. Correcting CFD walk

void cfdCalib(int x, int y, vector<int> runs){

	
	TChain *chain = new TChain("tHits");
	for (const auto& runNumber : runs){
		chain->Add(Form("/mnt/analysis/e23033/analysis/kevin/rootfiles/run%d_MonaCal.root",runNumber));
	}
	
	Long64_t nEvents = chain->GetEntries();

	chain->Draw("tAvg:qAvg>>hRaw(250,0,20,150,260,310)",Form("x==%d && y==%d && Sum$(y==15)>0",x,y),"COLZ");
	TH2F *hRaw = (TH2F*)gDirectory->Get("hRaw");
	TCanvas *cRaw = new TCanvas("cRaw","raw",800,600);
	hRaw->Draw("COLZ");
	TProfile *pf = hRaw->ProfileX();

	double meanY = pf->GetMean(2);
	
	TF1 *fit = new TF1("fit","[0]*exp(-[1]/x)+[2]");
	fit->SetParameters(285,0.5,2);
	double fitMin=0.6; double fitMax=15.5;
	
	pf->Fit("fit","RQ","",fitMin,fitMax);
	pf->SetMinimum(270);
	pf->SetMaximum(300);
	TCanvas *c2 = new TCanvas("c2","profile fit",800,600);
	pf->Draw();

	double height;
	double slope;
	double intercept;
	//double inverse;

	height = fit->GetParameter(0);
	slope = fit->GetParameter(1);
	intercept = fit->GetParameter(2);
	//inverse = fit->GetParameter(3);
	cout << "x," << "y," << "height," << "slope," << "intercept," << "meanY," << "fitMin," << "fitMax" << endl;
	cout << x << "," << y << "," << height << "," << slope << "," << intercept << "," << meanY << "," << fitMin << "," << fitMax << endl;
	//cout << inverse << endl;

	chain->Draw(Form("(tAvg + (%f - (%f*exp(-%f/qAvg)+%f))):qAvg>>h2(250,0,20,150,260,310)",meanY, height, slope, intercept),Form("x==%d && y==%d && qAvg>%f && Sum$(y==15)>0",x,y,fitMin),"COLZ");
	TH2F *hCal = (TH2F*)gDirectory->Get("h2");
	TProfile *pfC = hCal->ProfileX();
	TCanvas *c3 = new TCanvas("c3","cal profile",800,600);
	pfC->SetMinimum(270);
	pfC->SetMaximum(300);
	pfC->Draw();

}
