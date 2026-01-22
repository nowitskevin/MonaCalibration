//Kevin Eisenberg, 2025
//APPLY the position calibration parameters to VISUALIZE what one wall of mona looks like in a particular run. What you should see is a vertical square between -100 and +100 cm, with light in it. 
//2026: Currently, due to trigger anomalies, in-beam data will look zigzagged, while self-trigger cosmics look perfect.

void posCalibrateWall(int X, string runNumber){


	TTree *params = new TTree("params","params");
	params->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/positionParamsFull.csv","x/I:y/I:a/D:b/D:xL/D:xR/D");

	int x; int y; double a; double b; double xL; double xR;

	params->SetBranchAddress("x",&x);
	params->SetBranchAddress("y",&y);
	params->SetBranchAddress("a",&a);
	params->SetBranchAddress("b",&b);
	params->SetBranchAddress("xL",&xL);
	params->SetBranchAddress("xR",&xR);
	
	Long64_t nCalibrated = params->GetEntries();
	
	double fitA[16]; double fitB[16];
	for (int i = 0; i<16; i++){
		for (int j = 0; j<nCalibrated; j++){
			params->GetEntry(j);
			if (x==X && y==i){
				fitA[i] = a;
				fitB[i] = b;
				cout << a << " " << b << endl;
			}
		}
	}
	
	TFile *f = new TFile(Form("/mnt/analysis/e23033/analysis/kevin/rootfiles/run%s.root",runNumber.c_str()),"READ");

	TTree *t = (TTree*)f->Get("t");
	UShort_t time[2][18][16];

	if (runNumber == "5227" || runNumber == "5229"){
	t->SetBranchAddress("t[2][18][16]",time);
	}
	else {
	t->SetBranchAddress("time[2][18][16]",time);
	}
	TH2F *h2 = new TH2F("h2",Form("Wall %d position map run%s;Position (cm);Bar index",X,runNumber.c_str()),
                         200, -200, 200,
                         16, 0.0, 16);

	Long64_t nEvents = t->GetEntries();

	for (int i=0; i<nEvents; i++){
		t->GetEntry(i);

		for (int bar = 0; bar < 16; bar++){
		    int yy = 15 - bar;   

		    if (time[0][X][yy] == 0 || time[1][X][yy] == 0) continue;

		    double dt = time[0][X][yy] - time[1][X][yy];
		    double pos = fitA[yy]*dt + fitB[yy];
			if (dt > 4096 || dt < -4096) {throw runtime_error("Strange value encountered"); break;}

		    h2->Fill(pos, yy);
		}
	}

	TCanvas *C = new TCanvas("C","Wall Position",1200,800);
	gPad->SetLogz();
	h2->Draw("COLZ");
	
		
}
