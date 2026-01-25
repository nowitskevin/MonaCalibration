// Position calibration for all bars run 5200. 
//Kevin Eisenberg

#include <fstream>

void posCalibLoop(){
	gROOT->SetBatch();

	TFile *f = new TFile("/mnt/analysis/e23033/analysis/kevin/rootfiles/run5200.root");
    	TTree *t = (TTree*)f->Get("t");
	//output file
	ofstream csvFile("/mnt/analysis/e23033/analysis/kevin/outputs/posParamsAfterPulser.csv");
	csvFile << "x,y,a,b,xL,xR\n";
	
	//read file

	UShort_t time[2][18][16]; 
	t->SetBranchAddress("time[2][18][16]", time); 
	Long64_t nEntries = t->GetEntries();

	//read pulser params
	TTree *params = new TTree("params","params");
	params->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/pulserParams.csv","x/I:y/I:z/I:slope/D"); // PARAMETERS TO CONVERT FROM CHANNEL TO ns

	Long64_t nParams = params->GetEntries();

	int X; int Y; int Z; double slope;
	params->SetBranchAddress("x",&X);
	params->SetBranchAddress("y",&Y);
	params->SetBranchAddress("z",&Z);
	params->SetBranchAddress("slope",&slope);

	const int nWalls = 9;
	const int nVert  = 16;
	const int nSides = 2;
	const int nBars = nWalls*nVert;
	const int nPMTs  = nWalls * nVert * nSides;

	//MAKE DICTIONARY OF SLOPES
	vector<double> slopes(nPMTs);
	for (int i=0; i<nParams; i++) {
		params->GetEntry(i);
		int ind = X*nVert*nSides + Y*nSides + Z;
		slopes[ind] = slope;
		}
	
	//declare hists
	vector<TH1F*> hRaw(nBars);
	vector<TH1F*> hCal(nBars);
	cout << "Beginning position calibrations..." << endl;

	for (int x=0; x<nWalls; x++){
		for (int y=0; y<nVert; y++){
			int Index = x*nVert + y;
			hRaw[Index] = new TH1F(Form("hRaw[%d][%d]",x,y),"raw time diff",200,-40,40);
			hCal[Index] = new TH1F(Form("hCal[%d][%d]",x,y),"cal position:x (cm):counts",100,-200,200);
		}
	}

	cout << "Histograms declared" << endl;

	//raw data (don't need all entries)
	double dt;
	
	for (int i=0; i<5000000; i++){
		t->GetEntry(i); 
		for (int x=0; x<nWalls; x++){
			for (int y=0; y<nVert; y++){
				int Index = x*nVert + y;
				int leftSlopeIndex = x*nVert*nSides + y*nSides;
				int rightSlopeIndex = leftSlopeIndex+1;
				if (time[0][x][y]>0 && time[1][x][y]>0){
					dt = static_cast<double>(slopes[leftSlopeIndex]*time[0][x][y]) - static_cast<double>(slopes[rightSlopeIndex]*time[1][x][y]);
					hRaw[Index]->Fill( dt );
				}
			}
		}
	}
	cout << "Raw histograms filled" << endl;
	
	
	TF1 *fRaw = new TF1("fRaw","[p0]*(tanh((x-[p1])/[p2])-tanh((x-[p3])/[p4]))");
	TF1 *fCal = new TF1("fCal","[p0]*(tanh((x-[p1])/[p2])-tanh((x-[p3])/[p4]))");
		
	double fitmin = -40.0;
	double fitmax = 40.0;
	double amplitude,xL,tauL,xR,tauR;
	//fit parameters
	vector<double> a(nBars);
	vector<double> b(nBars);

	for (int i=0; i<nBars; i++){
		int x = i/nVert;
		int y = i%nVert;

		TH1F* h1 = hRaw[i];
		fRaw->SetParameters(h1->GetMaximum(),-11,10,14,10); //can change parameters 1 and 3 to troubleshoot. Usually works between (-)11 and (-)15
		cout << "fitting bar " << x << ", " << y << endl;
		h1->Fit("fRaw","R","",fitmin,fitmax);

		amplitude = fRaw->GetParameter(0);
		xL = fRaw->GetParameter(1);
		tauL = fRaw->GetParameter(2);
		xR = fRaw->GetParameter(3);
		tauR = fRaw->GetParameter(4);

		a[i] = 200/(xR-xL);
		b[i] = -100 - 200*xL/(xR-xL);
		
	}
	for (int i=0; i<1000000; i++){
		t->GetEntry(i); 
		for (int x=0; x<nWalls; x++){
			for (int y=0; y<nVert; y++){
				int Index = x*nVert + y;
				int leftSlopeIndex = x*nVert*nSides + y*nSides;
				int rightSlopeIndex = leftSlopeIndex+1;
				if (time[0][x][y]>0 && time[1][x][y]>0){
					dt = static_cast<double>(slopes[leftSlopeIndex]*time[0][x][y]) - static_cast<double>(slopes[rightSlopeIndex]*time[1][x][y]);
					hCal[Index]->Fill(a[Index]*dt+b[Index]);
				}
			}
		}
	}

	double calxL,calxR;
	for (int i=0; i<nBars; i++){
		int x = i/nVert;
		int y = i%nVert;

		TH1F* h2 = hCal[i];
		fCal->SetParameters(h2->GetMaximum(),-100,10,100,10);
		cout << "fitting bar " << x << ", " << y << endl;
		h2->Fit("fCal","RQ","",-200,200);

		calxL = fCal->GetParameter(1);
		calxR = fCal->GetParameter(3);
		
		csvFile << x << "," << y << "," << a[i] << "," << b[i] << "," << calxL << "," << calxR << "\n";
	}
	
	
csvFile.close();
	
}
