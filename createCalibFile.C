// This file creates a MoNA-calibrated ROOT file out of the raw data.
// Kevin Eisenberg, 2026

void createCalibFile(string runNumber){
	const int nWalls = 9;
	const int nVert  = 16;
	const int nSides = 2;
	//open raw file
	TFile *fRaw = TFile::Open(Form("/mnt/analysis/e23033/analysis/kevin/rootfiles/run%s.root",runNumber.c_str()));
	TTree *tRaw = (TTree*)fRaw->Get("t");
	cout << Form("Writing MoNA-calibrated file of run%s",runNumber.c_str()) << endl;
	
	//open pulser params
	TTree *pulserparams = new TTree("pulserparams","pulserparams");
	pulserparams->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/pulserParams.csv","x/I:y/I:z/I:slope/D");

	//open position params
	TTree *posparams = new TTree("posparams","posparams");
	posparams->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/posParamsAfterPulser.csv","x/I:y/I:a/D:b/D:xL/D:xR/D");

	//open charge params
	TTree *chargeparams = new TTree("chargeparams","chargeparams");
	chargeparams->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/chargeCalibParams.csv","x/I:y/I:z/I:a/D:b/D:peak/D");

	//assign branch addresses
	Long64_t nEvents = tRaw->GetEntries();
	
	Long64_t nPulsers = pulserparams->GetEntries();
	int xPuls; int yPuls; int zPuls; double slopePuls;
	pulserparams->SetBranchAddress("x",&xPuls);
	pulserparams->SetBranchAddress("y",&yPuls);
	pulserparams->SetBranchAddress("z",&zPuls);
	pulserparams->SetBranchAddress("slope",&slopePuls);

	Long64_t nPositions = posparams->GetEntries();
	int xPos; int yPos; double aPos; double bPos;
	posparams->SetBranchAddress("x",&xPos);
	posparams->SetBranchAddress("y",&yPos);
	posparams->SetBranchAddress("a",&aPos);
	posparams->SetBranchAddress("b",&bPos);

	Long64_t nCharged = chargeparams->GetEntries();
	int xQ; int yQ; int zQ; double aQ; double bQ;
	chargeparams->SetBranchAddress("x",&xQ);
	chargeparams->SetBranchAddress("y",&yQ);
	chargeparams->SetBranchAddress("z",&zQ);
	chargeparams->SetBranchAddress("a",&aQ);
	chargeparams->SetBranchAddress("b",&bQ);
	
	double pulserSlopes[nPulsers];
	double posSlopes[nPositions];
	double posIntercepts[nPositions];
	double chargeSlopes[nCharged];
	double chargeIntercepts[nCharged];
	//organize parameters into arrays
	for (int i=0; i<nPulsers; i++){
		pulserparams->GetEntry(i);
		int index = xPuls*nVert*nSides + yPuls*nSides + zPuls;
		pulserSlopes[index] = slopePuls;
	}
	for (int i=0; i<nPositions; i++){
		posparams->GetEntry(i);
		int index = xPos*nVert + yPos;
		posSlopes[index] = aPos;
		posIntercepts[index] = bPos;

	}
	for (int i=0; i<nCharged; i++){
		chargeparams->GetEntry(i);
		int index = xQ*nVert*nSides + yQ*nSides + zQ;
		chargeSlopes[index]=aQ;
		chargeIntercepts[index]=bQ;
	}

	//create new file
	cout << "creating output calibrated file" << endl;
	TFile *fCal = new TFile(Form("/mnt/analysis/e23033/analysis/kevin/rootfiles/run%s_monaCalibrated.root",runNumber.c_str()),"RECREATE");
	TTree *tCal = tRaw->CloneTree(-1,"fast");
	
	UShort_t time[2][18][16]; UShort_t light[2][18][16];
	tCal->SetBranchAddress("time[2][18][16]",time);
	tCal->SetBranchAddress("light[2][18][16]",light);

	//set new branch addresses
	double timecal[2][18][16]; //in ns
	double position[18][16]; //in cm
	double charge[2][18][16]; //in MeVee
	int eventNumber; //for tracking

	TBranch *b_Time = tCal->Branch("timecal",timecal,"timecal[2][18][16]/D");
	TBranch *b_Pos = tCal->Branch("position",position,"position[18][16]/D");
	TBranch *b_Charge = tCal->Branch("charge",charge,"charge[2][18][16]/D");
	TBranch *b_I = tCal->Branch("eventNumber",&eventNumber,"eventNumber/I");

	//main loop
	cout << "beginning calibrations. Events calibrated:" << endl;
	int fullIndex; int truncIndex;
	for (int i=0; i<nEvents; i++){
		tCal->GetEntry(i);
		if (i%10000 == 0){cout << "\r" << i << flush;}
		
		for (int xi=0; xi<nWalls; xi++){
			for (int yi=0; yi<nVert; yi++){
				for (int zi=0; zi<nSides; zi++){
					fullIndex = xi*nVert*nSides + yi*nSides + zi;
					truncIndex = xi*nVert + yi;
					
					timecal[zi][xi][yi] = pulserSlopes[fullIndex]*time[zi][xi][yi];
					if (zi==1){
						position[xi][yi] = posSlopes[truncIndex]*(timecal[zi-1][xi][yi]-timecal[zi][xi][yi])+posIntercepts[truncIndex];
					}

					charge[zi][xi][yi] = chargeSlopes[fullIndex]*light[zi][xi][yi]+chargeIntercepts[fullIndex];
					
				}
			}
		}
		eventNumber = i;

	b_Time->Fill();
	b_Pos->Fill();
	b_Charge->Fill();
	b_I->Fill();
	}
	cout << endl << "done. Closing files." << endl;
	//close files
	tCal->Write();
	fCal->Close();
	fRaw->Close();
	
}
