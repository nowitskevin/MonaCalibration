// This file creates a MoNA-calibrated ROOT file out of the raw data.
// Kevin Eisenberg, 2026
#include <vector>
#include <cmath>

void createCalibFile(int runNumber){
	const int nWalls = 9;
	const int nVert  = 16;
	const int nSides = 2;
	const double speedOfLight = 0.299792458; // m/ns
	const double barThickness = 0.10; // m
	const double closestMidDistance = 8.066 + (0.5 * barThickness); // m, middle of bar
	//open raw file
	TFile *fRaw = TFile::Open(Form("/mnt/analysis/e23033/analysis/kevin/uncalRootFiles/run%d_reglom.root",runNumber));
	TTree *tRaw = (TTree*)fRaw->Get("t");
	cout << Form("Writing MoNA-calibrated file of run%d",runNumber) << endl;
	
	//open pulser params
	TTree *pulserparams = new TTree("pulserparams","pulserparams");
	pulserparams->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/pulserParams.csv","x/I:y/I:z/I:slope/D");

	//open position params
	TTree *posparams = new TTree("posparams","posparams");
	posparams->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/posParamsAfterPulser.csv","x/I:y/I:a/D:b/D:xL/D:xR/D");

	//open charge params
	TTree *chargeparams = new TTree("chargeparams","chargeparams");
	chargeparams->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/chargeCalibParams.csv","x/I:y/I:z/I:a/D:b/D:peak/D");

	//open CFD params
	TTree *cfdparams = new TTree("cfdparams","cfdparams");
	cfdparams->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/cfdWalkParams.csv","x/I:y/I:amplitude/D:slope/D:intercept/D:meanY/D:fitMin/D:fitMax/D");
	
	//open wall offsets 
	TTree *wallOffsets = new TTree("wallOffsets","wallOffsets");
	wallOffsets->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/wallOffsetsMidBar.csv","x/I:y/I:offset/D");

	//open gamma peak locations
	TTree *gammaPeaks = new TTree("gammaPeaks","gammaPeaks");
	gammaPeaks->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/gammaPeaks.csv","x/I:y/I:peak/D");

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
	
	Long64_t nParams = cfdparams->GetEntries();
	int xP; int yP; double amp; double slope; double intercept; double mean; double fitMin; double fitMax;
	cfdparams->SetBranchAddress("x",&xP);
	cfdparams->SetBranchAddress("y",&yP);
	cfdparams->SetBranchAddress("amplitude",&amp);
	cfdparams->SetBranchAddress("slope",&slope);
	cfdparams->SetBranchAddress("intercept",&intercept);
	cfdparams->SetBranchAddress("meanY",&mean);
	cfdparams->SetBranchAddress("fitMin",&fitMin);
	cfdparams->SetBranchAddress("fitMax",&fitMax);

	Long64_t nOffsets = wallOffsets->GetEntries();
	int oX; int oY; double oVal;
	wallOffsets->SetBranchAddress("x", &oX);
	wallOffsets->SetBranchAddress("y", &oY);
	wallOffsets->SetBranchAddress("offset", &oVal);

	Long64_t nPeaks = gammaPeaks->GetEntries();
	int gX; int gY; double gPeak;
	gammaPeaks->SetBranchAddress("x", &gX);
	gammaPeaks->SetBranchAddress("y", &gY);
	gammaPeaks->SetBranchAddress("peak", &gPeak);

	//setup parameter arrays
	double pulserSlopes[nPulsers];
	double posSlopes[nPositions];
	double posIntercepts[nPositions];
	double chargeSlopes[nCharged];
	double chargeIntercepts[nCharged];
	double amps[nParams]; double slopes[nParams]; double ints[nParams];
	double means[nParams]; double fitMins[nParams]; double fitMaxs[nParams];
	vector<double> offsets(nWalls * nVert, 0.0);
	vector<double> peaks(nWalls, 0.0);

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
	for (int i=0; i<nParams; i++){
		cfdparams->GetEntry(i);
		int index = xP*nVert + yP;
		amps[index] = amp; slopes[index] = slope; ints[index] = intercept;
		means[index] = mean; fitMins[index] = fitMin; fitMaxs[index] = fitMax;
	}
	for (int i = 0; i<nOffsets; i++){
		wallOffsets->GetEntry(i);
		offsets[oX * nVert + oY] = oVal;
	}
	for (int i=0; i<nPeaks; i++){
		gammaPeaks->GetEntry(i);
		peaks[gX] = gPeak;
	}	
	//calculating offsets
	vector<double> absoluteOffset(nWalls * nVert, 0.0);
	for (int w = 0; w < nWalls; w++){
		int    gammaBar  = (w != 8) ? 8 : 9;
		double wallZ     = closestMidDistance + w * barThickness;
		double gammaBarY = (gammaBar - 8) * barThickness;
		double gammaToF0 = sqrt(wallZ*wallZ + gammaBarY*gammaBarY) / speedOfLight;
		for (int y=0; y<nVert; y++){
			int idx = w * nVert + y;
			absoluteOffset[idx] = peaks[w] + gammaToF0 + offsets[idx];
	    	}
	}

	//create new file
	cout << "creating output calibrated file" << endl;
	TFile *fCal = new TFile(Form("/mnt/analysis/e23033/analysis/kevin/rootfiles/run%d_MonaFull_FILENAME.root",runNumber),"RECREATE");
	TTree *tCal = tRaw->CloneTree(0);
	TTree *tHits = new TTree("tHits","Nonzero hits");

	UShort_t time[2][18][16]; UShort_t light[2][18][16];
	tRaw->SetBranchAddress("time[2][18][16]",time);
	tRaw->SetBranchAddress("light[2][18][16]",light);

	//set new branch addresses
	std::vector<int> vec_x, vec_y;
	std::vector<double> vec_tAvg, vec_tAvgCFDcorr, vec_qAvg, vec_pos;
	std::vector<double> vec_tL, vec_tR, vec_qL, vec_qR;
	std::vector<double> vec_tof, vec_beta;
	double normEvtNum;
	int evtNum;
	int runNum = runNumber;
	int nHits;  // multiplicity - handy to have

	tHits->Branch("x", &vec_x);
	tHits->Branch("y", &vec_y);
	tHits->Branch("tAvg", &vec_tAvg);
	tHits->Branch("tAvgCFDcorr",&vec_tAvgCFDcorr);
	tHits->Branch("tof",&vec_tof);
	tHits->Branch("beta",&vec_beta);
	tHits->Branch("qAvg", &vec_qAvg);
	tHits->Branch("pos", &vec_pos);
	tHits->Branch("tL", &vec_tL);
	tHits->Branch("tR", &vec_tR);
	tHits->Branch("qL", &vec_qL);
	tHits->Branch("qR", &vec_qR);
	tHits->Branch("nHits", &nHits, "nHits/I");
	tHits->Branch("normEvtNum", &normEvtNum, "normEvtNum/D");
	tHits->Branch("evtNum", &evtNum, "evtNum/I");
	tHits->Branch("runNum", &runNum, "runNum/I");
		
	//main loop
	cout << "beginning calibrations. Events calibrated:" << endl;
	int leftIndex, rightIndex; int truncIndex;
	cout << "Events: " << nEvents << endl;

	
	for (Long64_t i=0; i<nEvents; i++){//loop
		tRaw->GetEntry(i);

		vec_x.clear(); vec_y.clear();
    		vec_tAvg.clear(); vec_qAvg.clear(); vec_pos.clear(); vec_tAvgCFDcorr.clear();
    		vec_tL.clear(); vec_tR.clear();
    		vec_qL.clear(); vec_qR.clear();
		vec_tof.clear(); vec_beta.clear();

		if (i%10000 == 0){cout << "\r" << i << flush;}
		for (int xi=0; xi<nWalls; xi++){//x
			for (int yi=0; yi<nVert; yi++){//y
				if ((time[0][xi][yi] > 0) && (time[1][xi][yi] > 0) && (light[0][xi][yi] > 0) && (light[1][xi][yi] > 0)){//if
	
					leftIndex = xi*nVert*nSides + yi*nSides;
					rightIndex = leftIndex + 1;
					truncIndex = xi*nVert + yi;

					double tL = pulserSlopes[leftIndex]*time[0][xi][yi];
					double tR = pulserSlopes[rightIndex]*time[1][xi][yi];
					double qL = chargeSlopes[leftIndex]*light[0][xi][yi]+chargeIntercepts[leftIndex];
					double qR = chargeSlopes[rightIndex]*light[1][xi][yi]+chargeIntercepts[rightIndex];
					
					double position = posSlopes[truncIndex]*(tL - tR)+posIntercepts[truncIndex];
					double tAvg = (tL+tR)/2.0;
					
					if (qL>0 && qR>0){
						vec_x.push_back(xi);
						vec_y.push_back(yi);
						vec_tL.push_back(tL);
						vec_tR.push_back(tR);
						vec_qL.push_back(qL);
						vec_qR.push_back(qR);
						vec_tAvg.push_back(tAvg);
						double qAvg = sqrt(qL*qR);
						double tAvgCFDcorr;
						if (qAvg >= fitMins[truncIndex] && qAvg <= fitMaxs[truncIndex]){
							double CFDcorrection = means[truncIndex] - (amps[truncIndex] * TMath::Exp(-slopes[truncIndex] / qAvg) + ints[truncIndex]);
							tAvgCFDcorr = tAvg + CFDcorrection;
						}
						else {tAvgCFDcorr = tAvg;}
						vec_qAvg.push_back(qAvg);
						vec_tAvgCFDcorr.push_back(tAvgCFDcorr);
						vec_pos.push_back(position);
						double hitX = position /100.0;//cm to m
						double hitY = (yi-8) * barThickness; // m
						double wallZ = closestMidDistance + xi*barThickness; // m (could cache above)
						double hitDist = sqrt(wallZ*wallZ + hitY*hitY + hitX*hitX);
						double tof = -tAvgCFDcorr + absoluteOffset[xi*nVert+yi];
						double beta = (tof > 0) ? (hitDist/tof)/speedOfLight : -999.0;
						vec_tof.push_back(tof);
						vec_beta.push_back(beta);
											
					}//if q>0	
				}//if t&q>0
			}//y
		}//x
		nHits = max(vec_x.size(),vec_y.size());
		if (nHits>0){
			normEvtNum = (double(i)*1000000.0)/nEvents;
			evtNum = int(i);
			tHits->Fill();
		}
		tCal->Fill();
	}//loop
	cout << endl << Form("done with run%d. Closing files.",runNumber) << endl;
	//close files
	
	tCal->Write();
	tHits->Write();
	fCal->Close();
	fRaw->Close();
	
}
