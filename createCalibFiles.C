// This file creates a MoNA-calibrated ROOT file out of the raw data.
// Kevin Eisenberg, 2026
#include <vector>
#include <cmath>
#include <numeric>

void createCalibFile(int runNumber){
	const int nWalls = 9;
	const int nVert  = 16;
	const int nSides = 2;
	const double speedOfLight = 0.299792458; // m/ns
	const double barThickness = 0.10; // m
	const double closestMidDistance = 8.066 + (0.5 * barThickness); // m, middle of bar
	const double neutronMass = 939.56542; //MeV/c^2
	const double pi = 3.1415926535897932384626433;	
	
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
	TFile *fCal = new TFile(Form("/mnt/analysis/e23033/analysis/kevin/rootfiles/run%d_MonaFull.root",runNumber),"RECREATE");
	TTree *tMona = new TTree("tMona","Nonzero hits");

	UShort_t time[2][18][16]; UShort_t light[2][18][16];
	tRaw->SetBranchAddress("time[2][18][16]",time);
	tRaw->SetBranchAddress("light[2][18][16]",light);

	//set new branch addresses
	std::vector<int> vec_x, vec_y;
	//std::vector<double> vec_tAvg, 
	std::vector<double> vec_tAvgCFDcorr, vec_qAvg, vec_Xpos, vec_Ypos, vec_Zpos, vec_dist;
	//std::vector<double> vec_tL, vec_tR, vec_qL, vec_qR;
	std::vector<double> vec_tof, vec_beta, vec_gamma;
	std::vector<double> vec_E, vec_p, vec_px, vec_py, vec_pz, vec_th, vec_ph;
	double normEvtNum;
	int evtNum;
	int runNum = runNumber;
	int nHits;  // multiplicity - handy to have

	tMona->Branch("xInd", &vec_x); //index of bar wall
	tMona->Branch("yInd", &vec_y); //index of bar height
	//tMona->Branch("tAvg", &vec_tAvg);
	//tMona->Branch("tAvgCFDcorr",&vec_tAvgCFDcorr);
	tMona->Branch("tof",&vec_tof);
	tMona->Branch("beta",&vec_beta);
	tMona->Branch("gamma",&vec_gamma);
	tMona->Branch("qAvg", &vec_qAvg);
	tMona->Branch("Xpos", &vec_Xpos); //geometric x position along the bar (cm)
	tMona->Branch("Ypos",&vec_Ypos); //geometric y position vertically (cm) 
	tMona->Branch("Zpos",&vec_Zpos); //geometric z position along the beam axis (cm)
	tMona->Branch("dist",&vec_dist); //flight path
	tMona->Branch("E",&vec_E);
	tMona->Branch("p",&vec_p);
	tMona->Branch("px",&vec_px);
	tMona->Branch("py",&vec_py);
	tMona->Branch("pz",&vec_pz);
	tMona->Branch("theta",&vec_th);
	tMona->Branch("phi",&vec_ph);
	//tMona->Branch("tL", &vec_tL);
	//tMona->Branch("tR", &vec_tR);
	//tMona->Branch("qL", &vec_qL);
	//tMona->Branch("qR", &vec_qR);
	tMona->Branch("nHits", &nHits, "nHits/I");
	tMona->Branch("normEvtNum", &normEvtNum, "normEvtNum/D");
	tMona->Branch("evtNum", &evtNum, "evtNum/I");
	tMona->Branch("runNum", &runNum, "runNum/I");
		
	//main loop
	cout << "beginning calibrations. Events calibrated:" << endl;
	int leftIndex, rightIndex; int truncIndex;
	cout << "Events: " << nEvents << endl;

	
	for (Long64_t i=0; i<nEvents; i++){//loop
		tRaw->GetEntry(i);

		vec_x.clear(); vec_y.clear();
    		//vec_tAvg.clear(); 
		vec_qAvg.clear(); 
		vec_Xpos.clear(); vec_tAvgCFDcorr.clear(); vec_Ypos.clear(); vec_Zpos.clear(); 			vec_dist.clear();
    		//vec_tL.clear(); vec_tR.clear();
    		//vec_qL.clear(); vec_qR.clear();
		vec_tof.clear(); vec_beta.clear(); vec_gamma.clear();
		vec_E.clear(); vec_p.clear(); vec_px.clear(); vec_py.clear(); vec_pz.clear(); vec_th.clear(); vec_ph.clear();

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
						//vec_tL.push_back(tL);
						//vec_tR.push_back(tR);
						//vec_qL.push_back(qL);
						//vec_qR.push_back(qR);
						//vec_tAvg.push_back(tAvg);
						double qAvg = sqrt(qL*qR);
						double tAvgCFDcorr;
						if (qAvg >= fitMins[truncIndex] && qAvg <= fitMaxs[truncIndex]){
							double CFDcorrection = means[truncIndex] - (amps[truncIndex] * TMath::Exp(-slopes[truncIndex] / qAvg) + ints[truncIndex]);
							tAvgCFDcorr = tAvg + CFDcorrection;
						}
						else {tAvgCFDcorr = tAvg;}
						vec_qAvg.push_back(qAvg);
						vec_tAvgCFDcorr.push_back(tAvgCFDcorr);
						vec_Xpos.push_back(position); 
						double hitX = position /100.0;//cm to m
						double hitY = (yi-8) * barThickness; // m
						double wallZ = closestMidDistance + xi*barThickness; // m (could cache above)
						vec_Ypos.push_back(hitY*100.0);
						vec_Zpos.push_back(wallZ*100.0);
						double hitDist = sqrt(wallZ*wallZ + hitY*hitY + hitX*hitX);
						vec_dist.push_back(hitDist*100.0);
						double tof = -tAvgCFDcorr + absoluteOffset[xi*nVert+yi];
						double beta = (tof > 0) ? (hitDist/tof)/speedOfLight : -999.0;
						double gamma = (beta > 0 && beta < 0.9999) ? 1.0 / sqrt(1.0 - beta*beta) : -999.0;
						double E = -999.0; double p = -999.0; double px = -999.0; double py = -999.0; double pz = -999.0; double theta = -999.0; double phi = -999.0;
						if (gamma > 0){
							E = gamma * neutronMass;
							p = gamma * neutronMass * beta;
							theta = acos(wallZ/hitDist); //angle above beam axis
							phi = atan2(hitX,hitY); //azimuthal rotation around beam axis from -pi to pi, 0 pointing in +Y
							px = p*sin(theta)*sin(phi);
							py = p*sin(theta)*cos(phi);
							pz = p*cos(theta);
						}
						
						vec_tof.push_back(tof);
						vec_beta.push_back(beta);
						vec_gamma.push_back(gamma);
						vec_E.push_back(E);
						vec_p.push_back(p);
						vec_px.push_back(px);
						vec_py.push_back(py);
						vec_pz.push_back(pz);
						vec_th.push_back(theta*180/pi);
						vec_ph.push_back(phi*180/pi);
											
					}//if q>0	
				}//if t&q>0
			}//y
		}//x
		nHits = max(vec_x.size(),vec_y.size());
		if (nHits>0){
			normEvtNum = (double(i)*1000000.0)/nEvents;
			evtNum = int(i);
			if (nHits>1){ //SORT VECTORS BY INCREASING TOF... lambda funcs 
				vector<int> indVector(nHits);
				std::iota(indVector.begin(), indVector.end(), 0); //index vector

				sort(indVector.begin(), indVector.end(), [&](int a, int b){
					return vec_tof[a]<vec_tof[b];
				});

				auto sortByTof = [&](auto& vec){
					auto temp = vec; //create temporary copy and pass by ref
					for (int i=0; i<nHits; i++) vec[i] = temp[indVector[i]];
				};
				sortByTof(vec_x); sortByTof(vec_y); sortByTof(vec_tof); sortByTof(vec_qAvg); sortByTof(vec_Xpos); sortByTof(vec_Ypos); sortByTof(vec_Zpos); sortByTof(vec_dist); sortByTof(vec_beta); sortByTof(vec_gamma); sortByTof(vec_E); sortByTof(vec_p); sortByTof(vec_px); sortByTof(vec_py); sortByTof(vec_pz); sortByTof(vec_th); sortByTof(vec_ph);
 			}
			tMona->Fill();
		}
	}//loop
	cout << endl << Form("done with run%d. Closing files.",runNumber) << endl;
	//close files
	tMona->Write();
	fCal->Close();
	fRaw->Close();
	
}
