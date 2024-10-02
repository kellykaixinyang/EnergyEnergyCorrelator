#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TAttLine.h"
#include "THelix.h"
#include "TView.h"
#include "TRandom.h"
#include "TAttPad.h"
#include "TMath.h"
#include "TVector3.h"
#include "TView3D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

#include <iostream>
#include <string>
#include <vector>

//#include "SURF_FUNCTIONS.h"
//#include "coordinateTools.h"
//R__LOAD_LIBRARY(SURF_FUNCTIONS_C);

// Global variables
std::string title = "#mu^{+}#mu^{#minus} (3000 GeV)";
std::string fullPathDir = "/Users/kellyyanggoingharvard/Documents/MCdata/data/mupmum_3000GeV.root";
Float_t etaCutOff = 100;
Int_t multCutOff = 194;
TFile *f = new TFile(fullPathDir.c_str(), "read"); // Opening file

// Creating TTreeReader object and linking branches
TTreeReader* reader = new TTreeReader("trackTree");

// Setup branches for other particles
TTreeReaderValue<std::vector<vector<int>>> chgBranch(*reader,"genDau_chg");
TTreeReaderValue<std::vector<vector<float>>> etaBranch(*reader,"genDau_eta");
TTreeReaderValue<std::vector<vector<float>>> phiBranch(*reader,"genDau_phi"); 
TTreeReaderValue<std::vector<vector<float>>> ptBranch(*reader,"genDau_pt");
TTreeReaderValue<std::vector<float>> jetptBranch(*reader,"genJetPt");
TTreeReaderValue<std::vector<float>> jetetaBranch(*reader,"genJetEta");
TTreeReaderValue<std::vector<float>> jetphiBranch(*reader,"genJetPhi");
TTreeReaderValue<std::vector<int>> jetmultiplicityBranch(*reader,"genJetChargedMultiplicity");


// Function that returns 1D EEC distribution
TH1F createEnergyDistr() {
		
	std::string energyTitle = "Energy Energy Correlator for " + title;

	//Histogram for the EEC distribution 
	//const float EnergyBW = 0.3;
	TH1F hEnergy("hEnergy", energyTitle.c_str(),25,0,1);

	// **** EVENT LOOP ****
	reader->Restart(); // Restarting event loop
	//Int_t numTrigg = 0;
	//Int_t numSelectedEvents = 0;

	while (reader->Next())  { //looping through the events

		// Check to see if the event is in the multiplicity bin
		//Int_t eventIndex = reader->GetCurrentEntry();
		//if (multiplicityVector[eventIndex] < multCutOff) {continue;}
		//numSelectedEvents++;

		// ***** JET LOOP FOR GIVEN EVENT *****
		for (Int_t i = 0; i < etaBranch->size(); i++) { //etaBranch gives number of jets
		
			if (fabs((*jetetaBranch)[i]) >= 1.6) {continue;}
			//if (-0.1 < ((*jetphiBranch)[i]) < 1.2) {continue;}	
		
			// ***** PARTICLE LOOP FOR JET i *****
			//Particle 1 Loop
			for (Int_t j = 0; j < (*etaBranch)[i].size()-1; j++) {

            			// Checking if particle 1 is charged
            			if ((*chgBranch)[i][j] == 0) {continue;}
				
				// Checking if particle 1 charge reaches threshold
				if ((*chgBranch)[i][j] < 1) {continue;}
		
				Float_t eta1 = (*etaBranch)[i][j]; 
            			if (fabs(eta1) > etaCutOff) {continue;} // Checking eta range
				Float_t phi1 = (*phiBranch)[i][j];
			
				Float_t pt1 = (*ptBranch)[i][j];

				// Particle 2 Loop	
            			for (Int_t k = j + 1; k < (*etaBranch)[i].size(); k++) {
				
                			// Checking if particle 2 is charged
                			if ((*chgBranch)[i][k] == 0) {continue;}
					Float_t eta2 = (*etaBranch)[i][k];
					if (fabs(eta2) > etaCutOff) {continue;} // Checking eta range
					Float_t phi2 = (*phiBranch)[i][k];
					
					Float_t pt2 = (*ptBranch)[i][k];

					// Calculating delta eta, delta phi, and delta r
                			Float_t deltaEta = eta2 - eta1;
                			Float_t deltaPhi = TMath::ACos(TMath::Cos(phi2-phi1));	
					Float_t delta_r = sqrt(pow(deltaEta,2)+pow(deltaPhi,2)) / 0.8;

					//Filling 1D Histogram
					int n = 1;
 
					hEnergy.Fill(delta_r,pow((pt1*pt2), n));
				}//ending particle 2 loop
			}//ending particle 1 loop
		}//ending jet loop
	}//ending event loop
	//std::cout << "Number of selected events for Energy Distribution: " << numSelectedEvents << std::end1;

	// ***** HISTOGRAM CUSTOMIZATION *****
	hEnergy.GetXaxis()->SetTitle("#Delta r");
	
	hEnergy.SetYTitle("EEC");
	
	//hEnergy.GetXaxis()->SetTitleSize(0.05);
	//hEnergy.GetXaxis()->SetTitleFont(64);
	//hEnergy.GetYaxis()->SetTitleSize(0.05);
       // hEnergy.GetYaxis()->SetTitleFont(64);

	hEnergy.SetTitleFont(200, "T");
        hEnergy.SetTitleFont(200, "XYZ");
	
	
	
	hEnergy.SetLabelFont(132,"T");
	hEnergy.SetLabelFont(132,"XYZ");

	return hEnergy;
}

void EECCorr() {
	// Finding the multiplicities of each event and storing it in a vector
    	//std::vector<Int_t> multVec;
	//reader->Restart(); // Ensuring event loop starts from beginning
    	//while (reader->Next()) {
        //	Int_t multiplicity = 0; // Counter for N_ch
        //	for (Int_t i = 0; i < pxBranch->size(); i++) {
          //  		if ((*chgBranch)[i] != 0) {multiplicity++;}
        //	}
        //	multVec.push_back(multiplicity);
        
    //	}
    
	TFile *fout = new TFile("simpleEEC.root", "recreate");

	// Creating canvas for distribution
	TCanvas *cEnergy = new TCanvas("cEnergy", "Canvas for delta r Distribution", 800, 600);
	TH1F simpleEECHist = createEnergyDistr();

	cEnergy->cd();
	cEnergy->SetLogx();	
	cEnergy->SetLogy();

	simpleEECHist.Draw();
	cEnergy->Write();
	
	simpleEECHist.Write();

	delete cEnergy;
	f->Close();
	fout->Close();
}
    

	
