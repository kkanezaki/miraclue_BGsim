#include <vector>
#include <iostream>
using namespace std;

void intrinsic(){

	TString filename = "xe8.root";
	TFile* file = new TFile(filename, "read");
	TTree* tree = (TTree*)file->Get("tree");

	gROOT->SetStyle("ATLAS");

	TH1D* hist = new TH1D("hist","intrinsic BG;Energy[keV];counts/keV",100,0,100);
	TH1D* ion2 = new TH1D("ion2","intrinsic BG;Energy[keV];counts/keV",100,0,100);
	TH1D* fiducial = new TH1D("fiducial","intrinsic BG;Energy[keV];counts/keV",100,0,100);
	//hist->SetTitle("particle_id;")
	//tree->Draw("eDep>>Energy");
	//tree->Draw("eIn>>Energy", "eIn<1");
	//tree->Draw("particle>>h(10000000,1000000000,1000200000)");
	
	
	Int_t N = tree->GetEntries();
	vector<Int_t> *particle=0; 
	vector<Int_t> *TrackID=0; 
	vector<Double_t> *eIn=0; 
	vector<Double_t> *eDep=0;
	vector<Double_t> *x=0;
	vector<Double_t> *y=0;
	vector<Double_t> *z=0;

	tree->SetBranchAddress("particle", &particle);
	tree->SetBranchAddress("TrackID", &TrackID);
	tree->SetBranchAddress("eIn", &eIn);
	tree->SetBranchAddress("eDep", &eDep);
	tree->SetBranchAddress("x", &x);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("z", &z);

	Double_t x_0=0.;
	Double_t y_0=0.;
	Double_t z_0=0.;
	Int_t events=0;
	
	for( Int_t i=0; i<N; i++){
		tree->GetEntry(i);	

		if(eDep->size()!=0){
		//if(particle->size()!=0 && particle->size()!=1){
			//eDep_sum = std::accumulate(eDep->begin(),eDep->end(),0.);
			//cout << eDep->begin() << "\t" << eDep->end()<< endl;
			//hist->Fill(eDep_sum);
			events += 1;

			Double_t eDep_sum=0.;
			Int_t k=0;
			bool fiducial_draw=true;

			for(Int_t j=0; j<eDep->size(); ++j){
				eDep_sum += 1000*eDep->at(j); //keV
			
				if(TrackID->at(0)>1){
					x_0 = x->at(0);
					y_0 = y->at(0);
					z_0 = z->at(0);

					//Double_t r = pow( ( (x_0 - x->(j))**2 + (y_0 - y->(j))**2 + (z_0 - z->(j))**2),0.5); 

					if(particle->at(j)>1000000 && eDep!=0){
						k += 1;	
					}
				}

				if(particle->at(j)>1000000  && (x->at(j)>28 | y->at(j)>28 | z->at(j)>28) ){
					fiducial_draw=false;
				}
				
				/*
				cout << i <<"\t" << particle->at(j) << "\t" << eIn->at(j) 
					 << "\t" << eDep->at(j) << "\t" << TrackID->at(j) << endl;	
				*/
			}

			hist->Fill(eDep_sum);
			if(fiducial_draw){
				fiducial->Fill(eDep_sum);
			}if(fiducial_draw && k==2){
				ion2->Fill(eDep_sum);	
			}
					
		}
	}		

	hist->Draw();
	fiducial->SetLineColor(9);
	fiducial->Draw("same");
	ion2->SetLineColor(2);
	ion2->Draw("same");
	gPad->SetLogy();

	int A=130; //40 or 130
	double density=5.8971*8; //mg/cm3 1.7606

	cout << events << "/" << N << " events"<< endl;
	cout << "sigma = " << pow(10,24)*events/(N * 30. * 6.022*pow(10,23) * density * pow(10,-3)/A) << " barn" << endl;
}


