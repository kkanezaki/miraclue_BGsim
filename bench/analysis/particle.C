void particle(){

	TString filename = "out.root";
	TFile* file = new TFile(filename, "read");
	TTree* tree = (TTree*)file->Get("tree");

	double min = 1000000000;
	double max = 2000000000;
	double bins = min;

	TH1D* hist = new TH1D("particle_id","particle_id;PDG particle code;counts/bin",bins,min,max);
	//hist->SetTitle("particle_id;")
	tree->Draw("particle>>particle_id");


	//tree->Draw("particle>>h(10000000,1000000000,1000200000)");
	//Int_t N = tree->GetEntries();
	//for( int_t i=0, i<N, i++){
		
		 
}
