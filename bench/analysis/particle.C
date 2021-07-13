void particle(){

	TString filename = "out.root";
	TFile* file = new TFile(filename, "read");
	TTree* tree = (TTree*)file->Get("tree");

	TH1D* hist = new TH1D("particle_id","particle_id;PDG particle code;counts/bin",10000000,1000100000,1000200000);
	//hist->SetTitle("particle_id;")
	tree->Draw("particle>>particle_id");


	//tree->Draw("particle>>h(10000000,1000000000,1000200000)");
	//Int_t N = tree->GetEntries();
	//for( int_t i=0, i<N, i++){
		
		 
}
