void mat_test(){

	TString filename = "out.root";
	TFile* file = new TFile(filename, "read");
	TTree* tree = (TTree*)file->Get("tree");

	TH1D* hist = new TH1D("hist","Energy;Energy[MeV];counts/bin",5000,0,0.06);
	//hist->SetTitle("particle_id;")
	tree->Draw("eDep>>Energy");
	tree->Draw("eIn>>Energy");

	double e_max;
	double e_n = 0.565; //MeV
	double m_n = 939.6; //MeV
	double m_ar = 40.*m_n; //MeV

	e_max = e_n * (4*m_n*m_ar)/pow((m_n+m_ar),2);
	cout << e_max << endl;
	
	//tree->Draw("particle>>h(10000000,1000000000,1000200000)");
	//Int_t N = tree->GetEntries();
	//for( int_t i=0, i<N, i++){
		
		 
}
