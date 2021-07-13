{
  ifstream ifs("migdal_ibe_Ar.dat");
  TString str, dummy;
  int n, l;
  double factor = 511*511*1e3/2/TMath::Pi();
  double ene[6][3][300], prob[6][3][300];
  int g_index[6][3]={};
  while( ifs >> str ){
    if(str=="Principal"){
      ifs >> dummy >> dummy >> dummy >> dummy >> dummy >> n >> l >> dummy >> dummy >> dummy >> dummy;
      cout<<" header read n="<<n<<" l="<<l<<endl;
    }else{
      ene[n][l][g_index[n][l]] = atof(str)*0.001;//MeV->keV
      ifs >> prob[n][l][g_index[n][l]];
      prob[n][l][g_index[n][l]]*=factor;
      g_index[n][l]++;
    }
  }

  int i_color[6] = {0, 4, 1, 3, 2, 6};
  int i_style[3] = {1, 2, 3};
  TGraph *graph[6][3];
  for(int n=0; n<6; n++){
    for(int l=0; l<3; l++){
      if(g_index[n][l]>0){
	cout<<" n="<<n<<" l="<<l<<" i="<<g_index[n][l]<<endl;
	graph[n][l] = new TGraph(g_index[n][l], ene[n][l], prob[n][l]);
	graph[n][l]->SetLineColor(i_color[n]);
	graph[n][l]->SetLineStyle(i_style[l]);
	graph[n][l]->SetLineWidth(2*l+1);
      }
    }
  }

  TLegend *leg = new TLegend(0.75,0.35,0.9,0.9);
  for(int n=0; n<6; n++) for(int l=0; l<3; l++) if(g_index[n][l]>0) leg->AddEntry(graph[n][l],Form("n=%d, l=%d",n,l),"l");

  TCanvas *c = new TCanvas("c","c",1800,600);
  c->Divide(2,1,0.001,0.001);
  c->cd(1)->SetLogx();  c->cd(1)->SetLogy();
  c->cd(1)->DrawFrame(1e-3, 1e-7, 2e1, 1e2,"probability (q_{e}=511eV);E_{e} / keV;#frac{1 dp^{c}}{2#pidE_{e}} / keV^{-1}");
  for(int n=0; n<6; n++) for(int l=0; l<3; l++) if(g_index[n][l]>0) graph[n][l]->Draw("sameL");
  leg->Draw();

  c->cd(2)->DrawFrame(0, 0, 1e1, 7e-5,"probability (q_{e}=511eV);E_{e} / keV;#frac{1 dp^{c}}{2#pidE_{e}} / keV^{-1}");
  for(int n=0; n<6; n++) for(int l=0; l<3; l++) if(g_index[n][l]>0) graph[n][l]->Draw("sameL");

  c->Print("draw_migdal_ibe_Ar.png");
}
