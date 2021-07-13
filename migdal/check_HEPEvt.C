{

  TString nuclei;  int n, l;  double E_n;
  TPaveText *pt = new TPaveText(0,0,1,1);  pt->SetFillColor(0);
  ifstream ifs_info("event_list.HEPEvt_info");
  TString str0, str1;
  ifs_info >> dummy >> dummy >> dummy;
  while( ifs_info >> str0 >> str1 ){
    pt->AddText(str0+" :  "+str1);
    if(str0=="nuclei") nuclei = str1;
    if(str0=="n_select") n = str1.Atoi();
    if(str0=="l_select") l = str1.Atoi();
    if(str0=="E_n") E_n = str1.Atof();
  }

  // for check
  ifstream ifs_chk("event_list.HEPEvt");

  int n_particle, PDG_code[4];
  double px[4], py[4], pz[4], mass[4];
  TTree *tree = new TTree("tree","tree");
  tree->Branch("n_particle", &n_particle);
  tree->Branch("PDG_code",  PDG_code, "PDG_code[n_particle]/I");
  tree->Branch("px",  px, "px[n_particle]/D");
  tree->Branch("py",  py, "py[n_particle]/D");
  tree->Branch("pz",  pz, "pz[n_particle]/D");
  tree->Branch("mass",  mass, "mass[n_particle]/D");

  TString dummy;
  while( ifs_chk >> n_particle ){
    for( int i=0; i<n_particle; i++) ifs_chk >> dummy >> PDG_code[i] >> dummy >> dummy >> px[i] >> py[i] >> pz[i] >> mass[i];
    tree->Fill();
  }
  int n_ev = tree->GetEntries();

  double ene_min[4] = {0, 0, 0, 0}, ene_max[4] = {18, 75, 40, 40};
  if(nuclei=="Xe" && n==1 && l==0){ ene_max[0]=18;  ene_max[1]=75;  ene_max[2]=40;  ene_max[3]=40; }
  if(nuclei=="Xe" && n==1 && l==0 && E_n==911){ ene_max[0]=30;  ene_max[1]=75;  ene_max[2]=40;  ene_max[3]=40; }
  if(nuclei=="Xe" && n==2 && l==0){ ene_max[0]=18;  ene_max[1]=75;  ene_max[2]=40;  ene_max[3]=40; }
  if(nuclei=="Xe" && n==2 && l==1){ ene_max[0]=18;  ene_max[1]=75;  ene_max[2]=40;  ene_max[3]=40; }
  if(nuclei=="Ar" && n==1 && l==0){ ene_max[0]=60;  ene_max[1]=10;  ene_max[2]=10;  ene_max[3]=10; }

  TString title[4] = {"recoil nuclear", "auger e-", "de-excitation X-ray", "remaining de-excitation e-"};
  TString ene_unit[4] = {"; E_{NR} / keV", "; E_{e} / keV", "; E_{X} / keV", "; E_{dex}-E_{X} / keV"};
  TH1D *h1_ene[4], *h1_cos[4], *h1_phi[4], *h1_cos_CM[4];
  TH2D *h2_ene_cos[4], *h2_phi[4];
  TH2D *h2_ene = new TH2D("h2_ene","h2_ene"+ene_unit[0]+ene_unit[1], 40,ene_min[0],ene_max[0], 40,ene_min[1],ene_max[1]);
  for(int i=0; i<4; i++){
    h1_ene[i]     = new TH1D(Form("h1_ene_%d",    i), title[i]+ene_unit[i],        100, ene_min[i], ene_max[i]);
    h1_cos[i]     = new TH1D(Form("h1_cos_%d",    i), title[i]+";cos#theta_{LAB}", 100,         -1,          1);
    h1_phi[i]     = new TH1D(Form("h1_phi_%d",    i), title[i]+";#phi_{LAB}",      100,          0,          7);
    h1_cos_CM[i]  = new TH1D(Form("h1_cos_CM_%d", i), title[i]+";cos#theta_{CM}",  100,         -1,          1);
    h2_phi[i]     = new TH2D(Form("h2_phi_%d",    i), title[i]+";vx;vy",                        100,         -1,          1, 100, -1, 1);
    h2_ene_cos[i] = new TH2D(Form("h2_ene_cos_%d",i), title[i]+ene_unit[i]+";cos#theta_{LAB}",  100, ene_min[i], ene_max[i], 100, -1, 1);
    tree->Draw(Form("(TMath::Sqrt(px[%d]**2+py[%d]**2+pz[%d]**2+mass[%d]**2)-mass[%d])*1e6>>h1_ene_%d",i,i,i,i,i,i));
    tree->Draw(Form("pz[%d]/TMath::Sqrt(px[%d]**2+py[%d]**2+pz[%d]**2)>>h1_cos_%d",i,i,i,i,i));
    tree->Draw(Form("TMath::Cos(TMath::ACos(pz[%d]/TMath::Sqrt(px[%d]**2+py[%d]**2+pz[%d]**2))*2)>>h1_cos_CM_%d",i,i,i,i,i));
    tree->Draw(Form("TMath::ACos(px[%d]/TMath::Sqrt(px[%d]**2+py[%d]**2))>>h1_phi_%d",i,i,i,i));
    tree->Draw(Form("py[%d]/TMath::Sqrt(px[%d]**2+py[%d]**2+pz[%d]**2) : px[%d]/TMath::Sqrt(px[%d]**2+py[%d]**2+pz[%d]**2)>>h2_phi_%d",i,i,i,i,i,i,i,i,i));
    tree->Draw(Form("pz[%d]/TMath::Sqrt(px[%d]**2+py[%d]**2+pz[%d]**2) : (TMath::Sqrt(px[%d]**2+py[%d]**2+pz[%d]**2+mass[%d]**2)-mass[%d])*1e6>>h2_ene_cos_%d",i,i,i,i,i,i,i,i,i,i));
    h1_cos[i]->SetMinimum(0);
    h1_phi[i]->SetMinimum(0);
    h1_cos_CM[i]->SetMinimum(0);
  }
  tree->Draw("(TMath::Sqrt(px[1]**2+py[1]**2+pz[1]**2+mass[1]**2)-mass[1])*1e6 : (TMath::Sqrt(px[0]**2+py[0]**2+pz[0]**2+mass[0]**2)-mass[0])*1e6>>h2_ene");

  TCanvas *c = new TCanvas("c","c",2800,1900);
  c->Divide(6,5,0.001,0.001);
  for(int i=0; i<4; i++){
    c->cd(1+i*6);  h1_ene[i]->Draw();
    c->cd(2+i*6);  h1_cos[i]->Draw();
    c->cd(3+i*6);  h1_phi[i]->Draw();
    c->cd(4+i*6);  h2_phi[i]->Draw("COLZ");
    c->cd(5+i*6);  h1_cos_CM[i]->Draw();
    c->cd(6+i*6);  h2_ene_cos[i]->Draw("COLZ");
  }
  c->cd(25);  h2_ene->Draw("COLZ");
  c->cd(26);  h2_ene->Draw("surf");
  c->cd(27);  pt->Draw();

  c->Print("check_HEPEvt.png");

}
