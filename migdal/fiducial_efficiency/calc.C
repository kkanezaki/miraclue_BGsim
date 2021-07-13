{

  TH1D *h_eff = new TH1D("h_eff","h_eff", 90,0,30);
  TH2D *h_dxy = new TH2D("h_dxy","h_dxy",100,-1,1,100,-1,1);
  TH2D *h_dxz = new TH2D("h_dxz","h_dxz",100,-1,1,100,-1,1);
  TH2D *h_dyz = new TH2D("h_dyz","h_dyz",100,-1,1,100,-1,1);

  for(int ev=0; ev<100000; ev++){
    double x0 = gRandom->Uniform(-15.0,15.0);
    double y0 = gRandom->Uniform(-15.0,15.0);
    double z0 = gRandom->Uniform(  0.0,30.0);
    double costheta = gRandom->Uniform(-1.0,1.0);
    double phi      = gRandom->Uniform(0.0,2.0*TMath::Pi());
    double theta = TMath::ACos(costheta);
    double sintheta = TMath::Sin(theta);

    double dx = sintheta*TMath::Cos(phi);
    double dy = sintheta*TMath::Sin(phi);
    double dz = costheta;

    h_dxy->Fill(dx,dy);
    h_dxz->Fill(dx,dz);
    h_dyz->Fill(dy,dz);

    for(int i=0; i<90; i++){
      double dl = gRandom->Uniform(0.0,1.0);
      double length = 30.0/90.0*((double)i+dl);
      double x1 = x0 + dx*length;
      double y1 = y0 + dy*length;
      double z1 = z0 + dz*length;
      if(x1>-14 && x1<14 && y1>-14 && y1<14 && z1>1 && z1<29 &&
	 x0>-14 && x0<14 && y0>-14 && y0<14 && z0>1 && z0<29){
	h_eff->Fill(length);
      }
    }
  }

  h_eff->Scale(1/100000.0);

  TCanvas *c = new TCanvas("c","c",800,700);  c->Divide(2,2,0.001,0.001);
  c->cd(1);  h_eff->Draw("hist");
  c->cd(2);  h_dxy->Draw("COLZ");
  c->cd(3);  h_dxz->Draw("COLZ");
  c->cd(4);  h_dyz->Draw("COLZ");

  TFile *f = new TFile("eff.root","recreate");
  h_eff->Write();
  c->Write();



}
