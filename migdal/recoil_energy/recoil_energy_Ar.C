{

  double E_delta_max=100; //keV

  TGraph *g[3];
  for(int i=0; i<3; i++) g[i] = new TGraph();
  for(int i=0; i<3; i++) g[i]->SetLineColor(i+2);

  double costheta_CM[3] = {-1,0,1};
  for(int i=0; i<1000; i++){
    double M_n = 1e6;   //keV
    double M_N = 40e6; //keV
    double E_n = 565;   //keV
    double E_NR_max = 4*M_n*M_N/TMath::Power(M_n+M_N,2)*E_n;//keV
    double E_delta = E_delta_max*i/1000;//keV
    double root_factor = TMath::Sqrt(1.0-(M_N+M_n)*E_delta/M_N/E_n);
    for(int j=0; j<3; j++){
      double E_NR = E_NR_max*( TMath::Power(1.0-root_factor,2)/4.0 + (1-costheta_CM[j])*root_factor/2.0 );
      g[j]->SetPoint(i, E_delta, E_NR);
    }
  }

  TCanvas *c = new TCanvas("c","c",1200,800);
  c->cd(1)->DrawFrame(0,0,105,60,";#DeltaE / keV;E_{NR} / keV");
  for(int i=0; i<3; i++) g[i]->Draw("sameL");

  TLegend *leg = new TLegend(0.5,0.5,0.7,0.7);
  for(int i=0; i<3; i++) leg->AddEntry(g[i],Form("cos#theta_{CM}=%.1f",costheta_CM[i]),"l");
  leg->Draw();


  c->Print("recoil_energy_Ar.png");




}
