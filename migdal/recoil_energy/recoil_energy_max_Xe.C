{
  //E_NR_max = 4*M_n/M_N*E_n = 4/130*565
  TF1 *f = new TF1("f1","x*[0]",0,10000);
  f->SetParameter(0,4./130.);
  f->SetLineColor(2);

  TCanvas *c = new TCanvas("c","c",800,500);
  c->DrawFrame(0,0,10000,350,"E_n keV;E_NR_max keV");
  f->Draw("sameLP");

  TGraph *g = new TGraph();
  double px[6]={ 565, 911, 2500, 5000, 8000, 14800 };
  for(int i=0; i<6; i++) g->SetPoint(i, px[i], f->Eval(px[i]));
  g->SetMarkerStyle(4);
  g->Draw("sameP");

  TText *t[6];
  for(int i=0; i<6; i++){
    t[i] = new TText(500, 170+i*30, Form("E_n=%.0f  E_NR_max=%.1f",px[i],f->Eval(px[i])));
    t[i]->Draw("same");
  }


  c->Print("recoil_energy_max_Xe.png");
  system("send2kakunouko.sh recoil_energy_max_Xe.png");

}
