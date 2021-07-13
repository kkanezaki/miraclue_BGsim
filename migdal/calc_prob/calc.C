{
  double M_n=939.5693e3;//keV
  double M_u=931.4941024e3;
  double N_Xe=131.293;
  double N_Ar=39.948;
  double M_Xe=N_Xe*M_u;//keV
  double M_Ar=N_Ar*M_u;//keV
  double E_n=565.;//keV
  double M_e=511.;//keV
  double E_Xe=4*M_n*M_Xe/(M_n+M_Xe)/(M_n+M_Xe)*E_n;
  double E_Ar=4*M_n*M_Ar/(M_n+M_Ar)/(M_n+M_Ar)*E_n;
  double q2_Xe=2*M_e*M_e*E_Xe/M_Xe;
  double q2_Ar=2*M_e*M_e*E_Ar/M_Ar;
  double scaling_Xe=q2_Xe/0.511/0.511;
  double scaling_Ar=q2_Ar/0.511/0.511;
  cout<<"Ar: E_max="<<E_Ar<<" q2="<<q2_Ar<<" scaling="<<scaling_Ar<<endl;
  cout<<"Xe: E_max="<<E_Xe<<" q2="<<q2_Xe<<" scaling="<<scaling_Xe<<endl;
}



