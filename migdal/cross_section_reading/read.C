{
  ifstream ifs("Ar.txt");
  int ii=0;
  double ene, xsec, xsec_avg=0;
  while(ifs >> ene >> xsec){
    xsec_avg += xsec;
    cout<<" ii="<<ii<<" ene="<<ene<<" xsec="<<xsec<<" "<<xsec_avg<<endl;
    ii++;
  }
  xsec_avg /= (double)ii;
  cout<<"xsec_avg="<<xsec_avg<<endl;

}
