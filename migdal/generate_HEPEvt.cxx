#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TMath.h>
#include <TRandom.h>
using namespace std;

int main(int argc, char *argv[]){
  cout<<"option: [nuclei] [n] [l] [energy_of_incident_neutron(keV)] [anal_condition] [ev_total]"<<endl
      <<"  note: The \"anal_condition\" candidates are as follows: 000 010 001 100 110 101,"<<endl
      <<"        First  flag switch \"enable_precise_recoil\". 1: inelastic-like calculation, 0: normal elastic"<<endl
      <<"        Second flag switch \"enable_q_hosei\". Scaling of the probability of E_e by proportional to q (same to the Ibe paper)."<<endl
      <<"        Third  flag switch \"enable_q2_hosei\". Scaling of the probability of E_e by proportional to q^2 (for test)."<<endl;

  cout<<" ex) generate_HEPEvt Xe 1 0 565 110 10000"<<endl;
  cout<<" ex) generate_HEPEvt Xe 2 0 565 110 10000"<<endl;
  cout<<" ex) generate_HEPEvt Xe 2 1 565 110 10000"<<endl;
  cout<<" ex) generate_HEPEvt Ar 1 0 565 110 10000"<<endl<<endl;
  if(argc<5){ cout<<"ERROR: lack of option !"<<endl;  return 1; }

  TString nuclei = argv[1];
  int n_select   = atoi(argv[2]);
  int l_select   = atoi(argv[3]);
  double E_n     = atof(argv[4]);
  TString anal_condition = argv[5];
  double ev_total= atoi(argv[6]);

  // set E_nl
  TString filename;
  double E_dex, M_N, integ_prob;
  int PDG_code;
  if(nuclei=="Xe"){
    filename = "migdal_ibe_Xe.dat";
    M_N = 131e6; //keV
    PDG_code = 1000541310;
    if(n_select==1 && l_select==0){
      E_dex=34.566; //keV
      integ_prob=4.6e-6;
    }else if(n_select==2 && l_select==0){
      E_dex= 5.4; //keV
      integ_prob=2.9e-5;
    }else if(n_select==2 && l_select==1){
      E_dex= 4.9; //keV
      integ_prob=1.3e-4;
    }else{
      cout<<" n,l="<<n_select<<","<<l_select<<" is not supported in this program"<<endl;
      return 0;
    }
  }else if(nuclei=="Ar"){
    filename = "migdal_ibe_Ar.dat";
    M_N = 40e6; //keV
    PDG_code = 1000180400;
    if(n_select==1 && l_select==0){
      E_dex=3.206; //keV
      integ_prob=7.2e-5;
    }else{
      cout<<" n,l="<<n_select<<","<<l_select<<" is not supported in this program"<<endl;
      return 0;
    }
  }else{
    cout<<" nuclei="<<nuclei<<" is not supported in this program"<<endl;
    return 0;
  }

  // set X-ray data
  int num_X=0;
  double ene_X[50]={}, prob_X[50]={}, ddummy;
  ifstream ifs_X;
  if(     nuclei=="Xe") ifs_X.open("fl-tr-pr-54.dat");
  else if(nuclei=="Ar") ifs_X.open("fl-tr-pr-18.dat");
  if(! ifs_X){
    cout<<" ERROR: X-ray data fille is not found !!!"<<endl;
    return 0;
  }
  while(ifs_X >> ddummy){
    int to, from;
    if(ddummy==-1) ifs_X >> ddummy >> ddummy >> ddummy >> ddummy >> to;
    else if(ddummy==1) ifs_X >> ddummy >> to;
    else{
      //      from = (int)ddummy;
      double prob_tmp, ene_tmp;
      ifs_X >> prob_tmp >> ene_tmp;
      if( (n_select==1 && l_select==0 &&  1<=to && to<=2 ) || //1s
	  (n_select==2 && l_select==0 &&  3<=to && to<=4 ) || //2s
	  (n_select==2 && l_select==1 &&  5<=to && to<=7 ) || //2p
	  (n_select==3 && l_select==0 &&  8<=to && to<=9 ) || //3s
	  (n_select==3 && l_select==1 && 10<=to && to<=12) || //3p
	  (n_select==3 && l_select==2 && 13<=to && to<=15) || //3d
	  (n_select==4 && l_select==0 && 16<=to && to<=17) || //4s
	  (n_select==4 && l_select==1 && 18<=to && to<=20) || //4p
	  (n_select==4 && l_select==2 && 21<=to && to<=23) || //4d
	  (n_select==4 && l_select==3 && 24<=to && to<=26) || //4f
	  (n_select==5 && l_select==0 && 27<=to && to<=28) || //5s
	  (n_select==5 && l_select==1 && 29<=to && to<=31) || //5p
	  (n_select==5 && l_select==2 && 32<=to && to<=34) || //5d
	  (n_select==5 && l_select==3 && 35<=to && to<=37)    //5f
	  ){
	ene_X[num_X]  = ene_tmp*1000; // keV
	prob_X[num_X] = prob_tmp;
	cout<<" [X-ray] E="<<ene_X[num_X]<<", prob="<<prob_X[num_X]<<endl;
	num_X++;
      }
    }
  }

  // set analysis pattern
  bool enable_precise_recoil = false;
  bool enable_q_hosei = false;
  bool enable_q2_hosei = false;
  if(anal_condition=="000"){
    enable_precise_recoil = false;
    enable_q_hosei = false;
    enable_q2_hosei = false;
  }else if(anal_condition=="010"){
    enable_precise_recoil = false;
    enable_q_hosei = true;
    enable_q2_hosei = false;
  }else if(anal_condition=="001"){
    enable_precise_recoil = false;
    enable_q_hosei = false;
    enable_q2_hosei = true;
  }else if(anal_condition=="100"){
    enable_precise_recoil = true;
    enable_q_hosei = false;
    enable_q2_hosei = false;
  }else if(anal_condition=="110"){
    enable_precise_recoil = true;
    enable_q_hosei = true;
    enable_q2_hosei = false;
  }else if(anal_condition=="101"){
    enable_precise_recoil = true;
    enable_q_hosei = false;
    enable_q2_hosei = true;
  }else{
    cout<<" anal_condition="<<anal_condition<<" is not supported in this program"<<endl;
    return 0;
  }
  cout<<" enable_precise_recoil="<<enable_precise_recoil<<endl;
  cout<<" enable_q_hosei="<<enable_q_hosei<<endl;
  cout<<" enable_q2_hosei="<<enable_q2_hosei<<endl;

  // read migdal data for n,l
  ifstream ifs(filename.Data());
  TString str, dummy;
  int n_read, l_read;
  double ene[300], prob[300], prob_max=0;
  int g_index=0;
  while( ifs >> str ){
    if(str=="Principal"){
      ifs >> dummy >> dummy >> dummy >> dummy >> dummy >> n_read >> l_read >> dummy >> dummy >> dummy >> dummy;
      cout<<" header read n="<<n_read<<" l="<<l_read<<endl;
    }else{
      if(n_read==n_select && l_read==l_select){
	ene[g_index] = atof(str)*0.001;//MeV->keV
	ifs >> prob[g_index];
	if(prob[g_index]>prob_max) prob_max=prob[g_index];
	cout<<" g_index="<<g_index<<" ene="<<ene[g_index]<<" prob="<<prob[g_index]<<endl;
	g_index++;
      }
    }
  }

  // main loop
  double initial_neutron_num=0;
  ofstream ofs("event_list.HEPEvt");
  for(int ev =0; ev<ev_total; ev++){
    if(ev==9999) cout<<" "<<ev<<endl;
    if(ev%10000==0) cout<<ev<<endl;
    if(ev%1000==0) cout<<"|"<<flush;
    if(ev%100==0) cout<<"-"<<flush;

    // de-excitation X-ray
    double prob_X_sum=0;
    double E_X=0;
    for(int i=0; i<num_X; i++) prob_X_sum += prob_X[i];
    double ppp = gRandom->Uniform(0, prob_X_sum);
    for(int i=0; i<num_X; i++){
      if(ppp<prob_X[i]){
	E_X = ene_X[i];
	break;
      }else{
	ppp-=prob_X[i];
      }
    }
    if(E_X==0 || E_X>E_dex){
      cout<<"ERROR: energy of de-excitation X-ray ("<<E_X<<","<<E_dex<<") is wrong !!!"<<endl;
      return 1;
    }
    double costheta_X = gRandom->Uniform(-1,1);
    double phi_X = TMath::Pi()*2.0*gRandom->Uniform(0,1);//rad
    double dx_X = TMath::Cos(phi_X)*TMath::Sqrt(1-costheta_X*costheta_X);
    double dy_X = TMath::Sin(phi_X)*TMath::Sqrt(1-costheta_X*costheta_X);
    double dz_X = costheta_X;
    double p_X = E_X;
    double px_X = dx_X*p_X;
    double py_X = dy_X*p_X;
    double pz_X = dz_X*p_X;
    double costheta_X2 = gRandom->Uniform(-1,1);
    double phi_X2 = TMath::Pi()*2.0*gRandom->Uniform(0,1);//rad
    double dx_X2 = TMath::Cos(phi_X2)*TMath::Sqrt(1-costheta_X2*costheta_X2);
    double dy_X2 = TMath::Sin(phi_X2)*TMath::Sqrt(1-costheta_X2*costheta_X2);
    double dz_X2 = costheta_X2;
    double M_e = 511; //keV
    double E_X2 = E_dex-E_X;
    double p_X2 = sqrt(2*M_e*E_X2+E_X2*E_X2);
    double px_X2 = dx_X2*p_X2;
    double py_X2 = dy_X2*p_X2;
    double pz_X2 = dz_X2*p_X2;

    // auger electron
    double E_e;
    while(1){
      E_e = gRandom->Uniform(0.001,70);//keV
      double ppp = gRandom->Uniform(0,prob_max);
      bool ok_flag=false;
      for(int i=0; i<g_index-1; i++){
	if(ene[i]<E_e && E_e<ene[i+1]){
	  double p_compare = (prob[i]*(E_e-ene[i])+prob[i+1]*(ene[i+1]-E_e))/(ene[i+1]-ene[i]);
	  //cout<<" ene[i]="<<ene[i]<<" E_e="<<E_e<<" ene[i+1]="<<ene[i+1]<<" ppp="<<ppp<<"p_compare="<<p_compare<<endl;
	  if(ppp<p_compare) ok_flag=true;
	}
      }
      if(ok_flag) break;
    }
    double costheta_e = gRandom->Uniform(-1,1);
    double phi_e = TMath::Pi()*2.0*gRandom->Uniform(0,1);//rad
    double dx_e = TMath::Cos(phi_e)*TMath::Sqrt(1-costheta_e*costheta_e);
    double dy_e = TMath::Sin(phi_e)*TMath::Sqrt(1-costheta_e*costheta_e);
    double dz_e = costheta_e;
    double p_e = sqrt(2*M_e*E_e+E_e*E_e);
    double px_e = dx_e*p_e;
    double py_e = dy_e*p_e;
    double pz_e = dz_e*p_e;

    // nuclear
    double M_n = 1e6;   //keV
    //    double M_N = 40e6;  //keV
    //    double E_n = 565;   //keV
    double E_NR_max = 4*M_n*M_N/TMath::Power(M_n+M_N,2)*E_n; //keV
    //    cout<<" E_NR_max="<<E_NR_max<<endl;
    double E_delta = E_dex+E_e; //keV
    double root_factor = TMath::Sqrt(1.0-(M_N+M_n)*E_delta/M_N/E_n);
    double E_NR; //keV
    double costheta_CM;
    if(enable_precise_recoil==false) root_factor = 1.0;
    while(1){
      costheta_CM = gRandom->Uniform(-1.0, 1.0);
      E_NR = E_NR_max*( TMath::Power(1.0-root_factor,2)/4.0 + (1+costheta_CM)*root_factor/2.0 );
      if(enable_q_hosei==false && enable_q2_hosei==false) break;
      double ppp = gRandom->Uniform(0,1);
      double E_NR_max_actual = E_NR_max*(TMath::Power(1.0-root_factor,2)/4.0 + 1.0);
      if(enable_q_hosei==true){
	if( ppp < TMath::Power(E_NR/E_NR_max_actual, 1) ) break;
      }else if(enable_q2_hosei==true){
	if( ppp < TMath::Power(E_NR/E_NR_max_actual, 2) ) break;
      }
    }
    //    E_NR=0.511*0.511*M_N/2/M_e/M_e; //keV
    //    cout<<" E_NR="<<E_NR<<endl;

    double theta_CM = TMath::ACos(costheta_CM);//rad
    double theta_LAB = 0.5*theta_CM;//rad
    double costheta_LAB = TMath::Cos(theta_LAB);
    double phi_NR = TMath::Pi()*2.0*gRandom->Uniform(0,1);//rad
    double dx_NR = TMath::Cos(phi_NR)*TMath::Sqrt(1-costheta_LAB*costheta_LAB);
    double dy_NR = TMath::Sin(phi_NR)*TMath::Sqrt(1-costheta_LAB*costheta_LAB);
    double dz_NR = costheta_LAB;
    double p_NR = sqrt(2*M_N*E_NR+E_NR*E_NR);
    double px_NR = dx_NR*p_NR;
    double py_NR = dy_NR*p_NR;
    double pz_NR = dz_NR*p_NR;

    //HEPEvt
    // n_particle
    //state  PDG_code  first_daughter  last_daughter  px_in_GeV  py_in_GeV  pz_in_GeV  mass_in_GeV

    //state
    // 0:InitialState
    // 1:StableFinalState
    // 2: IntermediateState
    // ... etc.

    //PDG code
    // Xe131: 1000541310  Ar40: 1000180400
    // e-:    11
    // gamma: 22

    ofs << "4" <<endl;
    ofs << " 1 "<<PDG_code<<" 0 0 " <<px_NR*1e-6<<" "<<py_NR*1e-6<<" "<<pz_NR*1e-6<<" "<<M_N*1e-6<<endl;
    ofs << " 1    11      0 0 " <<px_e *1e-6<<" "<<py_e *1e-6<<" "<<pz_e *1e-6<<" "<<M_e*1e-6<<endl;
    ofs << " 1    22      0 0 " <<px_X *1e-6<<" "<<py_X *1e-6<<" "<<pz_X *1e-6<<" 0"<<endl;
    ofs << " 1    11      0 0 " <<px_X2*1e-6<<" "<<py_X2*1e-6<<" "<<pz_X2*1e-6<<" "<<M_e*1e-6<<endl;

    //p_total
    double q_e = sqrt(2*511*511*E_NR/M_N); //keV
    double prob_for_NR = integ_prob * q_e*q_e/0.511/0.511;
    //    cout<<" E_NR="<<E_NR<<" q_e="<<q_e<<" prob_for_NR="<<prob_for_NR<<endl;
    initial_neutron_num += 1.0/prob_for_NR;
  }

  cout<<" init_n_num = "<<initial_neutron_num<<endl;
  cout<<" ev_total = "<<ev_total<<endl;
  cout<<" init_n_num/ev_total = "<<initial_neutron_num/ev_total<<endl;
  cout<<" prob(=ev_total/init_n_num) = "<<ev_total/initial_neutron_num<<endl;

  ofstream ofs2("event_list.HEPEvt_info");
  ofs2 << "condition of event_list.HEPEvt" << endl << endl;
  ofs2 << "  nuclei         " << nuclei << endl;
  ofs2 << "  n_select       " << n_select << endl;
  ofs2 << "  l_select       " << l_select << endl;
  ofs2 << "  E_n            " << E_n << endl;
  ofs2 << "  anal_condition " << anal_condition << endl;
  ofs2 << "  ev_total       " << ev_total << endl << endl;

  ofs2 << "  filename       " << filename << endl;
  ofs2 << "  M_N            " << M_N << endl;
  ofs2 << "  PDG_code       " << PDG_code << endl;
  ofs2 << "  enable_precise_recoil  " << enable_precise_recoil << endl;
  ofs2 << "  enable_q_hosei         " << enable_q_hosei << endl;
  ofs2 << "  enable_q2_hosei        " << enable_q2_hosei << endl;

  ofs2 << "  number_of_neutron_per_migdal       " << initial_neutron_num/ev_total << endl;
  ofs2 << "  averaged_probability_of_migdal     " << ev_total/initial_neutron_num << endl;


}
