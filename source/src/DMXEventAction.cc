//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// History/Additions:
// 16 Jan 2002  Added analysis
//
//
// EventAction program
// --------------------------------------------------------------

#include "DMXEventAction.hh"
#include "DMXPrimaryGeneratorAction.hh"
#include "DMXDetectorConstruction.hh"
// note DMXTPCHit.hh and DMXLSHit.hh are included in DMXEventAction.hh
#include "DMXEventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include <fstream>
#include <iomanip>
#include "G4SystemOfUnits.hh"

DMXEventAction::DMXEventAction(G4String filename, G4String commandname, G4int event_num){

  eventMessenger = new DMXEventActionMessenger(this);

  char pcname[128], username[128], pwdname0[128], pwdname1[128], datename[128];
  FILE *fp;
  fp = popen("whoami | tr -d '\n'",                                                   "r");  fgets(username,sizeof(username),fp);
  fp = popen("hostname | tr -d '\n'",                                                 "r");  fgets(pcname,  sizeof(pcname),  fp);
  fp = popen("pwd | awk -F / '{$NF=\"\"; print $0}' | sed -e 's| |/|g' | tr -d '\n'", "r");  fgets(pwdname0,sizeof(pwdname0),fp);
  fp = popen("pwd | awk -F / '{print $NF}' | tr -d '\n'",                             "r");  fgets(pwdname1,sizeof(pwdname1),fp);
  fp = popen("date +'%Y/%m/%d-%H:%M:%S' | tr -d '\n'",                                "r");  fgets(datename,sizeof(datename),fp);

  outputfile = new std::ofstream;
  outputfile->open(filename);
  (*outputfile) << "BEGIN_OF_FILE_HEADER" << G4endl << G4endl
		<< "command  \t" << commandname << G4endl
		<< "simnum   \t" << event_num << G4endl
		<< "pwd      \t" << username << "@" << pcname << ":" << pwdname0 << pwdname1 << G4endl
		<< "date     \t" << datename << G4endl;

  //set default keta of precision
  lengthPrecision=4;
  energyPrecision=7;
  vectorPrecision=3;
  timingPrecision=4;
}

DMXEventAction::~DMXEventAction(){
  if(outputfile){ outputfile->close();  delete outputfile;  delete eventMessenger; }
}

void DMXEventAction::BeginOfEventAction(const G4Event* evt){
  G4int event_id = evt->GetEventID();
  
  if(event_id==0){
  /*
    const DMXDetectorConstruction* detConst = 
      dynamic_cast<const DMXDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    (*outputfile) << "detector \t" << detConst->GetDetectorName() << G4endl
		  << "parameter\t" << detConst->GetDetectorParamNum();
    for(int i=0; i<detConst->GetDetectorParamNum(); i++)
      (*outputfile) << " " << detConst->GetDetectorParam(i);
    (*outputfile) << G4endl
		  << "zoffset  \t" << detConst->GetDetectorZoffset() << "\tmm" << G4endl
		  << "zrange   \t" << detConst->GetDetectorZrange() << "\tmm" << G4endl
		  << "gas      \t" << detConst->GetGasName() << G4endl
		  << "pressure \t" << detConst->GetGasPressure() << "\tatm" << G4endl
		  << "stepsize \t" << detConst->GetDetectorStepSize() << "\tmm" << G4endl
		  << G4endl
		  << "G4Evt\tseed0\tseed1\tnParticle\tnDetector" << G4endl 
		  << " initParticle\tinitE(MeV)\tinitPos(mm)\tinitMom(vector)" << G4endl
		  << "  detectorName\tnhit\ttotalE(MeV)" << G4endl
		  << "   nStep\tparticle\tene\tde\tprepos\tpostpos\ttime\tprocess" << G4endl
                  << G4endl;
  */
    (*outputfile) << G4endl << G4endl << "END_OF_FILE_HEADER" << G4endl << G4endl << G4endl;
  }

    
  if(     event_id%10000==0) G4cout<<" "<<event_id<<" "<<G4endl;
  else if(event_id%1000 ==0) G4cout<<"|"<<std::flush;
  else if(event_id%100  ==0) G4cout<<"-"<<std::flush;
}

void DMXEventAction::EndOfEventAction(const G4Event* evt){

  // get event_id
  G4int event_id = evt->GetEventID();

  // get seeds
  const DMXPrimaryGeneratorAction* genAction = 
    dynamic_cast<const DMXPrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  const long* seeds = genAction->GetEventSeeds();

  // get initial condition from vertex
  G4int nparticle_pri=0;
  G4String particle_pri[10];
  G4double energy_pri[10];
  G4ThreeVector position_pri[10], direction_pri[10];
  G4int nVertex = evt->GetNumberOfPrimaryVertex();
  for(int i=0; i<nVertex; i++){
    G4PrimaryVertex *vertex = evt->GetPrimaryVertex(i);
    G4int nParticle = vertex->GetNumberOfParticle();
    for(int j=0; j<nParticle; j++){
      G4PrimaryParticle *particle = vertex->GetPrimary(j);
      particle_pri[ nparticle_pri] = particle->GetParticleDefinition()->GetParticleName();
      energy_pri[   nparticle_pri] = particle->GetKineticEnergy();
      position_pri[ nparticle_pri] = vertex->GetPosition();
      direction_pri[nparticle_pri] = particle->GetMomentumDirection();
      nparticle_pri++;
    }
  }

  // get collection ID
  //  const DMXDetectorConstruction* detConst = dynamic_cast<const DMXDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4bool tpc_flag=true, ls_flag=false;
  //  if(detConst->GetDetectorName()=="HP10L_1.4.0") ls_flag=true;
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4int tpcCollID, lsCollID;
  if(tpc_flag) tpcCollID = SDman->GetCollectionID("tpcCollection");
  if(ls_flag)  lsCollID  = SDman->GetCollectionID("lsCollection");

  // get ndetector
  G4int ndetector = 0;
  if(tpc_flag) ndetector++;
  if(ls_flag)  ndetector++;
  if(ndetector==0) return;

  // get nhit, totalEnergy
  G4int     tpcNhit=0, lsNhit=0;
  G4double  tpcTotalEnergy=0, lsTotalEnergy=0;
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(! HCE) return;
  DMXTPCHitsCollection* TPCHC;
  if(tpc_flag){
    TPCHC = (DMXTPCHitsCollection*)(HCE->GetHC(tpcCollID));
    if(TPCHC) tpcNhit = TPCHC->entries();
    for(G4int i=0; i<tpcNhit; i++) tpcTotalEnergy += (*TPCHC)[i]->GetEdep();
  }
  DMXLSHitsCollection*  LSHC;
  if(ls_flag ){
    LSHC =  (DMXLSHitsCollection* )(HCE->GetHC(lsCollID ));
    if(LSHC) lsNhit = LSHC->entries();
    for(G4int i=0; i<lsNhit; i++) lsTotalEnergy += (*LSHC)[i]->GetEdep();
  }
  if(tpcNhit+lsNhit==0) return;

  // save data
  // save: event header
  (*outputfile) << "g" << event_id << "\t" << seeds[0] << "\t" << seeds[1] << "\t" << nparticle_pri << "\t" << ndetector << G4endl;
  // save: particle header
  for(int i=0; i<nparticle_pri; i++)
    (*outputfile) << " " << particle_pri[i] << "\t" 
		  << std::setprecision(energyPrecision) << energy_pri[i]/MeV << "\t" 
		  << std::setprecision(lengthPrecision) << position_pri[i] << "\t" 
		  << std::setprecision(vectorPrecision) << direction_pri[i] << G4endl;
  // save: tpc data
  if(tpc_flag){
    (*outputfile) << "  TPC\t" << tpcNhit << "\t" << std::setprecision(energyPrecision) << tpcTotalEnergy/MeV << G4endl;
    for(int i=0; i<tpcNhit; i++){
      (*outputfile) <<"   n"<<i<<"  "
		    << ((*TPCHC)[i]->GetTrackID()) << "  " << ((*TPCHC)[i]->GetParticle()) << "  "
		    << std::setprecision(energyPrecision) << ((*TPCHC)[i]->GetParticleEnergy()) << "  " << ((*TPCHC)[i]->GetEdep())    << "  "
		    << std::setprecision(lengthPrecision) << ((*TPCHC)[i]->GetPrePos())         << "  " << ((*TPCHC)[i]->GetPostPos()) << "  "
		    << std::setprecision(timingPrecision) << ((*TPCHC)[i]->GetTime())/ns        << "  " << ((*TPCHC)[i]->GetProcess()) << G4endl;
    }
  }
  // save: ls data
  if(ls_flag){
    (*outputfile) << "  LS\t" << lsNhit << "\t" << std::setprecision(energyPrecision) << lsTotalEnergy/MeV << G4endl;
    for(int i=0; i<lsNhit; i++){
      (*outputfile) <<"   n"<<i<<"  "
		    << ((*LSHC)[i]->GetTrackID()) << "  " << ((*LSHC)[i]->GetParticle()) << "  "
		    << std::setprecision(energyPrecision) << ((*LSHC)[i]->GetParticleEnergy()) << "  " << ((*LSHC)[i]->GetEdep())    << "  "
		    << std::setprecision(lengthPrecision) << ((*LSHC)[i]->GetPrePos())         << "  " << ((*LSHC)[i]->GetPostPos()) << "  "
		    << std::setprecision(timingPrecision) << ((*LSHC)[i]->GetTime())/ns        << "  " << ((*LSHC)[i]->GetProcess()) << G4endl;
    }
  }
  // save: kaigyou
  (*outputfile) << G4endl;

}
