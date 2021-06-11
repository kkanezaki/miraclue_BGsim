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
// ScintSD (scintillator sensitive detector definition) program
// --------------------------------------------------------------

#include "DMXLSSD.hh"
#include "DMXLSHit.hh"
#include "DMXDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4ios.hh"

DMXLSSD::DMXLSSD(G4String name, G4String HCname) :G4VSensitiveDetector(name){
  collectionName.insert(HCname);
}
DMXLSSD::~DMXLSSD(){ }

void DMXLSSD::Initialize(G4HCofThisEvent*){
  lsCollection = new DMXLSHitsCollection(SensitiveDetectorName,collectionName[0]);
  HitID = -1;
}

G4bool DMXLSSD::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4ParticleDefinition* particleType = aStep->GetTrack()->GetDefinition();
  G4String particleName = particleType->GetParticleName();
  G4String processName = "none";
  if(aStep->GetTrack()->GetCreatorProcess() != 0) processName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();

  G4double stepl = 0.;
  if (particleType->GetPDGCharge() != 0.) stepl = aStep->GetStepLength();
  if ((edep==0.)&&(stepl==0.)) return false;      

  // fill in hit
  DMXLSHit* newHit = new DMXLSHit();
  newHit->SetEdep(edep);//edep
  newHit->SetPostPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetPrePos(aStep->GetPreStepPoint()->GetPosition());
  newHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
  newHit->SetParticle(particleName);
  newHit->SetParticleEnergy(aStep->GetPreStepPoint()->GetKineticEnergy() );
  newHit->SetProcess(processName);
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  HitID = lsCollection->insert(newHit);
  
  return true;
}

void DMXLSSD::EndOfEvent(G4HCofThisEvent* HCE){
  G4String HCname = collectionName[0];
  static G4int HCID = -1;
  if(HCID<0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
  HCE->AddHitsCollection(HCID, lsCollection);

  G4int nHits = lsCollection->entries();
  if (verboseLevel>=1) G4cout<<"     LS collection: "<<nHits<<" hits"<<G4endl;
  if (verboseLevel>=2) lsCollection->PrintAllHits();
}
