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
// ParticleSource header
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////

#ifndef DMXParticleSource_h
#define DMXParticleSource_h 1

#include <fstream>
#include "G4VPrimaryGenerator.hh"
#include "G4Navigator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"

#include "DMXParticleSourceMessenger.hh"


class DMXParticleSource : public G4VPrimaryGenerator {

 public:
  DMXParticleSource(G4String);
  ~DMXParticleSource ();
  void GeneratePrimaryVertex(G4Event *evt);

 public:
  // for ParticleSourceMessanger
  void SetVerbosity(G4int);
  void SetMode(G4String);
  void SetFile(G4String);
  void SetIfile(G4int);
  void SetPosDisType(G4String);
  void SetPosDisShape(G4String);
  void SetCentreCoords(G4ThreeVector);
  void SetHalfX(G4double);
  void SetHalfY(G4double);
  void SetHalfZ(G4double);
  void SetRadius(G4double);
  void SetConfine(G4String);
  void SetAngDistType(G4String);
  void SetParticleDirection(G4ThreeVector);
  void SetMinTheta(G4double);
  void SetMaxTheta(G4double);
  void SetMinPhi(G4double);
  void SetMaxPhi(G4double);
  void SetEnergyDisType(G4String);
  void SetMonoEnergy(G4double);
  void SetMinEnergy(G4double);
  void SetParticle(G4String);
  void SetIon(G4String);

  // for ParticleSource (this)
  G4ThreeVector GeneratePosition();
  G4ThreeVector GeneratePointSource();
  G4ThreeVector GeneratePointsInVolume();
  G4bool IsSourceConfined(G4ThreeVector);
  G4bool IsSourceInside(G4ThreeVector);
  G4double      GenerateTime();
  G4ParticleMomentum GenerateDirection();
  G4ParticleMomentum GenerateIsotropicFlux();
  G4double   GenerateMonoEnergy();
  G4double   GenerateFlatEnergy();
  G4double   GenerateNeutronEnergyOf252Cf();
  G4double   GenerateGammaEnergyOf133Ba();
  G4ParticleDefinition* GenerateParticle();
  G4double GenerateCharge(G4ParticleDefinition *pd);
 
  // divided
  void GeneratePrimaryVertex_Default(G4Event *evt);
  void GeneratePrimaryVertex_0nbb(G4Event *evt);
  void GeneratePrimaryVertex_File(G4Event *evt);
  void GeneratePrimaryVertex_nAIST565(G4Event *evt);
  void GeneratePrimaryVertex_gAIST565_1_0_1(G4Event *evt);
  void GeneratePrimaryVertex_gAIST565_1_0_3(G4Event *evt);
  //  void GeneratePrimaryVertex_MigdalXe1sAuger(G4Event *evt);
  //  void GeneratePrimaryVertex_WithXray(G4Event *evt);
  //  void GeneratePrimaryVertex_InterConv(G4Event *evt);

 private:
  G4String Mode, File;
  G4String SourcePosType, Shape;
  G4double Halfx, Halfy, Halfz;
  G4double Radius;
  G4ThreeVector CentreCoords;
  G4bool Confine;
  G4String VolName;
  G4double Time;
  G4String AngDistType;
  G4ThreeVector Direction;
  G4double MinTheta, MaxTheta, MinPhi, MaxPhi;
  G4double Theta, Phi;
  G4String EnergyDisType;
  G4double Energy, EnergyMin;
  G4String Particle;
  G4double IonParam[4];
  std::ifstream file_HEPEvt;

  // Verbose
  G4int verbosityLevel;

 private:
  DMXParticleSourceMessenger *theMessenger;
  G4Navigator *gNavigator;

};
#endif

