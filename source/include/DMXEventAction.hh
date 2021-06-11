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
// EventAction header
// --------------------------------------------------------------

#ifndef DMXEventAction_h
#define DMXEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "DMXTPCHit.hh"
#include "DMXLSHit.hh"

class DMXPrimaryGeneratorAction;
class DMXEventActionMessenger;

class DMXEventAction : public G4UserEventAction {

 public:
  DMXEventAction(G4String, G4String, G4int);
  virtual ~DMXEventAction();
  virtual void EndOfEventAction(const G4Event*);
  virtual void BeginOfEventAction(const G4Event*);
  void SetLengthPrecision(G4int val){ lengthPrecision = val; };
  void SetEnergyPrecision(G4int val){ energyPrecision = val; };
  void SetVectorPrecision(G4int val){ vectorPrecision = val; };
  void SetTimingPrecision(G4int val){ timingPrecision = val; };

 private:
  std::ofstream *outputfile;
  DMXEventActionMessenger* eventMessenger;
  G4double lengthPrecision, energyPrecision, vectorPrecision, timingPrecision;
};
#endif
