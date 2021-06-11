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
// ScintHit (scintillator sensitive detector definition) program
// --------------------------------------------------------------

#include "DMXTPCHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>

G4Allocator<DMXTPCHit> DMXTPCHitAllocator;

DMXTPCHit::DMXTPCHit(){
  edep           = 0.;
  prepos         = G4ThreeVector(0., 0., 0.);
  postpos        = G4ThreeVector(0., 0., 0.);
  time           = 0.;
  particleEnergy = 0.;
}
DMXTPCHit::~DMXTPCHit(){;}

DMXTPCHit::DMXTPCHit(const DMXTPCHit& right) : G4VHit(right){
  edep           = right.edep;
  prepos         = right.prepos;
  postpos        = right.postpos;
  time           = right.time;
  particleName   = right.particleName;
  particleEnergy = right.particleEnergy;
  processName    = right.processName;
  trackID        = right.trackID;
}

const DMXTPCHit& DMXTPCHit::operator=(const DMXTPCHit& right){
  edep           = right.edep;
  prepos         = right.prepos;
  postpos        = right.postpos;
  time           = right.time;
  particleName   = right.particleName;
  particleEnergy = right.particleEnergy;
  processName    = right.processName;
  trackID        = right.trackID;
  return *this;
}

int DMXTPCHit::operator==(const DMXTPCHit& right) const{
  return (this==&right) ? 1 : 0;
}
