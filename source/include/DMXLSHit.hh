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
// DetectorConstruction header
// --------------------------------------------------------------

#ifndef DMXLSHit_h
#define DMXLSHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class DMXLSHit : public G4VHit{
public:
  DMXLSHit();
  ~DMXLSHit();
  DMXLSHit(const DMXLSHit&);
  const DMXLSHit& operator=(const DMXLSHit&);
  int operator==(const DMXLSHit&) const;
  inline void* operator new(size_t);
  inline void  operator delete(void*);

public:
  void SetEdep           (G4double de)       { edep = de; };
  void SetPostPos        (G4ThreeVector xyz) { postpos = xyz; };
  void SetPrePos         (G4ThreeVector xyz) { prepos = xyz; };
  void SetParticle       (G4String name)     { particleName = name; };
  void SetParticleEnergy (G4double e1)       { particleEnergy = e1; };
  void SetTime           (G4double t2)       { time = t2; };
  void SetProcess        (G4String name)     { processName = name; };
  void SetTrackID        (G4int val)         { trackID = val; };

  G4double GetEdep()                         { return edep; };      
  G4ThreeVector GetPrePos()                  { return prepos; };
  G4ThreeVector GetPostPos()                 { return postpos; };
  G4String GetParticle()                     { return particleName;};
  G4double GetParticleEnergy()               { return particleEnergy;};
  G4double GetTime()                         { return time; };      
  G4String GetProcess()                      { return processName;};
  G4int    GetTrackID()                      { return trackID;};

private:
  G4double      edep;
  G4ThreeVector prepos, postpos;
  G4double      time;
  G4String      particleName, processName;
  G4double      particleEnergy;
  G4int         trackID;
};

typedef G4THitsCollection<DMXLSHit> DMXLSHitsCollection;
extern G4Allocator<DMXLSHit> DMXLSHitAllocator;

inline void* DMXLSHit::operator new(size_t){
  void* aHit;
  aHit = (void*) DMXLSHitAllocator.MallocSingle();
  return aHit;
}

inline void DMXLSHit::operator delete(void* aHit){
  DMXLSHitAllocator.FreeSingle((DMXLSHit*) aHit);
}

#endif
