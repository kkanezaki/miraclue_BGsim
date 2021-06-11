// DMXDetectorConstruction.hh
// 20201124 Kiseki Nakamura

#include "G4VUserDetectorConstruction.hh"
class DMXTPCSD;
class DMXLSSD;

class DMXDetectorConstruction : public G4VUserDetectorConstruction{

public:
  DMXDetectorConstruction(G4String);
  ~DMXDetectorConstruction();

public:
    G4String theDetectorName;
    G4int    theDetectorParamNum;
    G4double theDetectorParam[10];
    G4double theGasPressure;
    G4int    theGasTypeNum;
    G4int    theGasElementNum;
    G4String theGasType[10];
    G4double theGasRatio[10];
    G4String theGasName;
    G4double theDetectorZoffset;
    G4double theDetectorZrange;
  
private:
  virtual G4VPhysicalVolume* Construct();
  G4String config_filename;
};






/*

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

#ifndef DMXDetectorConstruction_h
#define DMXDetectorConstruction_h 1
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class DMXTPCSD;
class DMXLSSD;

class DMXDetectorConstruction : public G4VUserDetectorConstruction{
 public:
  DMXDetectorConstruction(G4String);
  ~DMXDetectorConstruction();

 public:
  G4VPhysicalVolume* Construct();
  G4String theDetectorName;
  G4int    theDetectorParamNum;
  G4double theDetectorParam[10];
  G4double theGasPressure;
  G4int    theGasTypeNum;
  G4int    theGasElementNum;
  G4String theGasType[10];
  G4double theGasRatio[10];
  G4String theGasName;
  G4double theDetectorZoffset;
  G4double theDetectorZrange;

  G4String GetDetectorName(void)    const{ return theDetectorName;        };
  G4double GetDetectorParamNum(void)const{ return theDetectorParamNum;    };
  G4double GetDetectorParam(int i)  const{ return theDetectorParam[i];    };
  G4String GetGasName(void)         const{ return theGasName;             };
  G4double GetGasPressure(void)     const{ return theGasPressure;         };
  G4double GetDetectorZoffset(void) const{ return theDetectorZoffset;     };
  G4double GetDetectorZrange(void)  const{ return theDetectorZrange;      };
  G4double GetDetectorStepSize(void)const{ return theStepSizeForDetector; };
 
private:
  void DefineMaterials();
  G4UserLimits* theUserLimitsForRoom;
  G4UserLimits* theUserLimitsForDetector;
  G4double theTimeCutForRoom;
  G4double theTimeCutForDetector;
  G4double theStepSizeForRoom;
  G4double theStepSizeForDetector;
  G4double theMinEkineForRoom;
  G4double theMinEkineForDetector;

  G4Cache<DMXTPCSD*> TPCSD;
  G4Cache<DMXLSSD*>  LSSD;
#include "DMXDetectorMaterial.ihh"

};
#endif

*/