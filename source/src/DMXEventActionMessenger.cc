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
// EventActionMessenger program
// --------------------------------------------------------------

#include "DMXEventActionMessenger.hh"
#include "DMXEventAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "globals.hh"

DMXEventActionMessenger::DMXEventActionMessenger(DMXEventAction* EvAct):eventAction(EvAct){
  dmxDirectory = new G4UIdirectory("/dmx/");
  dmxDirectory->SetGuidance("DM Example commands.");

  SetLengthPrecisionCmd = new G4UIcmdWithAnInteger("/dmx/lengthPrecision",this);
  SetLengthPrecisionCmd->SetGuidance("Set precision of length for output file. ex) 1 for 0.1mm, 4 for 0.1um");

  SetEnergyPrecisionCmd = new G4UIcmdWithAnInteger("/dmx/energyPrecision",this);
  SetEnergyPrecisionCmd->SetGuidance("Set precision of energy for output file. ex) 1 for 0.1MeV, 7 for 0.1eV");

  SetVectorPrecisionCmd = new G4UIcmdWithAnInteger("/dmx/vectorPrecision",this);
  SetVectorPrecisionCmd->SetGuidance("Set precision of vector for output file. ex) 1 for 0.1, 3 for 0.001");

  SetTimingPrecisionCmd = new G4UIcmdWithAnInteger("/dmx/timingPrecision",this);
  SetTimingPrecisionCmd->SetGuidance("Set precision of timing for output file. ex) 1 for 0.1us, 4 for 0.1ns");
}

DMXEventActionMessenger::~DMXEventActionMessenger(){
  delete SetLengthPrecisionCmd;
  delete SetEnergyPrecisionCmd;
  delete SetVectorPrecisionCmd;
  delete SetTimingPrecisionCmd;
}

void DMXEventActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
  if(command == SetLengthPrecisionCmd){ eventAction->SetLengthPrecision(SetLengthPrecisionCmd->GetNewIntValue(newValue));}
  if(command == SetEnergyPrecisionCmd){ eventAction->SetEnergyPrecision(SetEnergyPrecisionCmd->GetNewIntValue(newValue));}
  if(command == SetVectorPrecisionCmd){ eventAction->SetVectorPrecision(SetVectorPrecisionCmd->GetNewIntValue(newValue));}
  if(command == SetTimingPrecisionCmd){ eventAction->SetTimingPrecision(SetTimingPrecisionCmd->GetNewIntValue(newValue));}
}
