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
// ParticleSourceMessenger program
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>               

#include "DMXParticleSourceMessenger.hh"
#include "DMXParticleSource.hh"

#include "G4SystemOfUnits.hh"
#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"
#include "G4Tokenizer.hh"

DMXParticleSourceMessenger::DMXParticleSourceMessenger(DMXParticleSource *fPtclGun) : fParticleGun(fPtclGun){
  particleTable = G4ParticleTable::GetParticleTable();

  // create directory
  gunDirectory = new G4UIdirectory("/dmx/gun/");
  gunDirectory->SetGuidance("Particle Source control commands.");

  // list available particles
  listCmd = new G4UIcmdWithoutParameter("/dmx/gun/List",this);
  listCmd->SetGuidance("List available particles.");
  listCmd->SetGuidance(" Invoke G4ParticleTable.");

  // verbosity
  verbosityCmd = new G4UIcmdWithAnInteger("/dmx/gun/verbose",this);
  verbosityCmd->SetGuidance("Set Verbose level for gun");
  verbosityCmd->SetGuidance(" 0 : Silent");
  verbosityCmd->SetGuidance(" 1 : Limited information");
  verbosityCmd->SetGuidance(" 2 : Detailed information");
  verbosityCmd->SetParameterName("level",false);
  verbosityCmd->SetRange("level>=0 && level <=2");

  // source mode
  modeCmd = new G4UIcmdWithAString("/dmx/gun/mode",this);
  modeCmd->SetGuidance("Sets source mode.");
  modeCmd->SetGuidance("Default, 0nbb, ...");
  modeCmd->SetParameterName("Mode",true,true);
  modeCmd->SetDefaultValue("Default");
  //  modeCmd->SetCandidates("Default 0nbb File nAIST565 gAIST565_1.0.1 gAIST565_1.0.3");

  // source distribution type
  typeCmd = new G4UIcmdWithAString("/dmx/gun/type",this);
  typeCmd->SetGuidance("Sets source distribution type.");
  typeCmd->SetGuidance("Either Point or Volume");
  typeCmd->SetParameterName("DisType",true,true);
  typeCmd->SetDefaultValue("Point");
  typeCmd->SetCandidates("Point Volume");

  // source shape
  shapeCmd = new G4UIcmdWithAString("/dmx/gun/shape",this);
  shapeCmd->SetGuidance("Sets source shape type.");
  shapeCmd->SetParameterName("Shape",true,true);
  shapeCmd->SetDefaultValue("NULL");
  shapeCmd->SetCandidates("Sphere Cylinder Box");

  // centre coordinates
  centreCmd = new G4UIcmdWith3VectorAndUnit("/dmx/gun/centre",this);
  centreCmd->SetGuidance("Set centre coordinates of source.");
  centreCmd->SetParameterName("X","Y","Z",true,true);
  centreCmd->SetDefaultUnit("cm");
  centreCmd->SetUnitCandidates("nm um mm cm m km");

  // half x length of source
  halfxCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/halfx",this);
  halfxCmd->SetGuidance("Set x half length of source.");
  halfxCmd->SetParameterName("Halfx",true,true);
  halfxCmd->SetDefaultUnit("cm");
  halfxCmd->SetUnitCandidates("nm um mm cm m km");

  // half y length of source
  halfyCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/halfy",this);
  halfyCmd->SetGuidance("Set y half length of source.");
  halfyCmd->SetParameterName("Halfy",true,true);
  halfyCmd->SetDefaultUnit("cm");
  halfyCmd->SetUnitCandidates("nm um mm cm m km");

  // half height of source
  halfzCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/halfz",this);
  halfzCmd->SetGuidance("Set z half length of source.");
  halfzCmd->SetParameterName("Halfz",true,true);
  halfzCmd->SetDefaultUnit("cm");
  halfzCmd->SetUnitCandidates("nm um mm cm m km");

  // radius of source  
  radiusCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/radius",this);
  radiusCmd->SetGuidance("Set radius of source.");
  radiusCmd->SetParameterName("Radius",true,true);
  radiusCmd->SetDefaultUnit("cm");
  radiusCmd->SetUnitCandidates("nm um mm cm m km");

  // confine to volume
  confineCmd = new G4UIcmdWithAString("/dmx/gun/confine",this);
  confineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  confineCmd->SetGuidance("usage: confine VolName");
  confineCmd->SetParameterName("VolName",true,true);
  confineCmd->SetDefaultValue("NULL");

  // angular distribution
  angtypeCmd = new G4UIcmdWithAString("/dmx/gun/angtype",this);
  angtypeCmd->SetGuidance("Sets angular source distribution type");
  angtypeCmd->SetGuidance("Possible variables are: iso direction");
  angtypeCmd->SetParameterName("AngDis",true,true);
  angtypeCmd->SetDefaultValue("iso");
  angtypeCmd->SetCandidates("iso direction iso_part");

  // particle direction
  directionCmd = new G4UIcmdWith3Vector("/dmx/gun/direction",this);
  directionCmd->SetGuidance("Set momentum direction.");
  directionCmd->SetGuidance("Direction needs not to be a unit vector.");
  directionCmd->SetParameterName("Px","Py","Pz",true,true); 
  directionCmd->SetRange("Px != 0 || Py != 0 || Pz != 0");

  // particle min theta
  minthetaCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/mintheta",this);
  minthetaCmd->SetGuidance("Set minimum of theta of particle direction.");
  minthetaCmd->SetParameterName("MinTheta",true,true); 
  minthetaCmd->SetRange("0 <= MinTheta <= 180");
  minthetaCmd->SetUnitCandidates("deg rad");

  // particle max theta
  maxthetaCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/maxtheta",this);
  maxthetaCmd->SetGuidance("Set maximum of theta of particle direction.");
  maxthetaCmd->SetParameterName("MaxTheta",true,true); 
  maxthetaCmd->SetRange("0 <= MaxTheta <= 180");
  maxthetaCmd->SetUnitCandidates("deg rad");


  // particle min phi
  minphiCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/minphi",this);
  minphiCmd->SetGuidance("Set minimum of phi of particle direction.");
  minphiCmd->SetParameterName("MinPhi",true,true); 
  minphiCmd->SetRange("0 <= MinPhi <= 180");
  minphiCmd->SetUnitCandidates("deg rad");

  // particle max phi
  maxphiCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/maxphi",this);
  maxphiCmd->SetGuidance("Set maximum of phi of particle direction.");
  maxphiCmd->SetParameterName("MaxPhi",true,true); 
  maxphiCmd->SetRange("0 <= MaxPhi <= 360");
  maxphiCmd->SetUnitCandidates("deg rad");

  // energy distribution
  energytypeCmd = new G4UIcmdWithAString("/dmx/gun/energytype",this);
  energytypeCmd->SetGuidance("Sets energy distribution type");
  energytypeCmd->SetGuidance("Possible variables are: Mono, n_252Cf");
  energytypeCmd->SetParameterName("EnergyDis",true,true);
  energytypeCmd->SetDefaultValue("Mono");
  energytypeCmd->SetCandidates("Mono Flat n_252Cf g_133Ba");

  // particle energy
  energyCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/energy",this);
  energyCmd->SetGuidance("Set kinetic energy.");
  energyCmd->SetParameterName("Energy",true,true);
  energyCmd->SetDefaultUnit("GeV");
  energyCmd->SetUnitCategory("Energy");
  energyCmd->SetUnitCandidates("eV keV MeV GeV TeV");

  // particle energy_min
  energyminCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/energy_min",this);
  energyminCmd->SetGuidance("Set minimum kinetic energy for Flat distribution.");
  energyminCmd->SetParameterName("EnergyMin",true,true);
  energyminCmd->SetDefaultUnit("GeV");
  energyminCmd->SetUnitCategory("Energy");
  energyminCmd->SetUnitCandidates("eV keV MeV GeV TeV");

  // set particle  
  particleCmd = new G4UIcmdWithAString("/dmx/gun/particle",this);
  particleCmd->SetGuidance("Set particle to be generated.");
  particleCmd->SetGuidance(" (geantino is default)");
  particleCmd->SetGuidance(" (ion can be specified for shooting ions)");
  particleCmd->SetParameterName("particleName",true);
  particleCmd->SetDefaultValue("geantino");
  G4String candidateList; 
  for(G4int i=0; i<particleTable->entries(); i++){
    candidateList += particleTable->GetParticleName(i);
    candidateList += " ";
  }
  candidateList += "ion ";
  particleCmd->SetCandidates(candidateList);

  // ion 
  ionCmd = new G4UIcommand("/dmx/gun/ion",this);
  ionCmd->SetGuidance("Set properties of ion to be generated.");
  ionCmd->SetGuidance("[usage] /gun/ion Z A Q E");
  ionCmd->SetGuidance("        Z:(int) AtomicNumber");
  ionCmd->SetGuidance("        A:(int) AtomicMass");
  ionCmd->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
  ionCmd->SetGuidance("        E:(double) Excitation energy (in keV)");
  G4UIparameter* param;
  param = new G4UIparameter("Z",'i',false);
  param->SetDefaultValue("1");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("A",'i',false);
  param->SetDefaultValue("1");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("Q",'i',true);
  param->SetDefaultValue("0");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("E",'d',true);
  param->SetDefaultValue("0.0");
  ionCmd->SetParameter(param);
}

DMXParticleSourceMessenger::~DMXParticleSourceMessenger() {
  delete listCmd;  delete verbosityCmd;  delete modeCmd;
  delete typeCmd;  delete shapeCmd;  delete centreCmd;  delete halfzCmd;  delete radiusCmd;  delete confineCmd;
  delete angtypeCmd;  delete directionCmd;
  delete minthetaCmd;  delete maxthetaCmd;  delete minphiCmd;  delete maxphiCmd;
  delete energytypeCmd;  delete energyCmd;  delete energyminCmd;
  delete ionCmd;  delete particleCmd;
  delete gunDirectory;
}

void DMXParticleSourceMessenger::SetNewValue(G4UIcommand *command, G4String newValues){
  if(     command == listCmd )     particleTable->DumpTable();
  else if(command == verbosityCmd) fParticleGun->SetVerbosity(verbosityCmd->GetNewIntValue(newValues));
  else if(command == modeCmd)      fParticleGun->SetMode(newValues);
  else if(command == typeCmd)      fParticleGun->SetPosDisType(newValues);
  else if(command == shapeCmd)     fParticleGun->SetPosDisShape(newValues);
  else if(command == centreCmd)    fParticleGun->SetCentreCoords(centreCmd->GetNew3VectorValue(newValues));
  else if(command == halfxCmd)     fParticleGun->SetHalfX(halfxCmd->GetNewDoubleValue(newValues));
  else if(command == halfyCmd)     fParticleGun->SetHalfY(halfyCmd->GetNewDoubleValue(newValues));
  else if(command == halfzCmd)     fParticleGun->SetHalfZ(halfzCmd->GetNewDoubleValue(newValues));
  else if(command == radiusCmd)    fParticleGun->SetRadius(radiusCmd->GetNewDoubleValue(newValues));
  else if(command == confineCmd)   fParticleGun->SetConfine(newValues);
  else if(command == angtypeCmd)   fParticleGun->SetAngDistType(newValues);
  else if(command == directionCmd) fParticleGun->SetParticleDirection(directionCmd->GetNew3VectorValue(newValues));
  else if(command == minthetaCmd)  fParticleGun->SetMinTheta(minthetaCmd->GetNewDoubleValue(newValues));
  else if(command == maxthetaCmd)  fParticleGun->SetMaxTheta(maxthetaCmd->GetNewDoubleValue(newValues));
  else if(command == minphiCmd)    fParticleGun->SetMinPhi(minphiCmd->GetNewDoubleValue(newValues));
  else if(command == maxphiCmd)    fParticleGun->SetMaxPhi(maxphiCmd->GetNewDoubleValue(newValues));
  else if(command == energytypeCmd)fParticleGun->SetEnergyDisType(newValues);
  else if(command == energyCmd)    fParticleGun->SetMonoEnergy(energyCmd->GetNewDoubleValue(newValues));
  else if(command == energyminCmd) fParticleGun->SetMinEnergy(energyminCmd->GetNewDoubleValue(newValues));
  else if(command == particleCmd)  fParticleGun->SetParticle(newValues);
  else if(command == ionCmd)       fParticleGun->SetIon(newValues);
  else G4cout << "Error entering command" << G4endl;
}
