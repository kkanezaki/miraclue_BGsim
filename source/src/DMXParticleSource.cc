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
//
// ParticleSource program
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "DMXParticleSource.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

DMXParticleSource::DMXParticleSource(G4String filename){
  verbosityLevel = 2;
  Mode = "Default";
  SourcePosType = "Volume";
  Shape = "NULL";
  Halfx = 0.0;
  Halfy = 0.0;
  Halfz = 0.0;
  Radius = 0.0;
  CentreCoords = G4ThreeVector(0.0, 0.0, 0.0);
  Confine = false;
  VolName = "NULL";
  Time = 0.0;

  AngDistType = "iso";
  MinTheta = 0.0;
  MaxTheta = pi;
  MinPhi = 0.0;
  MaxPhi = twopi;
  Direction = G4ParticleMomentum(1.0, 0.0, 0.0);
  EnergyDisType = "none";
  Energy = 1*MeV;
  EnergyMin = 0*MeV;
  Particle = "geantino";

  file_HEPEvt.open(filename.data());

  theMessenger = new DMXParticleSourceMessenger(this);
  gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
}

DMXParticleSource::~DMXParticleSource(){ delete theMessenger; }

//// used in DMXParticleSourceMessanger.cc
void DMXParticleSource::SetVerbosity(int vL){ verbosityLevel = vL;  G4cout<<"Verbosity Set to: "<<verbosityLevel<<G4endl; }
void DMXParticleSource::SetMode(G4String mode){ Mode=mode; }
void DMXParticleSource::SetPosDisType(G4String PosType){ SourcePosType = PosType; }
void DMXParticleSource::SetPosDisShape(G4String shapeType){ Shape = shapeType; }
void DMXParticleSource::SetCentreCoords(G4ThreeVector coordsOfCentre){ CentreCoords = coordsOfCentre; }
void DMXParticleSource::SetHalfX(G4double xhalf){ Halfx = xhalf; }
void DMXParticleSource::SetHalfY(G4double yhalf){ Halfy = yhalf; }
void DMXParticleSource::SetHalfZ(G4double zhalf){ Halfz = zhalf; }
void DMXParticleSource::SetRadius(G4double radius){ Radius = radius; }
void DMXParticleSource::SetConfine(G4String VolName_in){
  G4cout<<" ConfineSourceToVolume : "<<VolName_in<<G4endl;
  VolName = VolName_in;
  G4bool found = false;
  G4PhysicalVolumeStore *PVStore = G4PhysicalVolumeStore::GetInstance();
  for(int i=0; i<(G4int)PVStore->size(); i++) if((*PVStore)[i]->GetName() == VolName){ found = true;  break; }
  if(     found  == true                ){ Confine = true;  G4cout<<"Volume "<<VolName<<" exists"<< G4endl; }
  else if(VolName=="NULL"               ){ Confine = false; }
  else if(VolName=="ConfineAllGasRegion"){ Confine = true;  }
  else{ G4cout<<" ### Error: Volume does not exist"<<G4endl;  exit(1); }
}
void DMXParticleSource::SetAngDistType(G4String atype){ AngDistType = atype; }
void DMXParticleSource::SetParticleDirection(G4ParticleMomentum aDirection){ Direction =  aDirection.unit(); }
void DMXParticleSource::SetMinTheta(G4double mintheta){ MinTheta = mintheta; }
void DMXParticleSource::SetMaxTheta(G4double maxtheta){ MaxTheta = maxtheta; }
void DMXParticleSource::SetMinPhi(G4double minphi){ MinPhi = minphi; }
void DMXParticleSource::SetMaxPhi(G4double maxphi){ MaxPhi = maxphi; }
void DMXParticleSource::SetEnergyDisType(G4String DisType){ EnergyDisType = DisType; }
void DMXParticleSource::SetMonoEnergy(G4double menergy){ Energy = menergy; }
void DMXParticleSource::SetMinEnergy(G4double energy){ EnergyMin = energy; }
void DMXParticleSource::SetParticle(G4String particle){ Particle = particle; }
void DMXParticleSource::SetIon(G4String ionparam){
  std::istringstream is(ionparam);
  is >> IonParam[0] >> IonParam[1] >> IonParam[2] >> IonParam[3];
}
/////////////////////////////////////////////////////////////////////////////////

G4ThreeVector DMXParticleSource::GeneratePosition(){
  G4ThreeVector pos;
  if(SourcePosType == "Point"){
    pos = GeneratePointSource();
    if(! IsSourceConfined(pos)){ G4cout<<" ERROR : particle_position is not confined"<<G4endl;  exit(1); }
    if(! IsSourceInside(pos)  ){ G4cout<<" ERROR : particle_position is not inside"  <<G4endl;  exit(1); }
  }else if(SourcePosType == "Volume"){
    G4int LoopCount = 0;
    while(1){
      pos = GeneratePointsInVolume();
      if(IsSourceConfined(pos) && IsSourceInside(pos)) break;
      LoopCount++;
      if(LoopCount==10000){ G4cout<<" ERROR : particle_position is not confined or inside"<<G4endl;  exit(1); }
    }
  }else{ G4cout<<"Error : SourcePosType undefined"<<G4endl;  exit(1); }
  return pos;
}

G4ThreeVector DMXParticleSource::GeneratePointSource(){
  if(verbosityLevel >= 1 && SourcePosType != "Point") G4cout << "Error SourcePosType is not set to Point" << G4endl;
  if(verbosityLevel >= 2) G4cout<<"GeneratePointSource"<<G4endl;
  return CentreCoords;
}

G4ThreeVector DMXParticleSource::GeneratePointsInVolume(){
    if(verbosityLevel >= 1 && SourcePosType != "Volume") G4cout << "Error SourcePosType not Volume" << G4endl;
    if(verbosityLevel >= 2) G4cout<<"GeneratePointsinVolume shape="<<Shape<<G4endl;
    G4ThreeVector RandPos;  G4double x=0., y=0., z=0.;
    if(Shape == "Sphere"){
        x = Radius*2.;  y = Radius*2.;  z = Radius*2.;
        while(((x*x)+(y*y)+(z*z)) > (Radius*Radius)){
            x = G4UniformRand();       y = G4UniformRand();       z = G4UniformRand();
            x = (x*2.*Radius)-Radius;  y = (y*2.*Radius)-Radius;  z = (z*2.*Radius)-Radius;
        }
    }else if(Shape == "Cylinder"){
        x = Radius*2.;  y = Radius*2.;
        while(1){
            x = G4UniformRand();       y = G4UniformRand();       z = G4UniformRand();
            x = (x*2.*Radius)-Radius;  y = (y*2.*Radius)-Radius;  z = (z*2.*Halfz)-Halfz;
            if( x*x+y*y < Radius*Radius ) break;
        }
    }else if(Shape == "Box"){
        x = G4UniformRand();     y = G4UniformRand();     z = G4UniformRand();
        x = (x*2.*Halfx)-Halfx;  y = (y*2.*Halfy)-Halfy;  z = (z*2.*Halfz)-Halfz;
    }else{ G4cout << "Error: Volume Shape (" << Shape << ") Does Not Exist" << G4endl;  exit(1); }
    RandPos.setX(x);  RandPos.setY(y);  RandPos.setZ(z);
    return CentreCoords + RandPos;
}

G4bool DMXParticleSource::IsSourceConfined(G4ThreeVector pos){
  if(Confine == false){ if(verbosityLevel >= 1) G4cout << "Confine is not used" << G4endl;  return(true); }
  G4ThreeVector null(0.,0.,0.);  G4ThreeVector *ptr;  ptr = &null;
  G4VPhysicalVolume *testVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,true);
  G4String testVolName = testVolume->GetName();
  if(verbosityLevel >= 1) G4cout<<" testVolName="<<testVolName<<" VolName="<<VolName<<" pos="<<pos<<G4endl;
  if(testVolName == VolName){
    if(verbosityLevel >= 1) G4cout << "Particle is in volume " << VolName << G4endl;
    return(true);
  }else if(VolName == "ConfineAllGasRegion" &&
	   ( testVolName == "GasRegion_phys"       || // GasRegion in the Chamber
	     testVolName == "DetectionRegion_phys" || // DetectionRegion in GasRegion
	     testVolName == "AnodePTFEHole_phys"   || // AnodePTFEHole in AnodePTFE
	     testVolName == "AnodeCuHole_phys"        // AnodeCuHole in AnodeCu
	     ) ){
    if(verbosityLevel >= 1) G4cout << "Particle is in gas-volume " << testVolName << G4endl;
    return(true);
  }else{
    return(false);
  }
}

G4bool DMXParticleSource::IsSourceInside(G4ThreeVector pos){
  G4ThreeVector null(0.,0.,0.);  G4ThreeVector *ptr;  ptr = &null;
  G4VPhysicalVolume *testVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,true);
  if(! testVolume){
    G4cout<<" ERROR : source position is out of the geometry !!! "<<G4endl;
    return false;
  }else{
    return true;
  }
}

G4double DMXParticleSource::GenerateTime(){ return Time; }

G4ThreeVector DMXParticleSource::GenerateDirection(){ return Direction; }

G4ThreeVector DMXParticleSource::GenerateIsotropicFlux(){
  G4double rndm = G4UniformRand();
  G4double costheta = std::cos(MinTheta) - rndm * (std::cos(MinTheta) - std::cos(MaxTheta));
  G4double sintheta = std::sqrt(1. - costheta*costheta);
  G4double rndm2 = G4UniformRand();
  Phi = MinPhi + (MaxPhi - MinPhi) * rndm2;
  Theta = std::acos(costheta);
  G4double sinphi = std::sin(Phi);  G4double cosphi = std::cos(Phi);
  G4double px = sintheta * cosphi;  G4double py = sintheta * sinphi;  G4double pz = costheta;
  G4double ResMag = std::sqrt((px*px) + (py*py) + (pz*pz));
  px = px/ResMag;  py = py/ResMag;  pz = pz/ResMag;
  // G4cout<<" Phi: min="<<MinPhi<<" Phi="<<Phi<<" max="<<MaxPhi<<G4endl;
  // G4cout<<" Theta: min="<<MinTheta<<" Theta="<<Theta<<" max="<<MaxTheta<<G4endl;
  // G4cout<<" p-xyz ="<<px<<", "<<py<<", "<<pz<<G4endl;

  // rotate (0,0,1) -> Direction
  G4double costheta_org = Direction.z();
  G4double sintheta_org = std::sqrt(1. - costheta_org*costheta_org);

  G4double phi_org = 0;
  if(Direction.x()!=0 || Direction.y()!=0) phi_org = std::acos(Direction.x()/std::sqrt(Direction.x()*Direction.x()+Direction.y()*Direction.y()));
  if(Direction.y()<0) phi_org = twopi - phi_org;
  G4double cosphi_org = std::cos(phi_org);
  G4double sinphi_org = std::sin(phi_org);
  G4double px2 = px *costheta_org - pz *sintheta_org;
  G4double pz2 = pz *costheta_org + px *sintheta_org;
  G4double py3 = py *cosphi_org   - px2*sinphi_org;
  G4double px3 = px2*cosphi_org   + py *sinphi_org;

  // G4cout<<" phi_org="<<phi_org<<" costheta_org="<<costheta_org<<G4endl;
  // G4cout<<" p3-xyz ="<<px3<<", "<<py3<<", "<<pz2<<G4endl;

  //output
  G4ThreeVector dir;  dir.setX(px3);  dir.setY(py3);  dir.setZ(pz2);
  if(verbosityLevel >= 2) G4cout<<"Generating isotropic vector: "<<dir<<G4endl;
  return dir;
}

G4double DMXParticleSource::GenerateMonoEnergy(){ return Energy; }

G4double DMXParticleSource::GenerateFlatEnergy(){ return EnergyMin + G4UniformRand() * (Energy - EnergyMin); }

G4double DMXParticleSource::GenerateNeutronEnergyOf252Cf(){
  G4double nSpectrum_ene[500];
  for(G4int i=0; i<500; i++) nSpectrum_ene[i] = 0.2*(G4double)i*MeV;
  G4double nSpectrum_rate[500] = { 393562857.889990807, 489810635.737437785, 527925875.019090474, 536465038.844141066, 527831784.613912761, 508844716.537064552, 483679601.528013587, 455043435.190826774, 424744948.254112601, 394008939.330489397, 363665016.815575004, 334268085.352385521, 306178770.431590617, 279618900.765984952, 254710731.751436442, 231505180.068130970, 210002415.520912707, 190167015.102001548, 171939177.720850199, 155243044.272337109, 139992867.003912628, 126097567.432055652, 113464079.440340102, 101999772.835808054, 91614179.367196888, 82220189.454486310, 73734847.946605504, 66079847.240171254, 59181793.373289771, 52972303.362770133, 47387978.725736059, 42370289.833692156, 37865397.765069112, 33823934.113031402, 30200754.364643447, 26954676.688944127, 24048215.020964820, 21447313.026108608, 19121083.735127624, 17041558.245631281, 15183445.806733830, 13523906.772438679, 12042339.274136277, 10720179.981663929, 9540718.962876616, 8488928.387543950, 7551304.631916756, 6715723.209050497, 5971305.863857417, 5308299.120500414, 4717963.544899073, 4192472.980263077, 3724823.023540826, 3308748.031408088, 2938645.972733068, 2609510.477853910, 2316869.471563305, 2056729.814931448, 1825527.419887668, 1620082.338979417, 1437558.370316030, 1275426.753941508, 1131433.570461270, 1003570.485476928, 890048.514148112, 789274.508952954, 699830.100459150, 620452.845681397, 550019.361450366, 487530.241233522, 432096.573121385, 382927.894323925, 339321.433614346, 300652.507816014, 266365.951755973, 235968.473206209, 209021.835298429, 185136.778821986, 163967.605785185, 145207.353719625, 128583.497512211, 113854.122131429, 100804.515539416, 89244.136410413, 79003.916065332, 69933.858332995, 61900.904908117, 54787.037237421, 48487.589067465, 42909.746566542, 37971.215420773, 33599.036530456, 29728.533923639, 26302.380283557, 23269.767076869, 20585.667690049, 18210.183249502, 16107.961932838, 14247.683588614, 12601.602382506, 11145.140990921, 9856.530578828, 8716.491436360, 7707.949716908, 6815.786225383, 6026.613655816, 5328.579078450, 4711.188833309, 4165.153304722, 3682.249333644, 3255.198275790, 2877.557936880, 2543.626814778, 2248.359254757, 1987.290280882, 1756.469005807, 1552.399644993, 1371.989271281, 1212.501543324, 1071.515728043, 946.890414197, 836.731382441, 739.363157841, 653.303824608, 577.242730492, 510.020750671, 450.612818436, 398.112463370, 351.718127168, 310.721053499, 274.494571507, 242.484613131, 214.201322690, 189.211633340, 167.132699357, 147.626085902, 130.392629179, 115.167889877, 101.718131622, 89.836763993, 79.341196596, 70.070056821, 61.880729369, 54.647180425, 48.258033640, 42.614868863, 37.630717893, 33.228734495, 29.341018546, 25.907576485, 22.875402314, 20.197665196, 17.832991315, 15.744829093, 13.900888096, 12.272643110, 10.834895820, 9.565387425, 8.444456275, 7.454735316, 6.580884712, 5.809355580, 5.128181201, 4.526792544, 3.995855252, 3.527125619, 3.113323334, 2.748019054, 2.425535073, 2.140857571, 1.889559088, 1.667730046, 1.471918247, 1.299075432, 1.146510080, 1.011845710, 0.892984050, 0.788072508, 0.695475442, 0.613748783, 0.541617623, 0.477956422, 0.421771529, 0.372185744, 0.328424687, 0.289804757, 0.255722505, 0.225645241, 0.199102752, 0.175679979, 0.155010551, 0.136771087, 0.120676150, 0.106473803, 0.093941678, 0.082883513, 0.073126087, 0.064516519, 0.056919880, 0.050217090, 0.044303056, 0.039085030, 0.034481162, 0.030419219, 0.026835458, 0.023673627, 0.020884086, 0.018423031, 0.016251807, 0.014336306, 0.012646430, 0.011155621, 0.009840445, 0.008680224, 0.007656713, 0.006753815, 0.005957324, 0.005254710, 0.004634914, 0.004088180, 0.003605902, 0.003180485, 0.002805230, 0.002474225, 0.002182255, 0.001924720, 0.001697561, 0.001497197, 0.001320470, 0.001164592, 0.001027105, 0.000905841, 0.000798886, 0.000704554, 0.000621354, 0.000547975, 0.000483257, 0.000426178, 0.000375838, 0.000331441, 0.000292286, 0.000257755, 0.000227301, 0.000200444, 0.000176758, 0.000155870, 0.000137449, 0.000121205, 0.000106879, 0.000094246, 0.000083105, 0.000073281, 0.000064617, 0.000056978, 0.000050241, 0.000044300, 0.000039062, 0.000034442, 0.000030369, 0.000026777, 0.000023610, 0.000020818, 0.000018355, 0.000016184, 0.000014269, 0.000012581, 0.000011093, 0.000009780, 0.000008623, 0.000007602, 0.000006703, 0.000005910, 0.000005210, 0.000004594, 0.000004050, 0.000003570, 0.000003148, 0.000002775, 0.000002447, 0.000002157, 0.000001902, 0.000001676, 0.000001478, 0.000001303, 0.000001149, 0.000001013, 0.000000893, 0.000000787, 0.000000694, 0.000000612, 0.000000539, 0.000000475, 0.000000419, 0.000000369, 0.000000326, 0.000000287, 0.000000253, 0.000000223, 0.000000197, 0.000000173, 0.000000153, 0.000000135, 0.000000119, 0.000000105, 0.000000092, 0.000000081, 0.000000072, 0.000000063, 0.000000056, 0.000000049, 0.000000043, 0.000000038, 0.000000034, 0.000000030, 0.000000026, 0.000000023, 0.000000020, 0.000000018, 0.000000016, 0.000000014, 0.000000012, 0.000000011, 0.000000010, 0.000000008, 0.000000007, 0.000000007, 0.000000006, 0.000000005, 0.000000004, 0.000000004, 0.000000003, 0.000000003, 0.000000003, 0.000000002, 0.000000002, 0.000000002, 0.000000002, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000000 };
  G4double total_rate=0;  for(G4int i=0; i<500; i++) total_rate += nSpectrum_rate[i];
  G4double div_rate = G4UniformRand()*total_rate;
  G4double nEnergy=0;
  for(G4int i=0; i<500; i++){
    if(div_rate<nSpectrum_rate[i]){
      nEnergy=nSpectrum_ene[i];
      break;
    }
    div_rate -= nSpectrum_rate[i];
  }
  return nEnergy +0.2*G4UniformRand() *MeV;
}

G4double DMXParticleSource::GenerateGammaEnergyOf133Ba(){
  G4double line_energy[9] = {53.161, 79.6139, 80.9971, 160.613, 223.234, 276.398, 302.853, 356.017, 383.851};
  G4double line_prob[9]   = { 2.199,  2.62,   34.06,     0.645,   0.450,   7.164,  18.33,   62.05,    8.94 };
  G4double prob_max=0;  for(G4int i=0; i<9; i++) prob_max+=line_prob[i];
  G4double rand = G4UniformRand()*prob_max;
  G4double gEnergy=0;
  for(G4int i=0; i<9; i++){
    if(rand<line_prob[i]){
      gEnergy = line_energy[i] *keV;
      break;
    }else{
      rand-=line_prob[i];
    }
  }
  return gEnergy;
}

G4ParticleDefinition* DMXParticleSource::GenerateParticle(){
  if(Particle !="ion"){//particle
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pd = particleTable->FindParticle(Particle);
    if(pd != NULL){
      return pd;
    }else{ G4cout<<" ERROR : Particle has unusual value"<<G4endl;  exit(1); }
  }else{//ion
    G4IonTable *ionTable = G4IonTable::GetIonTable();
    G4ParticleDefinition* ion =  ionTable->GetIon((G4int)IonParam[0], (G4int)IonParam[1], IonParam[3]);
    if(ion!=0){
      return ion;
    }else{//ion not found
      G4cout <<"Ion with Z="<<(G4int)IonParam[0]<<" A="<<(G4int)IonParam[1]<<" Charge="
	     <<(G4int)IonParam[2]<<" ExciteEnergy="<<IonParam[3]<<" is not be defined" << G4endl;
      exit(1);
    }
  }
  G4cout << "No particle has been defined!" << G4endl;
  exit(1);
}

G4double DMXParticleSource::GenerateCharge(G4ParticleDefinition *pd){
  if(pd->GetParticleName() != "ion") return pd->GetPDGCharge();
  else                               return IonParam[2]*eplus;
}


////////////////////////////////////////////////////////////////////////////////
void DMXParticleSource::GeneratePrimaryVertex(G4Event *evt){
  if(     Mode == "Default"        ) GeneratePrimaryVertex_Default(evt);
  else if(Mode == "0nbb"           ) GeneratePrimaryVertex_0nbb(evt);
  else if(Mode == "File"           ) GeneratePrimaryVertex_File(evt); //HEPEvt
  else if(Mode == "nAIST565"       ) GeneratePrimaryVertex_nAIST565(evt);
  else if(Mode == "gAIST565_1.0.1" ) GeneratePrimaryVertex_gAIST565_1_0_1(evt);
  else if(Mode == "gAIST565_1.0.3" ) GeneratePrimaryVertex_gAIST565_1_0_3(evt);
  else{ G4cout<<" ERROR : Mode has unusual value"<<G4endl;  return; }
  if(verbosityLevel > 1) G4cout << " Primary Vetex generated "<< G4endl;
}

// Default mode : generate single particle from single vertex.
void DMXParticleSource::GeneratePrimaryVertex_Default(G4Event *evt){
  G4ThreeVector pos = GeneratePosition();
  G4double time = GenerateTime();
  G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);

  G4ThreeVector dir;
  if(     AngDistType == "iso"      ) dir = GenerateIsotropicFlux();
  else if(AngDistType == "direction") dir = GenerateDirection();
  else{ G4cout<<"Error: AngDistType has unusual value"<<G4endl;  return; }
  G4double kin_ene;
  if(     EnergyDisType == "Mono"   ) kin_ene = GenerateMonoEnergy();
  else if(EnergyDisType == "Flat"   ) kin_ene = GenerateFlatEnergy();
  else if(EnergyDisType == "n_252Cf") kin_ene = GenerateNeutronEnergyOf252Cf();
  else if(EnergyDisType == "g_133Ba") kin_ene = GenerateGammaEnergyOf133Ba();
  else{ G4cout<<"Error: EnergyDisType has unusual value"<<G4endl;  return; }
  G4ThreeVector pol = G4ThreeVector(0.0, 0.0, 0.0);
  G4ParticleDefinition *pd = GenerateParticle();
  G4double charge = GenerateCharge(pd);
  G4double mass = pd->GetPDGMass();
  G4double tot_ene = kin_ene + mass;
  G4double pmom = std::sqrt(tot_ene*tot_ene-mass*mass);
  G4double px = pmom*dir.getX(), py = pmom*dir.getY(), pz = pmom*dir.getZ();
  if(verbosityLevel >= 1) G4cout<<" Particle: "<<pd->GetParticleName()<<G4endl
				<<"   Energy: "<<kin_ene<<G4endl
				<<" Position: "<<pos<<G4endl
				<<"Direction: "<<dir<<G4endl;
  G4PrimaryParticle* particle = new G4PrimaryParticle(pd,px,py,pz);
  particle->SetMass( mass );
  particle->SetCharge( charge );
  particle->SetPolarization(pol.x(), pol.y(), pol.z());

  vertex->SetPrimary( particle );
  evt->AddPrimaryVertex( vertex );
}

// File mode : generate particle from file (list of initial particle)
void DMXParticleSource::GeneratePrimaryVertex_File(G4Event *evt){
  G4ThreeVector pos = GeneratePosition();
  G4double time = GenerateTime();
  G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);

  G4int NHEP;  // number of entries
  file_HEPEvt >> NHEP;
  if( file_HEPEvt.eof() ){
    G4cout<<" end of HEPEvt file "<<G4endl;
    return;
  }
  for( G4int IHEP=0; IHEP<NHEP; IHEP++ ){
    G4int ISTHEP;   // status code
    G4int IDHEP;    // PDG code
    G4int JDAHEP1;  // first daughter
    G4int JDAHEP2;  // last daughter
    G4double PHEP1; // px in GeV
    G4double PHEP2; // py in GeV
    G4double PHEP3; // pz in GeV
    G4double PHEP5; // mass in GeV
    file_HEPEvt >> ISTHEP >> IDHEP >> JDAHEP1 >> JDAHEP2 >> PHEP1 >> PHEP2 >> PHEP3 >> PHEP5;
    if(verbosityLevel>=2) G4cout<<" IHEP/NHEP="<<IHEP<<"/"<<NHEP<<" ISTHEP="<<ISTHEP<<" IDHEP="<<IDHEP
				<<" P="<<PHEP1<<","<<PHEP2<<","<<PHEP3<<" M="<<PHEP5<<G4endl;
    G4PrimaryParticle* particle;
    if(IDHEP < 1000000000){ // not ion
      particle = new G4PrimaryParticle( IDHEP );
      particle->SetMomentum(PHEP1*GeV, PHEP2*GeV, PHEP3*GeV );
      particle->SetMass( PHEP5*GeV );
    }else{
      double Z = floor((IDHEP-1000000000)/10000);
      double A = (int)floor((IDHEP-1000000000)/10)%1000;
      double E = (int)(IDHEP-1000000000)%10;
      G4IonTable *ionTable = G4IonTable::GetIonTable();
      G4ParticleDefinition* ion =  ionTable->GetIon(Z, A, E);
      //      G4cout<<" PDG_code = "<<ion->GetPDGEncoding()<<G4endl;
      particle = new G4PrimaryParticle(ion, PHEP1*GeV, PHEP2*GeV, PHEP3*GeV);
    }
    vertex->SetPrimary( particle );
  }
  evt->AddPrimaryVertex( vertex );
}

// 0nbb mode : generate 2-electrons from single vertex
void DMXParticleSource::GeneratePrimaryVertex_0nbb(G4Event *evt){
  G4ThreeVector pos = GeneratePosition();
  G4double time = GenerateTime();
  G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);

  // particle (e-)
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pd = particleTable->FindParticle("e-");

  // energy angle
  G4double TotalEnergy = Energy;//2.458*MeV;
  G4double me = pd->GetPDGMass();
  if(TotalEnergy < 2.0*me){ G4cout<<" ERROR : Energy is too small for 0nbb generator"<<G4endl;  return; }
  G4double trial=1.0, pdf=0.0, cosTheta, eEnergy[2];
  G4ThreeVector eDirection_org[2], eDirection[2];
  while(trial>pdf){
    eEnergy[0] = G4UniformRand()*TotalEnergy;
    eEnergy[1] = TotalEnergy - eEnergy[0];
    cosTheta = G4UniformRand()*2.0-1.0;
    pdf = (1.-sqrt((1.-pow(me/(eEnergy[0]+me),2))*(1.-pow(me/(eEnergy[1]+me),2)))*cosTheta)/2.;
    trial = G4UniformRand();
  }
  G4double phi   = G4UniformRand()*2.0*M_PI;
  G4double theta = acos(G4UniformRand()*2.0-1.0);
  G4double psi   = G4UniformRand()*2.0*M_PI;
  eDirection_org[0] = G4ThreeVector(1.0, 0.0, 0.0);
  eDirection_org[1] = G4ThreeVector(cosTheta, sqrt(1.0-cosTheta*cosTheta), 0.0);
  for(int i=0; i<2; i++) eDirection[i] = eDirection_org[i].rotate(phi,theta,psi);

  // set vertex
  for(int i=0; i<2; i++){
    G4ThreeVector pol = G4ThreeVector(0.0, 0.0, 0.0);
    G4double charge = GenerateCharge(pd);
    G4double mass   = pd->GetPDGMass();
    G4double tot_ene = eEnergy[i] + mass;
    G4double pmom = std::sqrt(tot_ene*tot_ene-mass*mass);
    G4double px = pmom*eDirection[i].x();
    G4double py = pmom*eDirection[i].y();
    G4double pz = pmom*eDirection[i].z();
    if(verbosityLevel >= 1) G4cout<<" Particle: "<<pd->GetParticleName()<<G4endl
				  <<"   Energy: "<<eEnergy[i]<<G4endl
				  <<" Position: "<<pos<<G4endl
				  <<"Direction: "<<eDirection[i]<<G4endl;
    G4PrimaryParticle* particle = new G4PrimaryParticle(pd,px,py,pz);
    particle->SetMass( mass );
    particle->SetCharge( charge );
    particle->SetPolarization(pol.x(), pol.y(), pol.z());
    vertex->SetPrimary( particle );
  }
  evt->AddPrimaryVertex( vertex );
}

void DMXParticleSource::GeneratePrimaryVertex_nAIST565(G4Event *evt){
  G4ThreeVector pos = GeneratePosition();
  G4double time = GenerateTime();
  G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);

  double mass0=938.767, mass1=6535.254, mass2=939.5493, mass3=6536.116;//MeV
  double K0=2.3;//, Q=-1.63;//MeV
  double E0=K0+mass0;
  double p0=sqrt(E0*E0-mass0*mass0);
  double s_param=(E0+mass1)*(E0+mass1)-p0*p0;
  double beta=p0/(E0+mass1);
  double gama=1./sqrt(1.-beta*beta);
  double p2_cm = 1./2./sqrt(s_param)*sqrt(mass3*mass3*mass3*mass3 +mass2*mass2*mass2*mass2 +s_param*s_param
					  -2.*(mass2*mass2*mass3*mass3 +mass2*mass2*s_param +mass3*mass3*s_param));
  double E2_cm = sqrt(p2_cm*p2_cm + mass2*mass2);

  double costheta_cm = G4UniformRand()*2.0-1.0;
  double sintheta_cm = sqrt(1.-costheta_cm*costheta_cm);
  double p2_cm_z = p2_cm*costheta_cm;
  double p2_cm_x = p2_cm*sintheta_cm;
  double p2_x = p2_cm_x;
  double p2_z = gama*p2_cm_z + E2_cm*beta*gama;
  double k2 = sqrt(p2_z*p2_z+p2_x*p2_x + mass2*mass2)-mass2;
  
  double costheta = p2_z/sqrt(p2_z*p2_z+p2_x*p2_x);
  double sintheta = p2_x/sqrt(p2_z*p2_z+p2_x*p2_x);
  double phi = G4UniformRand()*2.0*M_PI;
  double dir_z = costheta;
  double dir_x = sintheta*cos(phi);
  double dir_y = sintheta*sin(phi);

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pd = particleTable->FindParticle("neutron");
  G4ThreeVector pol = G4ThreeVector(0.0, 0.0, 0.0);
  G4double charge = GenerateCharge(pd);
  G4double mass = pd->GetPDGMass();
  G4double tot_ene = k2 + mass;
  G4double pmom = std::sqrt(tot_ene*tot_ene-mass*mass);
  G4double px = pmom*dir_x, py = pmom*dir_y, pz = pmom*dir_z;
  if(verbosityLevel >= 1) G4cout<<" Particle: "<<pd->GetParticleName()<<G4endl
				<<"   Energy: "<<k2<<G4endl
				<<" Position: "<<pos<<G4endl
				<<"Direction: "<<dir_x<<","<<dir_y<<","<<dir_z<<G4endl;
  G4PrimaryParticle* particle = new G4PrimaryParticle(pd,px,py,pz);
  particle->SetMass( mass );
  particle->SetCharge( charge );
  particle->SetPolarization(pol.x(), pol.y(), pol.z());

  vertex->SetPrimary( particle );
  evt->AddPrimaryVertex( vertex );
}

void DMXParticleSource::GeneratePrimaryVertex_gAIST565_1_0_1(G4Event *evt){
  double r_outer_sphere = 2.0*m;

  double costheta = 2.0*G4UniformRand()-1.0;
  double theta = acos(costheta);
  double phi = G4UniformRand()*2.0*M_PI;
  double x_pos = r_outer_sphere*sin(theta)*cos(phi);
  double y_pos = r_outer_sphere*sin(theta)*sin(phi);
  double z_pos = r_outer_sphere*cos(theta);
  G4ThreeVector pos = G4ThreeVector(x_pos, y_pos, z_pos);
  costheta = 2.0*G4UniformRand()-1.0;
  theta = acos(costheta);
  phi = G4UniformRand()*2.0*M_PI;
  double dir_x = sin(theta)*cos(phi);
  double dir_y = sin(theta)*sin(phi);
  double dir_z = cos(theta);

  G4double time = GenerateTime();
  G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pd = particleTable->FindParticle("gamma");
  G4ThreeVector pol = G4ThreeVector(0.0, 0.0, 0.0);
  G4double charge = GenerateCharge(pd);
  G4double mass = pd->GetPDGMass();

  double ene_hist_wid = 50; // keV
  const int ene_hist_num = 200;
  double ene_hist[ene_hist_num]={31, 373, 363, 275, 217, 192, 128, 82, 74, 74, 158, 45, 33, 34, 38, 21, 42, 30, 19, 9, 13, 18, 7, 12, 8, 20, 8, 6, 17, 11, 10, 11, 25, 10, 20, 19, 6, 12, 9, 14, 7, 12, 11, 13, 636, 4, 1, 4, 5, 3, 7, 1, 6, 5, 2, 7, 5, 1, 0, 6, 3, 6, 3, 4, 2, 8, 1, 1, 7, 2, 30, 3, 0, 4, 2, 0, 2, 2, 4, 1, 0, 1, 0, 0, 11, 3, 4, 0, 7, 5, 0, 1, 1, 2, 2, 1, 5, 2, 29, 0, 3, 0, 1, 3, 1, 4, 1, 0, 0, 0, 1, 3, 2, 0, 1, 1, 0, 1, 29, 0, 19, 0, 1, 0, 0, 2, 1, 8, 0, 0, 1, 2, 5, 1, 3, 1, 12, 0, 1, 1, 0, 0, 5, 1, 0, 13, 0, 7, 0, 0, 3, 0, 247, 0, 1, 0, 16, 0, 27, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 24, 14, 0, 0, 0, 0, 0, 0, 27, 0, 35, 0, 0, 0, 0, 0, 23, 0, 0, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0};

  double ene_hist_sum[ene_hist_num]={};
  for(int i=0; i<ene_hist_num; i++){
    for(int j=0; j<=i; j++){
      ene_hist_sum[i] += ene_hist[j];
    }
  }
  double param = G4UniformRand()*ene_hist_sum[ene_hist_num-1];
  double ene = 0;
  for(int i=0; i<ene_hist_num; i++){
    if(ene_hist_sum[i]<=param && param <ene_hist_sum[i+1]){
      ene = ene_hist_wid*( (double)i+G4UniformRand());
      break;
    }
  }

  G4double tot_ene = ene*0.001 + mass; // MeV
  G4double pmom = std::sqrt(tot_ene*tot_ene-mass*mass);
  G4double px = pmom*dir_x, py = pmom*dir_y, pz = pmom*dir_z;
  if(verbosityLevel >= 1) G4cout<<" Particle: "<<pd->GetParticleName()<<G4endl
				<<"   Energy: "<<ene<<G4endl
				<<" Position: "<<pos<<G4endl
				<<"Direction: "<<dir_x<<","<<dir_y<<","<<dir_z<<G4endl;
  G4PrimaryParticle* particle = new G4PrimaryParticle(pd,px,py,pz);
  particle->SetMass( mass );
  particle->SetCharge( charge );
  particle->SetPolarization(pol.x(), pol.y(), pol.z());

  vertex->SetPrimary( particle );
  evt->AddPrimaryVertex( vertex );
}

void DMXParticleSource::GeneratePrimaryVertex_gAIST565_1_0_3(G4Event *evt){
  double r_outer_sphere = 2.0*m;

  double costheta = 2.0*G4UniformRand()-1.0;
  double theta = acos(costheta);
  double phi = G4UniformRand()*2.0*M_PI;
  double x_pos = r_outer_sphere*sin(theta)*cos(phi);
  double y_pos = r_outer_sphere*sin(theta)*sin(phi);
  double z_pos = r_outer_sphere*cos(theta);
  G4ThreeVector pos = G4ThreeVector(x_pos, y_pos, z_pos);
  costheta = 2.0*G4UniformRand()-1.0;
  theta = acos(costheta);
  phi = G4UniformRand()*2.0*M_PI;
  double dir_x = sin(theta)*cos(phi);
  double dir_y = sin(theta)*sin(phi);
  double dir_z = cos(theta);

  G4double time = GenerateTime();
  G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, time);

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pd = particleTable->FindParticle("gamma");
  G4ThreeVector pol = G4ThreeVector(0.0, 0.0, 0.0);
  G4double charge = GenerateCharge(pd);
  G4double mass = pd->GetPDGMass();

  double ene_hist_wid = 50; // keV
  const int ene_hist_num = 200;
  double ene_hist[ene_hist_num]={29, 285, 201, 118, 96, 66, 42, 27, 29, 22, 28, 12, 9, 13, 10, 12, 7, 5, 6, 2, 8, 3, 7, 5, 8, 9, 6, 3, 6, 1, 4, 1, 1, 3, 1, 32, 0, 3, 2, 5, 5, 7, 3, 2, 353, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 3, 3, 0, 0, 0, 0, 0, 0, 0, 3, 16, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 3, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 4, 1, 1, 0, 0, 14, 0, 0, 0, 2, 0, 1, 7, 1, 1, 1, 0, 1, 2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 3, 1, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  double ene_hist_sum[ene_hist_num]={};
  for(int i=0; i<ene_hist_num; i++){
    for(int j=0; j<=i; j++){
      ene_hist_sum[i] += ene_hist[j];
    }
  }
  double param = G4UniformRand()*ene_hist_sum[ene_hist_num-1];
  double ene = 0;
  for(int i=0; i<ene_hist_num; i++){
    if(ene_hist_sum[i]<=param && param <ene_hist_sum[i+1]){
      ene = ene_hist_wid*( (double)i+G4UniformRand());
      break;
    }
  }

  G4double tot_ene = ene*0.001 + mass; // MeV
  G4double pmom = std::sqrt(tot_ene*tot_ene-mass*mass);
  G4double px = pmom*dir_x, py = pmom*dir_y, pz = pmom*dir_z;
  if(verbosityLevel >= 1) G4cout<<" Particle: "<<pd->GetParticleName()<<G4endl
				<<"   Energy: "<<ene<<G4endl
				<<" Position: "<<pos<<G4endl
				<<"Direction: "<<dir_x<<","<<dir_y<<","<<dir_z<<G4endl;
  G4PrimaryParticle* particle = new G4PrimaryParticle(pd,px,py,pz);
  particle->SetMass( mass );
  particle->SetCharge( charge );
  particle->SetPolarization(pol.x(), pol.y(), pol.z());

  vertex->SetPrimary( particle );
  evt->AddPrimaryVertex( vertex );
}
