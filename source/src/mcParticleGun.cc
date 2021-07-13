#include "mcParticleGun.hh"
#include "mcParticleGunMessenger.hh"
#include "mcDetectorMessenger.hh"
#include "mcDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RotationMatrix.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>

using CLHEP::RandFlat;
/*
mcParticleGun::mcParticleGun()
:G4ParticleGun(1),
particleTable(G4ParticleTable::GetParticleTable()),
positionFlag(0),
particleFlag(0),
monoEnergy(1.0*keV),pMessenger(0)

{
    pMessenger = new mcParticleGunMessenger(this);
    
    
    tableFile.open("muonTable.dat");
    double f1, f2;
    while(tableFile >> f1 >> f2){
        muonE.push_back(f1);
        muonFlux.push_back(f2);
    }
    
    double totalFlux = 0;
    for(decltype(muonFlux.size()) i=0;i<muonFlux.size();i++){
        totalFlux += muonFlux.at(i);
    }
    double muonFactor = 1.0 / totalFlux;
    totalFlux = 0.0;
    
    for(decltype(muonFlux.size()) i=0;i<muonFlux.size();i++){
        totalFlux += muonFlux.at(i);
        muonPDF.push_back(muonFactor* totalFlux);
    }
    
    tableFileSp.open("total_table.dat");
    //double f1, f2;
    while(tableFileSp >> f1 >> f2){
        spE.push_back(f1);
        spFlux.push_back(f2);
    }
    
    totalFlux = 0;
    for(decltype(spFlux.size()) i=0;i<spFlux.size();i++){
        totalFlux += spFlux.at(i);
    }
    double spFactor = 1.0 / totalFlux;
    totalFlux = 0.0;
    
    for(decltype(spFlux.size()) i=0;i<spFlux.size();i++){
        totalFlux += spFlux.at(i);
        spPDF.push_back(spFactor* totalFlux);
    }
    
    tableFileFission.open("cf.dat");
    //double f1, f2;
    while(tableFileFission >> f1 >> f2){
        fissionE.push_back(f1);
        fissionFlux.push_back(f2);
    }
    
    totalFlux = 0;
    for(decltype(fissionFlux.size()) i=0;i<fissionFlux.size();i++){
        totalFlux += fissionFlux.at(i);
    }
    double fissionFactor = 1.0 / totalFlux;
    totalFlux = 0.0;
    
    for(decltype(fissionFlux.size()) i=0;i<fissionFlux.size();i++){
        totalFlux += fissionFlux.at(i);
        fissionPDF.push_back(fissionFactor* totalFlux);
    }
    
}
*/

mcParticleGun::mcParticleGun(G4String filename){
    verbosityLevel = 2;
    Mode = "nAIST565";
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

    pMessenger = new mcParticleGunMessenger(this);
}

mcParticleGun::~mcParticleGun()
{
    delete pMessenger;
}

/////////////////////////////////miraclue////////////////////////////////////////

G4ThreeVector mcParticleGun::GeneratePointSource(){
    if(verbosityLevel >= 1 && SourcePosType != "Point") G4cout << "Error SourcePosType is not set to Point" << G4endl;
    if(verbosityLevel >= 2) G4cout<<"GeneratePointSource"<<G4endl;
    return CentreCoords;
}

G4bool mcParticleGun::IsSourceInside(G4ThreeVector pos){
    G4ThreeVector null(0.,0.,0.);  G4ThreeVector *ptr;  ptr = &null;
    G4VPhysicalVolume *testVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,true);
    if(! testVolume){
        G4cout<<" ERROR : source position is out of the geometry !!! "<<G4endl;
        return false;
    }else{
        return true;
    }
}

G4bool mcParticleGun::IsSourceConfined(G4ThreeVector pos){
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

G4ThreeVector mcParticleGun::GeneratePointsInVolume(){
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

G4double mcParticleGun::GenerateTime(){ return Time; }

G4ThreeVector mcParticleGun::GeneratePosition(){
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



G4double mcParticleGun::GenerateCharge(G4ParticleDefinition *pd){
    if(pd->GetParticleName() != "ion") return pd->GetPDGCharge();
    else                               return IonParam[2]*eplus;
}




void mcParticleGun::GeneratePrimaryVertex(G4Event* anEvent){
    if(     Mode == "File"           ) GeneratePrimaryVertex_File(anEvent); //HEPEvt
    else if(Mode == "nAIST565"       ) GeneratePrimaryVertex_nAIST565(anEvent);
    //else if(Mode == "Default"        ) GeneratePrimaryVertex_Default(anEvent);
    //else if(Mode == "gAIST565_1.0.1" ) GeneratePrimaryVertex_gAIST565_1_0_1(anEvent);
    //else if(Mode == "gAIST565_1.0.3" ) GeneratePrimaryVertex_gAIST565_1_0_3(anEvent);
    else{ G4cout<<" ERROR : Mode has unusual value"<<G4endl;  return; }
}


void mcParticleGun::GeneratePrimaryVertex_nAIST565(G4Event *anEvent){
    //G4ThreeVector pos = GeneratePosition();
    G4ThreeVector pos = G4ThreeVector(-1*m,0,0);
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
    anEvent->AddPrimaryVertex( vertex );
}

void mcParticleGun::GeneratePrimaryVertex_File(G4Event *anEvent){
}


/*
void mcParticleGun::GeneratePrimaryVertex(G4Event* anEvent)
{
    G4ThreeVector vPos;
    if(positionFlag == Top)				vPos = PutTop();
    else if (positionFlag == Centre)	vPos = PutCentre();
    else if (positionFlag == flux)  	vPos = PutFlux();
    else 								vPos = particle_position;
    
    G4PrimaryVertex* vertex = new G4PrimaryVertex(vPos,particle_time);
    
    if (particleFlag == Muon){
        G4PrimaryParticle* particle[1]={0};
        GenerateMuon(particle);
        vertex->SetPrimary(particle[0]);
    }else if(particleFlag == Neutron){
        G4PrimaryParticle* particle[1]={0};
        GenerateNeutron(particle);
        vertex->SetPrimary(particle[0]);
    }else if(particleFlag == fluxNeutron){
        G4PrimaryParticle* particle[1]={0};
        GenerateFluxNeutron(particle,vPos);
        vertex->SetPrimary(particle[0]);
    }else if(particleFlag == fluxNeutronSp){
        G4PrimaryParticle* particle[1]={0};
        GenerateFluxNeutronSp(particle,vPos);
        vertex->SetPrimary(particle[0]);
    }else if(particleFlag == fission){
        G4PrimaryParticle* particle[1]={0};
        GenerateFission(particle);
        vertex->SetPrimary(particle[0]);
    }else{
        G4double mass = particle_definition->GetPDGMass();	
        for(G4int i=0; i<NumberOfParticlesToBeGenerated; i++){
            G4PrimaryParticle* particle = new G4PrimaryParticle(particle_definition);
            particle->SetKineticEnergy(particle_energy);
            particle->SetMass(mass);
            particle->SetMomentumDirection(particle_momentum_direction);
            particle->SetCharge(particle_charge);
            particle->SetPolarization(particle_polarization.x(),
                                      particle_polarization.y(),
                                      particle_polarization.z());
            vertex->SetPrimary(particle);
            
        }
    }
    anEvent->AddPrimaryVertex(vertex);
}
*/

void mcParticleGun::GenerateNeutron(G4PrimaryParticle* neutron[1])
{
    G4ParticleDefinition* pD = particleTable->FindParticle("neutron");
    neutron[0] = new G4PrimaryParticle(pD);
    neutron[0]->SetKineticEnergy(1.0 * GeV);
    neutron[0]->SetMomentumDirection(G4ThreeVector(0,0,-1.0));
}

void mcParticleGun::GenerateFluxNeutron(G4PrimaryParticle* neutron[1],G4ThreeVector vPos)
{
    G4ParticleDefinition* pD = particleTable->FindParticle("neutron");
    neutron[0] = new G4PrimaryParticle(pD);
    //Set energy
    neutron[0]->SetKineticEnergy(monoEnergy);
    
    G4double angleX = asin(2*RandFlat::shoot(0.0,1.0)-1.0);
    G4double delta = RandFlat::shoot(0.0,pi);
    
    G4ThreeVector momVec;
    momVec.setX(-1*vPos.getX()/vPos.getR());
    momVec.setY(-1*vPos.getY()/vPos.getR());
    momVec.setZ(-1*vPos.getZ()/vPos.getR());
    
    G4double theta = momVec.getTheta();
    G4double phi = momVec.getPhi();
    momVec.rotateZ(-1*phi);
    momVec.rotateY(-1*theta);
    momVec.rotateY(angleX);
    momVec.rotateZ(delta);
    momVec.rotateY(theta);
    momVec.rotateZ(phi);
    
    neutron[0]->SetMomentumDirection(momVec);
}

void mcParticleGun::GenerateFluxNeutronSp(G4PrimaryParticle* neutron[1],G4ThreeVector vPos)
{
    G4ParticleDefinition* pD = particleTable->FindParticle("neutron");
    neutron[0] = new G4PrimaryParticle(pD);
    //Set energy
    G4double prob = RandFlat::shoot(0.0,1.0);
    neutron[0]->SetKineticEnergy( LogLogInterpolatorCalculateSp(prob) * MeV); //real spectrum
    
    G4double angleX = asin(2*RandFlat::shoot(0.0,1.0)-1.0);
    G4double delta = RandFlat::shoot(0.0,pi);
    
    G4ThreeVector momVec;
    momVec.setX(-1*vPos.getX()/vPos.getR());
    momVec.setY(-1*vPos.getY()/vPos.getR());
    momVec.setZ(-1*vPos.getZ()/vPos.getR());
    
    G4double theta = momVec.getTheta();
    G4double phi = momVec.getPhi();
    momVec.rotateZ(-1*phi);
    momVec.rotateY(-1*theta);
    momVec.rotateY(angleX);
    momVec.rotateZ(delta);
    momVec.rotateY(theta);
    momVec.rotateZ(phi);
    
    neutron[0]->SetMomentumDirection(momVec);
}

void mcParticleGun::GenerateFission(G4PrimaryParticle* neutron[1])
{
    G4ParticleDefinition* pD = particleTable->FindParticle("neutron");
    neutron[0] = new G4PrimaryParticle(pD);
    //Set energy
    G4double prob = RandFlat::shoot(0.0,1.0);
    neutron[0]->SetKineticEnergy( LogLogInterpolatorCalculateFission(prob) * MeV); //real spectrum
    
    G4double px, py, pz; 
    G4double cs, sn, phi;
    cs    =  RandFlat::shoot(-1.0,1.0);
    sn    =  std::sqrt((1.0-cs)*(1.0+cs));   
    phi   =  RandFlat::shoot(0., CLHEP::twopi);   
    px    =  sn*std::cos(phi);
    py    =  sn*std::sin(phi);
    pz    =  cs; 
    
    neutron[0]->SetMomentumDirection(G4ThreeVector(px, py, pz));
}
void mcParticleGun::GenerateMuon(G4PrimaryParticle* muon[1])
{
    G4ParticleDefinition* pD = particleTable->FindParticle("mu-");
    muon[0] = new G4PrimaryParticle(pD);
    muon[0]->SetMomentumDirection(G4ThreeVector(0,0,-1.0));//real shoot
    G4double prob = RandFlat::shoot(0.0,1.0);
    muon[0]->SetKineticEnergy( LogLogInterpolatorCalculate(prob) * GeV); //real spectrum
}

G4double mcParticleGun::LogLogInterpolatorCalculate(G4double x){
    
    G4double value = 0;
    
    if (x > 0.9999999 || x < 0.0) return 0.0;
    
    if(x < muonPDF.at(1) || x == 0.0){
        value = 0.0;
    }else if( x > 1.0){
        value = 0.0;
    }else {
        size_t i = 0;
        for(i=0;i<muonPDF.size();i++){
            if(x < muonPDF.at(i)){ break; }
        }
        G4double e1 = muonPDF.at(i-1);
        G4double e2 = muonPDF.at(i);
        G4double d1 = muonE.at(i-1);
        G4double d2 = muonE.at(i);
        
        value = (std::log10(d1)*std::log10(e2/x) + std::log10(d2)*std::log10(x/e1)) / std::log10(e2/e1);
        value = std::pow(10.0,value);
    }
    return value;
}



G4double mcParticleGun::LogLogInterpolatorCalculateSp(G4double x){
    
    G4double value = 0;
    
    if (x > 0.9999999 || x < 0.0) return 0.0;
    
    if(x < spPDF.at(1) || x == 0.0){
        value = 0.0;
    }else if( x > 1.0){
        value = 0.0;
    }else {
        size_t i = 0;
        for(i=0;i<spPDF.size();i++){
            if(x < spPDF.at(i)){ break; }
        }
        G4double e1 = spPDF.at(i-1);
        G4double e2 = spPDF.at(i);
        G4double d1 = spE.at(i-1);
        G4double d2 = spE.at(i);
        
        value = (std::log10(d1)*std::log10(e2/x) + std::log10(d2)*std::log10(x/e1)) / std::log10(e2/e1);
        value = std::pow(10.0,value);
    }
    return value;
}

G4double mcParticleGun::LogLogInterpolatorCalculateFission(G4double x){
    
    G4double value = 0;
    
    if (x > 0.9999999 || x < 0.0) return 0.0;
    
    if(x < fissionPDF.at(1) || x == 0.0){
        value = 0.0;
    }else if( x > 1.0){
        value = 0.0;
    }else {
        size_t i = 0;
        for(i=0; i<fissionPDF.size();i++){
            if(x < fissionPDF.at(i)){ break; }
        }
        G4double e1 = fissionPDF.at(i-1);
        G4double e2 = fissionPDF.at(i);
        G4double d1 = fissionE.at(i-1);
        G4double d2 = fissionE.at(i);
        
        value = (std::log10(d1)*std::log10(e2/x) + std::log10(d2)*std::log10(x/e1)) / std::log10(e2/e1);
        value = std::pow(10.0,value);
    }
    return value;
}

const G4ThreeVector& mcParticleGun::PutCentre(){
    static G4ThreeVector vPos(0.0,0.0,0.0);
    return vPos;
}

const G4ThreeVector& mcParticleGun::PutTop(){
    static G4ThreeVector vPos(0.0,0.0,49.9*cm);
    return vPos;
}
const G4ThreeVector& mcParticleGun::PutFlux(){
    G4double radius = 40*cm;
    G4double phi = RandFlat::shoot(0.0,twopi);
    G4double theta = RandFlat::shoot(0.0,pi);
    G4double posX =radius*sin(theta)*cos(phi);
    G4double posY =radius*sin(theta)*sin(phi);
    G4double posZ =radius*cos(theta);
    static G4ThreeVector vPos(0,0,0);
    vPos.setX(posX);
    vPos.setY(posY);
    vPos.setZ(posZ);
    return vPos;
}
