#ifndef mcParticleGun_h
#define mcParticleGun_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Navigator.hh"
#include "G4ParticleGunMessenger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4Navigator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include <fstream>
#include <vector>

class G4ParticleTable;
class G4Event;
class mcDetectorConstraction;
class mcParticleGunMessenger;


class mcParticleGun : public G4ParticleGun
{
    friend class mcParticleGunMessenger;
public:
    mcParticleGun(G4String);
    ~mcParticleGun();
    void GeneratePrimaryVertex(G4Event* anEvent);
    
public:
    void GeneratePrimaryVertex_Default(G4Event *evt);
    //void GeneratePrimaryVertex_0nbb(G4Event *evt);
    void GeneratePrimaryVertex_File(G4Event *evt);
    void GeneratePrimaryVertex_nAIST565(G4Event *evt);
    //void GeneratePrimaryVertex_gAIST565_1_0_1(G4Event *evt);
    //void GeneratePrimaryVertex_gAIST565_1_0_3(G4Event *evt);

    void GenerateNeutron(G4PrimaryParticle* neutron[1]);
    void GenerateMuon(G4PrimaryParticle* muon[1]);
    void GenerateFluxNeutron(G4PrimaryParticle* neutron[1], G4ThreeVector vPos);
    void GenerateFluxNeutronSp(G4PrimaryParticle* neutron[1], G4ThreeVector vPos);
    void GenerateFission(G4PrimaryParticle* neutron[1]);
    
    const G4ThreeVector& PutCentre();
    const G4ThreeVector& PutTop();
    const G4ThreeVector& PutFlux();

    G4ThreeVector GeneratePosition();
    G4ThreeVector GeneratePointSource();
    G4double      GenerateTime();
    G4ThreeVector GeneratePointsInVolume();
    G4bool IsSourceConfined(G4ThreeVector);
    G4bool IsSourceInside(G4ThreeVector);
    G4double GenerateCharge(G4ParticleDefinition *pd);

    G4double LogLogInterpolatorCalculate(G4double);
    G4double LogLogInterpolatorCalculateSp(G4double);
    G4double LogLogInterpolatorCalculateFission(G4double);
    
private:
    std::ifstream file_HEPEvt;
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

    G4int verbosityLevel;


    std::ifstream tableFile;
    std::vector<long double> muonE;
    std::vector<long double> muonFlux;
    std::vector<long double> muonPDF;
    
    std::ifstream tableFileSp;
    std::vector<long double> spE;
    std::vector<long double> spFlux;
    std::vector<long double> spPDF;
    
    std::ifstream tableFileFission;
    std::vector<long double> fissionE;
    std::vector<long double> fissionFlux;
    std::vector<long double> fissionPDF;
    
    G4ParticleTable*					particleTable;
    const mcDetectorConstraction*		mcDC;
    
    G4int positionFlag;
    enum{ UserPos=0, Top, Centre,flux};
    G4int particleFlag;
    enum{ User=0, Muon, Neutron,fluxNeutron,fluxNeutronSp,fission};
    G4double monoEnergy;

private:
    mcParticleGunMessenger*	pMessenger;
    G4Navigator *gNavigator;
    
    
};

#endif
