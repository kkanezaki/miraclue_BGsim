#include "mcPrimaryGeneratorAction.hh"
#include "mcParticleGun.hh"
#include "mcDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4MuonMinus.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


mcPrimaryGeneratorAction::mcPrimaryGeneratorAction(const mcDetectorConstruction* mcDC)
:particleTable(G4ParticleTable::GetParticleTable())
//,mcDetector(mcDC)
{
    G4String filename; //とりあえず
    particleGun = new mcParticleGun(filename);
    G4ParticleDefinition* particle = particleTable->FindParticle("geantino");
    particleGun->SetParticleDefinition(particle);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(1.0,0.0,0.0));
    particleGun->SetParticleEnergy(electron_mass_c2);
    particleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));
    
    
}

/*
mcPrimaryGeneratorAction::mcPrimaryGeneratorAction(const mcDetectorConstruction* mcDC, G4int s0, G4int s1, G4String filename)
:particleTable(G4ParticleTable::GetParticleTable())
{
    particleGun = new mcParticleGun(filename);
    seeds[0] = s0;
    seeds[1] = s1;
    CLHEP::HepRandom::setTheSeeds(seeds);
}
*/


mcPrimaryGeneratorAction::~mcPrimaryGeneratorAction()
{
    delete particleGun;
}



void mcPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    particleGun->GeneratePrimaryVertex(anEvent);
    
}



