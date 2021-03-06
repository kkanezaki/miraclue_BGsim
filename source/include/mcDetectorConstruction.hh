#ifndef mcDetectorConstruction_h
#define mcDetectorConstruction_h 1

#include "mcAnalyzer.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
//class G4Orb;
//class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class G4UserLimits;
class mcDetectorMessenger;


class mcDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    
    mcDetectorConstruction();
    ~mcDetectorConstruction();
    
public:
    
    void SetSensorMaterial(G4String);
    void SetMaxStep(G4double);
    void SetMagField(G4double);
    void SetGasBoxSize(G4double);
    void SetNeutronShieldMaterial(G4String);
    void SetNeutronShieldSize(G4double);
    void SetNeutronShieldType(G4String);
    void SetGammaShield1Thickness(G4double);
    void SetGammaShield2Thickness(G4double);
    void SetScintiMaterial(G4String);
    
    G4VPhysicalVolume* Construct();
    
    void UpdateGeometry();
    
public:
    
    G4double GetWorldRadius()  const    {return WorldRadius;};
    
    const G4Material* GetSensorMaterial()  const {return sensorMaterial;};
    const G4Material* GetnShieldMaterial() const {return nShieldMaterial;};
    
    G4double    GetMaxStep()      const {return maxStep;};
    
    G4double    GetFieldValue()      const {return fieldValue;};
    
    const G4VPhysicalVolume* GetphysiWorld() const {return physWorld;};
    const G4VPhysicalVolume* GetSensor()     const {return physSensor;};
    
    void SetAnalyzer(mcAnalyzer*);
    
    
    
private:
    G4Material*        defaultMaterial;
    G4Material*        sensorMaterial;
    G4Material*        nShieldMaterial;
    G4String           ScintiMaterial;
    
    G4double           WorldRadius;
    
    //G4Orb*             solidWorld;
    G4Box*             solidWorld;
    G4LogicalVolume*   logicWorld;
    G4VPhysicalVolume* physWorld;

    G4LogicalVolume*   logicLab;
    G4LogicalVolume*   logicScinti;
    
    //G4Tubs*            solidSensor;
    G4Box*             solidSensor;
    G4LogicalVolume*   logicSensor;
    G4VPhysicalVolume* physSensor;
    
    G4UniformMagField* magField;      //pointer to the magnetic field
    G4double           fieldValue;

    G4double           nShieldSize;
    G4double           nShieldThetaMax;
    G4String           nShieldShape;

    G4double           l_gas;
    G4double           t_fiducial;
    G4double           d;
    G4double           theta_cone;
    G4double           l_cone;

    G4double           gShield1Thickness;
    G4double           gShield2Thickness;
    
    G4UserLimits*      pUserLimits;    //pointer to the UserLimits
    G4double           maxStep;          // max step length
    mcDetectorMessenger* detectorMessenger;  //pointer to the Messenger
    
    void DefineMaterials();
    mcAnalyzer* analyzer;

private:
    G4LogicalVolume* logicMother;
    void ConstructLaboratory();
    void ConstructBeamShield();
    void ConstructSphereBeamShield();
    void ConstructCubeBeamShield();
    void ConstructChamber();
    void ConstructGammaShield1();
    void ConstructTestShield();
    void ConstructScintillator();

};

#endif

