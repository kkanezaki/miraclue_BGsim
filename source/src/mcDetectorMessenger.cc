#include "mcDetectorMessenger.hh"
#include "mcDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"


mcDetectorMessenger::mcDetectorMessenger(mcDetectorConstruction* mcDet)
:mcDetector(mcDet)
{ 
    usrDir = new G4UIdirectory("/usr/");
    usrDir->SetGuidance("UI commands of this example");
    detDir = new G4UIdirectory("/usr/det/");
    detDir->SetGuidance("UI commands for detector setup");
    nShieldDir = new G4UIdirectory("/usr/nShield/");
    nShieldDir->SetGuidance("UI commands for neutron shield");

    
    MaterialCmd = new G4UIcmdWithAString("/usr/det/setMaterial",this);
    MaterialCmd->SetGuidance("Select Material of the sensor");
    MaterialCmd->SetParameterName("choice",false);
    MaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/usr/det/setField",this);
    MagFieldCmd->SetGuidance("Define magnetic field.");
    MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
    MagFieldCmd->SetParameterName("Bz",false);
    MagFieldCmd->SetUnitCategory("Magnetic flux density");
    MagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
    
    MaxStepCmd = new G4UIcmdWithADoubleAndUnit("/usr/det/setMaxStep",this);
    MaxStepCmd->SetGuidance("Set MaxStep ");
    MaxStepCmd->SetParameterName("MaxStep",false);
    MaxStepCmd->SetRange("MaxStep>0.");
    MaxStepCmd->SetUnitCategory("Length");    
    MaxStepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    ShieldRadiusCmd = new G4UIcmdWithADoubleAndUnit("/usr/nShield/radius",this);
    ShieldRadiusCmd->SetGuidance("neutron shield radius");
    ShieldRadiusCmd->SetParameterName("neutron shield radius",false);
    ShieldRadiusCmd->SetUnitCategory("Length");
    //ShieldRadiusCmd->SetDefaultValue(30*cm);
    ShieldRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    ShieldTypeCmd = new G4UIcmdWithAString("/usr/nShield/type",this);
    ShieldTypeCmd->SetGuidance("UI command for the shape of a neutron shield");
    //ShieldTypeCmd->SetCandidates("sphere hemisphere");
    ShieldTypeCmd->SetParameterName("neutron shield type",false);
    ShieldTypeCmd->SetDefaultValue("sphere");
    ShieldTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}


mcDetectorMessenger::~mcDetectorMessenger()
{
    delete MaterialCmd;
    delete MagFieldCmd;
    delete MaxStepCmd;
    delete ShieldRadiusCmd;
    delete ShieldTypeCmd;

    delete usrDir;
    delete nShieldDir;
    
}


void mcDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
    if( command == MaterialCmd ){
        mcDetector->SetSensorMaterial(newValue);
    } else if( command == MaxStepCmd ){
        mcDetector->SetMaxStep(MaxStepCmd->GetNewDoubleValue(newValue));
    } else if( command == MagFieldCmd ){
        mcDetector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));
    } else if( command == ShieldRadiusCmd ){
        mcDetector->SetNeutronShieldRadius(ShieldRadiusCmd->GetNewDoubleValue(newValue));
    } else if( command == ShieldTypeCmd ){
        mcDetector->SetNeutronShieldType(newValue);
    }
}

G4String mcDetectorMessenger::GetCurrentValue(G4UIcommand * command)
{
    G4String cv;
    
    if( command==MaterialCmd ){
        cv =  mcDetector->GetSensorMaterial()->GetName();
        
    } else if( command==MaxStepCmd ){
        cv =  MaxStepCmd->ConvertToString( mcDetector->GetMaxStep(),"mm");
        
    } else if( command==MagFieldCmd ){
        cv =  MagFieldCmd->ConvertToString( mcDetector->GetFieldValue(),"tesla");
        
    }
    
    return cv;
}

