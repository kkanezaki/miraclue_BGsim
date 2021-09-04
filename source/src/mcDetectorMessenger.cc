#include "mcDetectorMessenger.hh"
#include "mcDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithoutParameter.hh"


mcDetectorMessenger::mcDetectorMessenger(mcDetectorConstruction* mcDet)
:mcDetector(mcDet)
{ 
    usrDir = new G4UIdirectory("/usr/");
    usrDir->SetGuidance("UI commands of this example");
    detDir = new G4UIdirectory("/usr/det/");
    detDir->SetGuidance("UI commands for detector setup");
    nShieldDir = new G4UIdirectory("/usr/nShield/");
    nShieldDir->SetGuidance("UI commands for gamma shield");
    gShield1Dir = new G4UIdirectory("/usr/gShield1/");
    gShield1Dir->SetGuidance("UI commands for gamma ray shield on the neutron shield");
    gShield2Dir = new G4UIdirectory("/usr/gShield2/");
    gShield2Dir->SetGuidance("UI commands for gamma ray shield on the chamber");

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

    NeutronShieldSizeCmd = new G4UIcmdWithADoubleAndUnit("/usr/nShield/size",this);
    NeutronShieldSizeCmd->SetGuidance("neutron shield size");
    NeutronShieldSizeCmd->SetParameterName("neutron shield size",false);
    NeutronShieldSizeCmd->SetUnitCategory("Length");
    NeutronShieldSizeCmd->SetDefaultValue(25*cm);
    NeutronShieldSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    NeutronShieldTypeCmd = new G4UIcmdWithAString("/usr/nShield/type",this);
    NeutronShieldTypeCmd->SetGuidance("UI command for the shape of a neutron shield");
    NeutronShieldTypeCmd->SetCandidates("sphere hemisphere cube none");
    NeutronShieldTypeCmd->SetParameterName("neutron shield type",false);
    NeutronShieldTypeCmd->SetDefaultValue("none");
    NeutronShieldTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    NeutronShieldMaterialCmd = new G4UIcmdWithAString("/usr/nShield/setMaterial",this);
    NeutronShieldMaterialCmd->SetGuidance("UI command for the material of a neutron shield");
    NeutronShieldMaterialCmd->SetCandidates("polyethylene_boron10");
    NeutronShieldMaterialCmd->SetParameterName("neutron shield type",false);
    NeutronShieldMaterialCmd->SetDefaultValue("polyethylene_boron10");
    NeutronShieldMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    GammaShield1Cmd = new G4UIcmdWithABool("/usr/gShield1/putShield",this);
    GammaShield1Cmd->SetGuidance("UI command for putting gShield1 ");
    GammaShield1Cmd->SetParameterName("gamma shield type",false);
    GammaShield1Cmd->SetDefaultValue(false);
    GammaShield1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    GammaShield1ThicknessCmd = new G4UIcmdWithADoubleAndUnit("/usr/gShield1/thickness",this);
    GammaShield1ThicknessCmd->SetGuidance("UI command for the thickness of gamma shield 1");
    GammaShield1ThicknessCmd->SetParameterName("thickness of gamma shield 1",false);
    GammaShield1ThicknessCmd->SetUnitCategory("Length");
    GammaShield1ThicknessCmd->SetDefaultValue(5*cm);
    GammaShield1ThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}


mcDetectorMessenger::~mcDetectorMessenger()
{
    delete MaterialCmd;
    delete MagFieldCmd;
    delete MaxStepCmd;
    delete NeutronShieldSizeCmd;
    delete NeutronShieldTypeCmd;
    delete NeutronShieldMaterialCmd;
    delete GammaShield1Cmd;
    delete GammaShield1ThicknessCmd;

    delete usrDir;
    delete nShieldDir;
    delete gShield1Dir;
    delete gShield2Dir;
    
}


void mcDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
    if( command == MaterialCmd ){
        mcDetector->SetSensorMaterial(newValue);
    } else if( command == MaxStepCmd ){
        mcDetector->SetMaxStep(MaxStepCmd->GetNewDoubleValue(newValue));
    } else if( command == MagFieldCmd ){
        mcDetector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));
    } else if( command == NeutronShieldSizeCmd ){
        mcDetector->SetNeutronShieldSize(NeutronShieldSizeCmd->GetNewDoubleValue(newValue));
    } else if( command == NeutronShieldTypeCmd ){
        mcDetector->SetNeutronShieldType(newValue);
    } else if( command == GammaShield1Cmd){
        //mcDetector->;
    } else if( command == GammaShield1ThicknessCmd ){
        mcDetector->SetGammaShield1Thickness(GammaShield1ThicknessCmd->GetNewDoubleValue(newValue));
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

