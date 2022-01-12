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
    nShieldDir->SetGuidance("UI commands for neutron shield");
    nShield2Dir = new G4UIdirectory("/usr/nShield2/");
    nShield2Dir->SetGuidance("UI commands for layer neutron shield");
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

    GasBoxSizeCmd = new G4UIcmdWithADoubleAndUnit("/usr/det/setGasSize",this);
    GasBoxSizeCmd->SetGuidance("Set the size of gas box");
    GasBoxSizeCmd->SetParameterName("gas box size",false);
    GasBoxSizeCmd->SetUnitCategory("Length");
    GasBoxSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    NeutronShieldSizeCmd = new G4UIcmdWithADoubleAndUnit("/usr/nShield/size",this);
    NeutronShieldSizeCmd->SetGuidance("neutron shield size");
    NeutronShieldSizeCmd->SetParameterName("neutron shield size",false);
    NeutronShieldSizeCmd->SetUnitCategory("Length");
    NeutronShieldSizeCmd->SetDefaultValue(25*cm);
    NeutronShieldSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    NeutronShieldTypeCmd = new G4UIcmdWithAString("/usr/nShield/type",this);
    NeutronShieldTypeCmd->SetGuidance("UI command for the shape of a neutron shield");
    NeutronShieldTypeCmd->SetCandidates("sphere hemisphere cube layer none");
    NeutronShieldTypeCmd->SetParameterName("neutron shield type",false);
    NeutronShieldTypeCmd->SetDefaultValue("none");
    NeutronShieldTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    InnerShieldSizeCmd = new G4UIcmdWithADoubleAndUnit("/usr/nShield2/thickness",this);
    InnerShieldSizeCmd->SetGuidance("UI command for the length of inner neutron shield");
    InnerShieldSizeCmd->SetParameterName("size of inner shield",false);
    InnerShieldSizeCmd->SetUnitCategory("Length");
    InnerShieldSizeCmd->SetDefaultValue(25*cm);
    InnerShieldSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    OuterShieldThicknessCmd = new G4UIcmdWithADoubleAndUnit("/usr/nShield2/thickness",this);
    OuterShieldThicknessCmd->SetGuidance("UI command for the thickness of outer neutron shield");
    OuterShieldThicknessCmd->SetParameterName("thickness of outer shield",false);
    OuterShieldThicknessCmd->SetUnitCategory("Length");
    OuterShieldThicknessCmd->SetDefaultValue(1*cm);
    OuterShieldThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    NeutronShieldMaterialCmd = new G4UIcmdWithAString("/usr/nShield/setMaterial",this);
    NeutronShieldMaterialCmd->SetGuidance("UI command for the material of a neutron shield");
    NeutronShieldMaterialCmd->SetParameterName("neutron shield material",false);
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

    GammaShield2ThicknessCmd = new G4UIcmdWithADoubleAndUnit("/usr/gShield2/thickness",this);
    GammaShield2ThicknessCmd->SetGuidance("UI command for the thickness of gamma shield 1");
    GammaShield2ThicknessCmd->SetParameterName("thickness of gamma shield 2",false);
    GammaShield2ThicknessCmd->SetUnitCategory("Length");
    GammaShield2ThicknessCmd->SetDefaultValue(5*cm);
    GammaShield2ThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    ScintiMaterialCmd = new G4UIcmdWithAString("/usr/scinti/setMaterial",this);
    ScintiMaterialCmd->SetGuidance("UI command for the scintillator material");
    ScintiMaterialCmd->SetCandidates("BGO GAGG CsI NaI BC501A");
    ScintiMaterialCmd->SetParameterName("scintillator material", false);
    ScintiMaterialCmd->SetDefaultValue("BGO");
    ScintiMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}


mcDetectorMessenger::~mcDetectorMessenger()
{
    delete MaterialCmd;
    delete MagFieldCmd;
    delete MaxStepCmd;
    delete GasBoxSizeCmd;
    delete NeutronShieldSizeCmd;
    delete NeutronShieldTypeCmd;
    delete NeutronShieldMaterialCmd;
    delete GammaShield1Cmd;
    delete GammaShield1ThicknessCmd;
    delete GammaShield2ThicknessCmd;

    delete usrDir;
    delete nShieldDir;
    delete gShield1Dir;
    delete gShield2Dir;
    
}


void mcDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue) {
    if (command == MaterialCmd) {
        mcDetector->SetSensorMaterial(newValue);
    } else if (command == MaxStepCmd) {
        mcDetector->SetMaxStep(MaxStepCmd->GetNewDoubleValue(newValue));
    } else if (command == MagFieldCmd) {
        mcDetector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));
    } else if (command == GasBoxSizeCmd) {
        mcDetector->SetGasBoxSize(GasBoxSizeCmd->GetNewDoubleValue(newValue));
    } else if (command == NeutronShieldMaterialCmd) {
        mcDetector->SetNeutronShieldMaterial(newValue);
    } else if (command == NeutronShieldSizeCmd) {
        mcDetector->SetNeutronShieldSize(NeutronShieldSizeCmd->GetNewDoubleValue(newValue));
    } else if (command == NeutronShieldTypeCmd) {
        mcDetector->SetNeutronShieldType(newValue);
    } else if (command == GammaShield1ThicknessCmd) {
        mcDetector->SetGammaShield1Thickness(GammaShield1ThicknessCmd->GetNewDoubleValue(newValue));
    } else if (command == GammaShield2ThicknessCmd) {
        mcDetector->SetGammaShield2Thickness(GammaShield2ThicknessCmd->GetNewDoubleValue(newValue));
    } else if (command == ScintiMaterialCmd) {
        mcDetector->SetScintiMaterial(newValue);
    }
}

G4String mcDetectorMessenger::GetCurrentValue(G4UIcommand *command) {
    G4String cv;

    if (command == MaterialCmd) {
        cv = mcDetector->GetSensorMaterial()->GetName();

    } else if (command == MaxStepCmd) {
        cv = MaxStepCmd->ConvertToString(mcDetector->GetMaxStep(), "mm");

    } else if (command == MagFieldCmd) {
        cv = MagFieldCmd->ConvertToString(mcDetector->GetFieldValue(), "tesla");

    }

    return cv;
}

