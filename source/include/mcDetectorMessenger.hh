
#ifndef mcDetectorMessenger_h
#define mcDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class mcDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;


class mcDetectorMessenger: public G4UImessenger
{
public:
    mcDetectorMessenger(mcDetectorConstruction* );
    ~mcDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    G4String GetCurrentValue(G4UIcommand * command);
    
private:
    mcDetectorConstruction*    mcDetector;
    G4UIdirectory*             usrDir;
    G4UIdirectory*             detDir;
    G4UIdirectory*             nShieldDir;
    G4UIdirectory*             gShield1Dir;
    G4UIdirectory*             gShield2Dir;
    G4UIcmdWithAString*        MaterialCmd;
    G4UIcmdWithADoubleAndUnit* MaxStepCmd;
    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithADoubleAndUnit* NeutronShieldSizeCmd;
    G4UIcmdWithAString*        NeutronShieldTypeCmd;
    G4UIcmdWithAString*        NeutronShieldMaterialCmd;
    G4UIcmdWithABool *         GammaShield1Cmd;
    G4UIcmdWithADoubleAndUnit* GammaShield1ThicknessCmd;
    G4UIcmdWithABool *         GammaShield2Cmd;
    
};


#endif

