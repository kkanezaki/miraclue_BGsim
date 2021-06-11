/*

// DMXDetectorConstruction.cc
// 20201124 Kiseki Nakamura

#include "DMXDetectorConstruction.hh"
#include "DMXTPCSD.hh"
#include "DMXLSSD.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"


#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Ellipsoid.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4GenericPolycone.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ExtrudedSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"



DMXDetectorConstruction::DMXDetectorConstruction(G4String filename){
  config_filename = filename;

  //limit
  theUserLimitsForRoom     = 0;
  theUserLimitsForDetector = 0;
  theTimeCutForRoom      = 1000.*nanosecond;
  theTimeCutForDetector  = DBL_MAX;
  theStepSizeForRoom     = DBL_MAX;
  // theStepSizeForDetector = 0.1*mm;
  theMinEkineForRoom     = 250.0*eV;
  theMinEkineForDetector = 250.0*eV;

  theDetectorName = "none";
  theDetectorParamNum = 0;
  for(int i=0; i<10; i++) theDetectorParam[i] = 0.0;
  theDetectorZoffset = 0;
  theDetectorZrange = 0;
  theGasPressure = 0;
  theGasTypeNum = 0;
  theGasName = "gas";
  theStepSizeForDetector = 0.1*mm;
  std::ifstream file_param(filename,std::ios::in);
  G4cout<<" ### parameter read ###"<<G4endl;
  if(file_param.fail()==1){ G4cout<<" ### detector_param.dat not found !"<<G4endl;  exit(1); }
  else{
    G4String dummy;
    while( file_param >> dummy ){
      if(dummy=="detector")  file_param >> theDetectorName;
      if(dummy=="parameter"){
	file_param >> theDetectorParamNum;
	for(int i=0; i<theDetectorParamNum; i++) file_param >> theDetectorParam[i];
      }
      if(dummy=="pressure")  file_param >> theGasPressure;
      if(dummy=="gastype"){
	file_param >> theGasTypeNum;
	double totalRatio=0;
	for(int i=0; i<theGasTypeNum; i++){
	  file_param >> theGasType[i] >> theGasRatio[i];
	  totalRatio += theGasRatio[i];
	}
	for(int i=0; i<theGasTypeNum; i++){
	  theGasRatio[i] /= totalRatio;
	  theGasName = theGasName + "_" + theGasType[i] + "_" + std::to_string(theGasRatio[i]);
	}
      }
      if(dummy=="stepsize") file_param >> theStepSizeForDetector;
    }
  }
  if(theDetectorName=="none" || theGasPressure==0 || theGasTypeNum==0 || theStepSizeForDetector==0){
    G4cout<<" ### description in detector_param.dat is not sufficient !"<<G4endl;
    G4cout<<" theDetectorName="<<theDetectorName<<" theGasPressure="<<theGasPressure
	  <<" theGasTypeNum="<<theGasTypeNum<<" theStepSizeForDetector="<<theStepSizeForDetector<<G4endl;
    G4cout<<" theGasName="<<theGasName<<G4endl;
    exit(1);
  }

}




DMXDetectorConstruction::~DMXDetectorConstruction(){
}

G4VPhysicalVolume* DMXDetectorConstruction::Construct() {

  // read config file
  G4String detector_name = "none";
  int detector_param_num = 0;
  double detector_param[10] = {};
  double gas_pressure = 0;
  int gas_type_num = 0;
  G4String gas_name = "gas";
  G4String gas_type[10];
  double gas_ratio[10];
  double step_size = 0.1*mm;
  std::ifstream file_param(config_filename, std::ios::in);
  G4cout<<" ### parameter read ###"<<G4endl;
  if(file_param.fail()==1){ G4cout<<" ### detector_param.dat not found !"<<G4endl;  exit(1); }
  else{
    G4String dummy;
    while( file_param >> dummy ){
      if(dummy=="detector")  file_param >> detector_name;
      if(dummy=="parameter"){
	file_param >> detector_param_num;
	for(int i=0; i<detector_param_num; i++) file_param >> detector_param[i];
      }
      if(dummy=="pressure")  file_param >> gas_pressure;
      if(dummy=="gastype"){
	file_param >> gas_type_num;
	double total_ratio=0;
	for(int i=0; i<gas_type_num; i++){
	  file_param >> gas_type[i] >> gas_ratio[i];
	  total_ratio += gas_ratio[i];
	}
	for(int i=0; i<gas_type_num; i++){
	  gas_ratio[i] /= total_ratio;
	  gas_name = gas_name + "_" + gas_type[i] + "_" + std::to_string(gas_ratio[i]);
	}
      }
      if(dummy=="stepsize") file_param >> step_size;
    }
  }
  if(detector_name=="none" || gas_pressure<=0 || gas_type_num==0 || step_size==0){
    G4cout<<" ### description in detector_param.dat is not sufficient !"<<G4endl;
    G4cout<<" detector_name="<<detector_name<<" gas_pressure="<<gas_pressure
	  <<" gas_type_num="<<gas_type_num<<" step_size="<<step_size<<G4endl;
    G4cout<<" gas_name="<<gas_name<<G4endl;
    exit(1);
  }
  
  // make colours
  double alpha = 1.0;
  G4Colour white (1.0, 1.0, 1.0, alpha);  G4Colour grey  (0.5, 0.5, 0.5, alpha);
  G4Colour red   (1.0, 0.0, 0.0, alpha);  G4Colour blue  (0.0, 0.5, 1.0, alpha);
  G4Colour yellow(1.0, 1.0, 0.0, alpha);  G4Colour orange(.75, .55, 0.0, alpha);
  G4Colour green (0.0, 1.0, 0.0, alpha);
  
  // element definition
  G4Element* elH  = new G4Element("H",  "H",  1.,  1.007947*g/mole);
  G4Element* elHe = new G4Element("He", "He", 2.,  4.002602*g/mole);
  G4Element* elLi = new G4Element("Li", "Li", 3.,  6.9412*g/mole);
  G4Element* elBe = new G4Element("Be", "Be", 4.,  9.0121823*g/mole);
  G4Element* elB  = new G4Element("B",  "B",  5.,  10.8117*g/mole);
  G4Element* elC  = new G4Element("C",  "C",  6.,  12.0107*g/mole);
  G4Element* elN  = new G4Element("N",  "N",  7.,  14.00672*g/mole);
  G4Element* elO  = new G4Element("O",  "O",  8.,  15.99943*g/mole);
  G4Element* elF  = new G4Element("F",  "F",  9.,  18.9984031636*g/mole);
  G4Element* elNa = new G4Element("Na", "Na", 11., 22.989769282*g/mole);
  G4Element* elMg = new G4Element("Mg", "Mg", 12., 24.30506*g/mole);
  G4Element* elAl = new G4Element("Al", "Al", 13., 26.981538613*g/mole);
  G4Element* elSi = new G4Element("Si", "Si", 14., 28.08553*g/mole);
  G4Element* elP  = new G4Element("P",  "P",  15., 30.9737622*g/mole);
  G4Element* elS  = new G4Element("S",  "S",  16., 32.0655*g/mole);
  G4Element* elCl = new G4Element("Cl", "Cl", 17., 35.4532*g/mole);
  G4Element* elAr = new G4Element("Ar", "Ar", 18., 39.948*g/mole);
  G4Element* elK  = new G4Element("K",  "K",  19., 39.0983*g/mole);
  G4Element* elCa = new G4Element("Ca", "Ca", 20., 40.078*g/mole);
  G4Element* elSc = new G4Element("Sc", "Sc", 21., 44.9559126*g/mole);
  G4Element* elTi = new G4Element("Ti", "Ti", 22., 47.8671*g/mole);
  G4Element* elCr = new G4Element("Cr", "Cr", 24., 51.99616*g/mole);
  G4Element* elMn = new G4Element("Mn", "Mn", 25., 54.9380455*g/mole);
  G4Element* elFe = new G4Element("Fe", "Fe", 26., 55.8452*g/mole);
  G4Element* elCo = new G4Element("Co", "Co", 27., 58.9331955*g/mole);
  G4Element* elNi = new G4Element("Ni", "Ni", 28., 58.693442*g/mole);
  G4Element* elCu = new G4Element("Cu", "Cu", 29., 63.5463*g/mole);
  G4Element* elZn = new G4Element("Zn", "Zn", 30., 65.3824*g/mole);
  G4Element* elGa = new G4Element("Ga", "Ga", 31., 69.7231*g/mole);
  G4Element* elGe = new G4Element("Ge", "Ge", 32., 72.631*g/mole);
  G4Element* elBr = new G4Element("Br", "Br", 35., 79.904*g/mole);
  G4Element* elZr = new G4Element("Zr", "Zr", 40., 91.224*g/mole);
  G4Element* elMo = new G4Element("Mo", "Mo", 42., 95.961*g/mole);
  G4Element* elRu = new G4Element("Ru", "Ru", 44., 101.07*g/mole);
  G4Element* elAg = new G4Element("Ag", "Ag", 47., 107.8682*g/mole);
  G4Element* elSn = new G4Element("Sn", "Sn", 50., 118.710*g/mole);
  G4Element* elTe = new G4Element("Te", "Te", 52., 127.60*g/mole);
  G4Element* elI  = new G4Element("I",  "I",  53., 126.90447*g/mole);
  G4Element* elXe = new G4Element("Xe", "Xe", 54., 131.293*g/mole);
  G4Element* elCs = new G4Element("Cs", "Cs", 55., 132.90545192*g/mole);
  G4Element* elBa = new G4Element("Ba", "Ba", 56., 137.33*g/mole);
  G4Element* elLa = new G4Element("La", "La", 57., 138.90547*g/mole);
  G4Element* elGd = new G4Element("Gd", "Gd", 64., 157.25*g/mole);
  G4Element* elW  = new G4Element("W",  "W",  74., 183.84*g/mole);
  G4Element* elPt = new G4Element("Pt", "Pt", 78., 195.084*g/mole);
  G4Element* elAu = new G4Element("Au", "Au", 79., 196.9665694*g/mole);
  G4Element* elHg = new G4Element("Hg", "Hg", 80., 200.592*g/mole);
  G4Element* elPb = new G4Element("Pb", "Pb", 82., 207.2*g/mole);

  // make gas
  G4double density_Xe_1atm  = 5.8971 *mg/cm3;//noble gas detector(273K,1atm)
  G4double density_He_1atm  = 0.17850*mg/cm3;//noble gas detector(273K,1atm)
  G4double density_Ar_1atm  = 1.7606 *mg/cm3;//noble gas detector(273K,1atm)
  G4double density_CF4_1atm = 3.76   *mg/cm3;//Nishimura D-ron, Nakamura D-ron (15C,1atm)
  G4double density_1atm = 0.0;
  int gas_element_num = 0;
  for(int i=0; i<gas_type_num; i++){
    if(     gas_type[i]=="Xe"){  density_1atm+=density_Xe_1atm *gas_ratio[i];  gas_element_num+=1; }
    else if(gas_type[i]=="He"){  density_1atm+=density_He_1atm *gas_ratio[i];  gas_element_num+=1; }
    else if(gas_type[i]=="Ar"){  density_1atm+=density_Ar_1atm *gas_ratio[i];  gas_element_num+=1; }
    else if(gas_type[i]=="CF4"){ density_1atm+=density_CF4_1atm*gas_ratio[i];  gas_element_num+=2; }
    else{ G4cout<<"ERROR unknown gas_type[i] : "<<gas_type[i]<<G4endl; exit(1); }
  }
  G4cout<<" density="<<density_1atm<<"   "<<density_1atm/mg*cm3<<G4endl;
  G4Material* gas_mat = new G4Material("gas_mat", density_1atm*gas_pressure, gas_element_num,
				       kStateGas, 300.00*kelvin, gas_pressure*atmosphere);
  for(int i=0; i<gas_type_num; i++){
    if(gas_type[i]=="Xe") gas_mat->AddElement(elXe, gas_ratio[i]*density_Xe_1atm/density_1atm);
    if(gas_type[i]=="He") gas_mat->AddElement(elHe, gas_ratio[i]*density_He_1atm/density_1atm);
    if(gas_type[i]=="Ar") gas_mat->AddElement(elAr, gas_ratio[i]*density_Ar_1atm/density_1atm);
    if(gas_type[i]=="CF4"){
      gas_mat->AddElement(elC, gas_ratio[i]*density_CF4_1atm/density_1atm *12./88.);
      gas_mat->AddElement(elF, gas_ratio[i]*density_CF4_1atm/density_1atm *19.*4./88.);
    }
    G4cout<<" i="<<i<<" "<<gas_type[i]<<" "<<gas_ratio[i]<<G4endl;
  }

  // user limits   G4userLimits(step-length-max, track-length-max, time-cut, min-energy)
  G4UserLimits *user_limit = new G4UserLimits(step_size, DBL_MAX, DBL_MAX, 0, 0);

  // world
  G4Material* world_mat = new G4Material("world", 1.29*mg/cm3, 2); // air
  world_mat->AddElement(elN, 0.7);  world_mat->AddElement(elO, 0.3);
  G4double world_l = 100.0*m;
  G4Box* world_box = new G4Box("world_box", world_l, world_l, world_l);
  G4LogicalVolume* world_lv = new G4LogicalVolume( world_box, world_mat, "world", 0, 0, 0);
  G4VisAttributes* world_att = new G4VisAttributes( false, white );
  world_lv-> SetVisAttributes(world_att);
  G4PVPlacement* world_pv = new G4PVPlacement(0, G4ThreeVector(), "world_pv", world_lv, 0, false,0);

  // geometry

  if(      detector_name=="BOX_1.0.0"){
#include "GeometryBOX_1.0.0.icc"
  }else if(detector_name=="MigdalTestChamber_1.0.0"){
#include "GeometryMigdalTestChamber_1.0.0.icc"
  }else if(detector_name=="MigdalTestLab_1.0.3"){
#include "GeometryMigdalTestLab_1.0.3.icc"
  }else{
    G4cout<<"error detector="<<detector_name<<" is not defined"<<G4endl;
    exit(1);
  }

  // sensitive detector
  DMXTPCSD* tpcSD = new DMXTPCSD("/DMXDet/TPCSD","tpcCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector( tpcSD );
  SD_lv->SetSensitiveDetector( tpcSD );
  
  G4cout<<" Construct done"<<G4endl;
  return world_pv;
}
*/