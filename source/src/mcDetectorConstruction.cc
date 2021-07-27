#include "mcDetectorConstruction.hh"
#include "mcDetectorMessenger.hh"
#include "mcSensorSD.hh"
#include "mcAnalyzer.hh"

#include "DMXDetectorConstruction.hh"
#include "DMXDetectorMaterial.ihh"
#include "DMXTPCSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>

mcDetectorConstruction::mcDetectorConstruction()
        :defaultMaterial(0),sensorMaterial(0),
         WorldRadius(100*cm),
         solidWorld(0),logicWorld(0),physWorld(0),
         solidSensor(0),logicSensor(0),physSensor(0),
         magField(0),pUserLimits(0),maxStep(100.0*cm)
{

    // default parameter values of Sensor
    DefineMaterials();

    //SetSensorMaterial("NaI");
    SetSensorMaterial("Ar1atm");


    // create commands for interactive definition of the calorimeter
    detectorMessenger = new mcDetectorMessenger(this);
}

mcDetectorConstruction::~mcDetectorConstruction()
{
    delete detectorMessenger;
}

//------------------------------------------------------------------------//
// Begin of Construct()
G4VPhysicalVolume* mcDetectorConstruction::Construct()
{
    // Clean old geometry, if any
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    // make colours
    double alpha = 1.0;
    G4Colour white (1.0, 1.0, 1.0, alpha);  G4Colour grey  (0.5, 0.5, 0.5, alpha);
    G4Colour red   (1.0, 0.0, 0.0, alpha);  G4Colour blue  (0.0, 0.5, 1.0, alpha);
    G4Colour yellow(1.0, 1.0, 0.0, alpha);  G4Colour orange(.75, .55, 0.0, alpha);
    G4Colour green (0.0, 1.0, 0.0, alpha);

    // create geometry

    // World
    /*
    solidWorld = new G4Orb("World",WorldRadius);
    logicWorld = new G4LogicalVolume(solidWorld,     //its solid
                                     defaultMaterial,//its material
                                     "World");		 //its name

    physWorld = new G4PVPlacement(0,			     //no rotation
                                  G4ThreeVector(),	 //at (0,0,0)
                                  logicWorld,		 //its logical volume
                                  "World",		     //its name
                                  0,			     //its mother  volume
                                  false,			 //no boolean operation
                                  0);			     //copy number
    */

    //element
    G4Element* elH  = new G4Element("H",  "H",  1.,  1.007947*g/mole);
    G4Element* elB  = new G4Element("B",  "B",  5.,  10.8117*g/mole);
    G4Element* elC  = new G4Element("C",  "C",  6.,  12.0107*g/mole);
    G4Element* elN  = new G4Element("N",  "N",  7.,  14.00672*g/mole);
    G4Element* elO  = new G4Element("O",  "O",  8.,  15.99943*g/mole);
    G4Element* elAl = new G4Element("Al", "Al", 13., 26.981538613*g/mole);
    G4Element* elSi = new G4Element("Si", "Si", 14., 28.08553*g/mole);
    G4Element* elAr = new G4Element("Ar", "Ar", 18., 39.948*g/mole);
    G4Element* elCa = new G4Element("Ca", "Ca", 20., 40.078*g/mole);
    G4Element* elCr = new G4Element("Cr", "Cr", 24., 51.9961*g/mole);
    G4Element* elFe = new G4Element("Fe", "Fe", 26., 55.8452*g/mole);
    G4Element* elNi = new G4Element("Ni", "Ni", 28., 58.6934*g/mole);
    G4Element* elXe = new G4Element("Xe", "Xe", 54., 131.293*g/mole);
    G4Element* elPb = new G4Element("Pb", "Pb", 82., 207.2*g/mole);

    G4Material* vac_mat = new G4Material("vacuum", 1., 1.*g/mole, 1.e-20*g/cm3, kStateGas, 0.1*kelvin, 1.e-20*bar);
    world_mat = new G4Material("world", 1.29*mg/cm3, 3); // air
    world_mat->AddElement(elN, 0.76);  world_mat->AddElement(elO, 0.23); world_mat->AddElement(elAr, 0.01);

    G4double world_l = 100.0*m;
    G4Box* world_box = new G4Box("world_box", world_l, world_l, world_l);
    G4LogicalVolume* world_lv = new G4LogicalVolume( world_box, vac_mat, "world", 0, 0, 0);
    G4VisAttributes* world_att = new G4VisAttributes( false, white );
    world_lv-> SetVisAttributes(world_att);
    G4PVPlacement* world_pv = new G4PVPlacement(0, G4ThreeVector(), "world_pv", world_lv, 0, false,0);


    // geometry of chamber
    double step_size = 0.1*mm;

    // user limits   G4userLimits(step-length-max, track-length-max, time-cut, min-energy)
    G4UserLimits *user_limit = new G4UserLimits(step_size, DBL_MAX, DBL_MAX, 0, 0);

    G4double density_Ar_1atm = 1.7606 *mg/cm3; //noble gas detector(273K,1atm)
    G4double density_Xe_1atm = 5.8971 *mg/cm3; //noble gas detector(273K,1atm)

    G4Material* Ar1atm = new G4Material("Ar1atm", density_Ar_1atm, 1, kStateGas, 300.00*kelvin, 1*atmosphere);
    Ar1atm->AddElement(elAr, 1);
    G4Material* Xe8atm = new G4Material("Xe8atm", density_Xe_1atm*8, 1, kStateGas, 300.00*kelvin, 1*atmosphere);
    Xe8atm->AddElement(elXe, 1);

    G4double bulk_l = 4*mm;//100000*mm;
    G4double box_x = 1000*mm;
    G4double box_y = 300*mm;
    G4double box_z = 300*mm;
    G4double buffer_t = 300*mm;
    G4double chamber_t = 300*mm;

    G4double test_l = 300*mm;

    //////////////////////////////// devided gas box (for Intrinsic BG) ////////////////////////////////
    /*
    G4int Ndiv_x = 15;
    G4int Ndiv_y = 15;
    G4int Ndiv_z = 15;
    G4int Ncopy = Ndiv_x * Ndiv_y * Ndiv_z;

    G4Box* test_box = new G4Box("test", 0.5*test_l, 0.5*test_l, 0.5*test_l);
    G4LogicalVolume* test_lv = new G4LogicalVolume(test_box, Ar1atm, "test_lv");
    G4VisAttributes* test_att = new G4VisAttributes(1, blue);  test_lv->SetVisAttributes(test_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "test_pv", test_lv, world_pv, false, Ncopy);
    test_lv->SetUserLimits(user_limit);

    G4Box* rep1 = new G4Box("rep1", test_l*0.5, test_l*0.5, (test_l/Ndiv_z)*0.5);
    G4LogicalVolume* rep1_lv = new G4LogicalVolume(rep1, Ar1atm, "rep1_lv");
    G4VisAttributes* rep1_att = new G4VisAttributes(1, blue);  rep1_lv->SetVisAttributes(rep1_att);
    new G4PVReplica("rep1_pv", rep1_lv, test_lv, kZAxis, Ndiv_z, test_l/Ndiv_z);

    G4Box* rep2 = new G4Box("rep2", test_l*0.5, (test_l/Ndiv_y)*0.5, (test_l/Ndiv_z)*0.5);
    G4LogicalVolume* rep2_lv = new G4LogicalVolume(rep2, Ar1atm, "rep2_lv");
    G4VisAttributes* rep2_att = new G4VisAttributes(1, blue);  rep2_lv->SetVisAttributes(rep2_att);
    new G4PVReplica("rep2_pv", rep2_lv, rep1_lv, kYAxis, Ndiv_y, test_l/Ndiv_y);

    G4Box* rep3 = new G4Box("rep3", (test_l/Ndiv_x)*0.5, (test_l/Ndiv_y)*0.5, (test_l/Ndiv_z)*0.5);
    G4LogicalVolume* rep3_lv = new G4LogicalVolume(rep3, Ar1atm, "rep3_lv");
    G4VisAttributes* rep3_att = new G4VisAttributes(1, blue);  rep3_lv->SetVisAttributes(rep3_att);
    new G4PVReplica("rep3_pv", rep3_lv, rep2_lv, kXAxis, Ndiv_x, test_l/Ndiv_x);

    */
    /////////////////////////////////////////  for chamber BG //////////////////////////////////////////

    G4Material* SUS304 = new G4Material("SUS304", 8.03*g/cm3, 3);
    SUS304->AddElement(elFe, 0.7); SUS304->AddElement(elNi, 0.1); SUS304->AddElement(elCr, 0.2);
    G4Material* metalAl_mat = new G4Material("metalAl", 2.700*g/cm3, 1);
    metalAl_mat->AddElement(elAl, 1);

    G4double chamber_width = 3*mm; //SUS 3mm, Al 5mm, pula 1cm

    G4Box* chamber_box = new G4Box("chamber", 0.5*(chamber_t+chamber_width*2)*mm, 0.5*(chamber_t+chamber_width*2)*mm, 0.5*(chamber_t+chamber_width*2)*mm);
    G4LogicalVolume* chamber_lv = new G4LogicalVolume(chamber_box, SUS304, "chamber_lv");
    G4VisAttributes* chamber_att = new G4VisAttributes(1, blue);  chamber_lv->SetVisAttributes(chamber_att);
    //G4PVPlacement* chamber_pv = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "chamber_pv", chamber_lv, world_pv, false, 0);
    chamber_lv->SetUserLimits(user_limit);


    G4Box* vac_box = new G4Box("vac_box", 0.5*chamber_t*mm, 0.5*chamber_t, 0.5*chamber_t);
    G4LogicalVolume* vac_lv = new G4LogicalVolume(vac_box, vac_mat, "vac_lv");
    G4VisAttributes* vac_att = new G4VisAttributes(1, red);  vac_lv->SetVisAttributes(vac_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), vac_lv, "vac_pv", chamber_lv, false, 0);
    vac_lv->SetUserLimits(user_limit);

    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "chamber_pv", chamber_lv, world_pv, false, 0);


    ///////////////////////////////////////////// sandbox ///////////////////////////////////////////////
    /*
    G4Material* Pb_mat = new G4Material("Pb", 11.2*g/cm3, 1);
    Pb_mat->AddElement(elPb,1);
    G4Material* Polyeth = new G4Material("Polyeth",0.95*g/cm3,2);
    Polyeth->AddElement(elH, 0.66); Polyeth->AddElement(elC,0.34);
    G4Material* Polyeth_B2O3 = new G4Material("Polyeth_B2O3", 1.04*g/cm3, 4);
    Polyeth_B2O3->AddElement(elH,0.59); Polyeth_B2O3->AddElement(elC,0.33);
    Polyeth_B2O3->AddElement(elB,0.03); Polyeth_B2O3->AddElement(elO,0.05);

    G4double l_mat_test = 1000*mm;

    G4Box* mat_test_box = new G4Box("mat_test", 0.5*l_mat_test*mm, 0.5*l_mat_test*mm, 0.5*l_mat_test*mm);
    G4LogicalVolume* mat_test_lv = new G4LogicalVolume(mat_test_box, SUS304, "mat_test_lv");
    G4VisAttributes* mat_test_att = new G4VisAttributes(1, blue);  mat_test_lv->SetVisAttributes(mat_test_att);
    G4PVPlacement* mat_test_pv = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
                                                   "mat_test_pv", mat_test_lv, world_pv, false, 0);
    mat_test_lv->SetUserLimits(user_limit);
    */
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /*
    G4Box* sensor_box = new G4Box("sonsor", 10*mm, 0.5*test_l, 0.5*test_l);
    G4LogicalVolume* sensor_lv = new G4LogicalVolume(sensor_box, world_mat, "sensor_lv");
    G4VisAttributes* sensor_att = new G4VisAttributes(1, red);  sensor_lv->SetVisAttributes(sensor_att);
    new G4PVPlacement(0, G4ThreeVector(155*mm,0.,0.), "sensor_pv", sensor_lv, world_pv, false, 0);
    sensor_lv->SetUserLimits(user_limit);
    */

    /*
    G4Box* bulk_box = new G4Box("bulk_box", 0.5*bulk_l, 0.5*bulk_l, 0.5*bulk_l);
    G4LogicalVolume* bulk_lv = new G4LogicalVolume(bulk_box, vac_mat, "bulk_lv");
    G4VisAttributes* bulk_att = new G4VisAttributes(0, red);  bulk_lv->SetVisAttributes(bulk_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "bulk_pv", bulk_lv, world_pv, false, 0);
    bulk_lv->SetUserLimits(user_limit);

    G4Box* chamber_box = new G4Box("chamber_box", 0.5*box_x+buffer_t+chamber_t,
                                   0.5*box_y+buffer_t+chamber_t, 0.5*box_z+buffer_t+chamber_t);
    G4LogicalVolume* chamber_lv = new G4LogicalVolume(chamber_box, SUS304_mat, "chamber_lv");
    G4VisAttributes* chamber_att = new G4VisAttributes(1, green);
    chamber_lv->SetVisAttributes(chamber_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*box_z), chamber_lv, "chamber_pv", bulk_lv, false, 0);
    chamber_lv->SetUserLimits(user_limit);

    G4Box* gas_box = new G4Box("gas_box", 0.5*box_x+buffer_t, 0.5*box_y+buffer_t, 0.5*box_z+buffer_t);
    G4LogicalVolume* gas_lv = new G4LogicalVolume(gas_box, gas_mat, "gas_lv");
    G4VisAttributes* gas_att = new G4VisAttributes(1, blue);  gas_lv->SetVisAttributes(gas_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), gas_lv, "gas_pv", chamber_lv, false, 0);
    gas_lv->SetUserLimits(user_limit);

    G4Box* SD_box = new G4Box("SD_box", 0.5*box_x, 0.5*box_y, 0.5*box_z);
    G4LogicalVolume* SD_lv = new G4LogicalVolume(SD_box, gas_mat, "SD_lv");
    G4VisAttributes* SD_att = new G4VisAttributes(1, yellow);  SD_lv->SetVisAttributes(SD_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), SD_lv, "SD_pv", gas_lv, false, 0);
    SD_lv->SetUserLimits(user_limit);
    */

    ///////////////////////////////////////// geometry of Lab /////////////////////////////////////////
    /*
    //material
    G4Material* concrete_mat = new G4Material("concrete", 2.3*g/cm3, 6);
    concrete_mat->AddElement(elSi, 0.227915);  concrete_mat->AddElement(elO,  0.60541);
    concrete_mat->AddElement(elH,  0.09972 );  concrete_mat->AddElement(elCa, 0.04986);
    concrete_mat->AddElement(elAl, 0.014245);  concrete_mat->AddElement(elFe, 0.00285);
    G4Material* lab_mat = new G4Material("lab", 1.29*mg/cm3, 2); // air
    lab_mat->AddElement(elN, 0.7);  lab_mat->AddElement(elO, 0.3);
    G4Material* metalAl_mat = new G4Material("metalAl", 2.700*g/cm3, 1);
    metalAl_mat->AddElement(elAl, 1);

    // parameter
    G4double lab_w  = 11.5*m, lab_h = 11.5*m, concrete_t = 1.0*m, floor_h = 7.8*mm;
    G4double concrete_w = lab_w  + concrete_t, concrete_h = lab_h + concrete_t;


    // Concrete wall
    G4Box *concrete_box = new G4Box("concrete_box", 0.5*concrete_w, 0.5*concrete_h, 0.5*concrete_w);
    G4LogicalVolume* concrete_lv = new G4LogicalVolume(concrete_box, concrete_mat, "concrete_lv");
    G4VisAttributes* concrete_att = new G4VisAttributes(1, grey);
    concrete_lv->SetVisAttributes(concrete_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "concrete_pv", concrete_lv, world_pv, false, 0);
    concrete_lv->SetUserLimits(user_limit);

    // Lab air
    G4Box *lab_box = new G4Box("lab_box", 0.5*lab_w, 0.5*lab_h, 0.5*lab_w);
    G4LogicalVolume* lab_lv = new G4LogicalVolume(lab_box, lab_mat, "lab_lv");
    G4VisAttributes* lab_att = new G4VisAttributes(1, red); lab_lv->SetVisAttributes(lab_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), lab_lv, "lab_pv", concrete_lv, false, 0);
    lab_lv->SetUserLimits(user_limit);

    // Floor Al
    G4Box *floor_box = new G4Box("floor_box", 0.5*lab_w, 0.5*floor_h, 0.5*lab_w);
    G4LogicalVolume* floor_lv = new G4LogicalVolume(floor_box, metalAl_mat, "floor_lv");
    G4VisAttributes* floor_att = new G4VisAttributes(1, blue); floor_lv->SetVisAttributes(floor_att);
    new G4PVPlacement(0, G4ThreeVector(0.,-1.0*m,0.), floor_lv, "floor_pv", lab_lv, false, 0);
    floor_lv->SetUserLimits(user_limit);
    */
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    // Sensor
    /*
    solidSensor = new G4Tubs("Sensor",0.0*cm,2.54*cm,18.95*cm,0,CLHEP::twopi);
    logicSensor = new G4LogicalVolume(solidSensor,sensorMaterial,"Sensor");

    physSensor = new G4PVPlacement(0,G4ThreeVector(),logicSensor,"Sensor",logicWorld,false,1);
    physSensor = new G4PVPlacement(0,G4ThreeVector(20*cm,0,0),logicSensor,"Sensor",logicWorld,false,2);
    physSensor = new G4PVPlacement(0,G4ThreeVector(40*cm,0,0),logicSensor,"Sensor",logicWorld,false,3);
    */


    //------------------------------------------------
    // Sensitive detectors
    //------------------------------------------------


    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    mcSensorSD* aSensorSD = (mcSensorSD*)SDman->FindSensitiveDetector("mc/SensorSD");
    if ( aSensorSD == 0){
        aSensorSD = new mcSensorSD("mc/SensorSD");
        SDman->AddNewDetector( aSensorSD );
    }
    aSensorSD->SetAnalyzer(analyzer);

    //logicSensor->SetSensitiveDetector(aSensorSD);
    //test_lv->SetSensitiveDetector(aSensorSD);
    //rep3_lv->SetSensitiveDetector(aSensorSD);
    //mat_test_lv->SetSensitiveDetector(aSensorSD);
    vac_lv->SetSensitiveDetector(aSensorSD);

    /*
    // Set UserLimits
    G4double maxTrkLen = 10.0*WorldRadius;
    G4double maxTime   = 1000.0 * ns;
    pUserLimits = new G4UserLimits(maxStep, maxTrkLen, maxTime);
    logicWorld->SetUserLimits(pUserLimits);

    // Visualization attributes
    logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    logicSensor->SetVisAttributes(simpleBoxVisAtt);

    //return physWorld;
    */

    return world_pv;

}
// End of Construct()
//------------------------------------------------------------------------//


///////////////////////////////////////////////////////
void mcDetectorConstruction::DefineMaterials()
{
    //This function illustrates the possible ways to define materials

    G4String symbol, name;
    G4double a, z, density, temperature, pressure;
    // No use temperature & pressure in Geant4
    G4int iz, n, in;

    G4int ncomponents, natoms;
    G4double abundance, fractionmass;

    // a =mass of a mole;
    // z =mean number of protons;
    // iz=number of protons in an isotope;
    // n =number of nucleons in an isotope;


    //
    // define Elements
    //
    G4Element* H  = new G4Element("Hydrogen", "H",  1.,  1.00794   *g/mole);
    G4Element* He = new G4Element("Helium",   "He", 2.,  4.002602  *g/mole);
    G4Element* B  = new G4Element("Bolone",   "B",  5.,  10.81     *g/mole);
    G4Element* C  = new G4Element("Carbon",   "C",  6.,  12.011    *g/mole);
    G4Element* N  = new G4Element("Nitrogen", "N",  7.,  14.00674  *g/mole);
    G4Element* O  = new G4Element("Oxygen",   "O",  8.,  16.00     *g/mole);
    G4Element* F  = new G4Element("Fluorine", "F",  9.,  18.9984032*g/mole);
    G4Element* Na = new G4Element("Natrium",  "Na", 11., 22.99     *g/mole);
    G4Element* Al = new G4Element("Aluminium","Al", 13., 26.98     *g/mole);
    G4Element* Si = new G4Element("Silicon",  "Si", 14., 28.09     *g/mole);
    G4Element* Ar = new G4Element("Argon",    "Ar", 18., 39.948    *g/mole);
    G4Element* Ca = new G4Element("Calcium",  "Ca", 20., 40.078    *g/mole);
    G4Element* Cr = new G4Element("Chromium", "Cr", 24., 51.9961   *g/mole);
    G4Element* Fe = new G4Element("Iron",     "Fe", 26., 55.85     *g/mole);
    G4Element* Co = new G4Element("Cobalt",   "Co", 27., 58.9332   *g/mole);
    G4Element* Ni = new G4Element("Nickel",   "Ni", 28., 58.6934   *g/mole);
    G4Element* Cu = new G4Element("Copper",   "Cu", 29., 63.55     *g/mole);
    G4Element* I  = new G4Element("Iodine",   "I",  53., 126.90447 *g/mole);
    G4Element* Xe = new G4Element("Xenon",    "Xe", 54., 131.293   *g/mole);
    G4Element* Pb = new G4Element("Lead",     "Pb", 82., 207.2     *g/mole);

    //
    // define an Element from isotopes, by relative abundance
    //

    G4Isotope* U235 = new G4Isotope("Uranium235", 92, 235, 235.0*g/mole);
    G4Element* U  = new G4Element("Uranium", "U", 1);
    U->AddIsotope(U235, 1);

    // define simple materials

    //making Xenon
    G4Material* LXe = new G4Material("LXe",       3.02* g/cm3, 1, kStateLiquid, 173.15*kelvin,  1.0*atmosphere);
    G4Material* GXe = new G4Material("GXe", 10.0*5.894*mg/cm3, 1, kStateGas,    300.00*kelvin,  4.0*atmosphere);
    LXe->AddElement( Xe, 1);
    GXe->AddElement( Xe, 1);

    // define a material from elements for detector parts
    //G4Material* vacuum = new G4Material("vacuum", 1., 1.*g/mole, 1.e-20*g/cm3, kStateGas, 0.1*kelvin, 1.e-20*bar);

    G4Material* SUS304 = new G4Material("SUS304", 8.03*g/cm3, 3);
    SUS304->AddElement(Fe, 0.7); SUS304->AddElement(Ni, 0.1); SUS304->AddElement(Cr, 0.2);

    G4Material* GXe10atm = new G4Material("GXe10atm", 10.0*5.894*mg/cm3, 1, kStateGas, 300.00*kelvin, 10.0*atmosphere);
    GXe10atm->AddElement(Xe, 1);

    G4Material* PEEK = new G4Material("PEEK", 1.3*g/cm3, 3);
    PEEK->AddElement(C,20); PEEK->AddElement(O,3); PEEK->AddElement(H,12);

    G4Material* MPPC = new G4Material("MPPC", 1.0*g/cm3, 2);//??????
    MPPC->AddElement(Al,2); MPPC->AddElement(O,3);

    G4Material* polystyrene = new G4Material("polystyrene", 1.05*g/cm3, 2);
    polystyrene->AddElement(C,1); polystyrene->AddElement(H,1);

    G4Material* metalCu = new G4Material("metalCu", 8.960*g/cm3, 1);
    metalCu->AddElement(Cu,1);

    G4Material* PTFE = new G4Material("PTFE", 2.2*g/cm3, 2);
    PTFE->AddElement(C,1); PTFE->AddElement(F,2);

    G4Material* metalPb = new G4Material("metalPb", 11.340*g/cm3, 1);
    metalPb->AddElement(Pb, 1);

    G4Material* metalAl = new G4Material("metalAl", 2.700*g/cm3, 1);
    metalAl->AddElement(Al, 1);

    G4Material* quartz = new G4Material("quartz", 2.2*g/cm3, 2);
    quartz->AddElement(Si, 1); quartz->AddElement(O , 2);

    G4Material* photocathode = new G4Material("photocathode", 2.700*g/cm3, 1);
    photocathode->AddElement(Al, 1);


    //define gas pressure
    G4double density_Xe_1atm  = 5.8971 *mg/cm3;//noble gas detector(273K,1atm)
    G4double density_He_1atm  = 0.17850*mg/cm3;//noble gas detector(273K,1atm)
    G4double density_Ar_1atm  = 1.7606 *mg/cm3;//noble gas detector(273K,1atm)
    G4double density_CF4_1atm = 3.76   *mg/cm3;//Nishimura D-ron, Nakamura D-ron (15C,1atm)
    G4double density_1atm = 0.0;

    //define mixed gas
    G4int theGasElementNum = 1;
    G4double theGasPressure = 8.; //1[atom] Ar
    G4Material* GasMix = new G4Material("GasMix", density_1atm*theGasPressure, theGasElementNum, kStateGas, 300.00*kelvin, theGasPressure*atmosphere);
    GasMix->AddElement(Ar, theGasElementNum);
    //GasMix->AddElement(Ar, density_Ar_1atm/density_1atm);


    /*
    for(int i=0; i<theGasTypeNum; i++){
        if(     theGasType[i]=="Xe"){  density_1atm += density_Xe_1atm *theGasRatio[i];  theGasElementNum+=1; }
        else if(theGasType[i]=="He"){  density_1atm += density_He_1atm *theGasRatio[i];  theGasElementNum+=1; }
        else if(theGasType[i]=="Ar"){  density_1atm += density_Ar_1atm *theGasRatio[i];  theGasElementNum+=1; }
        else if(theGasType[i]=="CF4"){ density_1atm += density_CF4_1atm*theGasRatio[i];  theGasElementNum+=2; }
        else{ G4cout<<"ERROR unknown theGasType[i] : "<<theGasType[i]<<G4endl; exit(1); }
    }
    G4cout<<" density="<<density_1atm<<"   "<<density_1atm/mg*cm3<<G4endl;
    G4Material* GasMix;
    if(      theGasPressure > 0){ // set gas (normal use)
        GasMix = new G4Material("GasMix", density_1atm*theGasPressure, theGasElementNum, kStateGas, 300.00*kelvin, theGasPressure*atmosphere);
    }else if(theGasPressure==-1){ // set liquid
        if( theGasType[0]=="Xe" && theGasTypeNum==1 ){
            GasMix = new G4Material("GasMix", 2.96*g/cm3, theGasElementNum, kStateLiquid, 173.15*kelvin, 1.0*atmosphere); // LXe
        }
    }else{
        G4cout<<"ERROR: pressure or gastyle is wrong !"<<G4endl; exit(1);
    }
    for(int i=0; i<theGasTypeNum; i++){
        if(theGasType[i]=="Xe") GasMix->AddElement(Xe, theGasRatio[i]*density_Xe_1atm/density_1atm);
        if(theGasType[i]=="He") GasMix->AddElement(He, theGasRatio[i]*density_He_1atm/density_1atm);
        if(theGasType[i]=="Ar") GasMix->AddElement(Ar, theGasRatio[i]*density_Ar_1atm/density_1atm);
        if(theGasType[i]=="CF4"){
            GasMix->AddElement(C, theGasRatio[i]*density_CF4_1atm/density_1atm *12./88.);
            GasMix->AddElement(F, theGasRatio[i]*density_CF4_1atm/density_1atm *19.*4./88.);
        }
        G4cout<<" i="<<i<<" "<<theGasType[i]<<" "<<theGasRatio[i]<<G4endl;
    }
    */


    // other materials
    G4Material* NaI = new G4Material("NaI", 3.67*g/cm3, 2);
    NaI->AddElement(Na, 1);  NaI->AddElement(I, 1);

    G4Material* B4C = new G4Material("B4C", 2.52*g/cm3, 2);
    B4C->AddElement(B, 4);  B4C->AddElement(C, 1);

    G4Material* NiC = new G4Material("NiC", 2.0*g/cm3, 2); //density ?
    NiC->AddElement(Ni, 1);  NiC->AddElement(C, 1);

    G4Material* graphite = new G4Material("graphite", 2.26*g/cm3, 1);
    graphite->AddElement(C, 1);

    G4Material* metalSi = new G4Material("metalSi", 2.329*g/cm3, 1);
    metalSi->AddElement(Si, 1);

    G4Material* metalFe = new G4Material("MetalIron", 7.874*g/cm3, 1);
    metalFe->AddElement(Fe, 1);

    G4Material* ssteel = new G4Material("Steel", 7.7*g/cm3, 3);
    ssteel->AddElement(C, 0.04); ssteel->AddElement(Fe, 0.88); ssteel->AddElement(Co, 0.08);

    G4Material* Air = new G4Material("AIR", 1.2929*kg/m3, 2, kStateGas, 300.00*kelvin, 1.0*atmosphere);
    Air->AddElement(N, 0.8); Air->AddElement(O , 0.2);

    G4Material* LN2 = new G4Material("LN2", 0.8*g/cm3, 1, kStateLiquid, 77.*kelvin, 1.0*atmosphere);
    LN2->AddElement(N, 1);

    G4Material* concrete = new G4Material("Concrete", 2.3*g/cm3, 6);
    concrete->AddElement(Si, 0.227915); concrete->AddElement(O, 0.60541);   concrete->AddElement(H, 0.09972);
    concrete->AddElement(Ca, 0.04986);  concrete->AddElement(Al, 0.014245); concrete->AddElement(Fe, 0.00285);

    G4Material* water = new G4Material("water", 1.00*g/cm3, 2);
    water->AddElement(H, 2); water->AddElement(O, 1);

    G4Material* wood = new G4Material("wood", 0.9*g/cm3, 3);
    wood->AddElement(H, 4); wood->AddElement(O, 1); wood->AddElement(C, 2);


    G4Material* BC501A = new G4Material("BC501A", 0.874*g/cm3, 2);
    BC501A->AddElement(H, 482); BC501A->AddElement(C, 398);

    // assign materials
    world_mat = concrete;
    lab_mat = Air;
    air_mat = Air;
    SUS304_mat = SUS304;
    gas_mat = GasMix;
    PEEK_mat = PEEK;
    MPPC_mat = MPPC;
    polystyrene_mat = polystyrene;
    metalCu_mat = metalCu;
    PTFE_mat = PTFE;
    metalPb_mat = metalPb;
    metalAl_mat = metalAl;
    quartz_mat = quartz;
    photocathode_mat = metalAl;
    //vacuum_mat = vacuum;
    B4C_mat = B4C;
    NiC_mat = NiC;
    graphite_mat = graphite;
    metalSi_mat = metalSi;
    water_mat = water;
    BC501A_mat = BC501A;
    NaI_mat = NaI;

    // examples of vacuum
    /*
    G4Material* Vacuum =
    new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                   kStateGas, 2.73*kelvin, 3.e-18*pascal);

    G4Material* beam =
    new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,
                   kStateGas, STP_Temperature, 2.e-2*bar);
    beam->AddMaterial(Air, fractionmass=1.);

    G4cout << *(G4Material::GetMaterialTable()) << G4endl;


    //default materials of the World
    defaultMaterial  = Vacuum;
    */
}

///////////////////////////////////////////////////////
void mcDetectorConstruction::SetMaxStep(G4double value)
{
    //--------- example of User Limits -------------------------------
    // below is an example of how to set tracking constraints in a given
    // logical volume

    if (value >0.) {
        maxStep = value;
    } else {
        maxStep = DBL_MAX;
    }
    if (pUserLimits) {
        delete pUserLimits;
    }
    UpdateGeometry();
}
/////////////////////////////////////////////////////////////

void mcDetectorConstruction::SetSensorMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
    if (pttoMaterial) {
        sensorMaterial = pttoMaterial;
        G4cout << " mcDetectorConstruction::SetSensorMaterial:  ";
        G4cout << "Sensor material is " << materialChoice << G4endl;
        UpdateGeometry();
    } else {
        G4cout << " mcDetectorConstruction::SetSensorMaterial:  ";
        G4cout << materialChoice << " is not in the Material Table.";
        G4cout <<G4endl;
    }
}


#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void mcDetectorConstruction::SetMagField(G4double value)
{
    //apply a global uniform magnetic field along Z axis
    G4FieldManager* fieldMgr
            = G4TransportationManager::GetTransportationManager()->GetFieldManager();

    if(magField) delete magField;		//delete the existing magn field

    if(value!=0.){			// create a new one if non nul
        fieldValue = value;
        magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
        fieldMgr->SetDetectorField(magField);
        fieldMgr->CreateChordFinder(magField);
    } else {
        magField = 0;
        fieldMgr->SetDetectorField(magField);
    }
}

#include "G4RunManager.hh"

void mcDetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

void mcDetectorConstruction::SetAnalyzer(mcAnalyzer * analyzer_in){
    analyzer = analyzer_in;
}