#include "mcDetectorConstruction.hh"
#include "mcDetectorMessenger.hh"
#include "mcSensorSD.hh"
#include "mcAnalyzer.hh"

//#include "DMXDetectorConstruction.hh"
//#include "DMXDetectorMaterial.ihh"
//#include "DMXTPCSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
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
//#include "G4PhysicalConstants.hh"

#include <iostream>

mcDetectorConstruction::mcDetectorConstruction()
        :defaultMaterial(0),sensorMaterial(0),
         WorldRadius(100*cm),
         solidWorld(0),logicWorld(0),physWorld(0),
         solidSensor(0),logicSensor(0),physSensor(0),
         magField(0),pUserLimits(0),maxStep(100.0*cm),
         nShieldSize(0), nShieldThetaMax(0),
         gShield1Thickness(0)
{

    // default parameter values of Sensor
    DefineMaterials();

    SetSensorMaterial("Ar1atm");
    SetNeutronShieldMaterial("polyethylene_boron10");

    // user limits   G4userLimits(step-length-max, track-length-max, time-cut, min-energy)
    double step_size = 0.1*mm;
    pUserLimits = new G4UserLimits(step_size, DBL_MAX, DBL_MAX, 0, 0);

    // create commands for interactive definition of the calorimeter
    detectorMessenger = new mcDetectorMessenger(this);
    l_gas             = 30.*cm;
    nShieldSize       = 25.*cm;
    nShieldThetaMax   = 180.*deg;
    nShieldShape      = "cube";
    gShield1Thickness = 5*cm;
    gShield2Thickness = 5*cm;

    /*
    l_gas             = 30*cm;  // length of gas box
    t_fiducial        = 1*cm;   // thickness of fiducial cut area
    d                 = 100*cm; // distance of beam point to centere of chamber
    theta_cone        = atan( (l_gas/2. - t_fiducial )/ (d - (l_gas/2. - t_fiducial )) )*(180/M_PI)*deg;
    l_cone            = nShieldSize*tan(theta_cone);
    */

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


    // create geometry

    // World
    G4double world_l    = 100.0*m;
    G4Material* vac_mat = G4Material::GetMaterial("vacuum");

    solidWorld = new G4Box("world_box", world_l, world_l, world_l);
    logicWorld = new G4LogicalVolume( solidWorld, vac_mat, "world_lv", 0, 0, 0);
    G4VisAttributes* world_att = new G4VisAttributes( false );
    logicWorld->SetVisAttributes(world_att);
    physWorld = new G4PVPlacement(0, G4ThreeVector(), "world_pv", logicWorld, 0, false,0);

    /*
    G4double bulk_l = 4*mm; //100000*mm;
    G4double box_x = 1000*mm;
    G4double box_y = 300*mm;
    G4double box_z = 300*mm;
    G4double buffer_t = 300*mm;
    G4double chamber_t = 300*mm;

    G4double test_l = 300*mm;
    */


    ConstructLaboratory();
    ConstructBeamShield();
    //ConstructGammaShield1();
    //ConstructSphereBeamShield();
    //ConstructCubeBeamShield();
    ConstructChamber();


    //////////////////////////////// devided gas box (for Intrinsic BG) ////////////////////////////////
    /*
    G4int Ndiv_x = 15;
    G4int Ndiv_y = 15;
    G4int Ndiv_z = 15;
    G4int Ncopy = Ndiv_x * Ndiv_y * Ndiv_z;

    G4Colour ArColor (0,0.9,0.1,0.9);
    G4Material* gas_mat = G4Material::GetMaterial("Ar1atm");

    G4Box* test_box = new G4Box("test", 0.5*test_l, 0.5*test_l, 0.5*test_l);
    G4LogicalVolume* test_lv = new G4LogicalVolume(test_box, gas_mat, "test_lv");
    G4VisAttributes* test_att = new G4VisAttributes(1, ArColor);  test_lv->SetVisAttributes(test_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "test_pv", test_lv, physWorld, false, Ncopy);
    test_lv->SetUserLimits(pUserLimits);
    */

    /*
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
    /*
    G4Material* chamber_mat = G4Material::GetMaterial("Polyeth");
    G4double chamber_width = 10*mm; //SUS 3mm, Al 5mm, plastic 10mm

    G4Box* chamber_box = new G4Box("chamber", 0.5*(chamber_t+chamber_width*2)*mm, 0.5*(chamber_t+chamber_width*2)*mm, 0.5*(chamber_t+chamber_width*2)*mm);
    G4LogicalVolume* chamber_lv = new G4LogicalVolume(chamber_box, chamber_mat, "chamber_lv");
    G4VisAttributes* chamber_att = new G4VisAttributes(1, blue);  chamber_lv->SetVisAttributes(chamber_att);
    //G4PVPlacement* chamber_pv = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "chamber_pv", chamber_lv, world_pv, false, 0);
    chamber_lv->SetUserLimits(user_limit);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "chamber_pv", chamber_lv, world_pv, false, 0);


    G4Box* vac_box = new G4Box("vac_box", 0.5*chamber_t*mm, 0.5*chamber_t, 0.5*chamber_t);
    G4LogicalVolume* vac_lv = new G4LogicalVolume(vac_box, vac_mat, "vac_lv");
    G4VisAttributes* vac_att = new G4VisAttributes(1, red);  vac_lv->SetVisAttributes(vac_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), vac_lv, "vac_pv", chamber_lv, false, 0);
    vac_lv->SetUserLimits(user_limit);

    //new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "chamber_pv", chamber_lv, world_pv, false, 0);
    */
    ////////////////////////////////////////// material test ////////////////////////////////////////////
    /*
    G4Material* test_mat = G4Material::GetMaterial("LiF");
    G4double l_mat_test = 100*m;
    G4Colour testColor (0,0,1,0.9);

    G4Box* mat_test_box = new G4Box("mat_test", 0.5*l_mat_test*mm, 0.5*l_mat_test*mm, 0.5*l_mat_test*mm);
    G4LogicalVolume* mat_test_lv = new G4LogicalVolume(mat_test_box, test_mat, "mat_test_lv");
    G4VisAttributes* mat_test_att = new G4VisAttributes(1, testColor);  mat_test_lv->SetVisAttributes(mat_test_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),"mat_test_pv", mat_test_lv, physWorld, false, 0);
    mat_test_lv->SetUserLimits(pUserLimits);
    logicSensor = mat_test_lv;
    */


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

    /*
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    mcSensorSD* aSensorSD = (mcSensorSD*)SDman->FindSensitiveDetector("mc/SensorSD");
    if ( aSensorSD == 0){
        aSensorSD = new mcSensorSD("mc/SensorSD");
        SDman->AddNewDetector( aSensorSD );
    }
    aSensorSD->SetAnalyzer(analyzer);
    logicSensor->SetSensitiveDetector(aSensorSD);
    */

    //test_lv->SetSensitiveDetector(aSensorSD);

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
    */

    return physWorld;

}
// End of Construct()
//------------------------------------------------------------------------//

void mcDetectorConstruction::ConstructLaboratory()
{
    // parameters
    G4double lab_w  = 11.5*m, lab_h = 11.5*m, concrete_t = 1.0*m, floor_h = 7.8*mm;
    G4double concrete_w = lab_w  + concrete_t, concrete_h = lab_h + concrete_t;

    // colors
    G4Colour ConcreteColor (0.5, 0.5, 0.5, 1.0);
    G4Colour AirColor      (0.5, 0.5, 0.5, 1.0);
    G4Colour AlColor       (0.7, 0.7, 0.7, 1.0);

    // materials
    G4Material* concrete_mat = G4Material::GetMaterial("concrete");
    G4Material* air_mat      = G4Material::GetMaterial("Air");
    G4Material* floor_mat    = G4Material::GetMaterial("metalAl");

    // Concrete wall
    G4Box *concrete_box = new G4Box("concrete_box", 0.5*concrete_w, 0.5*concrete_h, 0.5*concrete_w);
    G4LogicalVolume* concrete_lv = new G4LogicalVolume(concrete_box, concrete_mat, "concrete_lv");
    G4VisAttributes* concrete_att = new G4VisAttributes(1, ConcreteColor);
    concrete_lv->SetVisAttributes(concrete_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "concrete_pv", concrete_lv, physWorld, false, 0);
    concrete_lv->SetUserLimits(pUserLimits);

    // Lab air
    G4Box *lab_box = new G4Box("lab_box", 0.5*lab_w, 0.5*lab_h, 0.5*lab_w);
    G4LogicalVolume* lab_lv = new G4LogicalVolume(lab_box, air_mat, "lab_lv");
    G4VisAttributes* lab_att = new G4VisAttributes(1, AirColor);
    lab_lv->SetVisAttributes(lab_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), lab_lv, "lab_pv", concrete_lv, false, 0);
    lab_lv->SetUserLimits(pUserLimits);
    logicLab = lab_lv;

    // Floor Al
    G4Box *floor_box = new G4Box("floor_box", 0.5*lab_w, 0.5*floor_h, 0.5*lab_w);
    G4LogicalVolume* floor_lv = new G4LogicalVolume(floor_box, floor_mat, "floor_lv");
    G4VisAttributes* floor_att = new G4VisAttributes(1, AlColor);
    floor_lv->SetVisAttributes(floor_att);
    new G4PVPlacement(0, G4ThreeVector(0.,-1.0*m,0.), floor_lv, "floor_pv", lab_lv, false, 0);
    floor_lv->SetUserLimits(pUserLimits);
}

void mcDetectorConstruction::ConstructBeamShield()
{
    if( (nShieldShape=="sphere") || (nShieldShape=="hemisphere") ){
        ConstructSphereBeamShield();
    } else if (nShieldShape=="cube"){
        ConstructCubeBeamShield();
    } else if (nShieldShape=="none"){
        // no shield
    } else {
        G4cout << "unkown nshield shape " << nShieldShape << G4endl;
    }
}

void mcDetectorConstruction::ConstructSphereBeamShield()
{
    G4double naikei = 0.*mm;
    G4double gaikei = nShieldSize;
    G4double startPhi = 0.*deg;
    G4double endPhi = 360.*deg;
    G4double startTheta = 10.*deg;
    G4double endTheta = nShieldThetaMax;

    G4RotationMatrix* rm_shield = new G4RotationMatrix();
    rm_shield->rotateY(180.*deg);

    //G4Material* shield_mat = G4Material::GetMaterial("polyethylene_boron10");
    G4Material* nshield_mat = nShieldMaterial;

    G4Sphere* nshield = new G4Sphere("shield", naikei, gaikei, startPhi, endPhi, startTheta, endTheta);
    G4LogicalVolume* nshield_lv = new G4LogicalVolume(nshield, nshield_mat, "nshield_lv");
    new G4PVPlacement(0, G4ThreeVector(0.,0.,-1.*m), nshield_lv, "nshield_pv", logicLab, false, 0);

}

void mcDetectorConstruction::ConstructCubeBeamShield()
{

    G4double l_gas_local             = 30*cm;  // length of gas box
    G4double t_fiducial_local        = 1*cm;   // thickness of fiducial cut area
    G4double d_local                 = 100*cm; // distance of beam point to centere of chamber
    G4double theta_cone_local        = atan( (l_gas_local/2. - t_fiducial_local )/ (d_local - (l_gas_local/2. - t_fiducial_local )) )*(180/M_PI)*deg;
    G4double l_cone_local            = nShieldSize*tan(theta_cone_local);
    l_cone = l_cone_local;

    //G4Material* cube_mat = G4Material::GetMaterial("polyethylene_boron10");
    //G4Material* cone_mat = G4Material::GetMaterial("Air");
    //G4Material* nshield_mat = G4Material::GetMaterial("polyethylene_boron10");
    G4Material* nshield_mat = nShieldMaterial;
    G4Colour ShieldColor = (0.1, 0.1, 0.1, 1.0);
    G4Material* gshield1_mat = G4Material::GetMaterial("metalPb");
    G4Colour gShield1Color = G4Colour(0.8,0.8,0.8, 1.0);
    G4Material* gshield2_mat = G4Material::GetMaterial("metalPb");
    G4Colour gShield2Color = G4Colour(0.3,0.3,0.3, 1.0);

    G4double halfedge_nshield    = nShieldSize;
    G4double halfedge_squarehole = (nShieldSize+gShield1Thickness)*tan(theta_cone_local);

    G4Box* cube = new G4Box("cube", halfedge_nshield, halfedge_nshield, halfedge_nshield);
    /*
    G4LogicalVolume* cube_lv = new G4LogicalVolume(cube, cube_mat, "cube_lv");
    G4VisAttributes* cube_att = new G4VisAttributes(1, ShieldColor);  cube_lv->SetVisAttributes(cube_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), cube_lv,"cube_pv", logicLab, false, 0);
    cube_lv->SetUserLimits(pUserLimits);
    */

    G4Trd* cone = new G4Trd("cone", 0., l_cone, 0, l_cone, nShieldSize/2.);
    /*
    G4LogicalVolume* cone_lv = new G4LogicalVolume(cone, cone_mat, "cone_lv");
    G4VisAttributes* cone_att = new G4VisAttributes(1, ShieldColor);  cone_lv->SetVisAttributes(cone_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.*m,0), cone_lv,"cone_pv", logicLab, false, 0);
    cone_lv->SetUserLimits(pUserLimits);
    */

    G4VSolid* nshield = new G4SubtractionSolid("nshield", cube, cone, 0, G4ThreeVector (0.,0.,nShieldSize/2.));
    G4LogicalVolume* nshield_lv = new G4LogicalVolume(nshield, nshield_mat, "nshield_lv");
    G4VisAttributes* nshield_att = new G4VisAttributes(1, ShieldColor);  nshield_lv->SetVisAttributes(nshield_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,-1*m), nshield_lv,"nshield_pv", logicLab, false, 0);
    nshield_lv->SetUserLimits(pUserLimits);

    // gamma shield 1
    if (gShield1Thickness != 0){
        G4Box* part1 = new G4Box("part1", nShieldSize, nShieldSize, gShield1Thickness/2.);
        G4Box* part2 = new G4Box("part2", halfedge_squarehole, halfedge_squarehole, gShield1Thickness/2.);
        G4VSolid* gshield1 = new G4SubtractionSolid("gshield1", part1, part2, 0, G4ThreeVector(0.,0.,0.));
        G4LogicalVolume* gshield1_lv = new G4LogicalVolume(gshield1, gshield1_mat, "gshield1_lv");
        G4VisAttributes* gshield1_att = new G4VisAttributes(1, gShield1Color);  gshield1_lv->SetVisAttributes(gshield1_att);
        new G4PVPlacement(0, G4ThreeVector(0.,0.,-100*cm+nShieldSize+gShield1Thickness/2.), gshield1_lv, "gshield1_pv", logicLab, false, 0);
        gshield1_lv->SetUserLimits(pUserLimits);
    }

    // gamma shield 2
    if (gShield2Thickness != 0) {
        G4Box *gshield2 = new G4Box("gshield2", l_gas / 2., l_gas / 2., gShield2Thickness / 2.);
        G4LogicalVolume *gshield2_lv = new G4LogicalVolume(gshield2, gshield2_mat, "gshield2_lv");
        new G4PVPlacement(0, G4ThreeVector(0., 0., -0.5 * (l_gas + gShield2Thickness)), gshield2_lv, "gshield2_lv", logicLab, false, 0);
    }
}

void mcDetectorConstruction::ConstructGammaShield1()
{

    G4Material* gshield1_mat = G4Material::GetMaterial("metalPb");
    G4Colour gShield1Color = G4Colour(0.8,0.8,0.8, 1.0);
    /*
    G4double l_gas             = 30*cm;  // length of gas box
    G4double t_fiducial        = 1*cm;   // thickness of fiducial cut area
    G4double d                 = 100*cm; // distance of beam point to centere of chamber
    G4double theta_cone        = atan( (l_gas/2. - t_fiducial )/ (d - (l_gas/2. - t_fiducial )) )*(180/M_PI)*deg;
    G4double l_cone            = nShieldSize*tan(theta_cone);
    */
    //gShield1Thickness          = 5*cm;


    G4Box* part1 = new G4Box("part1", nShieldSize, nShieldSize, gShield1Thickness);
    G4Box* part2 = new G4Box("part2", l_cone, l_cone, gShield1Thickness);

    G4VSolid* gshield1 = new G4SubtractionSolid("gshield1", part1, part2, 0, G4ThreeVector(0.,0.,0.));
    G4LogicalVolume* gshield1_lv = new G4LogicalVolume(gshield1, gshield1_mat, "gshield1_lv");
    G4VisAttributes* gshield1_att = new G4VisAttributes(1, gShield1Color);  gshield1_lv->SetVisAttributes(gshield1_att);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,-100*cm+nShieldSize+gShield1Thickness/2.), gshield1_lv, "cone_pv", logicLab, false, 0);
    gshield1_lv->SetUserLimits(pUserLimits);

}

void mcDetectorConstruction::ConstructChamber()
{

    G4int Ndiv_x = 30;
    G4int Ndiv_y = 30;
    G4int Ndiv_z = 30;
    //G4int Ncopy = Ndiv_x * Ndiv_y * Ndiv_z;

    l_gas = 30*cm;
    G4Colour ArColor (0,0.9,0.1,1.);
    G4Colour GasColor = ArColor;
    G4Material* gas_mat = G4Material::GetMaterial("vacuum");
    //G4Material* gas_mat = sensorMaterial;

    G4Box* gas_box = new G4Box("gas", 0.5*l_gas, 0.5*l_gas, 0.5*l_gas);
    G4LogicalVolume* gas_lv = new G4LogicalVolume(gas_box, gas_mat, "gas_lv");
    G4VisAttributes* gas_att = new G4VisAttributes(1, GasColor);  gas_lv->SetVisAttributes(gas_att);
    //new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "test_pv", gas_lv, physWorld, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), gas_lv,"gas_pv", logicLab, false, 0);
    gas_lv->SetUserLimits(pUserLimits);

    G4Box* rep1 = new G4Box("rep1", l_gas*0.5, l_gas*0.5, (l_gas/Ndiv_z)*0.5);
    G4LogicalVolume* rep1_lv = new G4LogicalVolume(rep1, gas_mat, "rep1_lv");
    G4VisAttributes* rep1_att = new G4VisAttributes(1, GasColor);  rep1_lv->SetVisAttributes(rep1_att);
    new G4PVReplica("rep1_pv", rep1_lv, gas_lv, kZAxis, Ndiv_z, l_gas/Ndiv_z);

    G4Box* rep2 = new G4Box("rep2", l_gas*0.5, (l_gas/Ndiv_y)*0.5, (l_gas/Ndiv_z)*0.5);
    G4LogicalVolume* rep2_lv = new G4LogicalVolume(rep2, gas_mat, "rep2_lv");
    G4VisAttributes* rep2_att = new G4VisAttributes(1, GasColor);  rep2_lv->SetVisAttributes(rep2_att);
    new G4PVReplica("rep2_pv", rep2_lv, rep1_lv, kYAxis, Ndiv_y, l_gas/Ndiv_y);

    G4Box* rep3 = new G4Box("rep3", (l_gas/Ndiv_x)*0.5, (l_gas/Ndiv_y)*0.5, (l_gas/Ndiv_z)*0.5);
    G4LogicalVolume* rep3_lv = new G4LogicalVolume(rep3, gas_mat, "rep3_lv");
    G4VisAttributes* rep3_att = new G4VisAttributes(1, GasColor);  rep3_lv->SetVisAttributes(rep3_att);
    new G4PVReplica("rep3_pv", rep3_lv, rep2_lv, kXAxis, Ndiv_x, l_gas/Ndiv_x);

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    mcSensorSD* aSensorSD = (mcSensorSD*)SDman->FindSensitiveDetector("mc/SensorSD");
    if ( aSensorSD == 0){
        aSensorSD = new mcSensorSD("mc/SensorSD");
        SDman->AddNewDetector( aSensorSD );
    }
    aSensorSD->SetAnalyzer(analyzer);
    rep3_lv->SetSensitiveDetector(aSensorSD);
}

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
    G4Element* Li = new G4Element("lithium",  "Li", 3.,  6.941     *g/mole);
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
    G4Material* vacuum = new G4Material("vacuum", 1., 1.*g/mole, 1.e-20*g/cm3, kStateGas, 0.1*kelvin, 1.e-20*bar);

    G4Material* SUS304 = new G4Material("SUS304", 8.03*g/cm3, 3);
    SUS304->AddElement(Fe, 0.7); SUS304->AddElement(Ni, 0.1); SUS304->AddElement(Cr, 0.2);

    G4Material* GXe10atm = new G4Material("GXe10atm", 10.0*5.894*mg/cm3, 1, kStateGas, 300.00*kelvin, 10.0*atmosphere);
    GXe10atm->AddElement(Xe, 1);

    G4Material* PEEK = new G4Material("PEEK", 1.3*g/cm3, 3);
    PEEK->AddElement(C,20); PEEK->AddElement(O,3); PEEK->AddElement(H,12);

    G4Material* MPPC = new G4Material("MPPC", 1.0*g/cm3, 2);//??????
    MPPC->AddElement(Al,2); MPPC->AddElement(O,3);

    G4Material* polyethylene = new G4Material("polyethylene",0.95*g/cm3,2);
    polyethylene->AddElement(H, 2); polyethylene->AddElement(C,1);

    G4Material* polyethylene_boron10 = new G4Material("polyethylene_boron10",0.98*g/cm3,4);
    polyethylene_boron10->AddElement(H, 0.129); polyethylene_boron10->AddElement(C,0.771);
    polyethylene_boron10->AddElement(B,0.031); polyethylene_boron10->AddElement(O, 0.069);

    G4Material* polyethylene_boron20 = new G4Material("polyethylene_boron20",1.04*g/cm3,4);
    polyethylene_boron20->AddElement(H, 0.114); polyethylene_boron20->AddElement(C,0.686);
    polyethylene_boron20->AddElement(B,0.062); polyethylene_boron20->AddElement(O, 0.138);

    G4Material* LiF = new G4Material("LiF", 2.5*g/cm3, 2);
    LiF->AddElement(Li, 1); LiF->AddElement(F, 1);

    G4Material* polyethylene_LiF01 = new G4Material("polyethylene_LiF01", 0.96*g/cm3, 2);
    polyethylene_LiF01->AddMaterial(polyethylene, 0.99); polyethylene_LiF01->AddMaterial(LiF,0.01);

    G4Material* polyethylene_LiF10 = new G4Material("polyethylene_LiF10", 1.01*g/cm3, 2);
    polyethylene_LiF10->AddMaterial(polyethylene, 0.9); polyethylene_LiF10->AddMaterial(LiF, 0.1);

    G4Material* polyethylene_LiF50 = new G4Material("polyethylene_LiF50", 1.38*g/cm3, 2);
    polyethylene_LiF50->AddMaterial(polyethylene, 0.5); polyethylene_LiF50->AddMaterial(LiF, 0.5);

    G4Material* acrylic = new G4Material("acrylic", 1.2*g/cm3, 3);
    acrylic->AddElement(H, 8); acrylic->AddElement(C, 5); acrylic->AddElement(O, 2);

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

    G4Material* Ar1atm = new G4Material("Ar1atm", density_Ar_1atm, 1, kStateGas, 300.00*kelvin, 1*atmosphere);
    Ar1atm->AddElement(Ar, 1);
    G4Material* Xe8atm = new G4Material("Xe8atm", density_Xe_1atm*8, 1, kStateGas, 300.00*kelvin, 1*atmosphere);
    Xe8atm->AddElement(Xe, 1);


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

    G4Material* Air = new G4Material("Air", 1.2929*kg/m3, 3, kStateGas, 300.00*kelvin, 1.0*atmosphere);
    //Air->AddElement(N, 0.8); Air->AddElement(O , 0.2);
    Air->AddElement(N,0.756); Air->AddElement(O,0.231); Air->AddElement(Ar,0.013);

    G4Material* LN2 = new G4Material("LN2", 0.8*g/cm3, 1, kStateLiquid, 77.*kelvin, 1.0*atmosphere);
    LN2->AddElement(N, 1);

    G4Material* concrete = new G4Material("concrete", 2.3*g/cm3, 6);
    concrete->AddElement(Si, 0.227915); concrete->AddElement(O, 0.60541);   concrete->AddElement(H, 0.09972);
    concrete->AddElement(Ca, 0.04986);  concrete->AddElement(Al, 0.014245); concrete->AddElement(Fe, 0.00285);

    G4Material* water = new G4Material("water", 1.00*g/cm3, 2);
    water->AddElement(H, 2); water->AddElement(O, 1);

    G4Material* wood = new G4Material("wood", 0.9*g/cm3, 3);
    wood->AddElement(H, 4); wood->AddElement(O, 1); wood->AddElement(C, 2);


    G4Material* BC501A = new G4Material("BC501A", 0.874*g/cm3, 2);
    BC501A->AddElement(H, 482); BC501A->AddElement(C, 398);

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

void mcDetectorConstruction::SetNeutronShieldMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
    if (pttoMaterial) {
        nShieldMaterial = pttoMaterial;
        G4cout << " mcDetectorConstruction::SetNeutronShieldMaterial:  ";
        G4cout << "Neutron Shield material is " << materialChoice << G4endl;
        UpdateGeometry();
    } else {
        G4cout << " mcDetectorConstruction::SetNeutronShieldMaterial:  ";
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

void mcDetectorConstruction::SetGasBoxSize(G4double value)
{
    l_gas = value;
    UpdateGeometry();
}

void mcDetectorConstruction::SetNeutronShieldSize(G4double value) {
    nShieldSize = value;
    UpdateGeometry();
}

void mcDetectorConstruction::SetNeutronShieldType(G4String nShieldType)
{
    if( nShieldType ==  "sphere" ){
        nShieldThetaMax = 180.*deg;
        nShieldShape = nShieldType;
        UpdateGeometry();
    }else if( nShieldType == "hemisphere" ) {
        nShieldThetaMax = 90. * deg;
        nShieldShape = nShieldType;
        UpdateGeometry();
    }else if( nShieldType == "cube" ) {
        nShieldShape = nShieldType;
        UpdateGeometry();
    }else if( nShieldType == "none" ){
        nShieldShape = nShieldType;
        UpdateGeometry();
    }else{
        G4cout << "unknown shield type: " << nShieldType << G4endl;
    }
}

void mcDetectorConstruction::SetGammaShield1Thickness(G4double value)
{
    gShield1Thickness = value;
    UpdateGeometry();
}

void mcDetectorConstruction::SetGammaShield2Thickness(G4double value)
{
    gShield2Thickness = value;
    UpdateGeometry();
}

#include "G4RunManager.hh"

void mcDetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

void mcDetectorConstruction::SetAnalyzer(mcAnalyzer * analyzer_in){
    analyzer = analyzer_in;
}