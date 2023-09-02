//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "calisteDetectorSD.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  // Envelope parameters
  G4double env_sizeXY = 20*cm, env_sizeZ = 20*cm;
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;
  // World
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  // Materials
  // G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");  // not used
  G4Material* Vacuum = new G4Material("Vacuum",
				      1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
				      kStateGas, 2.73*kelvin, 3.0E-18*pascal );
  G4Material* Be = nist->FindOrBuildMaterial("G4_Be");
  G4Material* Al = nist->FindOrBuildMaterial("G4_Al");
  G4Material* CdTe = nist->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
  G4Material* W = nist->FindOrBuildMaterial("G4_W");  // Tungsten
  //G4Material* Be = new G4Material("Beryllium", 4., 9.012*g/mole, 1.850*g/cm3);  // is this right?

  auto worldSolid = new G4Box("World",                           
    0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  
  auto worldLogicalVolume = new G4LogicalVolume(worldSolid, Vacuum, "World");
  auto worldPhysicalVolume = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    worldLogicalVolume,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

  // Envelope
  auto solidEnv = new G4Box("Envelope",                    
    0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.75 * env_sizeZ);
  auto logicEnv = new G4LogicalVolume(solidEnv, Vacuum, "Envelope");                                 
  auto physEnv = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    logicEnv,                 // its logical volume
    "Envelope",               // its name
    worldLogicalVolume,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  // Shape 1
  //

  G4ThreeVector pos1 = G4ThreeVector(0, 0*cm, 0*cm);
  G4double BeWindow_thick = 5 * mm;
  auto solidShape1 = new G4Box("BeWindow",                    // its name
    0.5 * env_sizeXY, 0.5 * env_sizeXY, BeWindow_thick);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
    Al,                                        // its material
    "BeWindow_L");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos1,                     // at position
    logicShape1,              // its logical volume
    "BeWindow_P",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  //
  // Shape 2
  //

  G4ThreeVector pos2 = G4ThreeVector(0, 0*cm, 7*cm);

  // Trapezoid shape
  //G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  //G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  //G4double shape2_dz  = 6*cm;
  //auto solidShape2 = new G4Trd("Shape2",  // its name
  //  0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
  //  0.5 * shape2_dz);  // its size

  //auto solidShape2 = new G4Box("Shape2",                    // its name
  //  0.5 * env_sizeXY, 0.5 * env_sizeXY, 1 * mm);

  //auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
  //  Be,                                        // its material
  //  "Shape2");                                         // its name

  //new G4PVPlacement(nullptr,  // no rotation
  //  pos2,                     // at position
  //  logicShape2,              // its logical volume
  //  "Shape2",                 // its name
  //  logicEnv,                 // its mother  volume
  //  false,                    // no boolean operation
  //  0,                        // copy number
  //  checkOverlaps);           // overlaps checking

  //-----------------------
  // - Sensitive detector -
  //-----------------------

  G4ThreeVector pos3 = G4ThreeVector(0, 0*cm, -5*cm);
  auto trackerSolid = new G4Box("tracker",  0.5 * 20 * cm, 0.5 * 20 * cm, 1 * mm);
  trackerLogicalVolume = new G4LogicalVolume(trackerSolid, CdTe, "Tracker", nullptr, nullptr, nullptr);
  auto trackerPhysicalVolume = new G4PVPlacement(nullptr, // no rotation
      pos3, // at (x, y, z)
      trackerLogicalVolume, // its logical volume
      "Tracker", // its name
      logicEnv, // its mother volume
      false,
      0, // no boolean operations
      checkOverlaps); // copy number

  //constexpr auto TRACKER_SENSITIVE_DETECTOR_NAME{"meddea/DetectorSD"};
  //auto *trackerSensitiveDetector = new calisteDetectorSD(TRACKER_SENSITIVE_DETECTOR_NAME);
  //SetSensitiveDetector(trackerLogicalVolume, trackerSensitiveDetector);
  trackerVisualizationStyle = new G4VisAttributes();
  trackerVisualizationStyle->SetColor(G4Color(1.0, 0.0, 0.0)); // red
  trackerLogicalVolume->SetVisAttributes(trackerVisualizationStyle);
  // Set Shape2 as scoring volume
  //
  fScoringVolume = trackerLogicalVolume;

  return worldPhysicalVolume;
}

void DetectorConstruction::ConstructSDandField()
{
  //G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  auto *sensitiveDetectorManager = G4SDManager::GetSDMpointer();

  constexpr auto TRACKER_SENSITIVE_DETECTOR_NAME{"meddea/DetectorSD"};
  auto *trackerSensitiveDetector = new calisteDetectorSD(TRACKER_SENSITIVE_DETECTOR_NAME);
  sensitiveDetectorManager->AddNewDetector(trackerSensitiveDetector);
  // the following is the problem line
  // G4cout << "MySensitiveDetector::I'M HERE!! " << G4endl;
  SetSensitiveDetector(trackerLogicalVolume, trackerSensitiveDetector);
  //G4cout << "MySensitiveDetector::I'M HERE 2!! " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
