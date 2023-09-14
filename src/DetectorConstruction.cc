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
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

namespace B1
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4VPhysicalVolume *DetectorConstruction::Construct()
  {
    // Get nist material manager
    G4NistManager *nist = G4NistManager::Instance();
    // Envelope parameters
    G4double env_sizeXY = 20 * cm, env_sizeZ = 20 * cm;
    // Option to switch on/off checking of volumes overlaps
    G4bool checkOverlaps = true;
    // World
    G4double world_sizeXY = 1.2 * env_sizeXY;
    G4double world_sizeZ = 1.2 * env_sizeZ;
    // Materials
    // G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");  // not used
    G4Material *Vacuum = new G4Material("Vacuum",
                                        1.0, 1.01 * g / mole, 1.0E-25 * g / cm3,
                                        kStateGas, 2.73 * kelvin, 3.0E-18 * pascal);
    G4Material *Be = nist->FindOrBuildMaterial("G4_Be");
    G4Material *Al = nist->FindOrBuildMaterial("G4_Al");
    G4Material *Si = nist->FindOrBuildMaterial("G4_Si");
    G4Material *CdTe = nist->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
    G4Material *W = nist->FindOrBuildMaterial("G4_W"); // Tungsten
                                                       // G4Material* Be = new G4Material("Beryllium", 4., 9.012*g/mole, 1.850*g/cm3);  // is this right?
    G4RotationMatrix *norot = new G4RotationMatrix;
    // geometry
    // All geometry inside of world and inside of envelope

    // Geometry values
    G4double FrontToRearBaffleSeperation = 21.8 * mm;
    G4double RearBaffleToDetectorSeparation = 5.0 * mm;
    G4double DetectorBoxDepth = 3 * cm;

    // Volume Dimensions
    G4double FrontBaffleThickness = 3.0 * mm;
    G4double FrontBaffleSideLength = 101.0 * mm;
    G4double FrontWindowSideLength = 16.0 * mm; // bigger than rear window to accomodate 6 deg field of view

    G4double RearBaffleThickness = 2.0 * mm;
    G4double RearBaffleSideLength = FrontBaffleSideLength;
    G4double RearWindowSideLength = 10.0 * mm;

    G4double BoxThickness = 0.5 * mm;
    G4double BoxSideLength = FrontBaffleSideLength + 2 * BoxThickness;
    G4double BoxDepth = 2 * (FrontToRearBaffleSeperation + RearBaffleToDetectorSeparation);
    G4double BoxInnerWallThickness = 3 * mm;

    // Detector Dimensions
    G4double DetectorThickness = 1.0 * mm;
    G4double DetectorSideLength = 10.0 * mm;
    G4int DetectorNum = 4;
    G4double DetectorSpacing = 2 * cm; // from Detector center to center!
    G4double DetectorStoolThickness = 14.2 * mm;
    G4double DetectorStoolSideLength = 12.0 * mm;
    G4double DetectorBoardThickness = 2.46 * mm;

    G4double DetectorFilterThickness = 1.5 * mm;
    G4double FrontApertureFilterThickness = 1.5 * mm;

    // Volume positions - origin is at detector center
    auto rearBafflePosition = G4ThreeVector(0, 0, DetectorThickness / 2.0 + RearBaffleToDetectorSeparation);
    auto frontBafflePosition = G4ThreeVector(0, 0, DetectorThickness / 2.0 + RearBaffleToDetectorSeparation + FrontToRearBaffleSeperation);
    auto boxPosition = G4ThreeVector(0, 0, 0);
    auto boxWallsPosition = G4ThreeVector(0, 0, rearBafflePosition[2] + FrontToRearBaffleSeperation / 2.0);
    auto detectorBoardPosition = G4ThreeVector(0, 0, -(DetectorThickness / 2.0 + DetectorStoolThickness + DetectorBoardThickness / 2.0));
    auto boxBottomPosition = G4ThreeVector(0, 0, detectorBoardPosition[2] - DetectorBoardThickness / 2.0 - BoxThickness / 2.0);

    G4ThreeVector pos1 = G4ThreeVector(0, 0 * cm, 0 * cm);

    auto worldSolid = new G4Box("World", 0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);
    auto worldLogicalVolume = new G4LogicalVolume(worldSolid, Vacuum, "World");
    auto worldPhysicalVolume = new G4PVPlacement(nullptr,            // no rotation
                                                 G4ThreeVector(),    // at (0,0,0)
                                                 worldLogicalVolume, // its logical volume
                                                 "World",            // its name
                                                 nullptr,            // its mother  volume
                                                 false,              // no boolean operation
                                                 0,                  // copy number
                                                 checkOverlaps);     // overlaps checking

    // Envelope
    auto solidEnv = new G4Box("Envelope", 0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);
    auto logicEnv = new G4LogicalVolume(solidEnv, Vacuum, "Envelope");
    auto physEnv = new G4PVPlacement(nullptr,            // no rotation
                                     G4ThreeVector(),    // at (0,0,0)
                                     logicEnv,           // its logical volume
                                     "Envelope",         // its name
                                     worldLogicalVolume, // its mother  volume
                                     false,              // no boolean operation
                                     0,                  // copy number
                                     checkOverlaps);     // overlaps checking

    // Front Baffle - square with 4 square hole
    // G4Box *frontouterBox = new G4Box("Front Window outer", FrontBaffleSideLength / 2.0, FrontBaffleSideLength / 2.0, FrontBaffleThickness / 2.0);
    // G4Box *frontinnerBox = new G4Box("Front Window inner", FrontWindowSideLength / 2.0, FrontWindowSideLength / 2.0, FrontBaffleThickness / 2.0 + 1 * cm); // add buffer to ensure cut occurs
    // G4SubtractionSolid *frontBaffle = new G4SubtractionSolid("FrontBaffle_S", frontouterBox, frontinnerBox);
    // G4LogicalVolume *frontBaffleLogicalVolume = new G4LogicalVolume(frontBaffle, W, "FrontBaffle_L", 0, 0, 0);
    // new G4PVPlacement(nullptr, frontBafflePosition, frontBaffleLogicalVolume, "FrontBaffle_P", logicEnv, false, 0, checkOverlaps);

    // Rear Baffle - square with square hole
    G4Box *rearHole = new G4Box("Rear Window inner", RearWindowSideLength / 2.0, RearWindowSideLength / 2.0, RearBaffleThickness / 2.0 + 1 * mm); // add buffer to ensure cut occurs
    G4Box *rearBaffleNoHoles = new G4Box("Rear Baffle No Holes", RearBaffleSideLength / 2.0, RearBaffleSideLength / 2.0, RearBaffleThickness / 2.0);
    G4SubtractionSolid *rearBaffle1hole = new G4SubtractionSolid("RearBaffle_S", rearBaffleNoHoles, rearHole, norot, G4ThreeVector(DetectorSpacing, DetectorSpacing, 0));
    G4SubtractionSolid *rearBaffle2hole = new G4SubtractionSolid("RearBaffle_S", rearBaffle1hole, rearHole, norot, G4ThreeVector(DetectorSpacing, -DetectorSpacing, 0));
    G4SubtractionSolid *rearBaffle3hole = new G4SubtractionSolid("RearBaffle_S", rearBaffle2hole, rearHole, norot, G4ThreeVector(-DetectorSpacing, -DetectorSpacing, 0));
    G4SubtractionSolid *rearBaffle4hole = new G4SubtractionSolid("RearBaffle_S", rearBaffle3hole, rearHole, norot, G4ThreeVector(-DetectorSpacing, DetectorSpacing, 0));

    G4LogicalVolume *rearBaffleLogicalVolume;
    rearBaffleLogicalVolume = new G4LogicalVolume(rearBaffle4hole, Al, "RearBaffle_L", 0, 0, 0);
    new G4PVPlacement(nullptr, rearBafflePosition, rearBaffleLogicalVolume, "RearBaffle_P", logicEnv, false, 0, checkOverlaps);

    // Front Baffle - square with square hole
    G4Box *frontHole = new G4Box("Front Window inner", FrontWindowSideLength / 2.0, FrontWindowSideLength / 2.0, FrontWindowSideLength / 2.0 + 1 * mm); // add buffer to ensure cut occurs
    G4Box *frontBaffleNoHoles = new G4Box("Front Baffle No Holes", FrontBaffleSideLength / 2.0, FrontBaffleSideLength / 2.0, FrontBaffleThickness / 2.0);
    G4SubtractionSolid *frontBaffle1hole = new G4SubtractionSolid("FrontBaffle_S", frontBaffleNoHoles, frontHole, norot, G4ThreeVector(DetectorSpacing, DetectorSpacing, 0));
    G4SubtractionSolid *frontBaffle2hole = new G4SubtractionSolid("FrontBaffle_S", frontBaffle1hole, frontHole, norot, G4ThreeVector(DetectorSpacing, -DetectorSpacing, 0));
    G4SubtractionSolid *frontBaffle3hole = new G4SubtractionSolid("FrontBaffle_S", frontBaffle2hole, frontHole, norot, G4ThreeVector(-DetectorSpacing, -DetectorSpacing, 0));
    G4SubtractionSolid *frontBaffle4hole = new G4SubtractionSolid("FrontBaffle_S", frontBaffle3hole, frontHole, norot, G4ThreeVector(-DetectorSpacing, DetectorSpacing, 0));

    G4LogicalVolume *frontBaffleLogicalVolume;
    frontBaffleLogicalVolume = new G4LogicalVolume(frontBaffle4hole, W, "FrontBaffle_L", 0, 0, 0);
    new G4PVPlacement(nullptr, frontBafflePosition, frontBaffleLogicalVolume, "FrontBaffle_P", logicEnv, false, 0, checkOverlaps);

    // Baffle box
    G4Box *outerBox = new G4Box("Box outer", BoxSideLength / 2., BoxSideLength / 2., BoxDepth / 2.0);
    G4Box *innerBox = new G4Box("Box inner", BoxSideLength / 2. - BoxThickness, BoxSideLength / 2. - BoxThickness, BoxDepth / 2.0 + 1 * cm);
    G4SubtractionSolid *box = new G4SubtractionSolid("Box_S", outerBox, innerBox);
    G4LogicalVolume *boxLogicalVolume;
    boxLogicalVolume = new G4LogicalVolume(box, Al, "Box_L", 0, 0, 0);
    new G4PVPlacement(nullptr, boxPosition, boxLogicalVolume, "Box_P", logicEnv, false, 0, checkOverlaps);

    auto boxVisStyle = new G4VisAttributes();
    boxVisStyle->SetColor(G4Color(0.3, 0.3, 0.3)); // red
    boxVisStyle->SetForceWireframe(true);
    boxLogicalVolume->SetVisAttributes(boxVisStyle);

    // Inter Walls detector blockers   - FrontBaffleThickness - RearBaffleThickness
    G4Box *Xwall = new G4Box("XWall", FrontBaffleSideLength / 2., BoxInnerWallThickness / 2., (FrontToRearBaffleSeperation - FrontBaffleThickness) / 2.0);
    G4Box *Ywall = new G4Box("YWall", BoxInnerWallThickness / 2., FrontBaffleSideLength / 2., (FrontToRearBaffleSeperation - FrontBaffleThickness) / 2.0);
    G4UnionSolid *boxInnerWalls = new G4UnionSolid("BoxInnerWalls", Xwall, Ywall, nullptr, G4ThreeVector(0, 0, 0));
    auto boxInnerWallsSolidVolume = new G4LogicalVolume(boxInnerWalls, Al, "BoxInnerWalls_L");
    new G4PVPlacement(nullptr, boxWallsPosition, boxInnerWallsSolidVolume, "BoxInnerWalls_P", logicEnv, false, 0, checkOverlaps);
    trackerVisualizationStyle = new G4VisAttributes();
    trackerVisualizationStyle->SetColor(G4Color(1.0, 0.0, 0.0));
    boxInnerWallsSolidVolume->SetVisAttributes(trackerVisualizationStyle);

    //-----------------------
    // - Sensitive detectors -
    //-----------------------
    for (G4int i = 0; i < DetectorNum; i++) // iterate over each detector
    {
      G4ThreeVector DetectorPos;
      G4ThreeVector detectorFilterPosition;
      G4ThreeVector frontFilterPosition;
      G4ThreeVector detectorStoolPosition;

      trackerVisualizationStyle = new G4VisAttributes();

      if (i == 0)
      {
        DetectorPos = G4ThreeVector(DetectorSpacing, DetectorSpacing, 0);
        trackerVisualizationStyle->SetColor(G4Color(1.0, 0.0, 0.0)); // red
      }
      else if (i == 1)
      {
        DetectorPos = G4ThreeVector(DetectorSpacing, -DetectorSpacing, 0);
        trackerVisualizationStyle->SetColor(G4Color(0.0, 1.0, 0.0)); // green
      }
      else if (i == 2)
      {
        DetectorPos = G4ThreeVector(-DetectorSpacing, -DetectorSpacing, 0);
        trackerVisualizationStyle->SetColor(G4Color(0.0, 0.0, 1.0)); // blue
      }
      else if (i == 3)
      {
        DetectorPos = G4ThreeVector(-DetectorSpacing, DetectorSpacing, 0);
        trackerVisualizationStyle->SetColor(G4Color(1.0, 1.0, 1.0)); // white
      }
      detectorFilterPosition = G4ThreeVector(DetectorPos[0], DetectorPos[1], DetectorThickness / 2.0 + RearBaffleToDetectorSeparation);
      frontFilterPosition = G4ThreeVector(DetectorPos[0], DetectorPos[1], DetectorThickness / 2.0 + RearBaffleToDetectorSeparation + FrontToRearBaffleSeperation);

      // G4RotationMatrix norot = G4RotationMatrix(0, 0, 0);
      // G4Transform3D tr = G4Transform3D(norot, frontFilterPosition);
      // frontBaffleHoles->AddNode(*aFrontHole, tr);

      // Detector Filter
      auto detFilterSolidShape = new G4Box("DetectorFilter_S", DetectorSideLength / 2.0, DetectorSideLength / 2.0, DetectorFilterThickness / 2.0);
      auto detFilterSolidVolume = new G4LogicalVolume(detFilterSolidShape, Be, "DetectorFilter_L");
      new G4PVPlacement(nullptr, detectorFilterPosition, detFilterSolidVolume, "DetectorFilter_P", logicEnv, false, 0, checkOverlaps);
      auto detFilterVisStyle = new G4VisAttributes();
      detFilterVisStyle->SetColor(G4Color(0.0, 0.0, 1.0)); // green
      detFilterSolidVolume->SetVisAttributes(detFilterVisStyle);

      // Front Aperture Filter
      auto frontFilterSolidShape = new G4Box("FrontFilter_S", FrontWindowSideLength / 2.0, FrontWindowSideLength / 2.0, DetectorFilterThickness / 2.0);
      auto frontFilterSolidVolume = new G4LogicalVolume(frontFilterSolidShape, Be, "FrontFilter_L");
      new G4PVPlacement(nullptr, frontFilterPosition, frontFilterSolidVolume, "FrontFilter_P", logicEnv, false, 0, checkOverlaps);
      auto frontFilterVisStyle = new G4VisAttributes();
      frontFilterVisStyle->SetColor(G4Color(0.0, 0.0, 1.0)); // green
      frontFilterSolidVolume->SetVisAttributes(frontFilterVisStyle);

      // Detector Stool or Electronics
      detectorStoolPosition = G4ThreeVector(DetectorPos[0], DetectorPos[1], -DetectorThickness / 2.0 - DetectorStoolThickness / 2.0);
      auto stoolSolid = new G4Box("DetectorStool_S", DetectorStoolSideLength / 2.0, DetectorStoolSideLength / 2.0, DetectorStoolThickness / 2.0);
      auto stoolSolidVolume = new G4LogicalVolume(stoolSolid, Si, "DetectorStool_L");
      new G4PVPlacement(nullptr, detectorStoolPosition, stoolSolidVolume, "DetectorStool_P", logicEnv, false, 0, checkOverlaps);

      // Sensitive Detector
      G4String detName = "tracker_" + std::to_string(i);
      auto trackerSolid = new G4Box(detName, DetectorSideLength / 2.0, DetectorSideLength / 2.0, DetectorThickness / 2.0);
      trackerLogicalVolume = new G4LogicalVolume(trackerSolid, CdTe, detName, nullptr, nullptr, nullptr);
      auto trackerPhysicalVolume = new G4PVPlacement(nullptr,
                                                     DetectorPos,
                                                     trackerLogicalVolume,
                                                     detName,
                                                     logicEnv,
                                                     false,
                                                     0,
                                                     checkOverlaps);

      trackerVisualizationStyle->SetColor(G4Color(1.0, 0.0, 0.0)); // red
      trackerLogicalVolume->SetVisAttributes(trackerVisualizationStyle);
      fScoringVolume = trackerLogicalVolume;
    }

    // Detector Board
    auto detectorBoard = new G4Box("DetectorBoard", BoxSideLength / 2., BoxSideLength / 2., DetectorBoardThickness / 2.);
    auto detectorBoardSolidVolume = new G4LogicalVolume(detectorBoard, Si, "DetectorBoard_L");
    new G4PVPlacement(nullptr, detectorBoardPosition, detectorBoardSolidVolume, "DetectorBoard_P", logicEnv, false, 0, checkOverlaps);
    auto VisualizationStyle = new G4VisAttributes();
    VisualizationStyle->SetColor(G4Color(0.0, 1.0, 0.0)); // green
    detectorBoardSolidVolume->SetVisAttributes(VisualizationStyle);

    // Baffle box - rear plate
    auto outerBoxRearSolidShape = new G4Box("BoxRear", BoxSideLength / 2., BoxSideLength / 2., BoxThickness / 2.0);
    auto outerBoxRearSolidVolume = new G4LogicalVolume(outerBoxRearSolidShape, Al, "BoxRear_L");
    new G4PVPlacement(nullptr, boxBottomPosition, outerBoxRearSolidVolume, "BoxRear_P", logicEnv, false, 0, checkOverlaps);

    // frontBaffleHoles->Voxelize();
    // G4SubtractionSolid *frontBaffle2 = new G4SubtractionSolid("FrontBaffle_S", frontBaffleNoHoles, frontBaffleHoles);
    // G4LogicalVolume *frontBaffleLogicalVolume2 = new G4LogicalVolume(frontBaffle2, W, "FrontBaffle_L2", 0, 0, 0);
    new G4PVPlacement(nullptr, frontBafflePosition, frontBaffleLogicalVolume, "FrontBaffle_P2", logicEnv, false, 0, checkOverlaps);
    // constexpr auto TRACKER_SENSITIVE_DETECTOR_NAME{"meddea/DetectorSD"};
    // auto *trackerSensitiveDetector = new calisteDetectorSD(TRACKER_SENSITIVE_DETECTOR_NAME);
    // SetSensitiveDetector(trackerLogicalVolume, trackerSensitiveDetector);

    // Set Shape2 as scoring volume
    //
    fScoringVolume = trackerLogicalVolume;

    return worldPhysicalVolume;
  }

  void DetectorConstruction::ConstructSDandField()
  {
    // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
    auto *sensitiveDetectorManager = G4SDManager::GetSDMpointer();

    constexpr auto TRACKER_SENSITIVE_DETECTOR_NAME{"meddea/DetectorSD"};
    auto *trackerSensitiveDetector = new calisteDetectorSD(TRACKER_SENSITIVE_DETECTOR_NAME);
    sensitiveDetectorManager->AddNewDetector(trackerSensitiveDetector);
    // the following is the problem line
    // G4cout << "MySensitiveDetector::I'M HERE!! " << G4endl;
    SetSensitiveDetector(trackerLogicalVolume, trackerSensitiveDetector);
    // G4cout << "MySensitiveDetector::I'M HERE 2!! " << G4endl;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
