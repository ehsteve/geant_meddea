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
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

namespace B1
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  SteppingAction::SteppingAction(EventAction *eventAction)
      : fEventAction(eventAction)
  {
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void SteppingAction::UserSteppingAction(const G4Step *step)
  {
    if (!fScoringVolume)
    {
      const auto detConstruction = static_cast<const DetectorConstruction *>(
          G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      fScoringVolume = detConstruction->GetScoringVolume();
    }

    // get volume of the current step
    G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4String volName;
    G4String nextVolName;

    volName = volume->GetName();

    // check if we are in scoring volume
    // if (volume != fScoringVolume) return;
    if (volName != "tracker_0" && volName != "tracker_1" && volName != "tracker_2" && volName != "tracker_3")
    {
      return;
    }
    // G4cout << volName << G4endl;

    AnalysisManager *analysis = AnalysisManager::Instance();

    G4int detectorNum;
    if (volName == "tracker_0")
    {
      detectorNum = 0;
    }
    else if (volName == "tracker_1")
    {
      detectorNum = 1;
    }
    else if (volName == "tracker_2")
    {
      detectorNum = 2;
    }
    else if (volName == "tracker_3")
    {
      detectorNum = 3;
    }

    // collect energy deposited in this step
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edepStep);

    // add to book
    G4bool entering = false;
    G4Track *track = step->GetTrack();

    if (track->GetVolume())
      volName = track->GetVolume()->GetName();
    if (track->GetNextVolume())
      nextVolName = track->GetNextVolume()->GetName();

    // Entering Detector
    if (volName != "tracker_1" && nextVolName == "tracker_1")
    {
      entering = true;
      // analysis->Update(track->GetKineticEnergy(),G4Threading::G4GetThreadId());
    }

    // Do the analysis related to this step
    analysis->analyseStepping(*track, entering, false, detectorNum);
    //analysis->Score(track->GetKineticEnergy() / keV, track->GetPosition() / mm);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
