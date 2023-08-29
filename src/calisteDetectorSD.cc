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

#include "AnalysisManager.hh"
#include "calisteDetectorSD.hh"
#include "DetectorHit.hh"

#include <fstream>
#include <iostream>

#include "G4HCofThisEvent.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

calisteDetectorSD::calisteDetectorSD(G4String name) : G4VSensitiveDetector(name)
{
    constexpr auto TRACKER_COLLECTION_NAME = "TrackerCollection";
    collectionName.insert(TRACKER_COLLECTION_NAME);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

calisteDetectorSD::~calisteDetectorSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void calisteDetectorSD::Initialize(G4HCofThisEvent* collection)
{
    trackerCollection = new DetectorHitsCollection(SensitiveDetectorName, collectionName[0]);
    static G4int collectionIdentifier = -1;

    if (collectionIdentifier < 0) {
        collectionIdentifier = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    collection->AddHitsCollection(collectionIdentifier, trackerCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto calisteDetectorSD::ProcessHits(G4Step* step, G4TouchableHistory*) -> G4bool
{
    if (step->GetTrack()->GetDefinition() != G4Gamma::GammaDefinition()) {
        return false;
    }

    // auto depositedEnergy = step->GetTotalEnergyDeposit();
    auto depositedEnergy = step->GetPreStepPoint()->GetKineticEnergy();
    if (depositedEnergy == 0.) {
        return false;
    }

    auto *newHit = new DetectorHit();
    newHit->SetTrackIdentifier(step->GetTrack()->GetTrackID());
    // newHit->SetChamberNumber(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
    newHit->SetDepositedEnergy(depositedEnergy);
    newHit->SetPosition(step->GetPostStepPoint()->GetPosition());
    trackerCollection->insert(newHit);

    // newHit->Print();
    // newHit->Draw();

    // std::ofstream dataFile("ASCII");
    // newHit->PrintToFile(dataFile);
    // dataFile.close();

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void calisteDetectorSD::EndOfEvent(G4HCofThisEvent*)
{
    G4int numberOfHits = trackerCollection->entries();

    if (verboseLevel > 0) {
        G4cout << G4endl
            << "Hits collection: in this event there are " << numberOfHits
            << " hits in the tracker chambers: " << G4endl;

        for (G4int i = 0; i < numberOfHits; i++) {
            (*trackerCollection)[i]->Print();
        }
    }

    for (G4int i = 0; i < numberOfHits; i++) {
        (*trackerCollection)[i]->PrintToFile();
    };
}
