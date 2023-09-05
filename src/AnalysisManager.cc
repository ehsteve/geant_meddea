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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "AnalysisManager.hh"
#include <sstream>
#include "G4AnalysisManager.hh"
#include "G4AutoLock.hh"

AnalysisManager *AnalysisManager::instance{nullptr};

namespace { 
  //Mutex to acquire access to singleton instance check/creation
  G4Mutex instanceMutex = G4MUTEX_INITIALIZER;
  //Mutex to acquire accss to histograms creation/access
  //It is also used to control all operations related to histos 
  //File writing and check analysis
  G4Mutex dataManipulationMutex = G4MUTEX_INITIALIZER;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager::AnalysisManager()
{
    histFileName = "DetectorHists";
    csvFile.open("DetectorHitsList.csv"); // open the file
    (csvFile) << "Energy (keV), x (mm), y (mm), z (mm)" << G4endl;  
//    dataFile2.open("TotalEnergy.out"); // open the file

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager::~AnalysisManager()
{
    csvFile.close(); // close the file
//    dataFile2.close(); // close the file
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto AnalysisManager::Instance() -> AnalysisManager*
{
    G4AutoLock l(&instanceMutex);
    // A new instance of AnalysisManager is created, if it does not exist:
    if (instance == nullptr) {
        instance = new AnalysisManager();
    }

    // The instance of AnalysisManager is returned:
    return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::Destroy()
{
    // The AnalysisManager instance is deleted, if it exists:
    if (instance != nullptr) {
        delete instance;
        instance = nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::Score(G4double depositedEnergy, G4ThreeVector position)
{
    // write data into the data file
    csvFile << depositedEnergy/keV << "," << position.x()/mm << "," << position.y()/mm << "," << position.z()/mm << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::ScoreTotalEnergy(G4double totalDepositedEnergy)
{
    //dataFile2 << totalDepositedEnergy << std::endl;
}

void AnalysisManager::book(G4bool isMaster)
{
    G4AutoLock l(&dataManipulationMutex);

    // Get/create analysis manager
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->SetDefaultFileType("root");

    // Open an output file: it is done in master and threads. The 
    // printout is done only by the master, for tidyness
    //if (isMaster)
    //    G4cout << "Opening output file " << histFileName << " ... ";
    man->OpenFile(histFileName);
    man->SetFirstHistoId(1);
    //if (isMaster)
    //    G4cout << " done" << G4endl;

   // Book 1D histograms
  man->CreateH1("h1","Energy, all /keV",  100,1.,100.);
  man->CreateH1("h2","Energy, entering detector /keV", 100,1.,100.);
  
  man->CreateH1("h3","Energy in detector, all /keV",  100,1.,100.);

  // Book 2D histograms (notice: the numbering is independent)
  man->CreateH2("d1","x-y, all /mm", 100,-100.,100.,100,-100.,100.); 
  man->CreateH2("d2","x-y, entering detector /mm", 200,-50.,50.,200,-50.,50.);
  man->CreateH2("d3","x-y, detector hit /mm", 200,-50.,50.,200,-50.,50.);

  // Book ntuples
  man->CreateNtuple("tree", "Track ntuple");
  man->CreateNtupleDColumn("energy");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleDColumn("dirx");
  man->CreateNtupleDColumn("diry");
  man->CreateNtupleDColumn("dirz");
  man->FinishNtuple();
}

void AnalysisManager::finish(G4bool isMaster)
{
  G4AutoLock l(&dataManipulationMutex);
  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  man->Clear();

  if (!isMaster)
    return;
}

void AnalysisManager::bookScore(G4double energy, G4ThreeVector position, G4int threadID){
    
}

void AnalysisManager::Update(G4double energy,G4int threadID){

}

void AnalysisManager::analyseStepping(const G4Track& track, G4bool entering, G4bool inDetector){
  G4AutoLock l(&dataManipulationMutex);
  eKin = track.GetKineticEnergy()/keV;
  G4ThreeVector pos = track.GetPosition()/mm;
  y = pos.y();
  z = pos.z();
  x = pos.x();
  G4ThreeVector dir = track.GetMomentumDirection();
  dirX = dir.x();
  dirY = dir.y();
  dirZ = dir.z();

  // Fill histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(1,eKin);
  man->FillH2(1,x,y);
  
  // Fill histograms and ntuple, tracks entering the detector
  if (entering) {
    // Fill and plot histograms
    man->FillH1(2,eKin);
    man->FillH2(2,x,y);

    man->FillNtupleDColumn(0,eKin);
    man->FillNtupleDColumn(1,x);
    man->FillNtupleDColumn(2,y);
    man->FillNtupleDColumn(3,z);
    man->FillNtupleDColumn(4,dirX);
    man->FillNtupleDColumn(5,dirY);
    man->FillNtupleDColumn(6,dirZ);
    man->AddNtupleRow();
  }
  if (inDetector) {
    // Fill and plot histograms
    man->FillH1(3,eKin);
    man->FillH2(3,x,y);
  }

}
