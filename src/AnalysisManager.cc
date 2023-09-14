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

namespace
{
  // Mutex to acquire access to singleton instance check/creation
  G4Mutex instanceMutex = G4MUTEX_INITIALIZER;
  // Mutex to acquire accss to histograms creation/access
  // It is also used to control all operations related to histos
  // File writing and check analysis
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

auto AnalysisManager::Instance() -> AnalysisManager *
{
  G4AutoLock l(&instanceMutex);
  // A new instance of AnalysisManager is created, if it does not exist:
  if (instance == nullptr)
  {
    instance = new AnalysisManager();
  }

  // The instance of AnalysisManager is returned:
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::Destroy()
{
  // The AnalysisManager instance is deleted, if it exists:
  if (instance != nullptr)
  {
    delete instance;
    instance = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::Score(G4double depositedEnergy, G4ThreeVector position)
{
  // write data into the data file
  csvFile << depositedEnergy / keV << "," << position.x() / mm << "," << position.y() / mm << "," << position.z() / mm << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::ScoreTotalEnergy(G4double totalDepositedEnergy)
{
  // dataFile2 << totalDepositedEnergy << std::endl;
}

void AnalysisManager::book(G4bool isMaster)
{
  G4AutoLock l(&dataManipulationMutex);

  // Get/create analysis manager
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  man->SetDefaultFileType("root");

  // Open an output file: it is done in master and threads. The
  // printout is done only by the master, for tidyness
  // if (isMaster)
  //    G4cout << "Opening output file " << histFileName << " ... ";
  man->OpenFile(histFileName);
  man->SetFirstHistoId(1);
  // if (isMaster)
  //     G4cout << " done" << G4endl;

  // Book 1D histograms
  man->CreateH1("det0", "Energy deposits Detector 0", 100, 1., 100., "keV", "cts"); // 0
  man->CreateH1("det1", "Energy deposits Detector 1", 100, 1., 100., "keV", "cts"); // 1
  man->CreateH1("det2", "Energy deposits Detector 2", 100, 1., 100., "keV", "cts"); // 2
  man->CreateH1("det4", "Energy deposits Detector 3", 100, 1., 100., "keV", "cts"); // 3
  man->CreateH1("alldet", "Energy deposits all dets", 100, 1., 100.);              // 4
  man->CreateH1("origE", "Original Energy", 100, 1., 100., "keV", "cts");          // 5

  man->CreateH1("h3", "Energy in detector, all /keV", 100, 1., 100.);

  // Book 2D histograms (notice: the numbering is independent)
  man->CreateH2("det0_xy", "Detector 0 x-y", 100, -40., 40., 100, -40., 40., "mm", "mm");      // 0
  man->CreateH2("det1_xy", "Detector 1 x-y", 100, -40., 40., 100, -40., 40., "mm", "mm");      // 1
  man->CreateH2("det2_xy", "Detector 2 x-y", 100, -40., 40., 100, -40., 40., "mm", "mm");      // 2
  man->CreateH2("det3_xy", "Detector 3 x-y", 100, -40., 40., 100, -40., 40., "mm", "mm");      // 3
  man->CreateH2("alldet_xy", "All Detectors x-y", 100, -40., 40., 100, -40., 40., "mm", "mm"); // 4
  man->CreateH2("orig_xy", "Original x-y", 100, -40., 40., 100, -40., 40., "mm", "mm"); // 5

  man->CreateH2("d2", "x-y, entering detector /mm", 200, -50., 50., 200, -50., 50.);
  man->CreateH2("d3", "x-y, detector hit /mm", 200, -50., 50., 200, -50., 50.);

  // Book ntuples
  man->CreateNtuple("tree", "Track ntuple");
  man->CreateNtupleDColumn("energy");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleDColumn("dirx");
  man->CreateNtupleDColumn("diry");
  man->CreateNtupleDColumn("dirz");
  man->CreateNtupleDColumn("detectorNum");
  man->FinishNtuple();
}

void AnalysisManager::finish(G4bool isMaster)
{
  G4AutoLock l(&dataManipulationMutex);
  // Save histograms
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  man->Clear();

  if (!isMaster)
    return;
}

void AnalysisManager::bookScore(G4double energy, G4ThreeVector position, G4int threadID)
{
}

void AnalysisManager::Update(G4double energy, G4int threadID)
{
}

void AnalysisManager::analyseStepping(const G4Track &track, G4bool entering, G4bool inDetector, G4int detectorNum)
{
  G4AutoLock l(&dataManipulationMutex);
  eKin = track.GetKineticEnergy() / keV;
  G4ThreeVector pos = track.GetPosition() / mm;
  y = pos.y();
  z = pos.z();
  x = pos.x();
  G4ThreeVector dir = track.GetMomentumDirection();
  dirX = dir.x();
  dirY = dir.y();
  dirZ = dir.z();

  // get original energy and position
  eKin0 = track.GetVertexKineticEnergy() / keV;
  G4ThreeVector pos0 = track.GetVertexPosition() / mm;
  y0 = pos0.y();
  z0 = pos0.z();
  x0 = pos0.x();

  // Fill histograms
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->FillH1(detectorNum, eKin);
  man->FillH2(detectorNum, x, y);

  man->FillH1(4, eKin); // for all detectors
  man->FillH2(4, x, y);

  man->FillH1(5, eKin0); // original energy
  man->FillH2(5, x0, y0); // original position

  man->FillNtupleDColumn(0, eKin);
  man->FillNtupleDColumn(1, x);
  man->FillNtupleDColumn(2, y);
  man->FillNtupleDColumn(3, z);
  man->FillNtupleDColumn(4, dirX);
  man->FillNtupleDColumn(5, dirY);
  man->FillNtupleDColumn(6, dirZ);
  man->FillNtupleDColumn(7, detectorNum);
  man->AddNtupleRow();

  // Fill histograms and ntuple, tracks entering the detector
  if (entering)
  {
    // Fill and plot histograms
    man->FillH1(6, eKin);
  }
  if (inDetector)
  {
    // Fill and plot histograms
    //man->FillH1(3, eKin);
    //man->FillH2(3, x, y);
  }
}
