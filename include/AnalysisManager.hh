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
/// \file B1/include/AnalysisManager.hh
/// \brief Definition of the B1::AnalysisManager class

#ifndef B1AnalysisManager_h
#define B1AnalysisManager_h 1

#include "globals.hh"
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include "G4ios.hh"
#include "G4AnalysisManager.hh"

class AnalysisManager {

public:
    // The analysis class is designed to be a singleton (i.e. only one instance can exist).
    // A member function called Instance is defined, which allows the user to get
    // a pointer to the existing instance or to create it, if it does not yet exist.
    static auto Instance() -> AnalysisManager*;

    // The analysis class instance can be deleted by calling the Destroy method.
    // (NOTE: The class destructor is protected, and can thus not be called directly)
    static void Destroy();
    void finish(G4bool isMaster);

    // Member function used to score the total energy deposit
    void ScoreTotalEnergy(G4double totalDepositedEnergy);

    // Member function used to dump hits into csv file
    void Score(G4double depositedEnergy, G4ThreeVector position);

    void analyseStepping(const G4Track& track, G4bool entering);
    void book(G4bool isMaster);
    void Update(G4double energy,G4int threadID);
    void bookScore(G4double energy, G4ThreeVector position, G4int threadID);

protected:
    // Constructor (protected)
    explicit AnalysisManager();

    // Destructor (protected)
    virtual ~AnalysisManager();

    // Prevent copying
    AnalysisManager(const AnalysisManager& only);
    
    auto operator=(const AnalysisManager& only) -> const AnalysisManager&;

private:    
    static AnalysisManager* instance; // The static instance of the AnalysisManager class
    G4String histFileName;
    G4String asciiFileName;

  // Quantities for the ntuple
  G4double eKin;
  G4double x;
  G4double y;
  G4double z;
  G4double dirX;
  G4double dirY;
  G4double dirZ;

    std::ofstream csvFile;

    //global counters: log separately for each thread (or sequential)
    std::map<G4int,G4int> *nEnteringTracks;
    std::map<G4int,G4double> *totEnteringEnergy;
};
#endif // ANALYSISMANAGER_HH
