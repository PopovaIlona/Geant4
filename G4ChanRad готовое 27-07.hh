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
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// base class for 'fast' parametrisation model describing Channeling Radiation
// created in some G4Envelope. 
// 
// History:
// 23.07.2021 A.A. Savchenko created it from XTR module

//

#ifndef G4ChanRad_h
#define G4ChanRad_h 1

#include <complex>
#include "globals.hh"
#include "Randomize.hh"

#include "G4LogicalVolume.hh"

#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Gamma.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VContinuousProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4Integrator.hh"
#include "G4ParticleChange.hh"


class G4VParticleChange;
class G4PhysicsFreeVector;
class G4PhysicsLinearVector;

class G4ChanRad : public G4VDiscreteProcess  // G4VContinuousProcess
{
public:

  explicit G4ChanRad (G4LogicalVolume *anEnvelope,G4Material*,
                    G4Material*, G4double,G4double,G4int,
                    const G4String & processName = "G4ChanRad",
                    G4ProcessType type = fElectromagnetic);
  virtual  ~G4ChanRad ();

 
  // These virtual has to be implemented in inherited particular crystals
  
 virtual  G4double SpectralCR(G4double energy);
  
 virtual G4double AngleCR(G4double varAngle);
 
  virtual G4bool IsApplicable(const G4ParticleDefinition&) override;

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
				   const G4Step&  aStep) override;

  virtual G4double GetMeanFreePath(const G4Track& aTrack,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;
  void BuildEnergyTable() ;
 // void BuildAngleForEnergyBank() ;

  void BuildTable(){} ;
 // void BuildAngleTable() ;
 // void BuildGlobalAngleTable() ;

  void GetNumberOfPhotons() ;  

  // Auxiliary functions 


  G4double  GetCRrandomEnergy( G4double scaledTkin, G4int iTkin );		//random choosing of CR photon energy
  G4double  GetCRenergy( G4int iPlace, G4double position, G4int iTransfer  );

  // G4double  GetRandomAngle( G4double energyXTR, G4int iTkin );
  // G4double  GetAngleXTR(G4int iCR,G4double position,G4int iAngle);

  G4double  GetGamma()   {return fGamma;}; 
  G4double  GetEnergy()  {return fEnergy;};                
  G4double  GetVarAngle(){return fVarAngle;};
               
  void SetGamma(G4double gamma)      {fGamma    = gamma;}; 
  void SetEnergy(G4double energy)    {fEnergy   = energy;};                
  void SetVarAngle(G4double varAngle){fVarAngle = varAngle;};               
  void SetAngleRadDistr(G4bool pAngleRadDistr){fAngleRadDistr=pAngleRadDistr;};               
                

  G4PhysicsLogVector* GetProtonVector(){ return fProtonEnergyVector;};
  G4int GetTotBin(){return fTotBin;};           
  G4PhysicsFreeVector* GetAngleVector(G4double energy, G4int n);

protected:

  G4ParticleDefinition* fPtrGamma ;    // pointer to TR photon

  G4double* fGammaCutInKineticEnergy ; // TR photon cut in energy array

  G4double         fGammaTkinCut ;     // Tkin cut of TR photon in current mat.
  G4LogicalVolume* fEnvelope ;
  G4PhysicsTable*  fAngleDistrTable ;
  G4PhysicsTable*  fEnergyDistrTable ;

  G4PhysicsLogVector* fProtonEnergyVector ;
  G4PhysicsLinearVector* fCREnergyVector ;

  G4double fTheMinEnCR;            //   min CR energy
  G4double fTheMaxEnCR;            //   max CR energy
  
  G4double fMinEnergyCR;               //  min CR energy in material
  G4double fMaxEnergyCR;               //  max CR energy in material
  
  // G4double fTheMaxAngle;               //  max theta of TR quanta
  // G4double fTheMinAngle;               //  max theta of TR quanta
  // G4double fMaxThetaTR;                //  max theta of TR quanta
  
  G4int    fBinCR;                     //  number of bins in CR vectors

  G4double fMinProtonTkin;             // min Tkin of proton in tables
  G4double fMaxProtonTkin;             // max Tkin of proton in tables
  G4int    fTotBin;                    // number of bins in log scale
  
  
  G4double fGamma;                     // current Lorentz factor
  G4double fEnergy;                    // energy and
  G4double fVarAngle;                  // angle squared
  G4double fLambda;



  G4bool   fExitFlux;					//flags for flux and angular distribution calculations
  G4bool   fAngleRadDistr;
  
  
  G4double fSigma1; 
  G4double fSigma2;                    //test values



  G4double fTotalDist;
 


  G4ParticleChange fParticleChange;

  G4PhysicsTable*                    fAngleForEnergyTable;
  std::vector<G4PhysicsTable*>       fAngleBank;

private:

  // copy constructor and hide assignment operator
  G4ChanRad(G4ChanRad &) = delete;
  G4ChanRad & operator=(const G4ChanRad &right) = delete;

};

#endif
