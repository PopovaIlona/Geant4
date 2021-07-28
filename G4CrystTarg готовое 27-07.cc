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

#include <complex>

#include "G4CrystTarg.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Integrator.hh"
#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4CrystTarg::G4CrystTarg(G4LogicalVolume *anEnvelope,
					 G4Material* foilMat,G4Material* gasMat, 
                                         G4double a, G4double b, G4int n,
                                         const G4String& processName) :
  G4ChanRad(anEnvelope,foilMat,gasMat,a,b,n,processName)
{
  G4cout<<"Channeling Radiation EM process is called"<<G4endl;
G4double t = SpectralCR(2);
  //  BuildTable();
}

///////////////////////////////////////////////////////////////////////////

G4CrystTarg::~G4CrystTarg()
{
  ;
}

///////////////////////////////////////////////////////////////////////////
//
//
 G4double G4CrystTarg::SpectralCR(G4double energy) 
 {
	G4double result=0.02;
	result+=energy-result;
	return result;
 }
  G4double G4CrystTarg::AngleCR(G4double varAngle)
  {
	  G4double result=0.05;
	result+=varAngle-result;
	
	return result;
  }



//
//
////////////////////////////////////////////////////////////////////////////








