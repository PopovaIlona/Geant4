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
//23.07.2021 Created by A.A. Savchenko

#include "G4ChanRad.hh"

#include "G4Timer.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Poisson.hh"
#include "G4MaterialTable.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "G4VSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"


#include "G4PhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4EmProcessSubType.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4ChanRad::G4ChanRad(G4LogicalVolume *anEnvelope,
				   G4Material* foilMat,G4Material* gasMat,
				   G4double a, G4double b,
				   G4int n,const G4String& processName,
				   G4ProcessType type) :
  G4VDiscreteProcess(processName, type),
  fGammaCutInKineticEnergy(nullptr),
  fGammaTkinCut(0.0),
 // fAngleDistrTable(nullptr),
  fEnergyDistrTable(nullptr)//,
 // fAngleForEnergyTable(nullptr)
{
  verboseLevel = 2;
  SetProcessSubType(fChanRad);
  fLambda = fTotalDist = 0.0;

  // Initialization of local constants
  fTheMinEnCR = 1.0*keV;  		//energy vector min value
  fTheMaxEnCR = 100.0*keV; 		//energy vector max value

  fBinCR          = 200; 		//energy vector bins

  fMinProtonTkin  = 1.0*GeV;
  fMaxProtonTkin  = 1.0*TeV;
  fTotBin         =  2;        	// Lorentz factor bins

  // Proton energy vector initialization
  fProtonEnergyVector = new G4PhysicsLogVector(fMinProtonTkin,    //vector in Log scale
					       fMaxProtonTkin,
					       fTotBin  );

  fCREnergyVector = new G4PhysicsLinearVector(fTheMinEnCR,		//vector in Linear scale, for angular distribution construction
					    fTheMaxEnCR,
					    fBinCR  );

    fEnvelope  = anEnvelope;
	
	
	fSigma1 = a; // test var
	fSigma2 = b; // test var
  
  
  //show here main properties of crystal
  
  // if(verboseLevel > 0)
    
  // if(fCrystThick == 0)
  // {
    // G4Exception("G4ChanRad::G4ChanRad()","VXTRELoss01",
    // FatalException,"No crystal");
  // }
  
  
  // default is CR dEdx, not flux after radiator
  
  fExitFlux      = false;
  fLambda = DBL_MAX;
  
  // Mean thicknesses of plates and gas gaps

  fTotalDist  = 1*mm; // test 1 mm //  size of crystal - thickness
  
  
  if(verboseLevel > 0)
    G4cout<<"total crystal thickness = "<<fTotalDist/cm<<" cm"<<G4endl;

  pParticleChange = &fParticleChange;
}

///////////////////////////////////////////////////////////////////////////

G4ChanRad::~G4ChanRad()
{
  delete fProtonEnergyVector;
  delete fCREnergyVector;
  if(fEnergyDistrTable) { 
    fEnergyDistrTable->clearAndDestroy(); 
    delete fEnergyDistrTable;
  }
  // if(fAngleRadDistr) {
    // fAngleDistrTable->clearAndDestroy();
    // delete fAngleDistrTable;
  // }
  // if(fAngleForEnergyTable) {
    // fAngleForEnergyTable->clearAndDestroy();
    // delete fAngleForEnergyTable;
  // }
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns condition for application of the model depending on particle type


G4bool G4ChanRad::IsApplicable(const G4ParticleDefinition& particle)
{
  return  ( particle.GetPDGCharge() != 0.0 );
}

/////////////////////////////////////////////////////////////////////////////////
//
// Calculate step size for XTR process inside raaditor

G4double G4ChanRad::GetMeanFreePath(const G4Track& aTrack,
					   G4double, // previousStepSize,
					   G4ForceCondition* condition)
{
  G4int iTkin, iPlace;
  G4double lambda, sigma, kinEnergy, mass, gamma;
  G4double charge, chargeSq, massRatio, TkinScaled;
  G4double E1,E2,W,W1,W2;

  *condition = NotForced;
  
  if( aTrack.GetVolume()->GetLogicalVolume() != fEnvelope ) lambda = DBL_MAX;
  else
  {
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    kinEnergy = aParticle->GetKineticEnergy();
    mass      = aParticle->GetDefinition()->GetPDGMass();
    gamma     = 1.0 + kinEnergy/mass;
    if(verboseLevel > 1)
    {
      G4cout<<" gamma = "<<gamma<<";   fGamma = "<<fGamma<<G4endl;
    }

    if ( std::fabs( gamma - fGamma ) < 0.05*gamma ) lambda = fLambda;
    else
    {
      charge = aParticle->GetDefinition()->GetPDGCharge();
      chargeSq  = charge*charge;
      massRatio = proton_mass_c2/mass;
      TkinScaled = kinEnergy*massRatio;

      for(iTkin = 0; iTkin < fTotBin; iTkin++)
      {
        if( TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))  break;    
      }
      iPlace = iTkin - 1;
	  
	  //G4cout << "( (*(*fEnergyDistrTable)(iPlace  ))(0)" <<  (*(*fEnergyDistrTable)(iPlace + 1 ))(0) << G4endl;
	  

      if(iTkin == 0) 
	  {
		  lambda = DBL_MAX;
		  G4cout << "Tkin is too small" << G4endl; // Tkin is too small, neglect of CR photon generation
	  }
		  
      else          // general case: Tkin between two vectors of the material
      {
        if(iTkin == fTotBin) 
        {
          sigma = (*(*fEnergyDistrTable)(iPlace))(0)*chargeSq;
        }
        else
        {
          E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1); 
          E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin);
          W = 1.0/(E2 - E1);
          W1 = (E2 - TkinScaled)*W;
          W2 = (TkinScaled - E1)*W;
          sigma = ( (*(*fEnergyDistrTable)(iPlace  ))(0)*W1 +
                (*(*fEnergyDistrTable)(iPlace+1))(0)*W2   )*chargeSq;
      
		
        }
        if (sigma < DBL_MIN) lambda = DBL_MAX;
        else                 lambda = 1./sigma; 
        fLambda = lambda;
        fGamma  = gamma;   
       // if(verboseLevel > 0)
        //{
	  G4cout<<" lambda = "<<lambda/mm<<" mm"<<G4endl;
       // }
      }
    }
  }  
  return lambda;
}

//////////////////////////////////////////////////////////////////////////
//
// Interface for build table from physics list

void G4ChanRad::BuildPhysicsTable(const G4ParticleDefinition& pd)
{
  if(pd.GetPDGCharge()  == 0.) 
  {
    G4Exception("G4ChanRad::BuildPhysicsTable", "Notification", JustWarning,
                 "CR initialisation for neutral particle ?!" );   
  }
  BuildEnergyTable();

  // if ( fAngleRadDistr ) 
  // {
    // if(verboseLevel > 0)
    // {
      // G4cout<<"Build angle for energy distribution according the current radiator"
	    // <<G4endl;
    // }
    // BuildAngleForEnergyBank();
  // }
}

 G4double G4ChanRad::SpectralCR(G4double energy) 
 {
	G4double result; //function overridden in CrystTarg
	return result;
 }
  G4double G4ChanRad::AngleCR(G4double varAngle)
  {
	  
	  G4double  result; //function overridden in CrystTarg
	return  result;
 }
//////////////////////////////////////////////////////////////////////////
//
// Build integral energy distribution of XTR photons

void G4ChanRad::BuildEnergyTable()
{
  G4int iTkin, iCR, iPlace;
  G4double energySum = 0.0;

  fEnergyDistrTable = new G4PhysicsTable(fTotBin);
 // if(fAngleRadDistr) fAngleDistrTable = new G4PhysicsTable(fTotBin);

  fGammaTkinCut = 0.0;
  
  // setting of min/max CR energies 
                               fMinEnergyCR = fTheMinEnCR;
	
                               fMaxEnergyCR = fTheMaxEnCR;
    
 // G4Integrator<G4ChanRad,G4double(G4ChanRad::*)(G4double)> integral;

  G4cout.precision(4);
  G4Timer timer;
  timer.Start();

 // if(verboseLevel > 0) 
 // {
    G4cout<<G4endl;
    G4cout<<"Lorentz Factor"<<"\t"<<"ChanRad photon number"<<G4endl;
    G4cout<<G4endl;
 // }
  for( iTkin = 0; iTkin < fTotBin; iTkin++ )      // Lorentz factor loop
  {
    G4PhysicsLinearVector* energyVector = new G4PhysicsLinearVector( fMinEnergyCR,
							       fMaxEnergyCR,
							       fBinCR  );

    fGamma = 1.0 + (fProtonEnergyVector->
		    GetLowEdgeEnergy(iTkin)/proton_mass_c2);
      
    energySum = 0.0;

    //energyVector->PutValue(fBinCR-1,energySum);

    for( iCR = fBinCR-1; iCR >=0; iCR-- )
    {
	// Legendre96 or Legendre10
	  G4double EE  = energyVector->GetLowEdgeEnergy(iCR)*1000;//+(energyVector->GetLowEdgeEnergy(iCR+1)-energyVector->GetLowEdgeEnergy(iCR))/2*1000;
      energySum += exp(-pow((EE-fSigma1),2)/fSigma2/fSigma2/2)/fSigma2/sqrt(2*pi);
	  // radiatorCof*fCofTR*integral.Legendre10(
       	     // this,&G4ChanRad::SpectralXTRdEdx,
	           // energyVector->GetLowEdgeEnergy(iCR),
	           // energyVector->GetLowEdgeEnergy(iCR+1) ); 
      	   
      G4cout<< "EE   " << EE <<
	   "    energySum   " << energySum <<G4endl;
	   
      energyVector->PutValue(iCR,energySum); ///fTotalDist
    }
    iPlace = iTkin;
    fEnergyDistrTable->insertAt(iPlace,energyVector);

    // if(verboseLevel > 0)
    // {
	G4cout
	  // <<iTkin<<"\t"
	  //   <<"fGamma = "
	  <<fGamma<<"\t"  //  <<"  fMaxThetaTR = "<<fMaxThetaTR
	  //  <<"sumN = "
	  <<energySum      // <<"; sumA = "<<angleSum
	  <<G4endl;
    // }
  }     
  timer.Stop();
  G4cout.precision(6);
  if(verboseLevel > 0) 
  {
    G4cout<<G4endl;
    G4cout<<"total time for build CR energy tables = "
	  <<timer.GetUserElapsed()<<" s"<<G4endl;
  }
  fGamma = 0.;
  return;
}

//////////////////////////////////////////////////////////////////////////
//
// Bank of angle distributions for given energies (slow!)

// void G4ChanRad::BuildAngleForEnergyBank()
// {
  // if( this->GetProcessName() == "TranspRegXTRadiator" || 
      // this->GetProcessName() == "TranspRegXTRmodel"   || 
      // this->GetProcessName() == "RegularXTRadiator"   || 
      // this->GetProcessName() == "RegularXTRmodel"       )
  // {
    // BuildAngleTable();
    // return;
  // }
  // G4int i, iTkin, iCR;
  // G4double angleSum  = 0.0;


  // fGammaTkinCut = 0.0;
  
  // // setting of min/max TR energies 
  
  // if(fGammaTkinCut > fTheMinEnCR)  fMinEnergyCR = fGammaTkinCut;
  // else                                 fMinEnergyCR = fTheMinEnCR;
	
  // if(fGammaTkinCut > fTheMaxEnCR) fMaxEnergyCR = 2.0*fGammaTkinCut;  
  // else                                fMaxEnergyCR = fTheMaxEnCR;

  // G4PhysicsLogVector* energyVector = new G4PhysicsLogVector( fMinEnergyCR,
							       // fMaxEnergyCR,
							       // fBinCR  );

  // G4Integrator<G4ChanRad,G4double(G4ChanRad::*)(G4double)> integral;

  // G4cout.precision(4);
  // G4Timer timer;
  // timer.Start();

  // for( iTkin = 0; iTkin < fTotBin; iTkin++ )      // Lorentz factor loop
  // {

    // fGamma = 1.0 + (fProtonEnergyVector->
		    // GetLowEdgeEnergy(iTkin)/proton_mass_c2);

    // fMaxThetaTR = 25*2500.0/(fGamma*fGamma) ;  // theta^2

    // fTheMinAngle = 1.0e-3; // was 5.e-6, e-6 !!!, e-5, e-4
 
    // if(      fMaxThetaTR > fTheMaxAngle )  fMaxThetaTR = fTheMaxAngle; 
    // else if( fMaxThetaTR < fTheMinAngle )  fMaxThetaTR = fTheMinAngle;
      
    // fAngleForEnergyTable = new G4PhysicsTable(fBinCR);

    // for( iCR = 0; iCR < fBinCR; iCR++ )
    // {
      // angleSum     = 0.0;
      // fEnergy      = energyVector->GetLowEdgeEnergy(iCR);    
      // G4PhysicsLinearVector* angleVector = new G4PhysicsLinearVector(0.0,
								   // fMaxThetaTR,
								   // fBinCR  );
    
      // angleVector ->PutValue(fBinCR - 1, angleSum);

      // for( i = fBinCR - 2; i >= 0; i-- )
      // {
	  // // Legendre96 or Legendre10

          // angleSum  += integral.Legendre10(
	           // this,&G4ChanRad::SpectralAngleXTRdEdx,
	           // angleVector->GetLowEdgeEnergy(i),
	           // angleVector->GetLowEdgeEnergy(i+1) );

          // angleVector ->PutValue(i, angleSum);
      // }
      // fAngleForEnergyTable->insertAt(iCR, angleVector);
    // }
    // fAngleBank.push_back(fAngleForEnergyTable); 
  // }    
  // timer.Stop();
  // G4cout.precision(6);
  // if(verboseLevel > 0) 
  // {
    // G4cout<<G4endl;
    // G4cout<<"total time for build X-ray TR angle for energy loss tables = "
	  // <<timer.GetUserElapsed()<<" s"<<G4endl;
  // }
  // fGamma = 0.;
  // delete energyVector;
// }

////////////////////////////////////////////////////////////////////////
//
// Build XTR angular distribution at given energy based on the model 
// of transparent regular radiator

// void G4ChanRad::BuildAngleTable()
// {
  // G4int iTkin, iCR;
  // G4double  energy;

  // fGammaTkinCut = 0.0;
                              
  // // setting of min/max TR energies 
  
  // if(fGammaTkinCut > fTheMinEnCR)  fMinEnergyCR = fGammaTkinCut;
  // else                                 fMinEnergyCR = fTheMinEnCR;
	
  // if(fGammaTkinCut > fTheMaxEnCR) fMaxEnergyCR = 2.0*fGammaTkinCut;  
  // else                                fMaxEnergyCR = fTheMaxEnCR;

  // G4cout.precision(4);
  // G4Timer timer;
  // timer.Start();
  // if(verboseLevel > 0) 
  // {
    // G4cout<<G4endl;
    // G4cout<<"Lorentz Factor"<<"\t"<<"XTR photon number"<<G4endl;
    // G4cout<<G4endl;
  // }
  // for( iTkin = 0; iTkin < fTotBin; iTkin++ )      // Lorentz factor loop
  // {
    
    // fGamma = 1.0 + (fProtonEnergyVector->
                            // GetLowEdgeEnergy(iTkin)/proton_mass_c2);

    // fMaxThetaTR = 25*2500.0/(fGamma*fGamma);  // theta^2

    // fTheMinAngle = 1.0e-3; // was 5.e-6, e-6 !!!, e-5, e-4
 
    // if( fMaxThetaTR > fTheMaxAngle )    fMaxThetaTR = fTheMaxAngle; 
    // else
    // {
       // if( fMaxThetaTR < fTheMinAngle )  fMaxThetaTR = fTheMinAngle;
    // }

    // fAngleForEnergyTable = new G4PhysicsTable(fBinCR);

    // for( iCR = 0; iCR < fBinCR; iCR++ )
    // {
      // // energy = fMinEnergyCR*(iCR+1);

      // energy = fCREnergyVector->GetLowEdgeEnergy(iCR);

      // G4PhysicsFreeVector* angleVector = GetAngleVector(energy,fBinCR);
      // // G4cout<<G4endl;

      // fAngleForEnergyTable->insertAt(iCR,angleVector);
    // }
    
    // fAngleBank.push_back(fAngleForEnergyTable); 
  // }     
  // timer.Stop();
  // G4cout.precision(6);
  // if(verboseLevel > 0) 
  // {
    // G4cout<<G4endl;
    // G4cout<<"total time for build XTR angle for given energy tables = "
	  // <<timer.GetUserElapsed()<<" s"<<G4endl;
  // }
  // fGamma = 0.;
  
  // return;
// } 

/////////////////////////////////////////////////////////////////////////
//
// Vector of angles and angle integral distributions

// G4PhysicsFreeVector* G4ChanRad::GetAngleVector(G4double energy, G4int n)
// {
  // G4double theta=0., result, tmp=0., cof1, cof2, cofMin, cofPHC, angleSum  = 0.;
  // G4int iTheta, k, /*kMax,*/ kMin;

  // G4PhysicsFreeVector* angleVector = new G4PhysicsFreeVector(n);
  
  // cofPHC  = 4.*pi*hbarc;
  // tmp     = (fSigma1 - fSigma2)/cofPHC/energy;
  // cof1    = fPlateThick*tmp;
  // cof2    = fGasThick*tmp;

  // cofMin  =  energy*(fPlateThick + fGasThick)/fGamma/fGamma;
  // cofMin += (fPlateThick*fSigma1 + fGasThick*fSigma2)/energy;
  // cofMin /= cofPHC;

  // kMin = G4int(cofMin);
  // if (cofMin > kMin) kMin++;

  // //kMax = kMin + fBinCR -1;

  // if(verboseLevel > 2)
  // {
    // G4cout<<"n-1 = "<<n-1<<"; theta = "
          // <<std::sqrt(fMaxThetaTR)*fGamma<<"; tmp = "
          // <<0. 
          // <<";    angleSum = "<<angleSum<<G4endl; 
  // }
  // //  angleVector->PutValue(n-1,fMaxThetaTR, angleSum);

  // for( iTheta = n - 1; iTheta >= 1; iTheta-- )
  // {

    // k = iTheta - 1 + kMin;

    // tmp    = pi*fPlateThick*(k + cof2)/(fPlateThick + fGasThick);

    // result = (k - cof1)*(k - cof1)*(k + cof2)*(k + cof2);

    // tmp = std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result;

    // if( k == kMin && kMin == G4int(cofMin) )
    // {
      // angleSum   += 0.5*tmp; // 0.5*std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result;
    // }
    // else if(iTheta == n-1);
    // else
    // {
      // angleSum   += tmp; // std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result;
    // }
    // theta = std::abs(k-cofMin)*cofPHC/energy/(fPlateThick + fGasThick);

    // if(verboseLevel > 2)
    // {
      // G4cout<<"iTheta = "<<iTheta<<"; k = "<<k<<"; theta = "
            // <<std::sqrt(theta)*fGamma<<"; tmp = "
            // <<tmp // std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result
            // <<";    angleSum = "<<angleSum<<G4endl;
    // }
    // angleVector->PutValue( iTheta, theta, angleSum );       
  // }
  // if (theta > 0.)
  // {
    // angleSum += 0.5*tmp;
    // theta = 0.;
  // }
  // if(verboseLevel > 2)
  // {
    // G4cout<<"iTheta = "<<iTheta<<"; theta = "
          // <<std::sqrt(theta)*fGamma<<"; tmp = "
          // <<tmp 
          // <<";    angleSum = "<<angleSum<<G4endl;
  // }
  // angleVector->PutValue( iTheta, theta, angleSum );

  // return angleVector;
// }

////////////////////////////////////////////////////////////////////////
//
// Build XTR angular distribution based on the model of transparent regular radiator

// void G4ChanRad::BuildGlobalAngleTable()
// {
  // G4int iTkin, iCR, iPlace;
  // G4double radiatorCof = 1.0;           // for tuning of XTR yield
  // G4double angleSum;
  // fAngleDistrTable = new G4PhysicsTable(fTotBin);

  // fGammaTkinCut = 0.0;
  
  // // setting of min/max TR energies 
  
  // if(fGammaTkinCut > fTheMinEnCR)  fMinEnergyCR = fGammaTkinCut;
  // else                                 fMinEnergyCR = fTheMinEnCR;
	
  // if(fGammaTkinCut > fTheMaxEnCR) fMaxEnergyCR = 2.0*fGammaTkinCut;  
  // else                                fMaxEnergyCR = fTheMaxEnCR;

  // G4cout.precision(4);
  // G4Timer timer;
  // timer.Start();
  // if(verboseLevel > 0) {
    // G4cout<<G4endl;
    // G4cout<<"Lorentz Factor"<<"\t"<<"XTR photon number"<<G4endl;
    // G4cout<<G4endl;
  // }
  // for( iTkin = 0; iTkin < fTotBin; iTkin++ )      // Lorentz factor loop
  // {
    
    // fGamma = 1.0 + (fProtonEnergyVector->
                            // GetLowEdgeEnergy(iTkin)/proton_mass_c2);

    // fMaxThetaTR = 25.0/(fGamma*fGamma);  // theta^2

    // fTheMinAngle = 1.0e-3; // was 5.e-6, e-6 !!!, e-5, e-4
 
    // if( fMaxThetaTR > fTheMaxAngle )    fMaxThetaTR = fTheMaxAngle; 
    // else
    // {
       // if( fMaxThetaTR < fTheMinAngle )  fMaxThetaTR = fTheMinAngle;
    // }
    // G4PhysicsLinearVector* angleVector = new G4PhysicsLinearVector(0.0,
                                                               // fMaxThetaTR,
                                                               // fBinCR      );

    // angleSum  = 0.0;

    // G4Integrator<G4ChanRad,G4double(G4ChanRad::*)(G4double)> integral;

   
    // angleVector->PutValue(fBinCR-1,angleSum);

    // for( iCR = fBinCR - 2; iCR >= 0; iCR-- )
    // {

      // angleSum  += radiatorCof*fCofTR*integral.Legendre96(
	           // this,&G4ChanRad::AngleXTRdEdx,
	           // angleVector->GetLowEdgeEnergy(iCR),
	           // angleVector->GetLowEdgeEnergy(iCR+1) );

      // angleVector ->PutValue(iCR,angleSum);
    // }
    // if(verboseLevel > 1) {
      // G4cout
	// // <<iTkin<<"\t"
	// //   <<"fGamma = "
	// <<fGamma<<"\t"  //  <<"  fMaxThetaTR = "<<fMaxThetaTR
	// //  <<"sumN = "<<energySum      // <<"; sumA = "
	// <<angleSum
	// <<G4endl;
    // }
    // iPlace = iTkin;
    // fAngleDistrTable->insertAt(iPlace,angleVector);
  // }     
  // timer.Stop();
  // G4cout.precision(6);
  // if(verboseLevel > 0) {
    // G4cout<<G4endl;
    // G4cout<<"total time for build X-ray TR angle tables = "
	  // <<timer.GetUserElapsed()<<" s"<<G4endl;
  // }
  // fGamma = 0.;
  
  // return;
// } 


//////////////////////////////////////////////////////////////////////////////
//
// The main function which is responsible for the treatment of a particle passage
// trough G4Envelope with discrete generation of G4Gamma

G4VParticleChange* G4ChanRad::PostStepDoIt( const G4Track& aTrack, 
		                                  const G4Step&  aStep   )
{
  G4int iTkin /*, iPlace*/;
  G4double energyCR, theta,theta2, phi, dirX, dirY, dirZ;
 

  fParticleChange.Initialize(aTrack);

  if(verboseLevel > 1)
  {
    G4cout<<"Start of G4ChanRad::PostStepDoIt "<<G4endl;
    G4cout<<"name of current material =  "
          <<aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()<<G4endl;
  }
  if( aTrack.GetVolume()->GetLogicalVolume() != fEnvelope ) 
  {
    if(verboseLevel > 0)
    {
      G4cout<<"Go out from G4ChanRad::PostStepDoIt: wrong volume "<<G4endl;
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  else
  {
    G4StepPoint* pPostStepPoint        = aStep.GetPostStepPoint();
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   
    // Now we are ready to Generate one CR photon

    G4double kinEnergy = aParticle->GetKineticEnergy();
    G4double mass      = aParticle->GetDefinition()->GetPDGMass();
    G4double gamma      = 1.0 + kinEnergy/mass;

    if(verboseLevel > 1 )
    {
      G4cout<<"gamma = "<<gamma<<G4endl;
    }
    G4double         massRatio   = proton_mass_c2/mass;
    G4double          TkinScaled = kinEnergy*massRatio;
    G4ThreeVector      position  = pPostStepPoint->GetPosition();
    G4ParticleMomentum direction = aParticle->GetMomentumDirection();
    G4double           startTime = pPostStepPoint->GetGlobalTime();

    for( iTkin = 0; iTkin < fTotBin; iTkin++ )
    {
      if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))  break;    
    }
    //iPlace = iTkin - 1;

    if(iTkin == 0) // Tkin is too small, neglect of TR photon generation
    {
      if( verboseLevel > 0)
      {
        G4cout<<"Go out from G4ChanRad::PostStepDoIt:iTkin = "<<iTkin<<G4endl;
      }
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    } 
    else          // general case: Tkin between two vectors of the material
    {
      fParticleChange.SetNumberOfSecondaries(1);

      energyCR = GetCRrandomEnergy(TkinScaled,iTkin);

      if( verboseLevel > 1)
      {
	G4cout<<"energyCR = "<<energyCR/keV<<" keV"<<G4endl;
      }
      // if ( fAngleRadDistr )
      // {
        // // theta = std::fabs(G4RandGauss::shoot(0.0,pi/gamma));

        // theta2 = GetRandomAngle(energyCR,iTkin);

        // if( theta2 > 0.) theta = std::sqrt(theta2);
        // else             theta = 0.;                // theta2;
      // }
      // else 
		  theta = std::fabs(G4RandGauss::shoot(0.0,pi/gamma));

      // theta = 0.;  // check no spread

      if( theta >= 0.1 ) theta = 0.1;

      // G4cout<<" : theta = "<<theta<<endl;

      phi = twopi*G4UniformRand();

      dirX = std::sin(theta)*std::cos(phi);
      dirY = std::sin(theta)*std::sin(phi);
      dirZ = std::cos(theta);

      G4ThreeVector directionTR(dirX,dirY,dirZ);
      directionTR.rotateUz(direction);
      directionTR.unit();

      G4DynamicParticle* aPhotonTR = new G4DynamicParticle(G4Gamma::Gamma(),
                                                           directionTR, energyCR);

      // A XTR photon is set on the particle track inside the radiator 
      // and is moved to the G4Envelope surface for standard X-ray TR models
      // only. The case of fExitFlux=true

      if( fExitFlux )
      {
        const G4RotationMatrix* rotM = pPostStepPoint->GetTouchable()->GetRotation();
        G4ThreeVector transl = pPostStepPoint->GetTouchable()->GetTranslation();
        G4AffineTransform transform = G4AffineTransform(rotM,transl);
        transform.Invert();
        G4ThreeVector localP = transform.TransformPoint(position);
        G4ThreeVector localV = transform.TransformAxis(directionTR);

        G4double distance = fEnvelope->GetSolid()->DistanceToOut(localP, localV);
        if(verboseLevel > 1)
        {
          G4cout<<"distance to exit = "<<distance/mm<<" mm"<<G4endl;
        }
        position         += distance*directionTR;
        startTime        += distance/c_light;
      }
      G4Track* aSecondaryTrack = new G4Track( aPhotonTR, 
		                                startTime, position );
      aSecondaryTrack->SetTouchableHandle(
                         aStep.GetPostStepPoint()->GetTouchableHandle());
      aSecondaryTrack->SetParentID( aTrack.GetTrackID() );

      fParticleChange.AddSecondary(aSecondaryTrack);
      fParticleChange.ProposeEnergy(kinEnergy);     
    }
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}






////////////////////////////////////////////////////////////////////////
//
// Computes matrix of Sandia photo absorption cross section coefficients for
// gas material

// void G4ChanRad::ComputeGasPhotoAbsCof() 
// {
  // const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  // const G4Material* mat = (*theMaterialTable)[fMatIndex2];
  // fGasPhotoAbsCof = mat->GetSandiaTable();
  // return;
// }

// //////////////////////////////////////////////////////////////////////
// //
// // Returns the value of linear photo absorption coefficient (in reciprocal 
// // length) for gas

// G4double G4ChanRad::GetGasLinearPhotoAbs(G4double omega) 
// {
  // G4double omega2, omega3, omega4; 

  // omega2 = omega*omega;
  // omega3 = omega2*omega;
  // omega4 = omega2*omega2;

  // const G4double* SandiaCof = fGasPhotoAbsCof->GetSandiaCofForMaterial(omega);
  // G4double cross = SandiaCof[0]/omega  + SandiaCof[1]/omega2 +
                   // SandiaCof[2]/omega3 + SandiaCof[3]/omega4;
  // return cross;

// }

// //////////////////////////////////////////////////////////////////////
// //
// // Calculates the product of linear cof by formation zone for plate. 
// // Omega is energy !!!

// G4double G4ChanRad::GetPlateZmuProduct( G4double omega ,
                                             // G4double gamma ,
                                             // G4double varAngle   ) 
// {
  // return GetPlateFormationZone(omega,gamma,varAngle)
    // * GetPlateLinearPhotoAbs(omega);
// }
// //////////////////////////////////////////////////////////////////////
// //
// // Calculates the product of linear cof by formation zone for plate. 
// // G4cout and output in file in some energy range.

// void G4ChanRad::GetPlateZmuProduct() 
// {
  // std::ofstream outPlate("plateZmu.dat", std::ios::out );
  // outPlate.setf( std::ios::scientific, std::ios::floatfield );

  // G4int i;
  // G4double omega, varAngle, gamma;
  // gamma = 10000.;
  // varAngle = 1/gamma/gamma;
  // if(verboseLevel > 0)
    // G4cout<<"energy, keV"<<"\t"<<"Zmu for plate"<<G4endl;
  // for(i=0;i<100;i++)
  // {
    // omega = (1.0 + i)*keV;
    // if(verboseLevel > 1)
      // G4cout<<omega/keV<<"\t"<<GetPlateZmuProduct(omega,gamma,varAngle)<<"\t";
    // if(verboseLevel > 0)
      // outPlate<<omega/keV<<"\t\t"<<GetPlateZmuProduct(omega,gamma,varAngle)<<G4endl;
  // }
  // return;
// }

//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof by formation zone for gas. 
// Omega is energy !!!

// G4double G4ChanRad::GetGasZmuProduct( G4double omega ,
                                             // G4double gamma ,
                                             // G4double varAngle   ) 
// {
  // return GetGasFormationZone(omega,gamma,varAngle)*GetGasLinearPhotoAbs(omega);
// }
// //////////////////////////////////////////////////////////////////////
// //
// // Calculates the product of linear cof byformation zone for gas. 
// // G4cout and output in file in some energy range.

// void G4ChanRad::GetGasZmuProduct() 
// {
  // std::ofstream outGas("gasZmu.dat", std::ios::out );
  // outGas.setf( std::ios::scientific, std::ios::floatfield );
  // G4int i;
  // G4double omega, varAngle, gamma;
  // gamma = 10000.;
  // varAngle = 1/gamma/gamma;
  // if(verboseLevel > 0)
    // G4cout<<"energy, keV"<<"\t"<<"Zmu for gas"<<G4endl;
  // for(i=0;i<100;i++)
  // {
    // omega = (1.0 + i)*keV;
    // if(verboseLevel > 1)
      // G4cout<<omega/keV<<"\t"<<GetGasZmuProduct(omega,gamma,varAngle)<<"\t";
    // if(verboseLevel > 0)
      // outGas<<omega/keV<<"\t\t"<<GetGasZmuProduct(omega,gamma,varAngle)<<G4endl;
  // }
  // return;
// }

// ////////////////////////////////////////////////////////////////////////
// //
// // Computes Compton cross section for plate material in 1/mm

// G4double G4ChanRad::GetPlateCompton(G4double omega) 
// {
  // G4int i, numberOfElements;
  // G4double xSection = 0., nowZ, sumZ = 0.;

  // const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  // numberOfElements = (*theMaterialTable)[fMatIndex1]->GetNumberOfElements();

  // for( i = 0; i < numberOfElements; i++ )
  // {
    // nowZ      = (*theMaterialTable)[fMatIndex1]->GetElement(i)->GetZ();
    // sumZ     += nowZ;
    // xSection += GetComptonPerAtom(omega,nowZ); // *nowZ;
  // }
  // xSection /= sumZ;
  // xSection *= (*theMaterialTable)[fMatIndex1]->GetElectronDensity();
  // return xSection;
// }


// ////////////////////////////////////////////////////////////////////////
// //
// // Computes Compton cross section for gas material in 1/mm

// G4double G4ChanRad::GetGasCompton(G4double omega) 
// {
  // G4int i, numberOfElements;
  // G4double xSection = 0., nowZ, sumZ = 0.;

  // const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  // numberOfElements = (*theMaterialTable)[fMatIndex2]->GetNumberOfElements();

  // for( i = 0; i < numberOfElements; i++ )
  // {
    // nowZ      = (*theMaterialTable)[fMatIndex2]->GetElement(i)->GetZ();
    // sumZ     += nowZ;
    // xSection += GetComptonPerAtom(omega,nowZ); // *nowZ;
  // }
  // xSection /= sumZ;
  // xSection *= (*theMaterialTable)[fMatIndex2]->GetElectronDensity();
  // return xSection;
// }

// ////////////////////////////////////////////////////////////////////////
// //
// // Computes Compton cross section per atom with Z electrons for gamma with
// // the energy GammaEnergy

// G4double G4ChanRad::GetComptonPerAtom(G4double GammaEnergy, G4double Z) 
// {
  // G4double CrossSection = 0.0;
  // if ( Z < 0.9999 )                 return CrossSection;
  // if ( GammaEnergy < 0.1*keV      ) return CrossSection;
  // if ( GammaEnergy > (100.*GeV/Z) ) return CrossSection;

  // static const G4double a = 20.0 , b = 230.0 , c = 440.0;

  // static const G4double
  // d1= 2.7965e-1*barn, d2=-1.8300e-1*barn, d3= 6.7527   *barn, d4=-1.9798e+1*barn,
  // e1= 1.9756e-5*barn, e2=-1.0205e-2*barn, e3=-7.3913e-2*barn, e4= 2.7079e-2*barn,
  // f1=-3.9178e-7*barn, f2= 6.8241e-5*barn, f3= 6.0480e-5*barn, f4= 3.0274e-4*barn;

  // G4double p1Z = Z*(d1 + e1*Z + f1*Z*Z), p2Z = Z*(d2 + e2*Z + f2*Z*Z),
           // p3Z = Z*(d3 + e3*Z + f3*Z*Z), p4Z = Z*(d4 + e4*Z + f4*Z*Z);

  // G4double T0  = 15.0*keV;
  // if (Z < 1.5) T0 = 40.0*keV;

  // G4double X   = std::max(GammaEnergy, T0) / electron_mass_c2;
  // CrossSection = p1Z*std::log(1.+2.*X)/X
               // + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);

  // //  modification for low energy. (special case for Hydrogen)

  // if (GammaEnergy < T0) 
  // {
    // G4double dT0 = 1.*keV;
    // X = (T0+dT0) / electron_mass_c2;
    // G4double sigma = p1Z*std::log(1.+2.*X)/X
                    // + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    // G4double   c1 = -T0*(sigma-CrossSection)/(CrossSection*dT0);
    // G4double   c2 = 0.150;
    // if (Z > 1.5) c2 = 0.375-0.0556*std::log(Z);
    // G4double    y = std::log(GammaEnergy/T0);
    // CrossSection *= std::exp(-y*(c1+c2*y));
  // }
  // //  G4cout << "e= " << GammaEnergy << " Z= " << Z << " cross= " << CrossSection << G4endl;
  // return CrossSection;  
// }

////////////////////////////////////////////////////////////////////////
//
// Check number of photons for a range of Lorentz factors from both energy 
// and angular tables

void G4ChanRad::GetNumberOfPhotons()
{
  G4int iTkin;
  G4double gamma, numberE;

  std::ofstream outEn("numberE.dat", std::ios::out );
  outEn.setf( std::ios::scientific, std::ios::floatfield );

  std::ofstream outAng("numberAng.dat", std::ios::out );
  outAng.setf( std::ios::scientific, std::ios::floatfield );

  for(iTkin=0;iTkin<fTotBin;iTkin++)      // Lorentz factor loop
  {
     gamma = 1.0 + (fProtonEnergyVector->
                            GetLowEdgeEnergy(iTkin)/proton_mass_c2);
     numberE = (*(*fEnergyDistrTable)(iTkin))(0);
     //  numberA = (*(*fAngleDistrTable)(iTkin))(0);
     if(verboseLevel > 1)
       G4cout<<gamma<<"\t\t"<<numberE<<"\t"    //  <<numberA
	     <<G4endl; 
     if(verboseLevel > 0)
       outEn<<gamma<<"\t\t"<<numberE<<G4endl; 
  }
  return;
}  

/////////////////////////////////////////////////////////////////////////
//
// Returns randon energy of a X-ray TR photon for given scaled kinetic energy
// of a charged particle

G4double G4ChanRad::GetCRrandomEnergy( G4double scaledTkin, G4int iTkin )
{
  G4int iTransfer, iPlace;
  G4double transfer = 0.0, position, E1, E2, W1, W2, W;

  iPlace = iTkin - 1;

  //  G4cout<<"iPlace = "<<iPlace<<endl;

  if(iTkin == fTotBin) // relativistic plato, try from left
  {
      position = (*(*fEnergyDistrTable)(iPlace))(0)*G4UniformRand();

      for(iTransfer=0;;iTransfer++)
      {
        if(position >= (*(*fEnergyDistrTable)(iPlace))(iTransfer)) break;
      }
      transfer = GetCRenergy(iPlace,position,iTransfer);
  }
  else
  {
    E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1); 
    E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin);
    W  = 1.0/(E2 - E1);
    W1 = (E2 - scaledTkin)*W;
    W2 = (scaledTkin - E1)*W;

    position =( (*(*fEnergyDistrTable)(iPlace))(0)*W1 + 
                    (*(*fEnergyDistrTable)(iPlace+1))(0)*W2 )*G4UniformRand();

        // G4cout<<position<<"\t";

    for(iTransfer=0;;iTransfer++)
    {
          if( position >=
          ( (*(*fEnergyDistrTable)(iPlace))(iTransfer)*W1 + 
            (*(*fEnergyDistrTable)(iPlace+1))(iTransfer)*W2) ) break;
    }
    transfer = GetCRenergy(iPlace,position,iTransfer);
    
  } 
  //  G4cout<<"XTR transfer = "<<transfer/keV<<" keV"<<endl; 
  if(transfer < 0.0 ) transfer = 0.0;
  return transfer;
}

////////////////////////////////////////////////////////////////////////
//
// Returns approximate position of X-ray photon energy during random sampling
// over integral energy distribution

G4double G4ChanRad::GetCRenergy( G4int    iPlace, 
                                       G4double  /* position */, 
                                       G4int    iTransfer )
{
  G4double x1, x2, y1, y2, result;

  if(iTransfer == 0)
  {
    result = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer);
  }  
  else
  {
    y1 = (*(*fEnergyDistrTable)(iPlace))(iTransfer-1);
    y2 = (*(*fEnergyDistrTable)(iPlace))(iTransfer);

    x1 = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1);
    x2 = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer);

    if ( x1 == x2 )    result = x2;
    else
    {
      if ( y1 == y2  ) result = x1 + (x2 - x1)*G4UniformRand();
      else
      {
        // result = x1 + (position - y1)*(x2 - x1)/(y2 - y1);
        result = x1 + (x2 - x1)*G4UniformRand();
      }
    }
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////
// //
// //  Get XTR photon angle at given energy and Tkin

// G4double G4ChanRad::GetRandomAngle( G4double energyXTR, G4int iTkin )
// {
  // G4int iCR, iAngle;
  // G4double position, angle;

  // if (iTkin == fTotBin) iTkin--;

  // fAngleForEnergyTable = fAngleBank[iTkin];

  // for( iCR = 0; iCR < fBinCR; iCR++ )
  // {
    // if( energyXTR < fCREnergyVector->GetLowEdgeEnergy(iCR) )  break;    
  // }
  // if (iCR == fBinCR) iCR--;
      
  // position = (*(*fAngleForEnergyTable)(iCR))(0)*G4UniformRand();

  // for( iAngle = 0;; iAngle++)
  // {
    // if( position >= (*(*fAngleForEnergyTable)(iCR))(iAngle) ) break;
  // }
  // angle = GetAngleXTR(iCR,position,iAngle);
  // return angle;
// }

// ////////////////////////////////////////////////////////////////////////
// //
// // Returns approximate position of X-ray photon angle at given energy during random sampling
// // over integral energy distribution

// G4double G4ChanRad::GetAngleXTR( G4int    iPlace, 
                                        // G4double position, 
                                        // G4int    iTransfer )
// {
  // G4double x1, x2, y1, y2, result;

  // if( iTransfer == 0 )
  // {
    // result = (*fAngleForEnergyTable)(iPlace)->GetLowEdgeEnergy(iTransfer);
  // }  
  // else
  // {
    // y1 = (*(*fAngleForEnergyTable)(iPlace))(iTransfer-1);
    // y2 = (*(*fAngleForEnergyTable)(iPlace))(iTransfer);

    // x1 = (*fAngleForEnergyTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1);
    // x2 = (*fAngleForEnergyTable)(iPlace)->GetLowEdgeEnergy(iTransfer);

    // if ( x1 == x2 )    result = x2;
    // else
    // {
      // if ( y1 == y2  ) result = x1 + (x2 - x1)*G4UniformRand();
      // else
      // {
        // result = x1 + (position - y1)*(x2 - x1)/(y2 - y1);
      // }
    // }
  // }
  // return result;
// }


//
//
///////////////////////////////////////////////////////////////////////

