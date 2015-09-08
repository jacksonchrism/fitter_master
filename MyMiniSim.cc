#include "MyMiniSim.hh"
#include <Randomize.hh>
#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <RAT/DB.hh>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom.h>
#include <RAT/Producer.hh>
#include <RAT/ProcBlock.hh>
#include <RAT/PMTCharge.hh>
#include <RAT/ChannelEfficiency.hh>
#include <RAT/PMTTransitTime.hh>
#include <RAT/PhotonTrajectory.hh>
#include <RAT/DS/MCPhoton.hh>
#include <RAT/TrackingAction.hh>

#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
using namespace CLHEP;

//class G4VPhysicalVolume;
//class G4VSolid;

double fCurrentp0;
double fCurrentp1;
double fPhotonWavelength;

namespace RAT {

MyMiniSim::MyMiniSim(double photonWavelength, double par0, double par1)
: MiniSim()
{


  std::cout << "making the minisim" << std::endl;
  //fix the random number seed
  CLHEP::HepRandom::setTheSeed(1234567);
  std::cout << "setting seed (again)" << std::endl;
//  fCurrentParam = parameterIterationValue;
  fCurrentp0 = par0;
  fCurrentp1 = par1;
  fPhotonWavelength = photonWavelength;
  fEventCtr = 0;
  const int nchannels = 90;
  const double maxangle = 90.;
  fIncidentPhotons = new TH1D("IncidentPhotons", "IncidentPhotons", nchannels, 0., maxangle);
  fSuccessfulPhotons = new TH1D("SuccessfulPhotons", "SuccessfulPhotons", nchannels, 0., maxangle);
  fReflectedPhotons = new TH1D("ReflectedPhotons", "ReflectedPhotons", nchannels, 0., maxangle);
  fSuccessfulZ = new TH1D("SuccessfulZ", "SuccessfulPhotons entry Z", 200, 0., 200);
  fSuccessfulRadius = new TH1D("SuccessfulRadius", "SuccessfulPhotons entry Radius", 200, 0., 200);

  fSuccessfulPhotons->SetDirectory(0);
  fIncidentPhotons->SetDirectory(0);
  //these are here for debugging
  fReflectedPhotons->SetDirectory(0);
  fSuccessfulZ->SetDirectory(0);
  fSuccessfulRadius->SetDirectory(0);

}

MyMiniSim::~MyMiniSim() {
  //delete fIncidentPhotons;
  // delete fSuccessfulPhotons;
  // delete fReflectedPhotons;
}

void MyMiniSim::ScaleProperly(TH1D* result)
{      
      const double normFactor = result->GetBinContent( 1 );
      const double normError = result->GetBinError( 1 );
      for( int iLoop = 1; iLoop <= result->GetNbinsX(); iLoop++ )
        {
          if( result->GetBinContent( iLoop ) == 0.0 )
            continue;
          const double binVal = result->GetBinContent( iLoop ) / normFactor;
          const double errVal = sqrt( pow( result->GetBinError( iLoop ) / result->GetBinContent( iLoop ), 2 ) + pow( normError / normFactor, 2 ) ) * binVal;

          result->SetBinContent( iLoop, binVal );
          result->SetBinError( iLoop, errVal );
        }
}


TH1D* MyMiniSim::GetSimAngResp(){

   //call this from the fitter to get minisim angular response
   fSuccessfulPhotons->Sumw2();
   fIncidentPhotons->Sumw2();
   TH1D* SimAngResponse = new TH1D("SimAngResponse", "SimAngResponse", 90, 0., 90.);
   SimAngResponse->SetDirectory(0);
   SimAngResponse->Divide(fSuccessfulPhotons,fIncidentPhotons);
   return SimAngResponse; 
}


TH1D* MyMiniSim::GetSuccessfulPhotons(){

   //call this from the fitter to get minisim angular response
   fSuccessfulPhotons->Sumw2();
   return fSuccessfulPhotons; 
}

TH1D* MyMiniSim::GetIncidentPhotons(){

   //call this from the fitter to get minisim angular response
   fIncidentPhotons->Sumw2();
   return fIncidentPhotons; 
}

TH1D* MyMiniSim::GetSuccessfulZ(){

   //call this from the fitter to get minisim angular response
   fSuccessfulZ->Sumw2();
   return fSuccessfulZ; 
}

TH1D* MyMiniSim::GetSuccessfulRadius(){

   //call this from the fitter to get minisim angular response
   fSuccessfulRadius->Sumw2();
   return fSuccessfulRadius; 
}

void MyMiniSim::GeneratePrimaries(G4Event* event)
{
  
  double dt = 0.0;

//  const int fNumPhotons = 1; //number of photons in each event
//  const int fNumPhotons = 90000; //number of photons in each event
  const double fEnergy = hbarc * twopi / (fPhotonWavelength * nm);
  G4ParticleDefinition* fOpticalPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

  for(int i=0; i < 180000; i++){
  //for(int i=0; i < 180; i++){
	//for(int i=0; i < fNumPhotons; i++){
        //double an = (G4UniformRand() * 180.)/2.;
        double an = (i % 180000)/2000.;
	//double an = (i % 180)/2.;
        if(an > 180.) std::cout << " an = " << an << std::endl;
        an = an * TMath::DegToRad();
        G4ThreeVector fPosition = G4ThreeVector(0., 1000.*sin(an), 1000.*cos(an));
        G4ThreeVector fNormal = -fPosition.unit();
        G4ThreeVector fX = fNormal.orthogonal().unit();
        G4ThreeVector fY = fNormal.cross( fX ).unit();

        const double theta = G4UniformRand() * 2.0 * pi;
        const double radius = G4UniformRand() * 300. * 300.;
        G4ThreeVector dx = ( cos( theta ) * fX + sin( theta ) * fY ) * sqrt( radius );
        dx += fPosition;

	G4PrimaryVertex* vertex = new G4PrimaryVertex(dx,dt);
        G4ThreeVector momentum = fNormal*fEnergy;
	G4PrimaryParticle* photon = new G4PrimaryParticle(fOpticalPhoton, momentum.x(), momentum.y(), momentum.z());

	//rand polarization
	double phi = (G4UniformRand()*2.0-1.0)*pi;
	G4ThreeVector e1 = fNormal.orthogonal().unit();
	G4ThreeVector e2 = fNormal.unit().cross(e1);
	G4ThreeVector pol = e1*cos(phi) + e2*sin(phi);
	photon->SetPolarization(pol.x(), pol.y(), pol.z());
	photon->SetMass(0.0);

	vertex->SetPrimary(photon);
	event->AddPrimaryVertex(vertex);

	} //loop over photons 
  

  
  /**
   //Chris' Version below

  double dt = 0.0;

//  const int fNumPhotons = 1; //number of photons in each event
//  const int fNumPhotons = 90000; //number of photons in each event
  const double fEnergy = hbarc * twopi / (fPhotonWavelength * nm);
  G4ParticleDefinition* fOpticalPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
  
  for(int i=0; i < 18000; i++){
  //for(int i=0; i < 1600; i++){
    //for(int i=0; i < fNumPhotons; i++){
    //double an = (G4UniformRand() * 180.)/2.;
    double an = (i % 18000)/200.;
    if(an > 180.) std::cout << " an = " << an << std::endl;
    an = an * TMath::DegToRad();

    //CJ 29/01/15
    //double an2 = G4UniformRand() * 2.0 * pi;

    //for(int j=-1; j<2; j++){
	
    G4ThreeVector fPosition = G4ThreeVector(0., 1000.*sin(an), 1000.*cos(an));
    //G4ThreeVector fPosition = G4ThreeVector(1000.*sin(an)*cos(an2), 1000.*sin(an)*sin(an2), 1000.*cos(an));
    //std::cout << "fPosition (X,Y,Z): " << fPosition.x() << ", " << fPosition.y() << ", " << fPosition.z() << ")" << std::endl;
    G4ThreeVector fNormal = -fPosition.unit();
    //std::cout << "fNormal (X,Y,Z): " << fNormal.x() << ", " << fNormal.y() << ", " << fNormal.z() << ")" << std::endl;
    G4ThreeVector fX = fNormal.orthogonal().unit();
    G4ThreeVector fY = fNormal.cross( fX ).unit();
    
    //G4ThreeVector faz = -fPosition.unit();
    //G4ThreeVector fax = fNormal.orthogonal().unit();
    //G4ThreeVector fay = fNormal.cross( fX ).unit();
    
    //for(int j=-45; j<46; j++){
    for(int j=-15; j<16; j++){

      for(int k=-15; k<16; k++){
    
	//std::cout << "j: " << j << std::endl;
	
	//Double_t j_ang = j * TMath::DegToRad();
	//std::cout << "j_ang: " << j_ang << std::endl;
	
	//G4ThreeVector fDirection = ((sin(j_ang)*fax) + (0*fay) + (cos(j_ang)*faz)); 
	//std::cout << "fDirection (X,Y,Z): " << fDirection.x() << ", " << fDirection.y() << ", " << fDirection.z() << ")" << std::endl;   
	
	
	//CJ 09/02/15   
	//const double theta = G4UniformRand() * 2.0 * pi;
	//const double radius = G4UniformRand() * 300. * 300.;
	//G4ThreeVector dx = ( cos( theta ) * fX + sin( theta ) * fY ) * sqrt( radius );
	G4ThreeVector dx = (( j*10 ) * fX) + (( k*10 ) * fY);
	dx += fPosition;
	//std::cout << "ijk: " << i << ", " << j << ", " << k << " and an: " << an << " and dx (X,Y,Z): " << dx.x() << ", " << dx.y() << ", " << dx.z() << ")" << std::endl;    
	
	//G4ThreeVector dx = fPosition;
	//std::cout << "dx (X,Y,Z): " << dx.x() << ", " << dx.y() << ", " << dx.z() << ")" << std::endl;    
	
	G4PrimaryVertex* vertex = new G4PrimaryVertex(dx,dt);
	G4ThreeVector momentum = fNormal*fEnergy;
	//G4ThreeVector momentumUnit=fDirection.unit();
	//G4ThreeVector momentum = fDirection.unit()*fEnergy;
	//G4ThreeVector momentum = momentumUnit*fEnergy;
	//std::cout << "momentum (X,Y,Z): " << momentum.x() << ", " << momentum.y() << ", " << momentum.z() << ")" << std::endl;
	
	G4PrimaryParticle* photon = new G4PrimaryParticle(fOpticalPhoton, momentum.x(), momentum.y(), momentum.z());
	
	//rand polarization
	double phi = (G4UniformRand()*2.0-1.0)*pi;
	//double phi = (0.5*2.0-1.0)*pi;
	G4ThreeVector e1 = fNormal.orthogonal().unit();
	G4ThreeVector e2 = fNormal.unit().cross(e1);
	G4ThreeVector pol = e1*cos(phi) + e2*sin(phi);
	photon->SetPolarization(pol.x(), pol.y(), pol.z());
	photon->SetMass(0.0);
	
	vertex->SetPrimary(photon);
	event->AddPrimaryVertex(vertex);
	//std::cout << " " << std::endl;
      
      }//k

    }//j
  
  } //loop over photons 
  **/
}



void MyMiniSim::BeginOfEventAction(const G4Event* /*event*/)
{

//  std::cout << "in BeginOfEventAction" << std::endl;
  fTrackingAction->SetTrackingLevel(TrackingAction::eCondensed);

}

void MyMiniSim::EndOfEventAction(const G4Event* g4ev)
{
  //this is where we see if the photon was successful
  // (i.e generated a photoelectron in the PMT bucket it entered)
  fEventCtr++;
//  std::cout << "in EndOfEventAction" << std::endl;
  TVector3 pmtDirection;
     pmtDirection.SetX(0.);
     pmtDirection.SetY(0.);
     pmtDirection.SetZ(-1.); //default is negative direction, correct for this in a sec

  // Trajectory Info
  G4TrajectoryContainer* trajectories = g4ev->GetTrajectoryContainer();
  if(trajectories != NULL) 
  {
  for( size_t iTrajectory = 0; iTrajectory < trajectories->size(); iTrajectory++ ) 
    {
      PhotonTrajectory* photonTrajectory = dynamic_cast< PhotonTrajectory* >( (*trajectories)[iTrajectory] );
      if( photonTrajectory != NULL ) // Add special MCPhoton information
        {
          std::vector<DS::MCPhoton> photons = photonTrajectory->GetMCPhotons();
          for( size_t iPhoton = 0; iPhoton < photons.size(); iPhoton++ )
            {

              TVector3 inPos = photons[iPhoton].GetInPosition();
              TVector3 inDir = photons[iPhoton].GetInDirection();
              const double angle = TMath::ACos( inDir.Dot( pmtDirection ) ) * TMath::RadToDeg();
              const double xyRadius = sqrt( inPos.X() * inPos.X() + inPos.Y() * inPos.Y() );

	      //std::cout << "inPos.Z() ( < 132.0): " << inPos.Z() << std::endl;
	      //std::cout << "xyRadius ( > 137.0): " << xyRadius << std::endl;
	      if( photons[iPhoton].GetFate() == RAT::DS::MCPhoton::ePhotoelectron){
		fSuccessfulZ->Fill(inPos.Z());
		fSuccessfulRadius->Fill(xyRadius);
	      }

              if( inPos.Z() < 132.0 ) continue; //ignore photons that entered PMT at 'invalid' z pos
	      if(xyRadius > 141.0) continue; //ignore photons that enter bucket outisde valid bucket radius

              fIncidentPhotons->Fill(angle);
              
              if( photons[iPhoton].GetFate() == RAT::DS::MCPhoton::eReflected) fReflectedPhotons->Fill(angle); 
              if( photons[iPhoton].GetFate() == RAT::DS::MCPhoton::ePhotoelectron) fSuccessfulPhotons->Fill(angle);
            }
        }
    }
   }//if not null trajectory

 if(false){
// if(fEventCtr%10000000 == 0)
 if(fEventCtr%1000000 == 0)
  {
  TH1D* SimAngResp = GetSimAngResp();
  char buffer[80];
  sprintf(buffer, "minisim_output/running_to_test_AR_no_aging_%.5f_%.2f.root", fCurrentp0, fCurrentp1);
  std::string str (buffer);
  TFile savef(str.c_str(), "recreate");
  MyMiniSim::ScaleProperly(SimAngResp);
  SimAngResp->Write();
  fIncidentPhotons->Write();
  fSuccessfulPhotons->Write();
  fReflectedPhotons->Write();
  savef.Close();
  }
  }

}



} // namespace RAT
