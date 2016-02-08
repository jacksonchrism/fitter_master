#include <RAT/DB.hh>
#include <RAT/G4Stream.hh>
#include <RAT/RunManager.hh>
#include <RAT/Log.hh>
#include <RAT/DS/Run.hh>
#include <RAT/Producer.hh>
#include <G4RunManager.hh>
#include <Randomize.hh>
#include <TFitter.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include "MyMiniSim.hh"
#include <RAT/ProcBlock.hh>
#include <RAT/PhotonThinning.hh>
#include <RAT/PMTCharge.hh>
#include <RAT/ChannelEfficiency.hh>
#include <RAT/PMTTransitTime.hh>
#include <RAT/HitPMTCollection.hh>
#include <RAT/PMTVariation.hh>
#include <RAT/DU/Utility.hh>

using namespace RAT;

std::string fname;
std::string fname2;
std::string dataset_info;
TFile* SNOdata;
TH1D* SNOpmtr;
TH1D* bestfit;

void ScaleProperly(TH1D* result)
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



void makeFile()
{ 
    //create file for saving parameter values
    std::ofstream outfile;
    outfile.open(fname.c_str(), std::ios::app);
    if(outfile.is_open()){
    std::cout << "Created file for saving stuff: " << fname.c_str()  << std::endl;
    outfile.close(); 
    }

    //create file for saving parameter values w/ bin values
    std::ofstream outfile2;
    outfile2.open(fname2.c_str(), std::ios::app);
    if(outfile2.is_open()){
    std::cout << "Created file for saving more stuff: " << fname2.c_str() << std::endl;
    outfile2 << 0.0 << '\t' << 1 << '\t' << 10 << '\t' << 20 << '\t' << 30 << '\t' << 40 << '\t' << 45 << '\t' << 50 << '\t' << 60 << std::endl;
    outfile2.close(); 
    }
}


void trackbinsvsParams(double p0, int nbins, double* bins)
{
    //this writes out a text file with parameter values
    std::cout << "trackbinsvsParams " << std::endl;
    std::ofstream thingy;
    thingy.open(fname2.c_str(), std::ios::app);
    if(thingy.is_open()){
    thingy << p0 << '\t';
        for(int i=0; i < nbins; i++){
        thingy << *(bins + i) << '\t';
	}
    thingy << std::endl;
    thingy.close();
    }
}

void trackChi2vsParams(double Chi2, double p0, double p1)
{
    //this writes out a text file with parameter values
    std::ofstream thingy;
    thingy.open(fname.c_str(), std::ios::app);
    if(thingy.is_open()){
    thingy << Chi2 << '\t' << p0 << '\t' << p1 << std::endl;
    thingy.close();
    }
}


int main(int argc, char* argv[])
{
   if ( argc != 5 ) // argc should be 4 for correct execution
   {   // We print argv[0]: it is the program name
      std::cout<<"usage: "<< argv[0]  << " <may02, apr03, oct03, sep01> <number of event iterations> <p0 value> >p1 value>\n";
      return -1; // exit
   }

   std::string dataset = argv[1];
   dataset_info = dataset;

   std::stringstream convert(argv[2]);
   int events_per_iteration;
   if(!(convert >> events_per_iteration))
     events_per_iteration=0;
   std::cout << "Number of iterations: " << events_per_iteration << std::endl;

   std::stringstream convert2(argv[3]);
   double p0_val;
   if(!(convert2 >> p0_val))
     p0_val=1;
   std::cout << "p0: " << p0_val << std::endl;

   std::stringstream convert3(argv[4]);
   double p1_val;
   if(!(convert3 >> p1_val))
     p1_val=1;
   std::cout << "p1: " << p1_val << std::endl;

   // start RAT logging otherwise horrible error
   RAT::Log::Init("wwfitter.log");

    //set the seed (again). Setting it here only did not set the seed in
    // minisim, so that's why it's set in both places.
    CLHEP::HepRandom::setTheSeed(1234567);
    std::cout << "seed set" << std::endl;

    // quiet down G4 output
    SetG4coutStream(G4Stream::DETAIL);
    SetG4cerrStream(G4Stream::WARN);

    // load ratdb tables
    std::cout << "Load DB..." << std::endl;
    RAT::DB::Get()->LoadDefaults();
    RAT::DB::Get()->Load("AGED_CONCENTRATOR_PARAMS.ratdb");
    RAT::DB::Get()->SetS("CONCENTRATOR", "cRAT", "model_type", "ConcOpticsAged");
    RAT::DB::Get()->SetS("DETECTOR", "","geo_file", "empty.geo");
    RAT::DB::Get()->SetS("DETECTOR","", "pmt_info_file", "singlepmt.ratdb");

    RAT::DB::Get()->SetI("PMTCALIB", "", "use_qhs_hhp", 0);
    //set starting values(?)
    RAT::DB::Get()->SetD("AGED_CONCENTRATOR_PARAMS", "", "p0", p0_val);
    RAT::DB::Get()->SetD("AGED_CONCENTRATOR_PARAMS", "", "p1", p1_val);

    //make a file name to keep params in
    //okay so the file name comes out kind of crazy/ugly.
    //I tried somethings to make it niceer but they didn't work so whatever
    time_t rawtime;
    struct tm* timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    //std::string temp = asctime(timeinfo);
    //char* temp2;
    //temp.copy(temp2, 19, 0);

    char buffer[120];
    char bufferA[120];
    char dataset_chars[5] = {dataset_info[0], dataset_info[1], dataset_info[2], dataset_info[3], dataset_info[4]};    
    sprintf(buffer, "/data/snoplus/home/jackson/SNO+_angular/fitter_exponential/fitOutput/params_");
    strcat(buffer,dataset_chars);
    strftime(bufferA,120, "_%F-%H-%M-%S.txt",timeinfo);
    strcat(buffer,bufferA);
    fname = buffer;

    char buffer2[120];
    char bufferB[120];
    sprintf(buffer2, "/data/snoplus/home/jackson/SNO+_angular/fitter_exponential/poppick/values_");
    strcat(buffer2,dataset_chars);
    strftime(bufferB,120, "_%F-%H-%M-%S.txt",timeinfo);
    strcat(buffer2,bufferB);
    fname2 = buffer2;
    //create the files
    makeFile();
    std::cout << asctime(timeinfo) << std::endl;


    // set up the geant4 run environment
    // RAT's RunManager does most of the work for us
    RAT::RunManager* run_manager = new RAT::RunManager;
    G4RunManager* g4_run_manager = G4RunManager::GetRunManager();
    g4_run_manager->Initialize();
    DU::Utility::Get()->BeginOfRun();
    PhotonThinning::BeginOfRun();
    PMTCharge::Get()->BeginOfRun(); // Must preceed ChannelEfficiency
    ChannelEfficiency::Get()->BeginOfRun();
    PMTTransitTime::Get()->BeginOfRun();
    HitPMTCollection::Get()->BeginOfRun();
    PMTVariation::Get()->BeginOfRun();

    //Set data set
    if(dataset.compare("sep00") == 0)
      {
	SNOdata = new TFile("/data/snoplus/home/jackson/SNO+_angular/snodata_allAngles/data_Sep00_386.root");
      }
    if(dataset.compare("sep01") == 0)
      {
	SNOdata = new TFile("/data/snoplus/home/jackson/SNO+_angular/snodata_allAngles/data_Sep01_386.root");
      }
    if(dataset.compare("may02") == 0)
      {
	SNOdata = new TFile("/data/snoplus/home/jackson/SNO+_angular/snodata_allAngles/data_May02_386.root");
      }
    if(dataset.compare("apr03") == 0)
      {
	SNOdata = new TFile("/data/snoplus/home/jackson/SNO+_angular/snodata_allAngles/data_May02_386.root");
      }

    SNOpmtr = (TH1D*) SNOdata->Get("AngularResponse");

    std::cout << "dataset " << dataset << " with " << events_per_iteration << " events" << std::endl;

    double fp0;
    double fp1;

    //get the degredation parameters
    DBLinkPtr agedConcentratorDB = DB::Get()->GetLink( "AGED_CONCENTRATOR_PARAMS"); 
    try
      {
        fp0 = agedConcentratorDB->GetD( "p0" );
        fp1 = agedConcentratorDB->GetD( "p1" );
      }
    catch( RAT::DBNotFoundError &e ) 
      { 
        RAT::Log::Die( "ConcentratorAgedOpticalModel::ConcentratorAgedOpticalModel: Cannot Find AGED_CONCENTRATOR Parameters" );
      }


    double fit_p0 = fp0;//par[0];
    double fit_p1 = fp1;//par[1];
    //std::cout << "WRITE THESE VLAUES TO DB" << std::endl;
    std::cout << "current fit p0: " << fit_p0 << std::endl;
    std::cout << "current fit p1: " << fit_p1 << std::endl;

    //write previoius iteration params to DB file for this iteration
    //writeParams(fit_p0, fit_p1);
    RAT::DB::Get()->SetD("AGED_CONCENTRATOR_PARAMS", "", "p0", fit_p0);
    RAT::DB::Get()->SetD("AGED_CONCENTRATOR_PARAMS", "", "p1", fit_p1);

    //std::cout << "CHECK THAT IT WROTE CORRECTLY" << std::endl;
    
//    RAT::DB::Get()->Load("AGED_CONCENTRATOR_PARAMS.ratdb");
    double p0;
    double p1;

    std::cout << "in fitter about to make a minisim" << std::endl;
    // run the MiniSim and extract angular response 
    RAT::MyMiniSim sim(386, fit_p0, fit_p1);
    std::cout << "a minisim is made, about to beamon" << std::endl;
    sim.BeamOn(events_per_iteration);
    TH1D* SIMpmtr = sim.GetSimAngResp();
    TH1D* SIMsuccP = sim.GetSuccessfulPhotons();
    TH1D* SIMincP = sim.GetIncidentPhotons();
    TH1D* SIMsuccZ = sim.GetSuccessfulZ();
    TH1D* SIMsuccRadius = sim.GetSuccessfulRadius();
    ScaleProperly(SIMpmtr);
    bestfit = (TH1D*)SIMpmtr->Clone("best_fit");

   //save bestfit
    char buffy[80];
    //char dataset_chars[5] = {dataset_info[0], dataset_info[1], dataset_info[2], dataset_info[3], dataset_info[4]};    
    //sprintf(buffy, "minisim_output/test_bestfit_reflectx1pt0_%s.root",dataset_chars);
    sprintf(buffy, "minisim_output/test_bestfit-%s-p0_%g-p1_%g.root",dataset_chars,p0_val,p1_val);
    std::string str(buffy);
    TFile bestf(str.c_str(), "recreate");
    //char buffy2[80];
    //sprintf(buffy2, "test_bestfit_%s_%.5f_%.5f.root",dataset_chars, p0, p1);
    //sprintf(buffy2, "test_bestfit_%s_p0%d_p1%d.root",dataset_chars,p0_val,p1_val);
    std::string str2(buffy);
    bestfit->SetName(str2.c_str());
    bestfit->SetTitle(str2.c_str());
    bestfit->Write();
    SNOpmtr->Write();
    SIMsuccP->Write();
    SIMincP->Write();
    SIMsuccZ->Write();
    SIMsuccRadius->Write();
    bestf.Close();


    //calculate chi2
    double chi2 = 0;
    for(int i = 0; i < 90; i++){
	double ex = SNOpmtr->GetBinContent(i+1);
	double ob = SIMpmtr->GetBinContent(i+1);
        //if(ex == 0.0) break;
	double err_ob = SIMpmtr->GetBinError(i+1);
	double err_ex = SNOpmtr->GetBinError(i+1);
	double err2 = err_ob*err_ob + err_ex*err_ex;
//	double err2 = err_ob*err_ob;

	double temp = (ob - ex) * (ob - ex) / err2;
	if(!std::isnan(temp)){
	  if(ex!=0.0 && ob!=0.0) {
	    chi2 += temp; 
	    std::cout << "i: " << i << "ch2 = " << chi2 << std::endl;
	  }
	}
    }

    trackChi2vsParams(chi2, p0, p1);
    std::cout << "ch2 = " << chi2 << std::endl;

    delete SIMpmtr;
    
    SNOdata->Close();
    
}
