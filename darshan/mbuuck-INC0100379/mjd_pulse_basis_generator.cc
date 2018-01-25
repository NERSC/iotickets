/*
* Micah Buuck (mbuuck@uw.edu)
* 2/15
* 
* Description: 
* 
* In: Run numbers and energy window for the basis generation
*
* Out: Root file containing waveform library
*
* TODO: Cleanup and tweaking
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cctype>
#include <ctime>
#include <vector>
#include <list>
#include "GATDataSet.hh"
#include "MGTSuperpulse.hh"
#include "MGTWaveformLibrary.hh"
#include "MGWFBaselineRemover.hh"
#include "MGWFFitter.hh"
#include "MGWFCalculateChiSquare.hh"
#include "MGWFExtremumFinder.hh"
#include "MGWFTimePointCalculator.hh"
#include "MGWFResampler.hh"
#include "MGTEvent.hh"
#include "TChain.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TMath.h" //TMath::Pi()
#include "TEventList.h"
#include "TObjString.h"
#include "TFitResult.h"

using namespace std;

bool isint(const char* anInput) {
  int it=0;
  int len=(int)strlen(anInput);
  while(isdigit(anInput[it]) && it<len) {
    it++;
  }
  return it==len;
}

bool isfloat(const char* anInput) {
  int deccount = 0;
  int it=0;
  int len=(int)strlen(anInput);
  while((isdigit(anInput[it]) || anInput[it]=='.' || anInput[it]=='e' ||
    anInput[it]=='E') && it<len && deccount<2) {
    if(isdigit(anInput[it])) it++;
    else if(anInput[it]=='.') {
      deccount++;
      it++;
    }
    else if(anInput[it]=='e' || anInput[it]=='E') {
      it++;
      if(it<len && (anInput[it]=='-' || anInput[it]=='+') &&
        isint(string(anInput).substr(it+1).c_str())) it=len;
      else if(it<len && isint(string(anInput).substr(it).c_str())) it=len;
    } 
  }
  return it==len;
}

double lorentzian(double* x, double* p) {
  return p[0]/(TMath::Pi()*p[1]*(1+pow((x[0]-p[2])/p[1],2)));
}

void PrintHelpMenu(double chidefault, string pathdefault,
  double fitwindowfractiondefault) {
  cout << "\nHelp Menu for mjd_pulse_basis_generator:" << endl;
  cout << "  Required Arguments:" << endl;
  cout << "    -d <detector name>" << endl;
  cout << "      Specify ORTEC or BEGe detector name." << endl;
  cout << "    -r <run range lower bound> (optional)<run range upper bound>" << endl;
  cout << "      or <path to file containing run numbers>" << endl;
  cout << "      Specify run or range of runs over which to generate pulse shape basis," << endl;
  cout << "      or a file containing a list of run numbers with a single number on each" << endl;
  cout << "      line." << endl;
  cout << "    -e <energy window lower bound> <energy window upper bound>" << endl;
  cout << "      Specify energy window over which to generate pulse shape basis." << endl;
  cout << "    -c <channel>" << endl;
  cout << "      Can specify channel instead of detector name. (Not recommended)" << endl;
  cout << "  Optional Arguments:" << endl;
  cout << "    --no-AoE" << endl;
  cout << "      Skip A over E cut before building pulse basis." << endl;
  cout << "    --lorentzian" << endl;
  cout << "      Use Lorentzian for AoE distribution fit function. Default is gaussian." << endl;
  cout << "    -x <cutoff chi-squared>" << endl;
  cout << "      Specify cutoff chi-squared value for creating new base pulse." << endl;
  cout << "      Default is " << chidefault << "." << endl;
  cout << "    -f <fit window fraction>" << endl;
  cout << "      Specify fraction of waveform to be fitted, centered around 50\% risetime." << endl;
  cout << "      Default is " << fitwindowfractiondefault << endl;
  cout << "    -s <(optional) path to directory>" << endl;
  cout << "      Specify directory to save output .root file. Default is " << endl;
  cout << "      " << pathdefault << ". (pwd)" << endl;
  cout << "    --no-granularity" << endl;
  cout << "      Skip granularity cut." << endl;
  cout << "    -v" << endl;
  cout << "      Set verbose output." << endl;
  cout << "    -h" << endl;
  cout << "      Print this help menu." << endl;
  cout << endl;
}

bool HandleCommandLineArgs(int argc, char** argv, string& selectedDetector,
                           int& selectedchannel, int& startrun, int& endrun,
                           stringstream& runList, double& startenergy,
                           double& endenergy, bool& doAoE, bool& useLorentzian,
                           bool& doGranularity, double& cutoffchi,
                           string& path, double& fitwindowfraction,
                           bool& verbose)
{
  bool fail = false;
  for(int iArg=1; iArg<argc; iArg++) {
    string arg=argv[iArg];
    if(arg[0]=='-') {
      switch (arg[1]) {
        case 'h':
          PrintHelpMenu(cutoffchi, path, fitwindowfraction);
          fail = true;
          break;
        case 'd':
          iArg++;
          if(iArg>=argc) {
            cerr << "Please input detector name!" << endl;
            fail = true;
            break;
          }
          selectedDetector = argv[iArg];
          break;
        case 'c':
          iArg++;
          if(iArg>=argc) {
            cerr << "Please input detector name!" << endl;
            fail = true;
            break;
          }
          if(!isint(argv[iArg])) {
            cerr << "Channel must be an integer!" << endl;
            fail = true;
            break;
          }
          selectedchannel = atoi(argv[iArg]);
          break;
        case 'r':
          iArg++;
          if(iArg>=argc) {
            cerr << "Please input run number(s)!" << endl;
            fail = true;
            break;
          }
          if(!isint(argv[iArg])) {
            runList << argv[iArg];
            startrun = -1;
            break;
          }
          if(startrun != -1) {
            startrun = atoi(argv[iArg]);
            if((iArg+1)<argc && isint(argv[iArg+1])) {
              iArg++;
              endrun = atoi(argv[iArg]);
              if((endrun-startrun)<=0) {
                cerr << "Run range must be positive!" << endl;
                fail = true;
                break;
              }
            }
          }
          break;
        case 'e':
          iArg++;
          if((iArg+1)>=argc) {
            cerr << "Please input calibration energy window!" << endl;
            fail = true;
            break;
          }
          if(!isfloat(argv[iArg])) {
            cerr << "Calibration energy lower bound must be a positive number!" << endl;
            fail = true;
            break;
          }
          startenergy = atof(argv[iArg]);
          iArg++;
          if(!isfloat(argv[iArg])) {
            cerr << "Calibration energy upper bound must be a positive number!" << endl;
            fail = true;
            break;
          }
          endenergy = atof(argv[iArg]);
          if((endenergy-startenergy)<=0) {
            cerr << "Calibration energy window must have positive width!" << endl;
            fail = true;
            break;
          }
          break;
        case '-':
          if(arg=="--no-AoE") doAoE = false;
          else if(arg=="--lorentzian") useLorentzian = true;
          else if(arg=="--no-granularity") doGranularity = false;
          else {
            cerr << "Unknown option: " << arg << endl;
            fail=true;
          }
          break;
        case 'x':
          iArg++;
          if(iArg>=argc) {
            cerr << "Please input cutoff chi-squared value!" << endl;
            fail = true;
            break;
          }
          if(!isfloat(argv[iArg]) || atof(argv[iArg]) <= 0) {
            cerr << "Cutoff chi-squared value must be a positive number!" << endl;
            fail = true;
            break;
          }
          cutoffchi = atof(argv[iArg]);
          break;
        case 'f':
          iArg++;
          if(iArg>=argc) {
            cerr << "Please specify fit window fraction!" << endl;
            fail = true;
            break;
          }
          if(!isfloat(argv[iArg]) || atof(argv[iArg]) <=0 || atof(argv[iArg]) > 1) {
            cerr << "Fit window fraction must be a positive number between 0 and 1!" << endl;
            fail = true;
            break;
          }
          fitwindowfraction = atof(argv[iArg]);
          break;
        case 's':
          if( (iArg+1 >= argc) || (argv[iArg+1][0] == '-') ) break;
          iArg++;
          path = argv[iArg];
          break;
        case 'v':
          verbose = true;
          break;
        default:
          cerr << "Unknown option: " << arg << endl;
          fail = true;
      }
    }
    else {
      cerr << "Unknown option: " << arg << endl;
      fail = true;
    }
    if(fail) break;
  }
  if(!fail) {
    if(selectedDetector.empty() && selectedchannel==0) {
      cerr << "Please input detector name (-d) or channel number (-c)!" << endl;
      fail = true;
    }
    if(!selectedDetector.empty() && selectedchannel!=0) {
      cerr << "Cannot input both detector name (-d) and channel number (-c)!" << endl;
      fail = true;
    }
    if(startrun==0) {
      cerr << "Please input run number(s)! (-r)" << endl;
      fail = true;
    }
    if(startenergy==0) {
      cerr << "Please input calibration energy window! (-e)" << endl;
      fail = true;
    }
    if(useLorentzian && !doAoE) {
      cout << "Can't use Lorentzian if fitting not enabled!" << endl;
    }
  }
  return fail;
}

inline bool great(MGTSuperpulse* wf1, MGTSuperpulse* wf2)
  {return (*wf1) > (*wf2);}

int main(int argc, char** argv) {
  string selectedDetector;
  int selectedchannel = 0;
  int startrun = 0;
  int endrun = 0;
  stringstream runList;
  double startenergy = 0;
  double endenergy = 0;
  bool doAoE = true;
  bool useLorentzian = false;
  bool doGranularity = true;
  double cutoffchi = 1.5;
  double fitwindowfraction = 0.5;
  bool verbose = false;
  FILE *fp = popen("pwd", "r");
  char buf[1024];
  fgets(buf, 1024, fp);
  buf[strlen(buf)-1] = '\0'; // strip newline char
  string path = buf;
  pclose(fp);
  if(argc==1) {
    PrintHelpMenu(cutoffchi, path, fitwindowfraction);
    return 0;
  }
  bool clerr = HandleCommandLineArgs(argc, argv, selectedDetector,
                                     selectedchannel, startrun, endrun,
                                     runList, startenergy, endenergy, doAoE,
                                     useLorentzian, doGranularity, cutoffchi,
                                     path, fitwindowfraction, verbose);
  if(clerr) return 0;

  /*=============================Main execution=============================*/

  // Print call info to stdout 
  for(int i=0; i<argc; i++) {
    cout << argv[i] << " ";
  }
  cout << "\n\n";

  // Start timer
  cout << "Saving files to " << path << "\n\n";
  time_t starttime, nowtime, lastupdate;
  time(&starttime);
  time(&lastupdate);

  //------------------------Set up friend output tree-------------------------
  TTree* friendTree = new TTree("friendTree","friendTree");
  std::vector<int> NMatchingPulses;
  friendTree->Branch("NMatchingPulses",&NMatchingPulses);

  //--------------------------------Load Data---------------------------------
  GATDataSet ds;
  if(startrun == -1) ds.AddRunNumbersFromFile(runList.str().c_str());
  else if(endrun > 0) ds.AddRunRange(startrun, endrun);
  else ds.AddRunNumber(startrun);
  TChain* chain = ds.GetChains();
  if(chain==NULL) {
    cerr << "Could not find data!" << endl;
    return 1;
  }
  MJTChannelMap* channelMap = ds.GetChannelMap();
  if(channelMap==NULL) {
    cerr << "Could not find channel map!" << endl;
    return 1;
  }
  if(!selectedDetector.empty()) {
    selectedchannel = channelMap->GetInt(selectedDetector,MJTChannelMap::kIDHiCol);
    if(selectedchannel==0) {
      cerr << "Could not find channel for detector: " << selectedDetector << ". Are you sure it is in this dataset?" << endl;
      return 1;
    }
  }
  bool isM1data = ds.GetRunNumber(0)<10000000;
  // Load and Friend trees using GATDataSet public methods
    // This is the same as using GATDataSet::GetChains(), but returns
    // a non-const TChain* instead of a const TChain*
  /*string bchainName = "MGTree";
  string chainName = "";
  bool isM1data = true;
  if(startrun != -1) {
    string file = GATDataSet::GetPathToRun(startrun, GATDataSet::kGatified);
    unsigned long start = file.rfind('/')+1;
    chainName = file.substr(start, file.rfind('_')-start);
    chainName += "Tree";
  }
  else {
    string line;
    ifstream runFile(runList.str().c_str());
    cout << "runList: " << runList.str() << endl;
    if(runFile.is_open()) {
      getline(runFile,line);
      int firstrun;
      if(isint(line.c_str())) firstrun = atoi(line.c_str());
      else {
        cerr << line << " is not a valid run number! Exiting..." << endl;
        return 1;
      }
      string file = GATDataSet::GetPathToRun(firstrun, GATDataSet::kGatified);
      unsigned long start = file.rfind('/')+1;
      chainName = file.substr(start, file.rfind('_')-start);
      chainName += "Tree";
      runFile.close();
    }
  }
  TChain* chain = new TChain(chainName.c_str(), chainName.c_str());
  TChain* bchain = new TChain(bchainName.c_str(), bchainName.c_str());
  if(startrun != -1) {
    int iRun=startrun;
    do {
      if(iRun>10000000) isM1data = false;
      cout << "Loading run " << iRun << endl;
      string file = GATDataSet::GetPathToRun(iRun, GATDataSet::kGatified);
      string bfile = GATDataSet::GetPathToRun(iRun, GATDataSet::kBuilt);
      if(file != "") chain->AddFile(file.c_str());
      if(bfile != "") bchain->AddFile(bfile.c_str());
      iRun++;
    } while(iRun <= endrun);
  }
  else {
    string line;
    ifstream runFile(runList.str().c_str());
    if(runFile.is_open()) {
      while(getline(runFile,line)) {
        int iRun;
        if(isint(line.c_str())) iRun = atoi(line.c_str());
        else {
          cerr << line << " not a valid run number!" << endl;
          continue;
        }
        if(iRun>10000000) isM1data = false;
        cout << "Loading run " << iRun << endl;
        string file = GATDataSet::GetPathToRun(iRun, GATDataSet::kGatified);
        string bfile = GATDataSet::GetPathToRun(iRun, GATDataSet::kBuilt);
        if(file != "") chain->AddFile(file.c_str());
        if(bfile != "") bchain->AddFile(bfile.c_str());
      }
      runFile.close();
    }
    else {
      cerr << "Unable to open file " << runList.str() << endl;
    }
  }
  if(chain->GetEntries() == 0) {
    cerr << "No valid runs were loaded." << endl;
    return 0;
  }
  chain->AddFriend(bchain);
  bchain->AddFriend(chain);*/


  //-------------------Declare data objects and transforms--------------------
  int flatTime = 1000;//ns
  MGWFBaselineRemover baseline;
  baseline.SetBaselineTime(flatTime);
  double baselineVar = 0;
  
  MGWFFitter fitter;
  if(verbose) fitter.SetVerbose();
  fitter.SetFitWindowFraction(fitwindowfraction);
  fitter.SetTolerance(1); //Tolerance of 1 now uses parabola interpolation to solve for min
  fitter.SetMaxTOffset(0.01);
  
  MGWFCalculateChiSquare calcChi;
  
  MGWFTimePointCalculator TPC;
  TPC.AddPoint(0.5);
  TPC.AddPoint(0.02);
  TPC.AddPoint(0.98);
  
  MGWFResampler resampler;
  resampler.SetReferenceTOffset();
  
  MGWFExtremumFinder exFinder;
  exFinder.SetFindMaximum();
  

  //------------------------Set up MGTWaveformLibrary-------------------------
  MGTWaveformLibrary wfLibrary;
  //This source designation is only appropriate if the correct energy range
  //  has been supplied, but won't automatically change to reflect a different
  //  energy range.
  wfLibrary.SetSourceType(MGWaveformLibrary::kComptonEdge);
  //The library type will be automatically corrected if it is changed during
  //  basis generation.
  wfLibrary.SetLibraryType(MGWaveformLibrary::kMGTSuperpulse);
  //Put assembly info in MGTWaveformLibrary
  stringstream assemblyInfo;
  assemblyInfo << "Detector: " << selectedDetector;
  if(startrun == -1) assemblyInfo << ", Runs: " << runList.str();
  else {
    assemblyInfo << ", Runs: " << startrun;
    if(endrun != 0) assemblyInfo << " to " << endrun;
  }
  assemblyInfo << ", Energies: " << startenergy << " to " << endenergy;
  assemblyInfo << ", Chi^2 cutoff: " << cutoffchi;
  assemblyInfo << ", Fit window fraction: " << fitwindowfraction;

  stringstream filename;
  filename << path << "/wfLib_det" << selectedDetector << "_energyWindow"
    << int(startenergy) << "keVto" << int(endenergy) << "keV.root";
  TFile f(filename.str().c_str(),"RECREATE");

  // Set branch pointers
  char energyBranchName[99];
  vector<double>* energy = 0;
  if(isM1data) sprintf(energyBranchName,"trapENFCal");
  else {
    cout << "Not M1 data, using energyCal!" << endl;
    sprintf(energyBranchName,"energyCal");
  }
  TBranch* energyBranch = chain->GetBranch(energyBranchName);
  energyBranch->SetAddress(&energy);

  vector<double>* channel=0;
  TBranch* channelBranch = chain->GetBranch("channel");
  channelBranch->SetAddress(&channel);

  MGTEvent* event=0;
  TBranch* eventBranch = chain->GetBranch("event");
  eventBranch->SetAddress(&event);

  double TSRisetime = 0;
  if(selectedchannel==112 || selectedchannel==114 || selectedchannel==120 ||
     selectedchannel==144 || selectedchannel==148 || selectedchannel==152 ||
     selectedchannel==644 || selectedchannel==628 || selectedchannel==626 ||
     selectedchannel==624 || selectedchannel==576 || selectedchannel==692 ||
     selectedchannel==690 || selectedchannel==688 || selectedchannel==640 ||
     selectedchannel==614 || selectedchannel==594) {
    TSRisetime = 50;
  }
  else if(selectedchannel==118 || selectedchannel==150 || selectedchannel==646
          || selectedchannel==642 || selectedchannel==662 ||
          selectedchannel==656 || selectedchannel==696 || selectedchannel==608
          || selectedchannel==592) {
    TSRisetime = 100;
  }
  else if(selectedchannel==146 || selectedchannel==674 || selectedchannel==610
          || selectedchannel==598) {
    TSRisetime = 200;
  }
  else if(!selectedDetector.empty()) {
    // Using detector name instead of channel number
    TSRisetime = 50;
  }
  else {
    cerr << "Not a known channel: " << selectedchannel << endl;
    cerr << "  setting TSRisetime to 50ns!" << endl;
    TSRisetime = 50;
  }

  TH1D* chi2Hist = new TH1D("chi2Hist","Histogram of ln of best chi2 fit",1000,-1,25);
  //-------------------------Calculate baseline noise-------------------------
  unsigned long NEntries = chain->GetEntries();
  list<pair<tuple<size_t,size_t,size_t>, MGTWaveform > > wfList;
  cout << NEntries << " entries to parse." << endl;
  int NChannelEntries = 0;
  int currentTree = -1; //Generates segfault later if set to 0
  char energyAndChannelCut[99];
  unsigned long NCut = 0;
  if(doGranularity) sprintf(energyAndChannelCut,"%s>%g && %s<%g && channel==%d && !wfDCBits && m<3",energyBranchName,startenergy,energyBranchName,endenergy,selectedchannel);
  else sprintf(energyAndChannelCut,"%s>%g && %s<%g && channel==%d && !wfDCBits",energyBranchName,startenergy,energyBranchName,endenergy,selectedchannel);
  cout << "energyAndChannelCut: " << energyAndChannelCut << endl;

  //-------------Calculate AoE cut and apply it-------------------------------
  int NMultisite = 0;
  TVectorD saveAoECutoff(1);
  saveAoECutoff = 0;
  if(doAoE) {
    char AoEName[99];
    cout << "TSRisetime: " << TSRisetime << endl;
    int NBins = 1000;
    sprintf(AoEName, "TSCurrent%gnsMax/%s >> AoEHist(%d,0.006,0.02)",
            TSRisetime, energyBranchName, NBins);
    cout << "AoEName: " << AoEName << endl;
    chain->Draw(AoEName,energyAndChannelCut,"GOFF"); //This line implicitly declares TH1F*
                                    // AoEHist, but we have to find it
                                    // with gROOT later.
    TH1F* AoEHist = dynamic_cast<TH1F*>(gROOT->FindObject("AoEHist"));
    if(AoEHist==NULL) {
      cerr << "AoEHist not found!" << endl;
      return 0;
    }
    
    //Compute mode of AoEHist (seems like this should be a TH1 method already)
    double mode = AoEHist->GetBinCenter(AoEHist->GetMaximumBin());
    TF1* f1 = 0;
    if(useLorentzian) {
      cout << "Lorentzian distribution used for A/E" << endl;
      f1 = new TF1("f1",lorentzian,mode,1.1*AoEHist->GetBinCenter(AoEHist->GetNbinsX()),3);
      f1->SetParameter(0,1);
      f1->SetParameter(1,0.0001);
      f1->SetParameter(2,mode);
    }
    else {
      cout << "Gaussian distribution used for A/E" << endl;
      f1 = new TF1("f1","gaus(0)",mode*0.99,mode*1.01);//1.1*AoEHist->GetBinCenter(AoEHist->GetNbinsX()));
      f1->SetParameter(0,1);
      f1->SetParameter(1,mode);
      f1->SetParameter(2,mode*0.01);
      f1->SetParLimits(0,0,NEntries);
      f1->SetParLimits(1,0.99*mode,1.01*mode);
      f1->SetParLimits(2,0,mode*0.1);
    }
    f1->SetNpx(NBins);
    cout << "\nAoE fit results. Histogram and fit will be saved." << endl;
    TFitResultPtr r = AoEHist->Fit(f1,"RS"); //R uses specified range and S saves to TFitResultPtr
    AoEHist->Write();
    r->Write("AoEFitPtr");
    double xmin, xmax;
    f1->GetRange(xmin,xmax);
    f1->SetRange(0,xmax);
    f1->Write("AoEFittedFunction");
    double AoEMiddle;
    double AoEScale;
    if(useLorentzian) {
      // Pull parameters from Lorentzian distribution
      AoEMiddle = f1->GetParameter(2);
      AoEScale = f1->GetParameter(1);
      AoEScale = TMath::Tan(0.9997*TMath::Pi()/2.0)*AoEScale;
    }
    else {
      AoEMiddle = f1->GetParameter(1);
      AoEScale = 3*f1->GetParameter(2);
    }
    char AoECut[99];
    sprintf(AoECut, "%s && TSCurrent%gnsMax/%s > %g", energyAndChannelCut, TSRisetime, energyBranchName, 0.0);//AoEMiddle-AoEScale);
    saveAoECutoff = AoEMiddle - AoEScale;
    cout << "AoECut: " << AoECut << endl;
    NCut = chain->Draw(">>goodWfs",AoECut);
    for(int iBin=AoEHist->GetXaxis()->FindBin(AoEMiddle-AoEScale); iBin<AoEHist->GetXaxis()->FindBin(AoEMiddle); iBin++) {
      cout << "iBin: " << iBin << ", bin center: " << AoEHist->GetBinCenter(iBin) << ", Bin Content: " << AoEHist->GetBinContent(iBin) << ", NMultisite: " << NMultisite << endl;
      NMultisite += AoEHist->GetBinContent(iBin) - f1->Eval(AoEHist->GetBinCenter(iBin));
    }
  }
  else NCut = chain->Draw(">>goodWfs",energyAndChannelCut);
  cout << NCut << " events passed all cuts." << endl;
  if(NMultisite) cout << "Expect roughly " << NMultisite << " multisite events." << endl; 
  TEventList* goodWfs = dynamic_cast<TEventList*>(gROOT->FindObject("goodWfs"));
  if(goodWfs==NULL) {
    cerr << "goodWfs not found!" << endl;
    return 0;
  }
  chain->SetEventList(goodWfs);
  //----------------Finished with AoE cut-------------------------------------

  for(unsigned long iEntry=0; iEntry<NCut; iEntry++) {
    // First loop over entries to calulate baseline noise
    time(&nowtime);
    if(difftime(nowtime, lastupdate) > 1) {
      double progress = (double)iEntry/(double)NCut;
      cout << 100*progress << "\% finished calculating baseline noise." << endl;
      time(&lastupdate);
    }
    //Only load data we need. Make sure correct tree is being pointed to, and
    //  load from energy, channel, and event
    int branchEntry = chain->LoadTree(goodWfs->GetEntry(iEntry));
    int newTree = chain->GetTreeNumber();
    size_t runNumber = ds.GetRunNumber(newTree);
    if(newTree!=currentTree) {
      currentTree = newTree;
      energyBranch = chain->GetBranch(energyBranchName);
      energyBranch->SetAddress(&energy);
      channelBranch = chain->GetBranch("channel");
      channelBranch->SetAddress(&channel);
      eventBranch = chain->GetBranch("event");
      eventBranch->SetAddress(&event);
    }

    energyBranch->GetEntry(branchEntry);
    channelBranch->GetEntry(branchEntry);
    bool eventLoaded = false;

    for(unsigned long iWfm=0; iWfm<channel->size(); iWfm++) {
      if(channel->at(iWfm)==selectedchannel && energy->at(iWfm)<endenergy && energy->at(iWfm)>startenergy) {
        //If event passes cuts, then load the waveforms.
        if(!eventLoaded) { eventBranch->GetEntry(branchEntry); eventLoaded = true; }
        MGTWaveform wf = *event->GetWaveform(iWfm);
        wf.SetLength(2016); //Chop off last two samples to account for multisampling errors
        baseline.TransformInPlace(wf);
        /*TPC.TransformInPlace(wf);
        if(TPC.GetFromStartRiseTime(2)-TPC.GetFromStartRiseTime(1)>wf->GetDuration()*fitwindowfraction
           || TPC.GetFromStartRiseTime(0)<wf->GetDuration()*0.45
           || TPC.GetFromStartRiseTime(0)>wf->GetDuration()*0.55) {
          cout << "Pulse " << iEntry << ", " << iWfm << " is bad. Skipping.\n";
        }
        else {*/
        wf.SetTOffset(0);
        wfList.push_back(make_pair(make_tuple(runNumber,branchEntry,iWfm),wf));
        baseline.CalculateBaselineAndRMS(wfList.back().second);
        baselineVar += baseline.GetBaselineRMS()*baseline.GetBaselineRMS();
        NChannelEntries++;
        //}
      }
    }
  }
  baselineVar /= double(NChannelEntries);
  cout << "baselineVar = " << baselineVar << endl;
  if(baselineVar <= 0) {
    cerr << "baselineVar not a positive number. Aborting!" << endl;
    return 0;
  }
  NCut = wfList.size();
  cout << '\n' << NCut << " waveforms passed all cuts.\n" << endl;
  //--------------------Finished calculating baseline noise-------------------

  //----------------------------Create superpulses----------------------------
  cout << "Creating superpulses" << endl;
  int i=0;
  auto wfit=wfList.begin();
  while (wfit!=wfList.end()) {
    // Second loop over waveforms to create superpulses
    time(&nowtime);
    if(difftime(nowtime, lastupdate) > 1) {
      double progress = (double)i/(double)wfList.size();
      cout << 100*progress << "\% finished creating superpulses." << endl;
      time(&lastupdate);
    }
    MGTWaveform* wf = &(wfit->second);
    double minchi = -1;
    //int minindex = -1;
    //double minTOffset = 0;
    //double minScale = 1;
    //double minVertOff = 0;
    //Perform a chisquare comparison to template pulses, and find minimum.
    //This loop gets skipped if NSuperpulses==0.
    for (unsigned long iSP=0; iSP<wfLibrary.GetSize(); iSP++) {
/*      if(i==18 && iSP==64) {
        verbose = true;
        fitter.SetVerbose();
      }
      else {
        verbose = false;
        fitter.SetVerbose(false);
      }
*/      if(verbose) cout << "Perform dynamic_cast" << endl;
      MGTSuperpulse* libSuper = dynamic_cast <MGTSuperpulse*> (wfLibrary[iSP]);
      if(verbose) cout << "dynamic_cast performed sucessfully" << endl;
      fitter.SetTemplate(&(libSuper->GetTemplate()));
      if(verbose) {
        cout << "fitter.SetTemplate(libSuper->GetTemplate()); successful" << endl;
        cout << "libSuper->GetTemplate().GetLength() = " << libSuper->GetTemplate().GetLength() << endl;
        cout << "libSuper->GetTemplate().GetVectorData()[0] = " << libSuper->GetTemplate().GetVectorData()[0] << endl;
        cout << "wf->GetLength() = " << wf->GetLength() << endl;
        cout << "wf->GetVectorData()[0] = " << wf->GetVectorData()[0] << endl;
      }
      fitter.CalculateParameters(*wf);
      if(verbose) cout << "Parameters calculated successfully" << endl;
      double chival = fitter.GetParameterValue(MGWFFitter::kChi2)/(baselineVar/((1+fitter.GetParameterValue(MGWFFitter::kNorm))/2.0));
      calcChi.SetTOffset(fitter.GetParameterValue(MGWFFitter::kTOffset));
      calcChi.SetNorm(fitter.GetParameterValue(MGWFFitter::kNorm));
      calcChi.SetVerticalOffset(fitter.GetParameterValue(MGWFFitter::kVerticalOffset));
      calcChi.TransformOutOfPlace(*wf, *fitter.GetTemplate());
      if(verbose) cout << "chival calculated successfully" << endl;
      if (minchi==-1 || minchi>chival) {
        minchi = chival;
        //minindex = iSP;
        //minTOffset = fitter.GetParameterValue(MGWFFitter::kTOffset);
        //minScale = fitter.GetParameterValue(MGWFFitter::kNorm);
        //minVertOff = fitter.GetParameterValue(MGWFFitter::kVerticalOffset);
        if(verbose) cout << "New minchi found" << endl;
      }
    }
    bool erased = false;
    if (minchi>cutoffchi || minchi==-1) { 
    //if minimum chisquared is above a preset value, or if no superpulses,
    //create new template
      chi2Hist->Fill(log(minchi));
      MGTSuperpulse* superP = new MGTSuperpulse;
      TPC.TransformInPlace(*wf);
      double tOff = TPC.GetFromStartRiseTime(0);
      tOff -= wf->GetDuration()*0.5;
      MGTWaveform tempWF = *wf;
      tempWF.SetTOffset(tOff);
      resampler.TransformOutOfPlace(*wf,tempWF);
      wf = &tempWF;
      wf->SetTOffset(0);
      superP->MakeSimilarTo(*wf);
      superP->SetTemplate(*wf);
      size_t runNumber = get<0>(wfit->first);
      size_t eventNumber = get<1>(wfit->first);
      size_t waveformNumber = get<2>(wfit->first);
      superP->SetTemplateID(runNumber, eventNumber, waveformNumber);
      superP->SetBuildMultiGraph();
      superP->AddPulse(*wf, runNumber, eventNumber, waveformNumber);
      wfit = wfList.erase(wfit);
      erased = true;
      wfLibrary.AddWaveformToEndOfLibrary(superP);
      if(verbose) cout << "New superpulse created." << endl;
    }
    if(!erased) {++wfit; ++i;}
  }


  //---------------------------Populate superpulses---------------------------
  cout << "\nPopulating superpulses." << endl;
  bool fail=false;
  bool failonce=false;
  wfit=wfList.begin();
  double wfListsize = wfList.size();
  i=0;
  while(wfit!=wfList.end()) {
    time(&nowtime);
    if(difftime(nowtime, lastupdate) > 1) {
      double progress = (double)i/wfListsize;
      cout << 100*progress << "\% finished populating superpulses." << endl;
      time(&lastupdate);
    }
    MGTWaveform* wf = &(wfit->second);
    double minchi = -1;
    int minindex = -1;
    for (unsigned long iSP=0; iSP<wfLibrary.GetSize(); iSP++) {
      MGTSuperpulse* libSuper = dynamic_cast <MGTSuperpulse*> (wfLibrary[iSP]);
      fitter.SetTemplate(&(libSuper->GetTemplate()));
      fitter.CalculateParameters(*wf);
      double chival = fitter.GetParameterValue(MGWFFitter::kChi2)/(baselineVar/((1+fitter.GetParameterValue(MGWFFitter::kNorm))/2.0));
      if (minchi==-1 || minchi>chival) {
        minchi = chival;
        minindex = iSP;
      }
    }
    if (minindex<0) {
      cerr << "No matching superpulse found!";
      break;
    }
    MGTSuperpulse* libSuper = dynamic_cast <MGTSuperpulse*> (wfLibrary[minindex]);
    fitter.SetTemplate(&(libSuper->GetTemplate()));
    fitter.Transform(wf);
    
    /*if(fabs(fitter.GetParameterValue(MGWFFitter::kTOffset))>1000) {
      cout << endl << "Error adding to iSP: " << minindex << endl;
      cout << "with NPulses: " << libSuper->GetNPulses() << endl;
      cout << "Waveform: " << i << endl;
      cout << "tOffset = " << fitter.GetParameterValue(MGWFFitter::kTOffset) << endl; 
      fail = true;
      failonce = true;
    }*/
      
    if(fitter.GetParameterValue(MGWFFitter::kNorm) <= 0) cout << "Waveform: " << i << ", libSuper: " << minindex << endl;
    exFinder.TransformInPlace( *libSuper );
    double max = exFinder.GetTheExtremumValue();
    if(baseline.GetBaseline(*libSuper) > 0.1*max) {
      cout << endl << "Baseline too high!" << endl;
      cout << "Baseline: " << baseline.GetBaselineMean() << endl;
      cout << "iSP: " << minindex << endl;
      cout << "NPulses: " << libSuper->GetNPulses() << endl;
      cout << "Waveform: " << i << endl;
    }
    if(!fail) {
      chi2Hist->Fill(log(minchi));
      size_t runNumber = get<0>(wfit->first);
      size_t eventNumber = get<1>(wfit->first);
      size_t waveformNumber = get<2>(wfit->first);
      libSuper->AddPulse(*wf,runNumber,eventNumber,waveformNumber);
    }
    wfit = wfList.erase(wfit);
    ++i;
    fail=false;
  }
  cout << "Finished populating superpulses." << endl << endl;


  //-----------------------------Normalize pulses-----------------------------
  // Baseline subtract and normalize max to 1
  cout << endl;
  for(unsigned long iSP=0; iSP < wfLibrary.GetSize(); iSP++) {
    MGTSuperpulse* libSuper = dynamic_cast<MGTSuperpulse*>(wfLibrary[iSP]);
    if(libSuper->GetNPulses()==0) continue;
    baseline.TransformInPlace( (*wfLibrary[iSP]) );
    exFinder.TransformInPlace( (*wfLibrary[iSP]) );
    (*wfLibrary[iSP]) /= exFinder.GetTheExtremumValue();
  }
  int deleted=0;
  for(unsigned long iSP=0; iSP < wfLibrary.GetSize(); iSP++) {
    MGTSuperpulse* libSuper = dynamic_cast<MGTSuperpulse*>(wfLibrary[iSP]);
    if(libSuper->GetNPulses()<2) {
      if(libSuper->GetNPulses()==1) deleted++;
      //delete libSuper;
      //wfLibrary.DeleteWaveformFromLibrary(iSP);
      //iSP--;
    }
  }

  if(!failonce) wfLibrary.SortDescending();
  for(unsigned long iSP=0; iSP < wfLibrary.GetSize(); iSP++) {
    MGTSuperpulse* libSuper = dynamic_cast<MGTSuperpulse*>(wfLibrary[iSP]);
    cout << "Superpulse " << iSP << " has " << libSuper->GetNPulses() << " pulses." << endl;
  }


  //-------------------------Save everything to disk--------------------------
  TObjString saveselectedDetector;
  saveselectedDetector.SetString(selectedDetector.c_str());
  TVectorD saveselectedchannel(1);
  saveselectedchannel = selectedchannel;
  TVectorD savecutoffchi(1);
  savecutoffchi = cutoffchi;
  TVectorD savebaselineVar(1);
  savebaselineVar = baselineVar;
  TString saveenergyBranchName;
  saveenergyBranchName = energyBranchName;
  f.cd();
  chi2Hist->Write("chi2Hist");
  saveselectedDetector.Write("saveselectedDetector");
  saveselectedchannel.Write("saveselectedchannel");
  savecutoffchi.Write("savecutoffchi");
  savebaselineVar.Write("savebaselineVar");
  saveAoECutoff.Write("saveAoECutoff");
  ds.Write("dataSet");
  wfLibrary.Write("wfLibrary");


  //----------------------Stop timer and delete objects-----------------------
  time(&nowtime);
  double seconds = difftime(nowtime,starttime);
  cout << "Code took " << seconds << " seconds to run." << endl;
  cout << "That's " << (double)NCut/seconds << " waveforms/second." << endl;
  while(wfLibrary.GetSize()>0) {
    MGTSuperpulse* libSuper = dynamic_cast <MGTSuperpulse*> (wfLibrary[wfLibrary.GetSize()-1]);
    delete libSuper;
    wfLibrary.DeleteWaveformFromEndOfLibrary();
  }
  cout /*<< "Deleted "*/<< deleted << " unique pulses from " << NCut << " total. ";
  cout << "(" << deleted/double(NCut)*100 << "\%)" << endl;
  double fiveOrMore = wfLibrary.GetPopularity(5);
  cout << fiveOrMore*100 << "\% of superpulses with popularity 5 or more." << endl;  
  cout << "A/E estimated " << NMultisite/double(NCut)*100 << "\%  multisite pulses." << endl;
  //delete chain;
  //delete bchain;
  return 0;
}
