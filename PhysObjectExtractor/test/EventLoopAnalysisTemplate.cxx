//////////////////////////////////////////////////////////////////////
// This template analysis code has been built with fragments from the 
// classes automatically obtained by the TTree MakeClass() method.
//
// The template shows the structre of a potential analysis code
// where more TTree friends can be added with more physics objects.
//
// Done with ROOT version 5.32/00
// from TTree Events/Events
// found on file: myoutput.root
//
//
// Compile me with:
// g++ -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnaalysisTemplate.cxx $(root-config --cflags --libs)
/////////////////////////////////////////////////////////////////////

//Include ROOT classes
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "TLatex.h"
#include "TStopwatch.h"
//Include C++ classes
#include <iostream>
#include <vector>
#include <map>
#include <string>

using namespace std;


/*
 * Base path to local filesystem or to EOS containing the datasets
 */
//const std::string samplesBasePath = "root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/";
const std::string samplesBasePath = "";


//book example histograms for specific variables
//copy them in the constructor if you add more
const int numberOfHistograms = 3;
TH1F* dataRunB_npv = new TH1F("dataRunB_npv","Number of primary vertices",25,5,30);
TH1F* dataRunC_npv = new TH1F("dataRunC_npv","Number of primary vertices",25,5,30);
TH1F* ZTT_npv = new TH1F("ZTT_npv","Number of primary vertices",25,5,30);




// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxtriggermap = 5;

class EventLoopAnalysisTemplate {
public :
	
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain	
   TTree          *tvertex;
   TTree          *ttrigger;
   //Add more trees for friendship

   Int_t           fCurrent; //!current Tree number in a TChain
  TString          labeltag;
   TString         filename;

  //array to keep histograms to be written and easily loop over them
   TH1F            *hists[numberOfHistograms];

   // Declaration of example leaf types
   Int_t           run;
   UInt_t          luminosityBlock;
   ULong64_t	   event;
   Int_t           PV_npvs;
   Int_t           triggermap_;
   string          triggermap_first[kMaxtriggermap];
   Int_t           triggermap_second[kMaxtriggermap];   //[triggermap_]


   // List of example branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_triggermap_;   //!
   TBranch        *b_triggermap_first;   //!
   TBranch        *b_triggermap_second;   //!



  EventLoopAnalysisTemplate(TString filename, TString labeltag);
  virtual ~EventLoopAnalysisTemplate();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  void analysis();
  bool MinimalSelection();
 
};

EventLoopAnalysisTemplate::EventLoopAnalysisTemplate(TString thefile, TString thelabel) : fChain(0)
{
  //Prepare some info for the object:
  filename = thefile;
  labeltag = thelabel;
  

  //Load histograms of interest to the object for automatic looping
  hists[0] = dataRunB_npv;
  hists[1] = dataRunC_npv;
  hists[2] = ZTT_npv;

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TTree* tree = 0;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f || !f->IsOpen()) {
         f = new TFile(filename);
      }
      //always, by default, use myevents as starting directory/tree
      TDirectory * dir = (TDirectory*)f->Get(filename+":/myevents");
      dir->GetObject("Events",tree);

      //Get trees for friendship
      tvertex = (TTree*)f->Get("mypvertex/Events");
      ttrigger = (TTree*)f->Get("mytriggers/Events");
      
      //Make friendship	
      tree->AddFriend(tvertex);
      tree->AddFriend(ttrigger);
	
   }
   Init(tree);
}

EventLoopAnalysisTemplate::~EventLoopAnalysisTemplate()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventLoopAnalysisTemplate::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}



Long64_t EventLoopAnalysisTemplate::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}


void EventLoopAnalysisTemplate::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("triggermap", &triggermap_, &b_triggermap_);
   fChain->SetBranchAddress("triggermap.first", triggermap_first, &b_triggermap_first);
   fChain->SetBranchAddress("triggermap.second", triggermap_second, &b_triggermap_second);
   Notify();
}


Bool_t EventLoopAnalysisTemplate::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventLoopAnalysisTemplate::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t EventLoopAnalysisTemplate::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void EventLoopAnalysisTemplate::Loop()
{
  if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	//Just an informative printout
	if(jentry%1000 == 0) {
	      cout<<"Processed "<<jentry<<" events out of "<<nentries<<endl;
	} 
       //Load the current event	
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       // if (Cut(ientry) < 0) continue;

       //Perform the analysis
       analysis();

    }
}


//-----------------------------------------------------------------
void EventLoopAnalysisTemplate::analysis()
{
//-----------------------------------------------------------------

  //cout<<"analysis() execution"<<endl;

  if (!MinimalSelection()) return;
  
  //fill histograms
  Int_t histsize = sizeof(hists)/sizeof(hists[0]);
  for (Int_t j=0;j<histsize;++j){

    TString histname = TString(hists[j]->GetName());
    TString thelabel = TString(histname.Tokenize("_")->At(0)->GetName());
    TString thevar = TString(histname.Tokenize("_")->At(1)->GetName());

    if (thelabel == labeltag && thevar == "npv"){
      hists[j]->Fill(PV_npvs);
    }

  }


}//------analysis()

/*
 * Perform a selection on the minimal requirements of an event
 */
//-----------------------------------------------------------------
bool EventLoopAnalysisTemplate::MinimalSelection()
{
//-----------------------------------------------------------------

  return true;

}//------MinimalSelection



//-----------------------------------------------------------------
int main()
{
//-----------------------------------------------------------------

  vector< pair<string,string> > sampleNames;
  //sampleNames.push_back(make_pair("GluGluToHToTauTau","ggH"));
  //sampleNames.push_back(make_pair("VBF_HToTauTau","qqH"));
  //sampleNames.push_back(make_pair("W1JetsToLNu","W1J"));
  //sampleNames.push_back(make_pair("W2JetsToLNu","W2J"));
  //sampleNames.push_back(make_pair("W3JetsToLNu","W3J"));
  //sampleNames.push_back(make_pair("TTbar","TT"));
  //sampleNames.push_back(make_pair("DYJetsToLL","ZLL"));
  sampleNames.push_back(make_pair("Run2012B_TauPlusX","dataRunB"));
  sampleNames.push_back(make_pair("Run2012C_TauPlusX","dataRunC"));
  sampleNames.push_back(make_pair("DYJetsToLL","ZTT"));


			
  //loop over sample files with names  defined above
  for(UInt_t j=0; j<sampleNames.size();++j){
    TString samplename = sampleNames.at(j).first;
    TString thelabel = sampleNames.at(j).second;
  
    cout << ">>> Processing sample " << samplename <<" with label "<<thelabel<<":" <<endl;
    TStopwatch time;
    time.Start();
 
    TString filename = samplesBasePath+samplename+".root";
    
    cout<<"Build the analysis object with file "<<filename<<endl;
    EventLoopAnalysisTemplate mytemplate(filename,thelabel);

    cout<<"Run the event loop"<<endl;
    mytemplate.Loop();

  }

  TFile* hfile = new TFile("histograms.root","RECREATE");
  dataRunB_npv->Write();
  dataRunC_npv->Write();
  ZTT_npv->Write();
  hfile->Close();
  return 1;

}


