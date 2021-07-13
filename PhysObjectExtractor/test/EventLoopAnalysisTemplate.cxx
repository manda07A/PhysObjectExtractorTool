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
// g++ -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnalysisTemplate.cxx $(root-config --cflags --libs)
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
#include "Math/Vector4D.h"
//Include C++ classes
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>

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

//Requiered trigger
string triggerRequest = "HLT_L2DoubleMu23_NoVertex";

// Fixed size dimensions of array or collections stored in the TTree if any.


class EventLoopAnalysisTemplate {
public :
	
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain	
  TTree          *tevents;
   TTree          *tvertex;
   TTree          *ttrigger;
  TTree           *tmuons;
  TTree           *ttaus;

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
   std::map<std::string, int> *triggermap;
   vector<float>   *muon_pt;
   vector<float>   *muon_eta;
   vector<float>   *muon_phi;
   vector<float>   *muon_tightid;
   vector<float>   *tau_pt;
   vector<float>   *tau_eta;
   vector<float>   *tau_phi;
   vector<float>   *tau_ch;
   vector<float>   *tau_iddecaymode;
   vector<float>   *tau_idisotight;
   vector<float>   *tau_idantieletight;
   vector<float>   *tau_idantimutight;
   vector<float>   *tau_reliso_all;
  

   // List of example branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_triggermap;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_tightid;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_ch;   //!
   TBranch        *b_tau_iddecaymode;   //!
   TBranch        *b_tau_idisotight;   //!
   TBranch        *b_tau_idantieletight;   //!
   TBranch        *b_tau_idantimutight;   //!
   TBranch        *b_tau_reliso_all;   //!

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
  bool isGoodMuon(Int_t idx);
  bool isGoodTau(Int_t idx);
  bool FindGoodMuons();
  bool FindGoodTaus();
  std::vector<int> FindMuonTauPair();
 
};

/*
 * Helper function to compute the difference in the azimuth coordinate taking
 * the boundary conditions at 2 * pi into account.
 */
//-------------------------------------------------------------------------
namespace Helper {
  template <typename T>
  float DeltaPhi(T v1, T v2, const T c = M_PI)
  {
//-------------------------------------------------------------------------
    float r = std::fmod(v2 - v1, 2.0 * c);
    if (r < -c) {
      r += 2.0 * c;
    }
    else if (r > c) {
      r -= 2.0 * c;
    }
    return r;
  }
}//-----------------Helper


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
      //always, by default, use mytrigger as starting directory/tree
      //because it is the most complex
      TDirectory * dir = (TDirectory*)f->Get(filename+":/mytriggers");
      dir->GetObject("Events",tree);

      //Get trees for friendship
      tevents = (TTree*)f->Get("myevents/Events");
      tvertex = (TTree*)f->Get("mypvertex/Events");
      tmuons = (TTree*)f->Get("mymuons/Events");
      ttaus = (TTree*)f->Get("mytaus/Events");

      //Make friendship	
      tree->AddFriend(tevents);
      tree->AddFriend(tvertex);
      tree->AddFriend(tmuons);
      tree->AddFriend(ttaus);

	
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
  //cout<<" Set the environment to read one entry"<<endl;
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

   // Set object pointer
   triggermap =0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_tightid = 0;
   tau_pt = 0;
   tau_eta = 0;
   tau_phi = 0;
   tau_ch = 0;
   tau_iddecaymode = 0;
   tau_idisotight = 0;
   tau_idantieletight = 0;
   tau_idantimutight = 0;
   tau_reliso_all = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //Comment out to be able to read map
   //https://root-forum.cern.ch/t/std-map-in-ttree-with-makeclass/14171
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("triggermap",&triggermap,&b_triggermap);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_tightid", &muon_tightid, &b_muon_tightid);
   fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_ch", &tau_ch, &b_tau_ch);
   fChain->SetBranchAddress("tau_iddecaymode", &tau_iddecaymode, &b_tau_iddecaymode);
   fChain->SetBranchAddress("tau_idisotight", &tau_idisotight, &b_tau_idisotight);
   fChain->SetBranchAddress("tau_idantieletight", &tau_idantieletight, &b_tau_idantieletight);
   fChain->SetBranchAddress("tau_idantimutight", &tau_idantimutight, &b_tau_idantimutight);
   fChain->SetBranchAddress("tau_reliso_all", &tau_reliso_all, &b_tau_reliso_all);
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
       //cout<<"Load the current event"<<endl;
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       // if (Cut(ientry) < 0) continue;

       analysis();

    }
}


//-----------------------------------------------------------------
void EventLoopAnalysisTemplate::analysis()
{
//-----------------------------------------------------------------

  //cout<<"analysis() execution"<<endl;
  if (!MinimalSelection()) return;
  if (!FindGoodMuons()) return;
  if (!FindGoodTaus()) return;
  vector<int> GoodMuonTauPair = FindMuonTauPair();
  if (GoodMuonTauPair[0]==-1 && GoodMuonTauPair[1]==-1) return;
  
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

  //cout<<"Applying minimal selection"<<endl;
  bool isTrigger = false;
  bool isMuons = false;
  bool isTaus = false;

  //Check trigger and acceptance bit
  for (map<string, int>::iterator it=triggermap->begin();it!=triggermap->end();it++){
    if(it->first.find(triggerRequest)!=string::npos &&
       it->second!=0){
	 //cout<<it->first<<"  "<<it->second<<endl;
      isTrigger = true;
    }
  }

  //Enforce presence of muons
  Int_t nmuons = muon_pt->size();
  if (nmuons>0){isMuons = true;}

  //Enforce presence of taus
  Int_t ntaus = tau_pt->size();
  if (ntaus>0){isTaus = true;}

  return (isTrigger && isMuons && isTaus);

}//------MinimalSelection


// Give index of muon and check if passes the good muon check
//-----------------------------------------------------------------
bool EventLoopAnalysisTemplate::isGoodMuon(Int_t idx) 
{
//-----------------------------------------------------------------
  bool isGoodMuon = false;

  float mu_eta_cut = 2.1;
  float mu_pt_cut = 17; //in GeV
  if (abs(muon_eta->at(idx))<mu_eta_cut && 
	muon_pt->at(idx)>mu_pt_cut && 
      bool(muon_tightid->at(idx)) ){
    isGoodMuon = true;
  }

  return isGoodMuon;
  
}//----------isGoodMuon




/*
 * Reduce the loop to the interesting events containing at least one interesting
 * muon candidate.
 */
//-----------------------------------------------------------------
bool EventLoopAnalysisTemplate::FindGoodMuons() 
{
//-----------------------------------------------------------------
  bool GoodMuonFound = false;
  
  Int_t nmuons = muon_pt->size();
  for (Int_t j=0; j<nmuons;++j){
    if (isGoodMuon(j)){
      GoodMuonFound = true;
      return GoodMuonFound;
    }
  }    

  return GoodMuonFound;

 }//-------FindGoodMuons


/*
 * Find the interesting taus in the tau collection
 *
 * The tau candidates in this collection represent hadronic decays of taus, which
 * means that the tau decays to combinations of pions and neutrinos in the final
 * state.
 */
//----------------------------------------------------------------- 
bool EventLoopAnalysisTemplate::isGoodTau(Int_t idx) 
{
//-----------------------------------------------------------------
 bool isGoodTau = false;
 
 float tau_eta_cut = 2.3;
 float tau_pt_cut = 20; //in GeV
 
 if (tau_ch->at(idx)!=0 &&
	abs(tau_eta->at(idx))<tau_eta_cut && 
	tau_pt->at(idx)>tau_pt_cut && 
	bool(tau_iddecaymode->at(idx)) &&
	bool(tau_idisotight->at(idx)) &&
	bool(tau_idantieletight->at(idx)) &&
	bool(tau_idantimutight->at(idx))){
      isGoodTau = true;
    }

 return isGoodTau;

}//-----------isGoodTau





/*
 * Reduce the loop to the interesting events containing at least one interesting
 * tau candidate.
 */
//-----------------------------------------------------------------
bool EventLoopAnalysisTemplate::FindGoodTaus() 
{
//-----------------------------------------------------------------
  bool GoodTauFound = false;

  Int_t ntaus = tau_pt->size();
  for (Int_t j=0; j<ntaus;++j){
    if (isGoodTau(j)){
      GoodTauFound = true;
      return GoodTauFound;
    }
  }    

  return GoodTauFound;

 }//-------FindGoodTaus


//-----------------------------------------------------------------
std::vector<int> EventLoopAnalysisTemplate::FindMuonTauPair() 
{
//-----------------------------------------------------------------

  //Find all possible pairs of muons and taus
  vector< pair<int,int> > comb;
  Int_t nmuons = muon_pt->size();
  Int_t ntaus = tau_pt->size();
  for(Int_t midx=0;midx<nmuons;++midx){
    for(Int_t tidx=0;tidx<ntaus;++tidx){
      comb.push_back(make_pair(midx,tidx));
    }
  }
  const size_t numComb= comb.size();

  //Find valid pairs based on delta r
  vector<int> validPair(numComb, 0);
  for(size_t i = 0; i < numComb; i++) {
    const int i1 = comb.at(i).first;
    const int i2 = comb.at(i).second;
    if(isGoodMuon(i1) && isGoodTau(i2)) {
      const float deltar = sqrt(
				pow(muon_eta->at(i1) - tau_eta->at(i2), 2) +
				pow(Helper::DeltaPhi(muon_phi->at(i1), tau_phi->at(i2)), 2));
      if (deltar > 0.5) {
	validPair[i] = 1;
      }
    }
  }

// Find best muon based on pt
  int idx_1 = -1;
  float maxPt = -1;
  for(size_t i = 0; i < numComb; i++) {
    if(validPair[i] == 0) continue;
    const int tmp = comb.at(i).first;
    if(maxPt < muon_pt->at(tmp)) {
      maxPt = muon_pt->at(tmp);
      idx_1 = tmp;
    }
  }

// Find best tau based on iso
  int idx_2 = -1;
  float minIso = 999;
  for(size_t i = 0; i < numComb; i++) {
    if(validPair[i] == 0) continue;
    if(int(comb.at(i).first) != idx_1) continue;
    const int tmp = comb.at(i).second;
    if(minIso > tau_reliso_all->at(tmp)) {
      minIso = tau_reliso_all->at(tmp);
      idx_2 = tmp;
    }
  }

  vector<int> thegoodidx;
  thegoodidx.push_back(idx_1);
  thegoodidx.push_back(idx_2);

  return thegoodidx;
  

}//---------FindMuonTauPair





//-----------------------------------------------------------------
int main()
{
//-----------------------------------------------------------------

  gROOT->ProcessLine("#include<map>");

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
  return 0;

}


