// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// Centrality QA
//
// ################################################################

#include "StSetupMaker.h"
#include "StRoot/StarRoot/StMemStat.h"

// ROOT includes
#include <set>
#include "TFile.h"
#include <sstream>
#include <fstream>

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"

// extra includes
#include "StJetPicoDefinitions.h"

ClassImp(StSetupMaker)

//______________________________________________________________________________
StSetupMaker::StSetupMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", bool mDoComments = kFALSE)
  : StJetFrameworkPicoBase(name)
{
  fDebugLevel = 0;
  fRunFlag = 0;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  mOutName = outName;
  fBadTowerListVers = 1;  // FIXME
  doRejectBadRuns = kFALSE;
  fBadRunListVers = 999;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fEmcTriggerEventType = 0; 
  fMBEventType = 2;
  doComments = mDoComments;
  fAnalysisMakerName = name;
}

//_____________________________________________________________________________
StSetupMaker::~StSetupMaker()
{ /*  */
  // destructor
}

//_____________________________________________________________________________
Int_t StSetupMaker::Init() {
  // initialize the histograms
  DeclareHistograms();

  // Add bad run lists
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp (200 GeV)
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_w_missing_HT)  AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2012_BadRuns_P12id_w_missing_HT.txt");
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_wo_missing_HT) AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2012_BadRuns_P12id_wo_missing_HT.txt");
        break;
  
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu (200 GeV)
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_w_missing_HT)  AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2014_BadRuns_P18ih_w_missing_HT.txt");
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_wo_missing_HT) AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2014_BadRuns_P18ih_wo_missing_HT.txt");
        break; 

    case StJetFrameworkPicoBase::Run14_AuAu200_MB : // Run14 AuAu (200 GeV) - HFT data set (MB trigger)
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_w_missing_HT)  AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2014_BadRuns_P18ih_w_missing_HT.txt");
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_wo_missing_HT) AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2014_BadRuns_P18ih_wo_missing_HT.txt");
        break;
  
    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu (200 GeV)
        AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2016_BadRuns_P16ij.txt");
        break; 
  
    default :
      AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Empty_BadRuns.txt");
  }

  // Add dead + bad tower lists
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp (200 GeV)
        if(fBadTowerListVers == 102) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2012_BadTowers_102.txt");
        if(fBadTowerListVers == 1)   AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2012_BadTowers_Rag.txt"); // Raghav's Zg list
        if(fBadTowerListVers == 169) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2012_AltBadTowers_155_ALT.txt"); // Alt list of 155, +14 = 169

        // tower threshold cuts
        if(fBadTowerListVers == 999)     AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2012_BadTowers_P12id.txt");
        if(fBadTowerListVers == 9990200) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2012_BadTowers_P12id_200MeV.txt");
        if(fBadTowerListVers == 9991000) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2012_BadTowers_P12id_1000MeV.txt");
        if(fBadTowerListVers == 9992000) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2012_BadTowers_P12id_2000MeV.txt");

        //AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_BadTowers.txt");
        AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2012_DeadTowers.txt");
        break;

    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu (200 GeV)
        if(fBadTowerListVers ==  1)  AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers.txt");   // original default
        if(fBadTowerListVers == 136) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_AltBadTowers_136.txt");
        if(fBadTowerListVers == 51)  AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers50_ALT.txt");// 50x + some manually added

        // P18ih - need new definitions from Nick (July 17, 2019)
        // tower threshold cuts
        if(fBadTowerListVers == 999)     AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers_P18ih.txt");
        if(fBadTowerListVers == 9990200) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers_P18ih_200MeV.txt");
        if(fBadTowerListVers == 9991000) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers_P18ih_1000MeV.txt");
        if(fBadTowerListVers == 9992000) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers_P18ih_2000MeV.txt");

        AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_DeadTowers.txt");
        break;

    case StJetFrameworkPicoBase::Run14_AuAu200_MB : // Run14 AuAu (200 GeV) - HFT dataset (MB trigger)
        if(fBadTowerListVers ==  1)  AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers.txt");   // original default
        if(fBadTowerListVers == 136) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_AltBadTowers_136.txt");
        if(fBadTowerListVers == 51)  AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers50_ALT.txt");// 50x + some manually added

        // P18ih - need new definitions from Nick (July 17, 2019)
        // tower threshold cuts
        if(fBadTowerListVers == 999)     AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers_P18ih.txt");
        if(fBadTowerListVers == 9990200) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers_P18ih_200MeV.txt");
        if(fBadTowerListVers == 9991000) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers_P18ih_1000MeV.txt");
        if(fBadTowerListVers == 9992000) AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers_P18ih_2000MeV.txt");

        AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_DeadTowers.txt");
        break;

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu (200 GeV)
        AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2016_BadTowers.txt");
        AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2016_DeadTowers.txt");
        break;

    default :
      AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_BadTowers.txt");
      AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_DeadTowers.txt");
  }

  return kStOK;
}

//____________________________________________________________________________
Int_t StSetupMaker::Finish() { 
  cout << "StSetupMaker::Finish()\n";

  return kStOK;
}

//_____________________________________________________________________________
void StSetupMaker::DeclareHistograms() {
  // declare histograms here
}

//
// write histograms
//_____________________________________________________________________________
void StSetupMaker::WriteHistograms() {
  // write histograms here
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StSetupMaker::Clear(Option_t *opt) {

}
// 
//  This method is called every event.
//_____________________________________________________________________________
Int_t StSetupMaker::Make() {
  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // get PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // get pointer to PicoEvent 
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  // this class does not inherit from base class: StJetFrameworkPicoBase, but we want to reduce redundancy
  //StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase();
  StJetFrameworkPicoBase *baseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!baseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get run number, check bad runs list if desired (kFALSE if bad)
  int fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !IsRunOK(fRunNumber) ) return kStOK;
  }

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get vertex 3-vector and z-vertex component
  TVector3 mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();
 
  // commented out for now, but centrality was configured with 30 cm z-vertex data 
  // Z-vertex cut - the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  return kStOK;
}

//
// Reset bad run list object
//____________________________________________________________________________
void StSetupMaker::ResetBadRunList( ){
  badRuns.clear();
}
//
// Add bad runs from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
//_________________________________________________________________________________
Bool_t StSetupMaker::AddBadRuns(TString csvfile){
  // open infile
  std::string line;
  std::ifstream inFile ( csvfile );

  __DEBUG(2, Form("Loading bad runs from %s", csvfile.Data()) );

  if( !inFile.good() ) {
    __WARNING(Form("Can't open %s", csvfile.Data()) );
    return kFALSE;
  }

  while(std::getline (inFile, line) ){
    if( line.size()==0 ) continue; // skip empty lines
    if( line[0] == '#' ) continue; // skip comments

    std::istringstream ss( line );
    while( ss ){
      std::string entry;
      std::getline( ss, entry, ',' );
      int ientry = atoi(entry.c_str());
      if(ientry) {
        badRuns.insert( ientry );
        __DEBUG(2, Form("Added bad run # %d", ientry));
      }
    }
  }

  return kTRUE;
}
//
// Function: check on if Run is OK or not
//____________________________________________________________________________________________
Bool_t StSetupMaker::IsRunOK( Int_t mRunId ){
  //if( badRuns.size()==0 ){
  if( badRuns.empty() ){
    __ERROR("StSetupMaker::IsRunOK: WARNING: You're trying to run without a bad run list. If you know what you're doing, deactivate this throw and recompile.");
    throw ( -1 );
  }
  if( badRuns.count( mRunId )>0 ){
    __DEBUG(9, Form("Reject. Run ID: %d", mRunId));
    return kFALSE;
  } else {
    __DEBUG(9, Form("Accept. Run ID: %d", mRunId));
    return kTRUE;
  }
}
//
// Function: reset bad tower list
//____________________________________________________________________________
void StSetupMaker::ResetBadTowerList( ){
  badTowers.clear();
}
//
// Function: check on if Tower is OK or not
//____________________________________________________________________________________________
Bool_t StSetupMaker::IsTowerOK( Int_t mTowId ){
  //if( badTowers.size()==0 ){
  if( badTowers.empty() ){
    __ERROR("StSetupMaker::IsTowerOK: WARNING: You're trying to run without a bad tower list. If you know what you're doing, deactivate this throw and recompile.");
    throw ( -1 );
  }
  if( badTowers.count( mTowId )>0 ){
    __DEBUG(9, Form("Reject. Tower ID: %d", mTowId));
    return kFALSE;
  } else {
    __DEBUG(9, Form("Accept. Tower ID: %d", mTowId));
    return kTRUE;
  }
}
//
// Add bad towers from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
//_________________________________________________________________________________
Bool_t StSetupMaker::AddBadTowers(TString csvfile){
  // open infile
  std::string line;
  std::ifstream inFile ( csvfile );

  __DEBUG(2, Form("Loading bad towers from %s", csvfile.Data()) );

  if( !inFile.good() ) {
    __WARNING(Form("Can't open %s", csvfile.Data()) );
    return kFALSE;
  }

  while(std::getline (inFile, line) ){
    if( line.size()==0 ) continue; // skip empty lines
    if( line[0] == '#' ) continue; // skip comments

    std::istringstream ss( line );
    while( ss ){
      std::string entry;
      std::getline( ss, entry, ',' );
      int ientry = atoi(entry.c_str());
      if(ientry) {
        badTowers.insert( ientry );
        __DEBUG(2, Form("Added bad tower # %d", ientry));
      }
    }
  }

  return kTRUE;
}
//
// Function: reset dead tower list
//____________________________________________________________________________
void StSetupMaker::ResetDeadTowerList( ){
  deadTowers.clear();
}
//
// Add dead towers from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
//_________________________________________________________________________________
Bool_t StSetupMaker::AddDeadTowers(TString csvfile){
  // open infile
  std::string line;
  std::ifstream inFile ( csvfile );

  __DEBUG(2, Form("Loading bad towers from %s", csvfile.Data()) );

  if( !inFile.good() ) {
    __WARNING(Form("Can't open %s", csvfile.Data()) );
    return kFALSE;
  }

  while(std::getline (inFile, line) ){
    if( line.size()==0 ) continue; // skip empty lines
    if( line[0] == '#' ) continue; // skip comments

    std::istringstream ss( line );
    while( ss ){
      std::string entry;
      std::getline( ss, entry, ',' );
      int ientry = atoi(entry.c_str());
      if(ientry) {
        deadTowers.insert( ientry );
        __DEBUG(2, Form("Added bad tower # %d", ientry));
      }
    }
  }

  return kTRUE;
}
//
// Function: check on if Tower is DEAD or not
//____________________________________________________________________________________________
Bool_t StSetupMaker::IsTowerDead( Int_t mTowId ){
  //if( deadTowers.size()==0 ){
  if( deadTowers.empty() ){
    __ERROR("StSetupMaker::IsTowerDead: WARNING: You're trying to run without a dead tower list. If you know what you're doing, deactivate this throw and recompile.");
    throw ( -1 );
  }
  if( deadTowers.count( mTowId )>0 ){
    __DEBUG(9, Form("Reject. Tower ID: %d", mTowId));
    return kTRUE;
  } else {
    __DEBUG(9, Form("Accept. Tower ID: %d", mTowId));
    return kFALSE;
  }
}
