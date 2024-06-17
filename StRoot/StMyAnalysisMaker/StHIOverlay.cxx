// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StHIOverlay.h"
#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"


// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"

// jet-framework includes
#include "StFJWrapper.h"
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StHIOverlay)

//________________________________________________________________________
StHIOverlay::StHIOverlay(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* filename = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StHIOverlay::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  fDoEffCorr = kFALSE;
  doRejectBadRuns = kFALSE;
  fCorrJetPt = kFALSE;
  fMinPtJet = 0.0;
  fTrackBias = 0.0;
  fTowerBias = 0.0;
  fJetRad = 0.4;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0;
  fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2; fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0; fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0;  fTowerPhiMaxCut = 2.0*TMath::Pi();
  fCentralityScaled = 0.; ref16 = -99; ref9 = -99; // FIXME - maybe not make global
  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0; fMBEventType = 2;
  fRho = 0x0;
  fRhoVal = 0;
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
  // fjw("HIMC", "HIMC");
  fJets = 0x0;
  // fJetsArr = 0x0;
  fJetsArr.clear();

  fMCFileListName = filename;

  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

  mBemcGeom = 0x0;

  fMinJetTrackPt = 0.2; fMaxJetTrackPt = 30.0;
  fJetTrackEtaMin = -1.0; fJetTrackEtaMax = 1.0;
  fJetTrackPhiMin = 0.0; fJetTrackPhiMax = 2.0*TMath::Pi();
  fJetTrackDCAcut = 3.0;
  fJetTracknHitsFit = 15;
  fJetTracknHitsRatio = 0.52;

  fJetTowerEMin = 0.2; fJetTowerEMax = 100.0;
  fJetTowerEtaMin = -1.0; fJetTowerEtaMax = 1.0;
  fJetTowerPhiMin = 0.0; fJetTowerPhiMax = 2.0*TMath::Pi();
  mTowerEnergyMin = 0.2;

  fSetNumberOfEvents = 0;
  fNumberofeventsoverlayed = 0;

  for (int i = 0; i < 100; i++){
    fMCEventTracks[i].clear();
    fMCEventTowers[i].clear();
    fRecoEventTracks[i].clear();
    fRecoEventTowers[i].clear();
    fMCD0Information[i] = {};
    fRecoD0Information[i] = {};
    fOrigin[i].SetXYZ(0,0,0);
  }

  randomevents = kTRUE;
  fPrintLevel = 0;

  if (!name) return;
  SetName(name);
}

//
//________________________________________________________________________
StHIOverlay::~StHIOverlay()
{ /*  */

}

//
//________________________________________________________________________
Int_t StHIOverlay::Init() {

  StJetFrameworkPicoBase::Init();

  // position object for Emc
  mBemcGeom = StEmcGeom::instance("bemc");
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  fJetsArr.clear();

  ifstream filelistforMCEvents(fMCFileListName.Data());

  if (!filelistforMCEvents.is_open()) {
    LOG_ERROR << "No MC File List! Exiting!" << endm; 
    return kStOk;
  }

  string line;

  // cout << "Files read in: " << endl;

  while(getline(filelistforMCEvents,line))
  {
    TString s(line);

    // cout << s << endl;
    filenamesforHIOverlay.push_back(s);  
  }

  TFile f("/star/u/droy1/Y2019/STAR/Momentum_resolution_SL16d.root");
  fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPion");
  fKaonMomResolution = (TF1*)f.Get("fKaon")->Clone("fKaon");
  fProtonMomResolution = (TF1*)f.Get("fProton")->Clone("fProton");

  TFile effweight("/star/u/droy1/Y2019/STAR/EffWeightsInCentralityBins.root");

  fPionWeight[0] = (TGraph *)effweight.Get("Pion_0_10");
  fKaonWeight[0] = (TGraph *)effweight.Get("Kaon_0_10");
  fProtonWeight[0] = (TGraph *)effweight.Get("Proton_0_10");
  fAProtonWeight[0] = (TGraph *)effweight.Get("AProton_0_10");

  fPionWeight[1] = (TGraph *)effweight.Get("Pion_10_40");
  fKaonWeight[1] = (TGraph *)effweight.Get("Kaon_10_40");
  fProtonWeight[1] = (TGraph *)effweight.Get("Proton_10_40");
  fAProtonWeight[1] = (TGraph *)effweight.Get("AProton_10_40");

  fPionWeight[2] = (TGraph *)effweight.Get("Pion_40_80");
  fKaonWeight[2] = (TGraph *)effweight.Get("Kaon_40_80");
  fProtonWeight[2] = (TGraph *)effweight.Get("Proton_40_80");
  fAProtonWeight[2] = (TGraph *)effweight.Get("AProton_40_80");

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StHIOverlay::Finish() { 
  cout << "StHIOverlay::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StHIOverlay::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StHIOverlay::Clear(Option_t *opt) {
  // fjw->Clear();
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StHIOverlay::Make() {
  // fjw->Clear();
  fJets->Delete();

  fJetsArr.clear();
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  for (int i = 0; i < 100; i++){
    fMCEventTracks[i].clear();
    fMCEventTowers[i].clear();
    fRecoEventTracks[i].clear();
    fRecoEventTowers[i].clear();
    fMCD0Information[i] = {};
    fRecoD0Information[i] = {};
    fOrigin[i].SetXYZ(0,0,0);
  }

  fNumberofeventsoverlayed = 0;

  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent 
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();
  deadTowers = mBaseMaker->GetDeadTowers();
  badTowers = mBaseMaker->GetBadTowers();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  //for(int i = 1; i<4801; i++) {  if(!mBaseMaker->IsTowerOK(i))  cout<<"tower: "<<i<<" is not good!!"<<endl;  }
  // can check for bad towers like this (if needed):
  // bool isTowOK = mBaseMaker->IsTowerOK(towerID);
  // ===========================================================================================

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  double zVtx_VPD = mPicoEvent->vzVpd();
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CenttMaker! Skip! " << endm;
    return kStWarn;
  }

  // centrality variables
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
  ref9 = mCentMaker->GetRef9();   // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16(); // binning from central -> peripheral
  int cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cout << "Centrality = " << grefMult << "\t" << ref16 << "\t" << cent16 << "\t" << fCentralityScaled << endl;


  if (fCentralityScaled >= 0 && fCentralityScaled < 10) centralitybinforefficiency = 0;
  else if (fCentralityScaled >= 10 && fCentralityScaled < 40) centralitybinforefficiency = 1;
  else if (fCentralityScaled >= 40 && fCentralityScaled < 80) centralitybinforefficiency = 2;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // cout << "1" << endl;

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // cout << "2" << endl;

  // ============================ end of CENTRALITY ============================== //

  if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5) return kStOK;

  // cout << "3" << endl;

  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  // cout << "4" << endl;

  int arrMB5_Run14[]  = {450005, 450015, 450025, 450050, 450060};

  bool matchMB = kFALSE;

  for(int i = 0; i < sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14); i++) {
    if(mPicoEvent->isTrigger(arrMB5_Run14[i])) matchMB = kTRUE;
    if(matchMB) break;
  }

  if (!matchMB) return kStOk;

  // cout << "5" << endl;

  if (abs(zVtx) > 6.) return kStOk;

  // cout << "6" << endl;

  if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > 2.) return kStOk;

  // cout << "7" << endl;

  if (abs(zVtx - zVtx_VPD) > 3) return kStOk;

  // cout << "8" << endl;

  // cout << "Centrality = " << fCentralityScaled << endl;
  // pick a random file
  // PickARandomFile();

  // run Towers:
  SampleMCEvents();

  // cout << "Jets Array Size (HIOverlay): " << fJetsArr.size() << endl;

  return kStOK;
}

void StHIOverlay::PickARandomFile(){


}


void StHIOverlay::SampleMCEvents(){

  fJetsArr.clear();

  TRandom3 *r1 = new TRandom3(0);
  int filenumber = r1->Integer(filenamesforHIOverlay.size());

  cout << filenumber << "\t" << filenamesforHIOverlay[filenumber] << endl;

  TFile *f = new TFile(filenamesforHIOverlay[filenumber].Data());

  fMCPico = (TTree *)f->Get("PicoDst");
  // fMCPico->SetDirectory(0);

  double pi0mass = Pico::mMass[0]; // GeV

  int nentries = fMCPico->GetEntriesFast();
  
  if (fPrintLevel) cout << "Entries in PicoDst = " << nentries << endl;

  // fMCPico->Print();

  TRandom3 *r2 = new TRandom3(0);

  vector<int> randomlisttoevents;

  randomlisttoevents.clear();

  // for (int i = 0; i < fSetNumberOfEvents; i++){ //This has been FIXED to be able to automate the number of MC events for overlay studies
  //   randomlisttoevents.push_back(r2->Integer(nentries));
  // }

  int numberofmceventstoconvulate = fSetNumberOfEvents;

  if (randomevents){
    if (centralitybinforefficiency == 0) numberofmceventstoconvulate = 4*fSetNumberOfEvents;
    else if (centralitybinforefficiency == 1) numberofmceventstoconvulate = fSetNumberOfEvents;
    else if (centralitybinforefficiency == 2) numberofmceventstoconvulate = 3*fSetNumberOfEvents;
  }
  else{
    numberofmceventstoconvulate = fSetNumberOfEvents;
  }

  // for (int i = 0; i < fSetNumberOfEvents; i++){
  //   randomlisttoevents.push_back(i);
  // }

  

    // for (int i = 0; i < fSetNumberOfEvents; i++){cout << randomlisttoevents[i] << endl;}

    // Structure from PicoDst that we need to mimic its functionalities.
    // I reckon this will be the most memory intensive process of the whole operation.
    // There is significant I/O happening here.

    // Investigate ways to cut I/O here.

  fMCPico->SetMakeClass(1);

  // Fixed size dimensions of array or collections stored in the TTree if any.
  static constexpr Int_t kMaxEvent = 1;
  static constexpr Int_t kMaxTrack = 300;
  static constexpr Int_t kMaxEmcTrigger = 89;
  static constexpr Int_t kMaxMtdTrigger = 1;
  static constexpr Int_t kMaxBTowHit = 4800;
  static constexpr Int_t kMaxBTofHit = 83;
  static constexpr Int_t kMaxMtdHit = 4;
  static constexpr Int_t kMaxBbcHit = 1;
  static constexpr Int_t kMaxEpdHit = 1;
  static constexpr Int_t kMaxFmsHit = 1;
  static constexpr Int_t kMaxEmcPidTraits = 7;
  static constexpr Int_t kMaxBTofPidTraits = 31;
  static constexpr Int_t kMaxMtdPidTraits = 1;
  static constexpr Int_t kMaxTrackCovMatrix = 107;
  static constexpr Int_t kMaxBEmcSmdEHit = 1;
  static constexpr Int_t kMaxBEmcSmdPHit = 1;
  static constexpr Int_t kMaxETofHit = 1;
  static constexpr Int_t kMaxETofPidTraits = 1;
  static constexpr Int_t kMaxMcVertex = 300;
  static constexpr Int_t kMaxMcTrack = 800;

  // Declaration of leaf types
  Int_t           Event_;
  Int_t           Event_mRunId[kMaxEvent];   //[Event_]
  Int_t           Event_mEventId[kMaxEvent];   //[Event_]
  Int_t           Event_mPrimaryVertexX[kMaxEvent];
  Int_t           Event_mPrimaryVertexY[kMaxEvent];
  Int_t           Event_mPrimaryVertexZ[kMaxEvent];

  Int_t           Track_;
  Float_t         Track_mGMomentumX[kMaxTrack];   //[Track_]
  Float_t         Track_mGMomentumY[kMaxTrack];   //[Track_]
  Float_t         Track_mGMomentumZ[kMaxTrack];   //[Track_]
  Float_t         Track_mOriginX[kMaxTrack];   //[Track_]
  Float_t         Track_mOriginY[kMaxTrack];   //[Track_]
  Float_t         Track_mOriginZ[kMaxTrack];   //[Track_]
 
  Char_t          Track_mNHitsFit[kMaxTrack];   //[Track_]
  UChar_t         Track_mNHitsMax[kMaxTrack];   //[Track_]
  
  Short_t         Track_mNSigmaPion[kMaxTrack];   //[Track_]
  Short_t         Track_mNSigmaKaon[kMaxTrack];   //[Track_]
  Short_t         Track_mNSigmaProton[kMaxTrack];   //[Track_]
  Short_t         Track_mNSigmaElectron[kMaxTrack];   //[Track_]
  Short_t         Track_mBEmcMatchedTowerIndex[kMaxTrack];   //[Track_]
  UShort_t        Track_mIdTruth[kMaxTrack];   //[Track_]

  Int_t           BTowHit_;
  Short_t         BTowHit_mE[kMaxBTowHit];   //[BTowHit_]
  
  Int_t           McVertex_;
  Int_t           McVertex_mId[kMaxMcVertex];   //[McVertex_]
  UShort_t        McVertex_mNoDaughters[kMaxMcVertex];   //[McVertex_]
  Int_t           McVertex_mIdParTrk[kMaxMcVertex];   //[McVertex_]
  Int_t           McVertex_mIsInterm[kMaxMcVertex];   //[McVertex_]
  Float_t         McVertex_mTime[kMaxMcVertex];   //[McVertex_]
  Float_t         McVertex_mVx[kMaxMcVertex];   //[McVertex_]
  Float_t         McVertex_mVy[kMaxMcVertex];   //[McVertex_]
  Float_t         McVertex_mVz[kMaxMcVertex];   //[McVertex_]
  Int_t           McTrack_;
  UShort_t        McTrack_mId[kMaxMcTrack];   //[McTrack_]
  Int_t           McTrack_mGePid[kMaxMcTrack];   //[McTrack_]
  Char_t          McTrack_mCharge[kMaxMcTrack];   //[McTrack_]
  UChar_t         McTrack_mHits[kMaxMcTrack][22];   //[McTrack_]
  Float_t         McTrack_mPx[kMaxMcTrack];   //[McTrack_]
  Float_t         McTrack_mPy[kMaxMcTrack];   //[McTrack_]
  Float_t         McTrack_mPz[kMaxMcTrack];   //[McTrack_]
  Float_t         McTrack_mE[kMaxMcTrack];   //[McTrack_]
  Bool_t          McTrack_mIsFromShower[kMaxMcTrack];   //[McTrack_]
  Short_t         McTrack_mIdVtxStart[kMaxMcTrack];   //[McTrack_]
  Short_t         McTrack_mIdVtxStop[kMaxMcTrack];   //[McTrack_]
  Short_t         McTrack_mIdVtxItrmd[kMaxMcTrack];   //[McTrack_]

  // List of branches
  TBranch        *b_Event_;   //!
  TBranch        *b_Event_mRunId;   //!
  TBranch        *b_Event_mEventId;   //!
  TBranch        *b_Event_mPrimaryVertexX;   //!
  TBranch        *b_Event_mPrimaryVertexY;   //!
  TBranch        *b_Event_mPrimaryVertexZ;   //!
  
  TBranch        *b_Track_;   //!
  TBranch        *b_Track_mGMomentumX;   //!
  TBranch        *b_Track_mGMomentumY;   //!
  TBranch        *b_Track_mGMomentumZ;   //!
  TBranch        *b_Track_mOriginX;   //!
  TBranch        *b_Track_mOriginY;   //!
  TBranch        *b_Track_mOriginZ;   //!
  
  TBranch        *b_Track_mNHitsFit;   //!
  TBranch        *b_Track_mNHitsMax;   //!
  TBranch        *b_Track_mNSigmaPion;   //!
  TBranch        *b_Track_mNSigmaKaon;   //!
  TBranch        *b_Track_mNSigmaProton;   //!
  TBranch        *b_Track_mNSigmaElectron;   //!
  TBranch        *b_Track_mTopologyMap;   //!
  TBranch        *b_Track_mBEmcMatchedTowerIndex;   //!
  TBranch        *b_Track_mIdTruth;   //!

  TBranch        *b_BTowHit_;   //!
  TBranch        *b_BTowHit_mE;   //!

  TBranch        *b_McVertex_;   //!
  TBranch        *b_McVertex_mId;   //!
  TBranch        *b_McVertex_mNoDaughters;   //!
  TBranch        *b_McVertex_mIdParTrk;   //!
  TBranch        *b_McVertex_mIsInterm;   //!
  TBranch        *b_McVertex_mTime;   //!
  TBranch        *b_McVertex_mVx;   //!
  TBranch        *b_McVertex_mVy;   //!
  TBranch        *b_McVertex_mVz;   //!
  TBranch        *b_McTrack_;   //!
  TBranch        *b_McTrack_mId;   //!
  TBranch        *b_McTrack_mGePid;   //!
  TBranch        *b_McTrack_mCharge;   //!
  TBranch        *b_McTrack_mHits;   //!
  TBranch        *b_McTrack_mPx;   //!
  TBranch        *b_McTrack_mPy;   //!
  TBranch        *b_McTrack_mPz;   //!
  TBranch        *b_McTrack_mE;   //!
  TBranch        *b_McTrack_mIsFromShower;   //!
  TBranch        *b_McTrack_mIdVtxStart;   //!
  TBranch        *b_McTrack_mIdVtxStop;   //!
  TBranch        *b_McTrack_mIdVtxItrmd;   //!

  fMCPico->SetBranchStatus("*", false);
  fMCPico->SetBranchStatus("Event", true);
  fMCPico->SetBranchStatus("Event.mRunId", true);
  fMCPico->SetBranchStatus("Event.mEventId", true);

  fMCPico->SetBranchStatus("Event.mPrimaryVertexX", true);
  fMCPico->SetBranchStatus("Event.mPrimaryVertexY", true);
  fMCPico->SetBranchStatus("Event.mPrimaryVertexZ", true);

  fMCPico->SetBranchStatus("Track", true);
  fMCPico->SetBranchStatus("Track.mGMomentumX", true);
  fMCPico->SetBranchStatus("Track.mGMomentumY", true);
  fMCPico->SetBranchStatus("Track.mGMomentumZ", true);
  fMCPico->SetBranchStatus("Track.mOriginX", true);
  fMCPico->SetBranchStatus("Track.mOriginY", true);
  fMCPico->SetBranchStatus("Track.mOriginZ", true);
  fMCPico->SetBranchStatus("Track.mNHitsFit", true);
  fMCPico->SetBranchStatus("Track.mNHitsMax", true);
  fMCPico->SetBranchStatus("Track.mNSigmaPion", true);
  fMCPico->SetBranchStatus("Track.mNSigmaKaon", true);
  fMCPico->SetBranchStatus("Track.mNSigmaProton", true);
  fMCPico->SetBranchStatus("Track.mNSigmaElectron", true);
  fMCPico->SetBranchStatus("Track.mBEmcMatchedTowerIndex", true);
  fMCPico->SetBranchStatus("Track.mIdTruth", true);

  fMCPico->SetBranchStatus("BTowHit", true);
  fMCPico->SetBranchStatus("BTowHit", true);

  fMCPico->SetBranchStatus("McVertex", true);
  fMCPico->SetBranchStatus("McVertex.mId", true);
  fMCPico->SetBranchStatus("McVertex.mNoDaughters", true);
  fMCPico->SetBranchStatus("McVertex.mIdParTrk", true);
  fMCPico->SetBranchStatus("McVertex.mIsInterm", true);
  fMCPico->SetBranchStatus("McVertex.mTime", true);
  fMCPico->SetBranchStatus("McVertex.mVx", true);
  fMCPico->SetBranchStatus("McVertex.mVy", true);
  fMCPico->SetBranchStatus("McVertex.mVz", true);
  fMCPico->SetBranchStatus("McTrack", true);
  fMCPico->SetBranchStatus("McTrack.mId", true);
  fMCPico->SetBranchStatus("McTrack.mGePid", true);
  fMCPico->SetBranchStatus("McTrack.mCharge", true);
  fMCPico->SetBranchStatus("McTrack.mHits[22]", true);
  fMCPico->SetBranchStatus("McTrack.mPx", true);
  fMCPico->SetBranchStatus("McTrack.mPy", true);
  fMCPico->SetBranchStatus("McTrack.mPz", true);
  fMCPico->SetBranchStatus("McTrack.mE", true);
  fMCPico->SetBranchStatus("McTrack.mIsFromShower", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxStart", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxStop", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxItrmd", true);


  // fMCPico->SetBranchStatus("Event", true);

  fMCPico->SetBranchAddress("Event", &Event_, &b_Event_);
  fMCPico->SetBranchAddress("Event.mRunId", Event_mRunId, &b_Event_mRunId);
  fMCPico->SetBranchAddress("Event.mEventId", Event_mEventId, &b_Event_mEventId);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexX", Event_mPrimaryVertexX, &b_Event_mPrimaryVertexX);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexY", Event_mPrimaryVertexY, &b_Event_mPrimaryVertexY);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexZ", Event_mPrimaryVertexZ, &b_Event_mPrimaryVertexZ);
  
  fMCPico->SetBranchAddress("Track", &Track_, &b_Track_);
  fMCPico->SetBranchAddress("Track.mGMomentumX", Track_mGMomentumX, &b_Track_mGMomentumX);
  fMCPico->SetBranchAddress("Track.mGMomentumY", Track_mGMomentumY, &b_Track_mGMomentumY);
  fMCPico->SetBranchAddress("Track.mGMomentumZ", Track_mGMomentumZ, &b_Track_mGMomentumZ);
  fMCPico->SetBranchAddress("Track.mOriginX", Track_mOriginX, &b_Track_mOriginX);
  fMCPico->SetBranchAddress("Track.mOriginY", Track_mOriginY, &b_Track_mOriginY);
  fMCPico->SetBranchAddress("Track.mOriginZ", Track_mOriginZ, &b_Track_mOriginZ);
  fMCPico->SetBranchAddress("Track.mNHitsFit", Track_mNHitsFit, &b_Track_mNHitsFit);
  fMCPico->SetBranchAddress("Track.mNHitsMax", Track_mNHitsMax, &b_Track_mNHitsMax);
  fMCPico->SetBranchAddress("Track.mNSigmaPion", Track_mNSigmaPion, &b_Track_mNSigmaPion);
  fMCPico->SetBranchAddress("Track.mNSigmaKaon", Track_mNSigmaKaon, &b_Track_mNSigmaKaon);
  fMCPico->SetBranchAddress("Track.mNSigmaProton", Track_mNSigmaProton, &b_Track_mNSigmaProton);
  fMCPico->SetBranchAddress("Track.mNSigmaElectron", Track_mNSigmaElectron, &b_Track_mNSigmaElectron);
  fMCPico->SetBranchAddress("Track.mBEmcMatchedTowerIndex", Track_mBEmcMatchedTowerIndex, &b_Track_mBEmcMatchedTowerIndex);
  fMCPico->SetBranchAddress("Track.mIdTruth", Track_mIdTruth, &b_Track_mIdTruth);
  
  fMCPico->SetBranchAddress("BTowHit", &BTowHit_, &b_BTowHit_);
  fMCPico->SetBranchAddress("BTowHit.mE", BTowHit_mE, &b_BTowHit_mE);
  
  fMCPico->SetBranchAddress("McVertex", &McVertex_, &b_McVertex_);
  fMCPico->SetBranchAddress("McVertex.mId", McVertex_mId, &b_McVertex_mId);
  fMCPico->SetBranchAddress("McVertex.mNoDaughters", McVertex_mNoDaughters, &b_McVertex_mNoDaughters);
  fMCPico->SetBranchAddress("McVertex.mIdParTrk", McVertex_mIdParTrk, &b_McVertex_mIdParTrk);
  fMCPico->SetBranchAddress("McVertex.mIsInterm", McVertex_mIsInterm, &b_McVertex_mIsInterm);
  fMCPico->SetBranchAddress("McVertex.mTime", McVertex_mTime, &b_McVertex_mTime);
  fMCPico->SetBranchAddress("McVertex.mVx", McVertex_mVx, &b_McVertex_mVx);
  fMCPico->SetBranchAddress("McVertex.mVy", McVertex_mVy, &b_McVertex_mVy);
  fMCPico->SetBranchAddress("McVertex.mVz", McVertex_mVz, &b_McVertex_mVz);
  
  fMCPico->SetBranchAddress("McTrack", &McTrack_, &b_McTrack_);
  fMCPico->SetBranchAddress("McTrack.mId", McTrack_mId, &b_McTrack_mId);
  fMCPico->SetBranchAddress("McTrack.mGePid", McTrack_mGePid, &b_McTrack_mGePid);
  fMCPico->SetBranchAddress("McTrack.mCharge", McTrack_mCharge, &b_McTrack_mCharge);
  fMCPico->SetBranchAddress("McTrack.mHits[22]", McTrack_mHits, &b_McTrack_mHits);
  fMCPico->SetBranchAddress("McTrack.mPx", McTrack_mPx, &b_McTrack_mPx);
  fMCPico->SetBranchAddress("McTrack.mPy", McTrack_mPy, &b_McTrack_mPy);
  fMCPico->SetBranchAddress("McTrack.mPz", McTrack_mPz, &b_McTrack_mPz);
  fMCPico->SetBranchAddress("McTrack.mE", McTrack_mE, &b_McTrack_mE);
  fMCPico->SetBranchAddress("McTrack.mIsFromShower", McTrack_mIsFromShower, &b_McTrack_mIsFromShower);
  fMCPico->SetBranchAddress("McTrack.mIdVtxStart", McTrack_mIdVtxStart, &b_McTrack_mIdVtxStart);
  fMCPico->SetBranchAddress("McTrack.mIdVtxStop", McTrack_mIdVtxStop, &b_McTrack_mIdVtxStop);
  fMCPico->SetBranchAddress("McTrack.mIdVtxItrmd", McTrack_mIdVtxItrmd, &b_McTrack_mIdVtxItrmd);


  int counter = 0;
  int eventsconsidered = 0;

  vector<int> vertexids;
  vector<int> pionids;
  vector<int> kaonids;

  vertexids.clear();
  pionids.clear();
  kaonids.clear();

  while (counter < numberofmceventstoconvulate){

    int eventnum;
    if (randomevents){
      eventnum = r2->Integer(nentries);
      if (std::find(randomlisttoevents.begin(), randomlisttoevents.end(), eventnum) != randomlisttoevents.end()) continue;
    }

    else{
      eventnum = counter;
    }

    randomlisttoevents.push_back(eventnum);
    // cout << "Size of randomlisttoevents = " << randomlisttoevents.size() << endl;

    if (randomlisttoevents.size() < nentries){

      // fMCPico->LoadTree(eventnum);
      fMCPico->GetEntry(eventnum);
      // fMCPico->GetEntry(randomlisttoevents[eventsconsidered]);

      if (fPrintLevel) cout << "Event Number = " << Event_mEventId[0] << "\t" << counter << endl;

      // MC Entries here
      
      vertexids.clear();
      pionids.clear();
      kaonids.clear();

      // cout << " =================== " << counter << "\t" << Event_mEventId[0] << "\t" << McTrack_ << "===================" << endl;

      // Saving the D0 vertices
      for (int mc = 0; mc < McTrack_; mc++){
        if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38){
          TLorentzVector v;
          v.SetPxPyPzE(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], McTrack_mE[mc]);
          // track variables
          double pt = v.Perp();
          double phi = v.Phi();
          if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
          if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
          double eta = v.PseudoRapidity();
          if (pt < 1.0 || pt > 10.0) continue; // Only D0s we care about need to have pT > 1 GeV
          if ((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax) || (phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax)) continue;

          bool goodvertex = kTRUE;

          int pionid = -99;
          int kaonid = -99;

          vector <int> tmp;
          tmp.clear();

          // cout << "D0 pT = " << pt << "\t" << eta << "\t" << phi << endl;

          for(int daug = 0; daug < McTrack_; daug++){
            if (McTrack_mIdVtxStart[daug] != McTrack_mIdVtxStop[mc]) continue; // We only want the kaon and pion that originated from the D0 we are interested in.

            TLorentzVector v2;
            v2.SetPxPyPzE(McTrack_mPx[daug], McTrack_mPy[daug], McTrack_mPz[daug], McTrack_mE[daug]);
            double pt2 = v2.Perp();
            double phi2 = v2.Phi();
            if(phi2 < 0.0)    phi2 += 2.0*pi;  // force from 0-2pi
            if(phi2 > 2.0*pi) phi2 -= 2.0*pi;  // force from 0-2pi
            double eta2 = v2.PseudoRapidity();

            // cout << "Daug = " << McTrack_mGePid[daug] << "\t" << pt2 << "\t" << eta2 << "\t" << phi2 << endl;

            // cout << "Daug daughters = " << McVertex_mNoDaughters[McTrack_mIdVtxStop[daug]] << endl;

            // if (McVertex_mNoDaughters[McTrack_mIdVtxStop[daug]] < 5){
            //   for(int daug2 = 0; daug2 < McTrack_; daug2++){
            //     if (McTrack_mIdVtxStart[daug2] != McTrack_mIdVtxStop[daug]) continue;

            //     TLorentzVector v3;
            //     v3.SetPxPyPzE(McTrack_mPx[daug2], McTrack_mPy[daug2], McTrack_mPz[daug2], McTrack_mE[daug2]);
            //     double pt3 = v3.Perp();
            //     double phi3 = v3.Phi();
            //     if(phi3 < 0.0)    phi3 += 2.0*pi;  // force from 0-2pi
            //     if(phi3 > 2.0*pi) phi3 -= 2.0*pi;  // force from 0-2pi
            //     double eta3 = v3.PseudoRapidity();

            //     cout << "Daug of daug = " << McTrack_mGePid[daug2] << "\t" << pt3 << "\t" << eta3 << "\t" << phi3 << endl;
            //   }
            // }

            if (McTrack_mGePid[daug] != 8 && McTrack_mGePid[daug] != 9 && McTrack_mGePid[daug] != 11 && McTrack_mGePid[daug] != 12)
            { 
              goodvertex = kFALSE; break;
            }

            if (McTrack_mIdVtxStop[daug] != 0) {goodvertex = kFALSE; break;} // If these are not final state kaon pions, then they are of no use to us. 

            // Push the pion and kaon ids into the vector
            if (McTrack_mGePid[daug] == 8 || McTrack_mGePid[daug] == 9){
              pionid = daug;
            }
            
            else if (McTrack_mGePid[daug] == 11 || McTrack_mGePid[daug] == 12){
              kaonid = daug;
            }
          
          }

          // cout << "Vertex Stop = " << McTrack_mIdVtxStop[mc] << endl;
          if (goodvertex){
            if (pionid == -99 || kaonid == -99) { cout << "Something wrong with D0. Exiting." << endl; continue; }
            vertexids.push_back(McTrack_mIdVtxStop[mc]);
            pionids.push_back(pionid);
            kaonids.push_back(kaonid);
            // cout << "=============== D0 ================ " << pt << "\t" << eta << "\t" << phi << endl;
          } 
        }
      }

      // cout << "=============== Number of D0s =============== " << vertexids.size() << "\t" << pionids.size() << endl;

      if (vertexids.size() == 0) continue; //While loop continues.

      int D0Counter = 0;


      // Loop over each D0
      for (int D0 = 0; D0 < vertexids.size(); D0++){

        // cout << "This is D0 # " << D0 << "***********************" << endl;

        // I only want to use TLorentzVector to propagate the information ahead. All the processing is done within this class.
        // This means I need to find a way to make sure I can identify the D0 track when I save out the information.
        // Since only D0 needs to be identified, I am saving the mass information for the D0 as 1.865.
        // All the other tracks are saved out with pi0mass, because ultimately, jet constituents are assumed to have that mass.

        // How do I propagate charge information though?
        // Do I need it? Mostly, nope! In fact, once I separate track and tower, it should be enough.

        // This loop fills the input vector for the MC Side

        // cout << "Number of MC tracks = " << McTrack_ << endl;

        for (int mc = 0; mc < McTrack_; mc++){
          TLorentzVector v;
          if ((McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38)) v.SetXYZM(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], 1.865); // This is a D0
          else v.SetXYZM(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], pi0mass);

          // track variables
          double pt = v.Perp();
          double phi = v.Phi();
          if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
          if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
          double eta = v.PseudoRapidity();

          // if (McTrack_mGePid[mc] == 8 || McTrack_mGePid[mc] == 9 || McTrack_mGePid[mc] == 11 || McTrack_mGePid[mc] == 12){
          //   cout << "Vertex Start Stop = " << McTrack_mIdVtxStart[mc] << "\t" << McTrack_mIdVtxStop[mc] << endl;
          // }

          if (McTrack_mIdVtxStop[mc] != 0 && (McTrack_mGePid[mc] != 37 && McTrack_mGePid[mc] != 38)) continue; // Only Final State Particles are included in JetMaker. D0s are the only exception

          // if (ignoretracks[D0].size()!=0){ if (std::find(ignoretracks[D0].begin(), ignoretracks[D0].end(), mc) != ignoretracks[D0].end()) continue; } // Ignoring the particles which decayed from K & Pi which came from the D0

          // if ((McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38) && pt < 1.0) continue; // We don't need D0s < 1 GeV
          // The previous check is redundant and hence done away with.
          if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38) {
            if (McTrack_mIdVtxStop[mc] != vertexids[D0]) continue; // Only consider the current D0
          }

          // if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38){D0Counter++;}
          // Unstable particles which shouldn't make it to the end are discarded by hand. The list provisionally includes:
          /*
            Lambda, Eta, Sigma0, Xi0, Muon, Neutrino, KS0, KL0
          */

          if (McTrack_mGePid[mc] == 4 || McTrack_mGePid[mc] == 5 || McTrack_mGePid[mc] == 6 || McTrack_mGePid[mc] == 10 || McTrack_mGePid[mc] == 16 || McTrack_mGePid[mc]== 17 || McTrack_mGePid[mc] == 18 || McTrack_mGePid[mc] == 20 || McTrack_mGePid[mc] == 22 ) continue;

          // Have to discard Kaons and Pions which come from the Current D0

          if (McTrack_mIdVtxStart[mc] == vertexids[D0]){
            continue;
          }

          if ((pt < fMinJetTrackPt) || (pt > fMaxJetTrackPt) || (eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax) || (phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax)) continue;
          if (McTrack_mCharge[mc] != 0 || McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38) {
            if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38) {
              D0Counter++; 
              // cout << "D0 = " << pt << "\t" << eta << "\t" << phi << endl;
            }
            fMCEventTracks[counter].push_back(v); //Neutral particles which are not D0 are saved in towers.
            // cout << "Constituents = " << McTrack_mGePid[mc] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
          }
          else {
            fMCEventTowers[counter].push_back(v);
            // cout << "Towers = " << McTrack_mGePid[mc] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
          } 
        }

        // Make the jets here. If they are not what we want them to be, discard the whole event and keep going.

        StFJWrapper             *fjw = new StFJWrapper("HIMC", "HIMC"); //!fastjet wrapper
        fjw->Clear();
        fJets->Delete();

        // setup fj wrapper
        fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
        fastjet::JetAlgorithm           algorithm = fastjet::antikt_algorithm;
        fastjet::Strategy               strategy = fastjet::Best;

        fjw->SetAreaType(fastjet::active_area_explicit_ghosts);
        fjw->SetStrategy(strategy);
        fjw->SetGhostArea(0.005);
        fjw->SetR(0.4);
        fjw->SetAlgorithm(algorithm);        //fJetAlgo);
        fjw->SetRecombScheme(recombScheme);  //fRecombScheme);
        fjw->SetMaxRap(1.2);

        
        // Loop over all saved particles for MC in StHIOverlay vectors and enter them into fastjet wrapper 
        for(unsigned short iTracks = 0; iTracks < fMCEventTracks[counter].size(); iTracks++){
           
          TLorentzVector v;
          v = fMCEventTracks[counter][iTracks];
          // track variables
          
          double px = v.X();
          double py = v.Y();
          double pz = v.Z();
          double p = v.P();
          double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
          double mass = v.M();

          fjw->AddInputVector(px, py, pz, energy, iTracks ); // includes E

          double pt = v.Pt();
          double phi = v.Phi();
          if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
          if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
          double eta = v.Eta();

          // if (int(mass)==1){  cout << "D0 is in Input loop = " << mass << endl; }

          // cout << "Input = " << pt << "\t" << eta << "\t" << phi << endl;

        } // track loop

        for (unsigned short iTowers = 0; iTowers < fMCEventTowers[counter].size(); iTowers++){

          TLorentzVector v;
          v = fMCEventTowers[counter][iTowers];
          // tower variables
          
          double px = v.X();
          double py = v.Y();
          double pz = v.Z();
          double p = v.P();
          double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
          double mass = v.M();

          fjw->AddInputVector(px, py, pz, energy, iTowers+10000); // includes E

          double pt = v.Pt();
          double phi = v.Phi();
          if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
          if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
          double eta = v.Eta();

          // cout << "Input = " << pt << "\t" << eta << "\t" << phi << endl;

        } // tower loop

        // Jet Running

        fjw->Run();

        std::vector<fastjet::PseudoJet> jets_incl = fjw->GetInclusiveJets();

        // sort jets according to jet pt
        static Int_t indexes[9999] = {-1};
        GetSortedArray(indexes, jets_incl);

        int jetCount = 0; // This is the index of the D0 jet in the jet tree.
        bool acceptedevent = kFALSE;

        for(UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet) {

          Int_t ij = indexes[ijet];

          // PERFORM CUTS ON inclusive JETS before saving
          // cut on min jet pt
          if(jets_incl[ij].perp() < 5.0) continue;
          // cut on eta acceptance
          if((jets_incl[ij].eta() < -0.6) || (jets_incl[ij].eta() > 0.6)) continue;

          // fill jet constituents
          vector<fastjet::PseudoJet> constituents = fjw->GetJetConstituents(ij);

          bool D0Jet = kFALSE;

          for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
            // get user defined index
            Int_t uid = constituents[ic].user_index();

            if(int(uid) >= 0 && int(uid) < 10000 ){
              TLorentzVector v;
              v = fMCEventTracks[counter][uid];

              if (int(v.M())==1) {
                D0Jet = kTRUE; 
                break;
              }
            }
          }

          if (!D0Jet) continue;

          StJet *jet = new ((*fJets)[jetCount])
            StJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());

          jet->SetLabel(ij);

          // area vector and components
          fastjet::PseudoJet area(fjw->GetJetAreaVector(ij));
          jet->SetArea(area.perp());  // same as fjw->GetJetArea(ij)
          jet->SetAreaEta(area.eta());
          jet->SetAreaPhi(area.phi());
          jet->SetAreaE(area.E());

          jet->SetJetConstituents(constituents);

          int nt = 0;
          int nc = 0;
          int ng = 0;

          jet->SetNumberOfTracks(constituents.size());
          jet->SetNumberOfTowers(constituents.size());

          for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
            Int_t uid = constituents[ic].user_index();

            if(int(uid) >= 0 && int(uid) < 10000 ){
              jet->AddTrackAt(abs(uid), nt);
              // TLorentzVector v;
              // v = fMCEventTracks[counter][uid];
              // if (int(v.M())==1) {
              //   double mcjetpx = jet->Px();
              //   double mcjetpy = jet->Py();
              //   double mcjetpt = jet->Pt();
              //   double mcD0px = v.Px();
              //   double mcD0py = v.Py();
              //   double mcD0pt = v.Pt();
              //   double z = (mcjetpx*mcD0px + mcjetpy*mcD0py)/(pow(mcjetpt, 2));
              //   cout << "Found jet with pT = " << mcjetpt << " D0 pT \t" << mcD0pt << " Z = " << z << endl;
              // }
              nt++;
            }
            else if ( int(uid)>=10000 ){
              // cout << uid  << "\t" << nc << endl;
              jet->AddTowerAt(abs(uid), nc);
              nc++;
            }
            else if (uid < 0){
              ng++;
            }
          }
          jet->SetNumberOfTracks(nt);
          jet->SetNumberOfTowers(nc);

          TClonesArray *tempjets = (TClonesArray *)fJets->Clone();
          fJetsArr.push_back(tempjets);

          acceptedevent = kTRUE;

          if (fPrintLevel) cout << "Jet Found with pT eta phi = " << jet->Pt() << "\t" << jet->Eta() << "\t" << jet->Phi() << endl;

          break; //Once you encountered a D0 jet to your liking, break.
        }

        delete fjw;

        if (!acceptedevent) {
          fMCEventTracks[counter].clear();
          fMCEventTowers[counter].clear();
          continue; // This is the continuation for the D0 loop.
        }

        // cout << "MC Jet Array Size = " << fJetsArr.size() << endl;

        // cout << "MC Events Size for counter # " << counter << "\t = "  << fMCEventTracks[counter].size() << "\t" << fMCEventTowers[counter].size() << endl;

        const int numberoftowers = BTowHit_;
        double towerenergy[numberoftowers];
        for (int i = 0; i < numberoftowers; i++){towerenergy[i] = 0.;}

        // This loop fills the input vector for the RECO side
        for (int reco = 0; reco < Track_; reco++){

          TVector3 o;
          o.SetXYZ(Track_mOriginX[reco], Track_mOriginY[reco], Track_mOriginZ[reco]);
          TVector3 g;
          g.SetXYZ(Track_mGMomentumX[reco], Track_mGMomentumY[reco], Track_mGMomentumZ[reco]);

          // track variables
          double pt = g.Perp();
          double phi = g.Phi();
          if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
          if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
          double eta = g.PseudoRapidity();
          double px = g.x();
          double py = g.y();
          double pz = g.z();
          double p = g.Mag();
          double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
          short charge = (Track_mNHitsFit[reco] > 0) ? 1 : -1;

          double dca = (mVertex - o).Mag();
          
          //// Variables For FastSim

          double pt_new = pt;
          double phi_new = phi;
          double eta_new = eta;
          double px_new = px;
          double py_new = py;
          double pz_new = pz;
          double p_new = p;
          double energy_new = energy;
          short charge_new = charge;

          bool mctrackavailable = kTRUE;

          int mcid = Track_mIdTruth[reco];

          if (mcid <= 0) mctrackavailable = kFALSE;

          bool isatrackfromD0 = kFALSE;

          // Here, we have two paths to take. If the track needs replacement, we replace it with the fastsim method that is standardised.
          // Else, the pt, eta, phi are sent as is to the final vector.

          if (mctrackavailable){
            if (McTrack_mIdVtxStart[mcid] == vertexids[D0]) isatrackfromD0 = kTRUE; // Kaons and Pions that come from the current D0 need to be tossed, and replaced by the fast sim version 

            TVector3 mg(McTrack_mPx[mcid], McTrack_mPy[mcid], McTrack_mPz[mcid]);

            double relativesmearing = TMath::Sqrt(pow(mg.Px() - px, 2) + pow(mg.Py() - py, 2))/(mg.Pt());

            double fastsimsmearing;

            int pid = McTrack_mGePid[mcid];

            if     (pid ==8  || pid == 9) fastsimsmearing =  fPionMomResolution->Eval(mg.Pt());// Pion
            else if(pid ==11 || pid == 12)fastsimsmearing =  fKaonMomResolution->Eval(mg.Pt());// Kaon
            else if(pid ==15 || pid == 14)fastsimsmearing =  fProtonMomResolution->Eval(mg.Pt());// Proton
            else fastsimsmearing = fPionMomResolution->Eval(mg.Pt());// Catch all: pions

            if (relativesmearing > 3*fastsimsmearing){
              TVector3 fastsimsmearedmom = FastSimMom(mg, pid);
              pt_new = fastsimsmearedmom.Perp();
              phi_new = fastsimsmearedmom.Phi();
              if(phi_new < 0.0)    phi_new += 2.0*pi;  // force from 0-2pi
              if(phi_new > 2.0*pi) phi_new -= 2.0*pi;  // force from 0-2pi
              eta_new = fastsimsmearedmom.PseudoRapidity();
              px_new = fastsimsmearedmom.x();
              py_new = fastsimsmearedmom.y();
              pz_new = fastsimsmearedmom.z();
              p_new = fastsimsmearedmom.Mag();
              energy_new = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
            }
          }

          // jet track acceptance cuts now
          if(pt_new < fMinJetTrackPt) continue;
          if(pt_new > fMaxJetTrackPt) continue; // 20.0 STAR, 100.0 ALICE
          if((eta_new < fJetTrackEtaMin) || (eta_new > fJetTrackEtaMax)) continue;
          if(phi_new < 0.0)    phi_new += 2.0*pi;  // force from 0-2pi
          if(phi_new > 2.0*pi) phi_new -= 2.0*pi;  // force from 0-2pi
          if((phi_new < fJetTrackPhiMin) || (phi_new > fJetTrackPhiMax)) continue;
              
          // additional quality cuts for tracks
          if(dca > fJetTrackDCAcut)            continue;
          if(abs(Track_mNHitsFit[reco]) < fJetTracknHitsFit)     continue;
          if(abs(Track_mNHitsFit[reco])/Track_mNHitsMax[reco] < fJetTracknHitsRatio) continue;

          // cout << "TRACK = " << px_new << "\t" << py_new << "\t" << pz_new << endl;

          int particleid = -99;

          if      (abs(Track_mNSigmaPion[reco]/1000.) < 2) particleid = 1;
          else if (abs(Track_mNSigmaKaon[reco]/1000.) < 2.) particleid = 2;
          else if (abs(Track_mNSigmaProton[reco]/1000.) < 2.) particleid = 3*charge;

          if (!KeepTrack(particleid, centralitybinforefficiency, pt_new) && !isatrackfromD0) continue; // To match the efficiency, we start tossing random tracks.

          TLorentzVector v;
          v.SetXYZM(px_new, py_new, pz_new, pi0mass);

          if (!isatrackfromD0) fRecoEventTracks[counter].push_back(v);
          //// Fill the Tower Array with Energy Depositions here

          int matchedTowerIndex = abs(GetMatchedBtowID(Track_mBEmcMatchedTowerIndex[reco], g, o, charge_new)) - 1; // towerIndex = towerID - 1
          
          if (matchedTowerIndex >= 0){towerenergy[matchedTowerIndex]+=energy_new;} //This place takes care of all energy depositions due to tracks accepted for jet reco. It also includes the kaon and pion from D0.        
        }

        TVector3 mcKaon(McTrack_mPx[kaonids[D0]], McTrack_mPy[kaonids[D0]], McTrack_mPz[kaonids[D0]]);
        TVector3 mcPion(McTrack_mPx[pionids[D0]], McTrack_mPy[pionids[D0]], McTrack_mPz[pionids[D0]]);

        fMCD0Information[counter] = {mcPion, mcKaon}; 

        TVector3 recoKaon = FastSimMom(mcKaon, 11);
        TVector3 recoPion = FastSimMom(mcPion, 8);

        fRecoD0Information[counter] = {recoPion, recoKaon};

        fOrigin[counter].SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

        // fOrigin[counter].push_back(eventorigin); 

        TVector3 recoD0;
        recoD0 = recoKaon + recoPion;

        TLorentzVector v;
        v.SetXYZM(recoD0.x(), recoD0.y(), recoD0.z(), 1.865); // This is a D0

        fRecoEventTracks[counter].push_back(v);

        // This loop fills the input vector for the towers

        for (int tower = 0; tower < numberoftowers; tower++){
          int towerID = tower + 1;
          if(towerID < 0) continue; // double check these aren't still in the event list

          TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
          double towerPhi = towerPosition.Phi();
          if(towerPhi < 0.0)    towerPhi += 2.0*pi;  // force from 0-2pi
          if(towerPhi > 2.0*pi) towerPhi -= 2.0*pi;  // force from 0-2pi
          double towerEta = towerPosition.PseudoRapidity();

          // check for bad (and dead) towers
          bool TowerOK = mBaseMaker->IsTowerOK(towerID);      // kTRUE means GOOD
          bool TowerDead = mBaseMaker->IsTowerDead(towerID);  // kTRUE means BAD
          if(!TowerOK)  { continue; }
          if(TowerDead) { continue; }

          // jet track acceptance cuts njow
          if((towerEta < fJetTowerEtaMin) || (towerEta > fJetTowerEtaMax)) continue;
          if((towerPhi < fJetTowerPhiMin) || (towerPhi > fJetTowerPhiMax)) continue;

          double towerEunCorr = BTowHit_mE[tower]/1000.;  // uncorrected energy
          double towerE = BTowHit_mE[tower]/1000.;        // corrected energy (hadronically - done below)
          double towEtunCorr = towerE / (1.0*TMath::CosH(towerEta));

          // cut on min tower energy after filling histos
          if(towerEunCorr < mTowerEnergyMin) continue; // if we don't have enough E to start with, why mess around

          // =======================================================================
          // HADRONIC CORRECTION
          
          double sumEt = (towerEunCorr - towerenergy[tower])/(1.0*TMath::CosH(towerEta));
          double towerEt = sumEt;
          if(towerEt < mTowerEnergyMin) continue;
          towerE = towerEt  * 1.0*TMath::CosH(towerEta);

          Double_t p = 1.0*TMath::Sqrt(towerE*towerE - pi0mass*pi0mass);

          double posX = towerPosition.x();
          double posY = towerPosition.y();
          double posZ = towerPosition.z();

          Double_t r = TMath::Sqrt(posX*posX + posY*posY + posZ*posZ) ;

          TLorentzVector v;
          v.SetXYZM(p*posX/r, p*posY/r, p*posZ/r, pi0mass);

          // cout << "TOWER = " << v.X() << "\t" << v.Y() << "\t" << v.Z() << endl;

          fRecoEventTowers[counter].push_back(v);
        }

        // cout << "Reco Events Size for counter # " << counter << "\t = " << fRecoEventTracks[counter].size() << "\t" << fRecoEventTowers[counter].size() << endl;

        counter++; // Since one MC event can be invoked multiple times (due to having multiple D0s, this step is necessary.)
      }
      eventsconsidered++;
    // }
    }

    else break;
  }

  // cout << "Counter = " << counter << endl;
  fNumberofeventsoverlayed = counter;

  fMCPico->Reset();
  f->Close();

}

// vector <int> StHIOverlay::TracksToToss(int stopvertex, int *)

int StHIOverlay::GetMatchedBtowID(int trkbemcid, TVector3 gMom, TVector3 org, int charge){
  Double_t bemc_radius = mBemcGeom->Radius();
  // Magnetic field in Tesla 
  Double_t mBField_tesla = Bfield / 10.0; //Check this definition. Magnetic fields are minefields of error in STAR

  // Needed for projection of the track onto the barrel radius
  TVector3 bemc_pos, bemc_mom;

  // BEMC hardware indices 
  Int_t h_m, h_e, h_s = 0;

  // tower index: if no tower can be matched, assign 0
  // picoTrk->setBEmcMatchedTowerIndex(0);
  Int_t tow_id = 0;
  Bool_t close_match = false;

  // int trkbemcid = trk->bemcTowerIndex();

  // Check if the track can be projected onto the current radius
  // if not, track can't be matched.
  // By JetCorr request the global track projection to BEMC is used.
  if ( mEmcPosition->projTrack(&bemc_pos, &bemc_mom, gMom, org, Bfield, charge, bemc_radius) ) {
    // First, examine track eta. If it falls in two regions:
    // 0 < |eta| < etaMin()
    // etaMax() < |eta| < 1.0
    // then shift the eta for the projection slightly into the neighboring tower
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(org, trkbemcid + 1);
    // cout << bemc_pos.Phi() << "\t" << towerPosition.Phi() << "\t" << bemc_pos.PseudoRapidity() << "\t" << towerPosition.PseudoRapidity() << endl;

    if ( fabs(bemc_pos.PseudoRapidity()) < mBemcGeom->EtaMin() ) {
      Double_t unsigned_eta = mBemcGeom->EtaMin() + 0.001;
      // Double_t unsigned_eta = mBemcGeom->EtaMin() + 0.000001;
      Double_t unsigned_theta = 2.0 * atan(exp(-1.0 * unsigned_eta));
      Double_t signed_theta = (bemc_pos.PseudoRapidity() >= 0 ? 1.0 : -1.0) * unsigned_theta;
      bemc_pos.SetTheta(signed_theta);
      close_match = true;
    } 
    else if ( fabs(bemc_pos.PseudoRapidity()) > mBemcGeom->EtaMax() &&
      fabs(bemc_pos.PseudoRapidity()) < 1.0 ) {
      Double_t unsigned_eta = mBemcGeom->EtaMax() - 0.001;
      // Double_t unsigned_eta = mBemcGeom->EtaMax() - 0.000001;
      Double_t unsigned_theta = 2.0 * atan(exp(-1.0 * unsigned_eta));
      Double_t signed_theta = (bemc_pos.PseudoRapidity() >= 0 ? 1.0 : -1.0) * unsigned_theta;
      bemc_pos.SetTheta(signed_theta);
      close_match = true;
    }


    // Get the BEMC hardware location in (m, e, s) and translate to id
    // If StEmcGeom::getBin() != 0: track was not matched to a tower.
    // Its outside of the BEMC eta range (> 1.0).

    if ( mBemcGeom->getBin(bemc_pos.Phi(),bemc_pos.PseudoRapidity(),h_m,h_e,h_s) == 0 ) {
      // If StEmcGeom::getId() == 0: the track was matched successfully. Otherwise, 
      // the track was not matched to a tower at this radius, the track was projected
      // into the gap between modules in phi. 
      if ( h_s != -1 ) {
        mBemcGeom->getId(h_m,h_e,h_s,tow_id);
        if (close_match) {
          return -1*tow_id;
        }

        else{
          return tow_id;
        }
      }

      // Track fell in between modules in phi. We will find which module it is closer
      // to by shifting phi slightly.
      else {
        // Value of the "dead space" per module in phi:
        // 2*pi/60 (amount of azimuth covered per module)
        // 2*0.0495324 (active size of module)

        Double_t dphi = (TMath::Pi() / 60.0) - 0.0495324;
        // Shift the projected phi by dphi in positive and negative directions
        // if we look for the projection for both of these, only one should give
        // a tower id, and the other should still be in the inter-tower space

        TVector3 bemc_pos_shift_pos(bemc_pos); 
        bemc_pos_shift_pos.SetPhi(bemc_pos_shift_pos.Phi() + dphi);
        TVector3 bemc_pos_shift_neg(bemc_pos); 
        bemc_pos_shift_neg.SetPhi(bemc_pos_shift_neg.Phi() - dphi);

        if ( mBemcGeom->getBin(bemc_pos_shift_pos.Phi(),bemc_pos_shift_pos.PseudoRapidity(),h_m,h_e,h_s) == 0 && h_s != -1 ){
          mBemcGeom->getId(h_m,h_e,h_s,tow_id);        
          return -1*tow_id;
        }

        else if( mBemcGeom->getBin(bemc_pos_shift_neg.Phi(),bemc_pos_shift_neg.PseudoRapidity(),h_m,h_e,h_s) == 0 && h_s != -1 ){
          mBemcGeom->getId(h_m,h_e,h_s,tow_id);
          return -1*tow_id;
        } 
      }
    }
  }

  return tow_id;

}


TVector3 StHIOverlay::FastSimMom(TVector3 p, int pid)
{
    float pt = p.Perp();
    float pt1 = pt;
    // if(pt1>2) pt1 = 2;//Used for high pt-hat bin smearing test
    if(pt1>10) pt1 = 10;//Used for high pt-hat bin smearing test
    float sPt = -1;
    
    if(pid ==8 || pid == 9)sPt = gRandom->Gaus(pt, pt * fPionMomResolution->Eval(pt1));// Pion
    else if(pid ==11 || pid == 12)sPt = gRandom->Gaus(pt, pt * fKaonMomResolution->Eval(pt1));// Kaon
    else if(pid ==15 || pid == 14)sPt = gRandom->Gaus(pt, pt * fProtonMomResolution->Eval(pt1));// Proton
    else sPt = gRandom->Gaus(pt, pt * fPionMomResolution->Eval(pt1));// Catch all: pions

    TVector3 smearedmom(sPt * cos(p.Phi()), sPt * sin(p.Phi()), sPt * sinh(p.PseudoRapidity()));
    
    return smearedmom;
}


Bool_t StHIOverlay::KeepTrack(int particleid, int centralitybin, double pt)
{
  bool keeptrack = kTRUE;

  TRandom3 *r = new TRandom3(0);
  double rando = r->Rndm();

  // cout << "Random = " << rando << "\t" << fPionWeight->Eval(pt) << "\t" << fKaonWeight->Eval(pt) << "\t" << fProtonWeight->Eval(pt) << "\t" << fAProtonWeight->Eval(pt) << endl;

  if (particleid == 1 || particleid == -99) keeptrack = (rando > fPionWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE; //Either a pion
  else if (particleid == 2) keeptrack = (rando > fKaonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == 3) keeptrack = (rando > fProtonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == -3) keeptrack = (rando > fAProtonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;

  return keeptrack;
}

Bool_t StHIOverlay::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
{
  static Float_t pt[9999] = {0};
  const Int_t n = (Int_t)array.size();
  if(n < 1) return kFALSE;

  for(Int_t i = 0; i < n; i++)
    pt[i] = array[i].perp();

  TMath::Sort(n, pt, indexes);

  return kTRUE;
}
