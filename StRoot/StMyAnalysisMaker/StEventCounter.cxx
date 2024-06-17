// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StEventCounter.h"
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
#include "TLorentzVector.h"


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
// #include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"

// // KFParticle includes
// #include "StRoot/StEvent/StDcaGeometry.h"
// #include "StRoot/StarRoot/KFParticle.h"

// Bichsel includes
#include "StBichsel/Bichsel.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StEventCounter)

//________________________________________________________________________
StEventCounter::StEventCounter(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{ 
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StEventCounter::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fMCEventsWithoutCent = kFALSE;
  fCentralitySelectionCut = -99;
  fDoEffCorr = kFALSE;
  doRejectBadRuns = kFALSE;
  fTopoLevel = 1;
  fCorrJetPt = kFALSE;
  fMinPtJet = 0.0;
  fTrackBias = 0.0;
  fTowerBias = 0.0;
  fJetRad = 0.4;
  fEventZVtxMinCut = -60.0; fEventZVtxMaxCut = 60.0;
  fMaxEventTrackPt = 50.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 50.0;
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
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

  for(int i=0; i<10; i++) {numberofevents[i] = 0; numberoftracks[i] = 0;}
  //Mass cut
  fInvMassSignal1 = 1.80;
  fInvMassSignal2 = 1.92;

  fInvMassULBg1 = 1.70;
  fInvMassULBg2 = 2.08;

  fInvMassLSBg1 = 1.70;
  fInvMassLSBg2 = 2.08;

  fd0 = kFALSE;
  fd0BgUS = kFALSE;
  fd0BgLS = kFALSE;

  // fd0TrackIndices.clear();
  // fd0BgUSTrackIndices.clear();
  // fd0BgLSTrackIndices.clear();
}

//
//________________________________________________________________________
StEventCounter::~StEventCounter()
{ 

  if (hEventZvertex_whole) delete hEventZvertex_whole;
  if (hEventZvertex_VPD) delete hEventZvertex_VPD;
  if (hEventZvertex_diff) delete hEventZvertex_diff;

  if (hCentrality) delete hCentrality;
  if (hMultiplicity) delete hMultiplicity;
  if (cuthistogram_event) delete cuthistogram_event;

  for (int i = 0; i < 10; i++){
    if (hAcceptedPt[i]) delete hAcceptedPt[i];
    if (hAcceptedEta[i]) delete hAcceptedEta[i];
  }

  if (z_pi) delete z_pi;
  if (z_ka) delete z_ka;
  if (normalised_invbetavpT_tof_pi) delete normalised_invbetavpT_tof_pi;
  if (normalised_invbetavpT_tof_ka) delete normalised_invbetavpT_tof_ka;
  if (dEdXvp) delete dEdXvp;
  if (invbetavp) delete invbetavp;

  if (decaylengthd0US) delete decaylengthd0US;
  if (distancepikUS) delete distancepikUS;
  if (distanced0PVUS) delete distanced0PVUS;
  if (dcakPVUS) delete dcakPVUS;
  if (dcapiPVUS) delete dcapiPVUS;

  if (decaylengthd0LS) delete decaylengthd0LS;
  if (distancepikLS) delete distancepikLS;
  if (distanced0PVLS) delete distanced0PVLS;
  if (dcakPVLS) delete dcakPVLS;
  if (dcapiPVLS) delete dcapiPVLS;

  if (kaonpt) delete kaonpt;
  if (pionpt) delete pionpt;
  if (d0pt) delete d0pt;
  if (kaonpionpt) delete kaonpionpt;
  if (kaoneta) delete kaoneta;
  if (pioneta) delete pioneta;
  if (d0eta) delete d0eta;
  if (kaonphi) delete kaonphi;
  if (pionphi) delete pionphi;
  if (d0phi) delete d0phi;

  if (kaonbgUSpt) delete kaonbgUSpt;
  if (pionbgUSpt) delete pionbgUSpt;
  if (d0bgUSpt) delete d0bgUSpt;
  if (kaonpionbgUSpt) delete kaonpionbgUSpt;
  if (kaonbgUSeta) delete kaonbgUSeta;
  if (pionbgUSeta) delete pionbgUSeta;
  if (d0bgUSeta) delete d0bgUSeta;
  if (kaonbgUSphi) delete kaonbgUSphi;
  if (pionbgUSphi) delete pionbgUSphi;
  if (d0bgUSphi) delete d0bgUSphi;

  if (kaonbgLSpt) delete kaonbgLSpt;
  if (pionbgLSpt) delete pionbgLSpt;
  if (d0bgLSpt) delete d0bgLSpt;
  if (kaonpionbgLSpt) delete kaonpionbgLSpt;
  if (kaonbgLSeta) delete kaonbgLSeta;
  if (pionbgLSeta) delete pionbgLSeta;
  if (d0bgLSeta) delete d0bgLSeta;
  if (kaonbgLSphi) delete kaonbgLSphi;
  if (pionbgLSphi) delete pionbgLSphi;
  if (d0bgLSphi) delete d0bgLSphi;

  if (invmass) delete invmass;
  if (invmassbg) delete invmassbg;

  if (hNumberOfD0s) delete hNumberOfD0s;
  if (hNumberOfD0BgUS) delete hNumberOfD0BgUS;
  if (hNumberOfD0BgLS) delete hNumberOfD0BgLS;


  if(mEmcPosition) delete mEmcPosition;

}

//________________________________________________________________________
Int_t StEventCounter::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StEventCounter::Finish() { 

  cout << "StEventCounter::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),  "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StEventCounter::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StEventCounter::DeclareHistograms() {

  //gPtvpPt = new TH2F("gPtvpPt", "gPtvpPt", 1000, 0, 20, 1000, 0, 20);

  // binning for cent histograms
  int nHistCentBins = 20;

  // binning for mult histograms
  double kHistMultMax = 800.;
  int kHistMultBins = 400;

  // pp specific settings
  if(doppAnalysis) {
    kHistMultMax = 100.;
    kHistMultBins = 100.;
  }

  hEventZvertex_whole = new TH1F("hEventZvertex_whole", "z-vertex distribution for all events", 200, -100., 100.);
  hEventZvertex_VPD = new TH1F("hEventZvertex_VPD", "z-vertex distribution from VPD", 200, -100., 100.);
  hEventZvertex_diff = new TH1F("hEventZvertex_diff", "z-vertex distribution Difference", 200, -100., 100.);

  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);

  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  cuthistogram_event = new TH1F("cuthistogram_event", "cuthistogram_event", 10, 0, 10);

  // hAcceptedEvents = new TH1F("hAcceptedEvents", "hAcceptedEvents", 10000, 0, 10000);

  // Pt After Track Cuts

  for (int i = 0; i < 10; i++){
    hAcceptedPt[i] = new TH1F(Form("hAcceptedPt_%i", i), Form("hAcceptedPt_%i", i), 1000, 0, 30);
    hAcceptedEta[i] = new TH1F(Form("hAcceptedEta%i", i), Form("hAcceptedEta%i", i), 1000, -2, 2);
  }
  

  z_pi = new TH2F("z_pi", "z_pi", 3000, 0, 30, 2000, -10, 10);
  z_ka = new TH2F("z_ka", "z_ka", 3000, 0, 30, 2000, -10, 10);
  normalised_invbetavpT_tof_pi = new TH2F("normalised_invbetavpT_tof_pi", "normalised_invbetavpT_tof_pi", 3000, 0, 30, 2000, -10, 10);
  normalised_invbetavpT_tof_ka = new TH2F("normalised_invbetavpT_tof_ka", "normalised_invbetavpT_tof_ka", 3000, 0, 30, 2000, -10, 10);

  dEdXvp = new TH2F("dEdXvp", "dEdXvp", 3000, 0, 30, 4000, -20, 20);
  invbetavp = new TH2F("invbetavp", "invbetavp", 3000, 0, 30, 4000, -20, 20);


  // Topological histograms

  decaylengthd0US = new TH1F("decaylengthd0US", "decaylengthd0US", 10000, 0, 1);
  distancepikUS = new TH1F("distancepikUS", "distancepikUS", 10000, 0, 1);
  distanced0PVUS = new TH1F("distanced0PVUS", "distanced0PVUS", 10000, 0, 1);
  dcakPVUS = new TH1F("dcakPVUS", "dcakPVUS", 10000, 0, 1);
  dcapiPVUS = new TH1F("dcapiPVUS", "dcapiPVUS", 10000, 0, 1);

  decaylengthd0LS = new TH1F("decaylengthd0LS", "decaylengthd0LS", 10000, 0, 1);
  distancepikLS = new TH1F("distancepikLS", "distancepikLS", 10000, 0, 1);
  distanced0PVLS = new TH1F("distanced0PVLS", "distanced0PVLS", 10000, 0, 1);
  dcakPVLS = new TH1F("dcakPVLS", "dcakPVLS", 10000, 0, 1);
  dcapiPVLS = new TH1F("dcapiPVLS", "dcapiPVLS", 10000, 0, 1);

  // hD0US_pT = new TH1F("hD0US_pT", "hD0US_pT", 1000, 0, 50);
  // hD0LS_pT = new TH1F("hD0LS_pT", "hD0LS_pT", 1000, 0, 50);


  //Invariant Mass Plots

  kaonpt = new TH1F("kaonpt", "kaonpt", 1000, 0, 30);
  pionpt = new TH1F("pionpt", "pionpt", 1000, 0, 30);
  d0pt = new TH1F("d0pt", "d0pt", 1000, 0, 30);
  kaonpionpt = new TH2F("kaonpionpt", "kaonpionpt", 1000, 0, 30, 1000, 0, 30);
  kaoneta = new TH1F("kaoneta", "kaoneta", 1000, -1, 1);
  pioneta = new TH1F("pioneta", "pioneta", 1000, -1, 1);
  d0eta = new TH1F("d0eta", "d0eta", 1000, -5, 5);
  kaonphi = new TH1F("kaonphi", "kaonphi", 1000, -10, 10);
  pionphi = new TH1F("pionphi", "pionphi", 1000, -10, 10);
  d0phi = new TH1F("d0phi", "d0phi", 1000, -10, 10);


  kaonbgUSpt = new TH1F("kaonbgUSpt", "kaonbgUSpt", 1000, 0, 30);
  pionbgUSpt = new TH1F("pionbgUSpt", "pionbgUSpt", 1000, 0, 30);
  d0bgUSpt = new TH1F("d0bgUSpt", "d0bgUSpt", 1000, 0, 30);
  kaonpionbgUSpt = new TH2F("kaonpionbgUSpt", "kaonpionbgUSpt", 1000, 0, 30, 1000, 0, 30);
  kaonbgUSeta = new TH1F("kaonbgUSeta", "kaonbgUSeta", 1000, -1, 1);
  pionbgUSeta = new TH1F("pionbgUSeta", "pionbgUSeta", 1000, -1, 1);
  d0bgUSeta = new TH1F("d0bgUSeta", "d0bgUSeta", 1000, -5, 5);
  kaonbgUSphi = new TH1F("kaonbgUSphi", "kaonbgUSphi", 1000, -10, 10);
  pionbgUSphi = new TH1F("pionbgUSphi", "pionbgUSphi", 1000, -10, 10);
  d0bgUSphi = new TH1F("d0bgUSphi", "d0bgUSphi", 1000, -10, 10);


  kaonbgLSpt = new TH1F("kaonbgLSpt", "kaonbgLSpt", 1000, 0, 30);
  pionbgLSpt = new TH1F("pionbgLSpt", "pionbgLSpt", 1000, 0, 30);
  d0bgLSpt = new TH1F("d0bgLSpt", "d0bgLSpt", 1000, 0, 30);
  kaonpionbgLSpt = new TH2F("kaonpionbgLSpt", "kaonpionbgLSpt", 1000, 0, 30, 1000, 0, 30);
  kaonbgLSeta = new TH1F("kaonbgLSeta", "kaonbgLSeta", 1000, -1, 1);
  pionbgLSeta = new TH1F("pionbgLSeta", "pionbgLSeta", 1000, -1, 1);
  d0bgLSeta = new TH1F("d0bgLSeta", "d0bgLSeta", 1000, -5, 5);
  kaonbgLSphi = new TH1F("kaonbgLSphi", "kaonbgLSphi", 1000, -10, 10);
  pionbgLSphi = new TH1F("pionbgLSphi", "pionbgLSphi", 1000, -10, 10);
  d0bgLSphi = new TH1F("d0bgLSphi", "d0bgLSphi", 1000, -10, 10);

  invmass = new TH1F("invmass", "invmass", 1000, 0.5, 2.5);
  invmassbg = new TH1F("invmassbg", "invmassbg", 1000, 0.5, 2.5);

  hNumberOfD0s = new TH1F("hNumberOfD0s", "hNumberOfD0s", 21, -0.5, 20.5);
  hNumberOfD0BgUS = new TH1F("hNumberOfD0BgUS", "hNumberOfD0BgUS", 21, -0.5, 20.5);
  hNumberOfD0BgLS = new TH1F("hNumberOfD0BgLS", "hNumberOfD0BgLS", 21, -0.5, 20.5);

}
//
// write histograms
//_____________________________________________________________________________
void StEventCounter::WriteHistograms() {

  // Centrality Histograms

  hEventZvertex_whole->Write();
  hEventZvertex_VPD->Write();
  hEventZvertex_diff->Write();

  hCentrality->Write();
  hMultiplicity->Write();

  // Event Cuts Histograms
  const char *event_cuts[10] = {"Total Events", "MinBias","|V_{z}| < 6 cm", "|V_{z} - V_{z(VPD)}| < 3 cm", "nBEMCMatch > 0", "nBTOFMatch > 0", "HT1", "HT2", "HT3", ""};
  const char *track_cuts[10] = {"Total Tracks", "nHitsDedx > 20", "#frac{nHitsDedx}{nHitsMax} > 0.52", "DCA < 3 cm", "BEMC Hit", "BTOF Hit", "Charged Tracks", "", "", ""};

  for (int i=1; i <= 10; i++){
    cuthistogram_event->SetBinContent(i, numberofevents[i-1]);
    cuthistogram_event->GetXaxis()->SetBinLabel(i, event_cuts[i-1]);
    // cuthistogram_track->SetBinContent(i, numberoftracks[i-1]);
    // cuthistogram_track->GetXaxis()->SetBinLabel(i, track_cuts[i-1]);
  }
  cuthistogram_event->Write();
  // cuthistogram_track->Write();

  // Track Cut Histograms
  for (int i = 0; i < 8; i++){
    hAcceptedPt[i]->Write();
    hAcceptedEta[i]->Write();
  }

  // Particle Identification Histograms
  z_pi->Write();
  z_ka->Write();
  normalised_invbetavpT_tof_pi->Write();
  normalised_invbetavpT_tof_ka->Write();

  dEdXvp->Write();
  invbetavp->Write();

  // Topological Histograms

  decaylengthd0US->Write();
  distancepikUS->Write();
  distanced0PVUS->Write();
  dcakPVUS->Write();
  dcapiPVUS->Write();

  decaylengthd0LS->Write();
  distancepikLS->Write();
  distanced0PVLS->Write();
  dcakPVLS->Write();
  dcapiPVLS->Write();

  // Pt and Invmass Histograms for KPi pairs after topological cuts

  kaonpt->Write();
  pionpt->Write();
  d0pt->Write();
  kaonpionpt->Write();
  kaoneta->Write();
  pioneta->Write();
  d0eta->Write();
  kaonphi->Write();
  pionphi->Write();
  d0phi->Write();

  kaonbgUSpt->Write();
  pionbgUSpt->Write();
  d0bgUSpt->Write();
  kaonpionbgUSpt->Write();
  kaonbgUSeta->Write();
  pionbgUSeta->Write();
  d0bgUSeta->Write();
  kaonbgUSphi->Write();
  pionbgUSphi->Write();
  d0bgUSphi->Write();

  kaonbgLSpt->Write();
  pionbgLSpt->Write();
  d0bgLSpt->Write();
  kaonpionbgLSpt->Write();
  kaonbgLSeta->Write();
  pionbgLSeta->Write();
  d0bgLSeta->Write();
  kaonbgLSphi->Write();
  pionbgLSphi->Write();
  d0bgLSphi->Write();

  invmass->Write();
  invmassbg->Write();

  hNumberOfD0s->Write();
  hNumberOfD0BgUS->Write();
  hNumberOfD0BgLS->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StEventCounter::Clear(Option_t *opt) {
  // fJets->Clear();
  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StEventCounter::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  fd0 = kFALSE;
  fd0BgUS = kFALSE;
  fd0BgLS = kFALSE;

  //zero out global vectors
  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();

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

  // cout << "B = " << Bfield << endl;
  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  double zVtx_VPD = mPicoEvent->vzVpd();

  hEventZvertex_whole->Fill(zVtx);
  hEventZvertex_VPD->Fill(zVtx_VPD);
  hEventZvertex_diff->Fill(zVtx - zVtx_VPD);

  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


  // HFT Analysis Specific Cuts For Event

    // ============================ CENTRALITY ============================== //
    // get CentMaker pointer
    mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
    if(!mCentMaker) {
      LOG_WARN << " No CentMaker! Skip! " << endm;
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

    //cout << "The centrality bin is " << fCentralityScaled << endl;
    //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
    // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

    // cut on unset centrality, > 80%

  if (!fMCEventsWithoutCent)
  {
    if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

    cout << "Here" << endl;

    // fill histograms
    hCentrality->Fill(fCentralityScaled);
    hMultiplicity->Fill(refCorr2);

  

    // cut on centrality for analysis before doing anything
    if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

    // ============================ end of CENTRALITY ============================== //

    // ========================= Trigger Info =============================== //
    // looking at the EMCal triggers - used for QA and deciding on HT triggers
    //  FillEmcTriggers();

    // get trigger IDs from PicoEvent class and loop over them
    vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
    if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"EventTriggers: ";
    for(unsigned int i=0; i<mytriggers.size(); i++) {
      if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
    }
    if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<endl;

    numberofevents[0]++;

    // check for MB/HT event
    bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
    bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
    bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
    bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);

    bool fRunForMB = kFALSE;
    if(doppAnalysis) fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
    if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

    int arrMB5_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
    int arrHT1_Run14[] = {450201, 450211, 460201};
    int arrHT2_Run14[] = {450202, 450212, 460202, 460212};
    int arrHT3_Run14[] = {450203, 450213, 460203};

    if (!doppAnalysis){

      bool matchMB = kFALSE;

      for(int i = 0; i < sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14); i++) {
        if(mPicoEvent->isTrigger(arrMB5_Run14[i])) matchMB = kTRUE;
        if(matchMB) break;
      }

      if (!matchMB) return kStOk;

    }

    numberofevents[1]++;
    // ======================== end of Triggers ============================= //

    

    if (abs(zVtx) > 6) return kStOk;
    numberofevents[2]++;

    if (abs(zVtx - zVtx_VPD) > 3) return kStOk;
    numberofevents[3]++;
     
    if (mPicoEvent->nBEMCMatch() == 0) return kStOk;
    numberofevents[4]++;

    if (mPicoEvent->nBTOFMatch() == 0) return kStOk;
    numberofevents[5]++;

    if (!doppAnalysis){

      bool matchHT1 = kFALSE;
      bool matchHT2 = kFALSE;
      bool matchHT3 = kFALSE;

      for(int i = 0; i < sizeof(arrHT1_Run14)/sizeof(*arrHT1_Run14); i++) {
        if(mPicoEvent->isTrigger(arrHT1_Run14[i])) matchHT1 = kTRUE;
        if(matchHT1) break;
      }

      for(int i = 0; i < sizeof(arrHT2_Run14)/sizeof(*arrHT2_Run14); i++) {
        if(mPicoEvent->isTrigger(arrHT2_Run14[i])) matchHT2 = kTRUE;
        if(matchHT2) break;
      }

      for(int i = 0; i < sizeof(arrHT3_Run14)/sizeof(*arrHT3_Run14); i++) {
        if(mPicoEvent->isTrigger(arrHT3_Run14[i])) matchHT3 = kTRUE;
        if(matchHT3) break;
      }

      if (matchHT1 || matchHT2 || matchHT3) numberofevents[6]++;


      if (matchHT2 || matchHT3) numberofevents[7]++;

     
      if (matchHT3) numberofevents[8]++;
    // cout << "MinBias Event Accepted." << endl;
    }
  }

  return kStOK;
}

