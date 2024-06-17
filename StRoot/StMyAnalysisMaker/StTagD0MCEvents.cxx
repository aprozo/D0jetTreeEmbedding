// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StTagD0MCEvents.h"
#include "StMemStat.h"

#include <algorithm>

// ROOT includes
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLorentzVector.h"


// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoMcTrack.h"
#include "StRoot/StPicoEvent/StPicoMcVertex.h"
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

ClassImp(StTagD0MCEvents)

//________________________________________________________________________
StTagD0MCEvents::StTagD0MCEvents(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{ 
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StTagD0MCEvents::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  fDoEffCorr = kFALSE;
  doRejectBadRuns = kFALSE;
  fTopoLevel = 1;
  fAnaCentBin = 2;
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

  fReplaceBadTracks = kFALSE;
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

  fPrintLevel = 0;

  fd0MCTrackIndices.clear();
  fd0RecoTrackIndices.clear();

  fpionReco4Momenta.clear();
  fkaonReco4Momenta.clear();
  fd0Reco4Momenta.clear();

  fpionRecoMomenta.clear();
  fkaonRecoMomenta.clear();
  fd0RecoMomenta.clear();

  fDroppedTracks.clear();
  fDroppedMCTracks.clear();

  fTracksSmearedWell.clear();

  fKaonMomResolution = NULL;
  fPionMomResolution = NULL;
  fProtonMomResolution = NULL;

  fPionWeight = NULL;
  fKaonWeight = NULL;
  fProtonWeight = NULL;
  fAProtonWeight = NULL;

  d1 = NULL;
  d2 = NULL;
}

//
//________________________________________________________________________
StTagD0MCEvents::~StTagD0MCEvents()
{ 
  if(mEmcPosition) delete mEmcPosition;
}

//________________________________________________________________________
Int_t StTagD0MCEvents::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  TFile f("/star/u/droy1/Y2019/STAR/Momentum_resolution_SL16d.root");
  fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPion");
  fKaonMomResolution = (TF1*)f.Get("fKaon")->Clone("fKaon");
  fProtonMomResolution = (TF1*)f.Get("fProton")->Clone("fProton");

  TFile effweight("/star/u/droy1/Y2019/STAR/EffWeightsInCentralityBins.root");

  if (fAnaCentBin == 0){
    fPionWeight = (TGraph *)effweight.Get("Pion_0_10");
    fKaonWeight = (TGraph *)effweight.Get("Kaon_0_10");
    fProtonWeight = (TGraph *)effweight.Get("Proton_0_10");
    fAProtonWeight = (TGraph *)effweight.Get("AProton_0_10");
  }

  else if (fAnaCentBin == 1){
    fPionWeight = (TGraph *)effweight.Get("Pion_10_40");
    fKaonWeight = (TGraph *)effweight.Get("Kaon_10_40");
    fProtonWeight = (TGraph *)effweight.Get("Proton_10_40");
    fAProtonWeight = (TGraph *)effweight.Get("AProton_10_40");
  }

  else if (fAnaCentBin == 2){
    fPionWeight = (TGraph *)effweight.Get("Pion_40_80");
    fKaonWeight = (TGraph *)effweight.Get("Kaon_40_80");
    fProtonWeight = (TGraph *)effweight.Get("Proton_40_80");
    fAProtonWeight = (TGraph *)effweight.Get("AProton_40_80");
  }

  else{
    LOG_WARN << " Analysis Centrality Class Not Recognised. Try again! " << endm;
    return kStWarn;
  }

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StTagD0MCEvents::Finish() { 

  cout << "StTagD0MCEvents::Finish()\n";

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

  cout<<"End of StTagD0MCEvents::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StTagD0MCEvents::DeclareHistograms() {

  invmass = new TH1F("invmass", "invmass", 1000, 0.5, 2.5);

  hNumberOfD0s = new TH1F("hNumberOfD0s", "hNumberOfD0s", 21, -0.5, 20.5);
  hNumberOfD0BgUS = new TH1F("hNumberOfD0BgUS", "hNumberOfD0BgUS", 21, -0.5, 20.5);
  hNumberOfD0BgLS = new TH1F("hNumberOfD0BgLS", "hNumberOfD0BgLS", 21, -0.5, 20.5);

  hRecoVMCKaonPt = new TH2F("hRecoVMCKaonPt", "RecovMC Kaon Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);
  hRecoVMCPionPt = new TH2F("hRecoVMCPionPt", "RecovMC Pion Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);

  hRecoVMCKaonPhi = new TH2F("hRecoVMCKaonPhi", "RecovMC Kaon Phi", 80, -10, 10, 80, -10, 10);
  hRecoVMCPionPhi = new TH2F("hRecoVMCPionPhi", "RecovMC Pion Phi", 80, -10, 10, 80, -10, 10);

  hRecoVMCKaonEta = new TH2F("hRecoVMCKaonEta", "RecovMC Kaon Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);
  hRecoVMCPionEta = new TH2F("hRecoVMCPionEta", "RecovMC Pion Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);

  hRecoVMCKaonPt_Smeared = new TH2F("hRecoVMCKaonPt_Smeared", "RecovMC_Smeared Kaon Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);
  hRecoVMCPionPt_Smeared = new TH2F("hRecoVMCPionPt_Smeared", "RecovMC_Smeared Pion Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);

  hRecoVMCKaonPhi_Smeared = new TH2F("hRecoVMCKaonPhi_Smeared", "RecovMC_Smeared Kaon Phi", 80, -10, 10, 80, -10, 10);
  hRecoVMCPionPhi_Smeared = new TH2F("hRecoVMCPionPhi_Smeared", "RecovMC_Smeared Pion Phi", 80, -10, 10, 80, -10, 10);

  hRecoVMCKaonEta_Smeared = new TH2F("hRecoVMCKaonEta_Smeared", "RecovMC_Smeared Kaon Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);
  hRecoVMCPionEta_Smeared = new TH2F("hRecoVMCPionEta_Smeared", "RecovMC_Smeared Pion Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);

  hMCD0Pt = new TH1F("hMCD0Pt", "hMCD0Pt", 100, 0, 10);
  hMCD0Eta = new TH1F("hMCD0Eta", "hMCD0Eta", 200, -2, 2);

  hTracksPtQM = new TH1F("hTracksPtQM", "hTracksPtQM", 60, 0, 30);
  hTracksPt   = new TH1F("hTracksPt", "hTracksPt", 60, 0, 30);

  hMCPtvDCA = new TH2F("hMCPtvDCA", "hMCPtvDCA", 60, 0, 30, 2000, 0, 1000);
  hMCPtvNMatch = new TH2F("hMCPtvNMatch", "hMCPtvNMatch", 60, 0, 30, 10, 0.5, 10.5);
  hMinMaxDCA = new TH2F("hMinMaxDCA", "hMinMaxDCA", 500, 0, 500, 500, 0, 500);
}
//
// write histograms
//_____________________________________________________________________________
void StTagD0MCEvents::WriteHistograms() {

  invmass->Write();
  hMCD0Pt->Write();
  hMCD0Eta->Write();

  hRecoVMCKaonPt->Write();
  hRecoVMCPionPt->Write();

  hRecoVMCKaonPhi->Write();
  hRecoVMCPionPhi->Write();

  hRecoVMCKaonEta->Write();
  hRecoVMCPionEta->Write();

  hRecoVMCKaonPt_Smeared->Write();
  hRecoVMCPionPt_Smeared->Write();

  hRecoVMCKaonPhi_Smeared->Write();
  hRecoVMCPionPhi_Smeared->Write();

  hRecoVMCKaonEta_Smeared->Write();
  hRecoVMCPionEta_Smeared->Write();

  hTracksPtQM->Write();
  hTracksPt->Write();
  hMCPtvDCA->Write();
  hMCPtvNMatch->Write();
  hMinMaxDCA->Write();

}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StTagD0MCEvents::Clear(Option_t *opt) {
  // fJets->Clear();
  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();
  fd0RecoTrackIndices.clear();
  fd0MCTrackIndices.clear();
  fpionRecoMomenta.clear();
  fkaonRecoMomenta.clear();
  fd0RecoMomenta.clear();
  fpionReco4Momenta.clear();
  fkaonReco4Momenta.clear();
  fd0Reco4Momenta.clear();
  fTracksSmearedWell.clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StTagD0MCEvents::Make() {
  if (fPrintLevel) cout << "Called StTagD0MCEvents" << endl;
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  fd0 = kFALSE;
  fd0BgUS = kFALSE;
  fd0BgLS = kFALSE;

  //zero out global vectors
  fd0MCTrackIndices.clear();
  fd0RecoTrackIndices.clear();
  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();

  fpionRecoMomenta.clear();
  fkaonRecoMomenta.clear();
  fd0RecoMomenta.clear();
  fpionReco4Momenta.clear();
  fkaonReco4Momenta.clear();
  fd0Reco4Momenta.clear();

  fDroppedTracks.clear();
  fDroppedMCTracks.clear();

  fTracksSmearedWell.clear();

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

  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


  // HFT Analysis Specific Cuts For Event


  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    // LOG_WARN << " No CentMaker! Skip! " << endm;
    // return kStWarn;
  }

  int grefMult; // see StPicoEvent
  int refMult;  // see StPicoEvent
  int cent16; // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin;
  double refCorr2;

  if (!mCentMaker){
    grefMult = 0;
    refMult = 0;
    ref9 = 0;
    ref16 = 0;
    cent16 = 0;
    centbin = 0;
    refCorr2 = 0.0;
    fCentralityScaled = 0;
  }

  else{
    // centrality variables
    grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
    refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
    ref9 = mCentMaker->GetRef9();   // binning from central -> peripheral
    ref16 = mCentMaker->GetRef16(); // binning from central -> peripheral
    cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
    centbin = mCentMaker->GetRef16();
    refCorr2 = mCentMaker->GetRefCorr2();
    fCentralityScaled = mCentMaker->GetCentScaled();
  }

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

  numberofevents[1]++;
  // ======================== end of Triggers ============================= //

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

  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // cout << "Event Number = " << mPicoEvent->eventId() << endl;

  fvertextotrack = GetMCTrackListForVertex();

  GetD0TracksPair();
  ReplacePoorlySmearedTracks();
  
  delete[] fvertextotrack;

  DoQAForDifferencesWithQM();

  return kStOK;
}

void StTagD0MCEvents::DoQAForDifferencesWithQM(){
    const Int_t ntracks = mPicoDst->numberOfTracks();
    double pi0mass = Pico::mMass[0]; // GeV

    for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){

        double pt = fTracksSmearedWell[iTracks][0];
        double phi = fTracksSmearedWell[iTracks][1];
        if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
        if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
        double eta = fTracksSmearedWell[iTracks][2];
        double px = fTracksSmearedWell[iTracks][3];
        double py = fTracksSmearedWell[iTracks][4];
        double pz = fTracksSmearedWell[iTracks][5];
        double p = fTracksSmearedWell[iTracks][6];
        double energy = fTracksSmearedWell[iTracks][7];
        short charge = fTracksSmearedWell[iTracks][8];
        int matchedTowerIndex = fTracksSmearedWell[iTracks][9]; // towerIndex = towerID - 1

        if(pt < 0.2) continue;
        if(pt > 30.) continue;
        if((eta < -1.) || (eta > 1.)) continue;

        hTracksPtQM->Fill(pt);
        // hTracksPtvDCAQM->Fill(pt, dca);

        StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
        if(!trk){ continue; }

        double dca = trk->gDCA(mVertex).Mag();
        int nHitsFit = trk->nHitsFit();
        int nHitsMax = trk->nHitsMax();
        double nHitsRatio = 1.0*nHitsFit/nHitsMax;

        // additional quality cuts for tracks
        if(dca > fTrackDCAcut)            continue;
        if(nHitsFit < fTracknHitsFit)     continue;
        if(nHitsRatio < fTracknHitsRatio) continue;

        hTracksPt->Fill(pt);
    }

    const Int_t nmctracks = mPicoDst->numberOfMcTracks();

    for (unsigned short track = 0; track < nmctracks; track++){
        StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track));
        mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track));

        if (mctrk->idVtxStop() != 0) continue;

        int numberofmatchedtracks = 0;
        vector<int> matchedreco;
        matchedreco.clear();
        vector<double> matchedrecodca;
        matchedrecodca.clear();

        double mcpt  = mctrk->p().Perp();
        double mceta = mctrk->p().PseudoRapidity();

        if(mcpt < 0.2) continue;
        if(mcpt > 30.) continue;
        if((mceta < -1) || (mceta > 1)) continue;

        for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
            StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
            if(!trk){ continue; }
            if (trk->idTruth() - 1 != track) continue;

            double pt = fTracksSmearedWell[iTracks][0];
            double eta = fTracksSmearedWell[iTracks][2]; 
            
            if(pt < 0.2) continue;
            if(pt > 30.) continue;
            if((eta < -1.) || (eta > 1.)) continue;

            double dca = trk->gDCA(mVertex).Mag();
            int nHitsFit = trk->nHitsFit();
            int nHitsMax = trk->nHitsMax();
            double nHitsRatio = 1.0*nHitsFit/nHitsMax;

            numberofmatchedtracks++;
            matchedreco.push_back(iTracks);
            matchedrecodca.push_back(dca);
        }

        if (matchedreco.size() <= 1) continue;
        hMCPtvNMatch->Fill(mcpt, matchedreco.size());

        double maxElement = *std::max_element(matchedrecodca.begin(), matchedrecodca.end());
        double minElement = *std::min_element(matchedrecodca.begin(), matchedrecodca.end());

        hMinMaxDCA->Fill(minElement, maxElement);

        for(unsigned short iTracks = 0; iTracks < matchedreco.size(); iTracks++){
            StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(matchedreco[iTracks]));
            if(!trk){ continue; }

            double pt = fTracksSmearedWell[matchedreco[iTracks]][0];
            double eta = fTracksSmearedWell[matchedreco[iTracks]][2]; 
            
            if(pt < 0.2) continue;
            if(pt > 30.) continue;
            if((eta < -1.) || (eta > 1.)) continue;

            double dca = trk->gDCA(mVertex).Mag();
            int nHitsFit = trk->nHitsFit();
            int nHitsMax = trk->nHitsMax();
            double nHitsRatio = 1.0*nHitsFit/nHitsMax;

            hMCPtvDCA->Fill(mcpt, dca);
        }
    }

}


vector <int> *StTagD0MCEvents::GetMCTrackListForVertex(){
  const int numvertices = 1000;

  vector<int> *arr = new vector<int> [numvertices];
  for (int i = 0; i < numvertices; i++){
    arr[i].clear();
  }

  for (int mc = 0; mc< mPicoDst->numberOfMcTracks(); mc++){
    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(mc));
    int idvxstart = mctrk->idVtxStart() -  1;

    arr[idvxstart].push_back(mc);
  }

  if (fPrintLevel){
    for (int i = 0; i < numvertices; i++){
      if (arr[i].size() != 0) cout << "Tracks in Vertex = " << i << ": " << "\t" ;
      for (int t = 0; t < arr[i].size(); t++){
        cout << arr[i][t] << "\t";
      }
      if (arr[i].size() != 0) cout << endl;
    }
  }

  return arr;
}

/* Get the reconstructed D0 Kaon and Pion Track IDs and Vtx ID of production.
If either track is not matched, we have -999 as the trackid */


void StTagD0MCEvents::GetD0TracksPair(){

  unsigned int nmctracks = mPicoDst->numberOfMcTracks();

  for (unsigned short track1 = 0; track1 < nmctracks; track1++){
    StPicoMcTrack *mctrk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track1));
    if (mctrk1->geantId() == 37 || mctrk1->geantId() == 38) {
      int idvxstop = mctrk1->idVtxStop() -  1;

      if (mctrk1->p().Perp() < 1 || mctrk1->p().Perp() > 10) continue;
      if (mctrk1->p().PseudoRapidity() < -1 || mctrk1->p().PseudoRapidity() > 1) continue;
      
      StPicoMcVertex *mcvx1 = static_cast<StPicoMcVertex*>(mPicoDst->mcVertex(idvxstop));
      if (fPrintLevel) 
        cout << "MC D0s in the event and daughters = " << mctrk1->geantId() << "\t" << mctrk1->p().Perp() << "\t" <<  mctrk1->p().PseudoRapidity() << "\t" << mcvx1->numberOfDaughters() << endl;

      hMCD0Eta->Fill(mctrk1->p().PseudoRapidity());
      hMCD0Pt->Fill(mctrk1->p().Perp());
    }
  }

  for (unsigned short track1 = 0; track1 < nmctracks; track1++){
    StPicoMcTrack *mctrk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track1));

    for (unsigned short track2 = 0; track2 < nmctracks; track2++){
      StPicoMcTrack *mctrk2 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track2));
      
      int firsttrack = mctrk1->geantId();
      int secondtrack = mctrk2->geantId();

      if (!((firsttrack == 8 && secondtrack == 12) || (firsttrack == 9 && secondtrack == 11))) continue; // Makes sure it's always pions + kaons (not the other way around)

      if (!IsD0MCTrack(track1) || !IsD0MCTrack(track2)) continue; // If both tracks are not from D0s, then discard

      int idvxstart1 = mctrk1->idVtxStart();
      int idvxstart2 = mctrk2->idVtxStart();

      if (idvxstart1 != idvxstart2) continue; // If both tracks are not from the same D0, then discard

      int recotrack1 = GetMatchedRecoTrackFromMCTrack(track1);
      int recotrack2 = GetMatchedRecoTrackFromMCTrack(track2);

      int idvx = idvxstart1;
      fd0MCTrackIndices.push_back({track1, track2, idvx});
      fd0RecoTrackIndices.push_back({recotrack1, recotrack2, idvx});

      TVector3 p1 = FastSimMom(mctrk1->p(), mctrk1->geantId());
      TVector3 p2 = FastSimMom(mctrk2->p(), mctrk2->geantId());

      TLorentzVector v1;
      TLorentzVector v2;

      v1.SetXYZM(p1.X(), p1.Y(), p1.Z(), Mpion);
      v2.SetXYZM(p2.X(), p2.Y(), p2.Z(), Mkaon);

      fpionReco4Momenta.push_back(v1);
      fkaonReco4Momenta.push_back(v2);
      fd0Reco4Momenta.push_back(v1+v2);

      fpionRecoMomenta.push_back(p1);
      fkaonRecoMomenta.push_back(p2);
      fd0RecoMomenta.push_back(p1+p2);

      invmass->Fill((v1+v2).M());
    }
  }

  if (fPrintLevel){
    for (int D0 = 0; D0 < fd0RecoMomenta.size(); D0++){
      cout << "Reco D0s in the event = " << fd0RecoMomenta[D0].Perp() << "\t" << fd0RecoMomenta[D0].PseudoRapidity() << "\t" << fd0RecoMomenta[D0].Phi() << endl;
    }
  }

  MCTracksToDiscard();
  FindTracksToDrop();
}

void StTagD0MCEvents::MCTracksToDiscard(){ // This is for MC tracks
  for (int d0 = 0; d0 < fd0MCTrackIndices.size(); d0++){
    StPicoMcTrack *daugtrk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(fd0MCTrackIndices[d0][0]));
    StPicoMcTrack *daugtrk2 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(fd0MCTrackIndices[d0][1]));

    int daugstopvx1 = daugtrk1->idVtxStop() - 1;
    int daugstopvx2 = daugtrk2->idVtxStop() - 1;

    vector<int> tmp;
    tmp.clear();

    if (fPrintLevel > 1) cout << "Pion" << endl;
    GetAllTracksFromVertex(daugstopvx1, tmp);
    if (fPrintLevel > 1) cout << "Kaon" << endl;
    GetAllTracksFromVertex(daugstopvx2, tmp);

    if (fPrintLevel){
      cout << "Dropped Track List for D0 # = " << d0 << endl;  
      for (int i = 0; i < tmp.size(); i++){
        StPicoMcTrack *t = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(tmp[i]));
        cout << tmp[i] << "\t" << t->geantId() << "\t" << t->p().Perp() << "\t" << t->p().PseudoRapidity() << "\t" << endl;
      }
      cout << endl;
    }
    fDroppedMCTracks.push_back(tmp);
  }

  if (fPrintLevel)cout << "Dropped Track Array Size = " << fDroppedMCTracks.size() << endl;
}

void StTagD0MCEvents::FindTracksToDrop()
{
  unsigned int ntracks = mPicoDst->numberOfTracks();

  for (int d0 = 0; d0 < fd0RecoTrackIndices.size(); d0++){

    std::vector<int> temptrackstodrop;
    temptrackstodrop.clear();

    for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
      if(!trk){ continue; }

      if (iTracks == fd0RecoTrackIndices[d0][0] || iTracks == fd0RecoTrackIndices[d0][1]) continue;

      if (!AcceptTrack(trk, Bfield, mVertex)) continue;

      if (!KeepTrack(iTracks)) {
        temptrackstodrop.push_back(iTracks);
        continue;
        // cout << "Track Dropped = " << iTracks << endl;
      }

      if (IsTrackADescendantOfD0Daughters(iTracks, d0)){
        temptrackstodrop.push_back(iTracks);
      }
    }

    fDroppedTracks.push_back(temptrackstodrop);

  }

}

void StTagD0MCEvents::ReplacePoorlySmearedTracks(){

  const Int_t ntracks = mPicoDst->numberOfTracks();
  double pi0mass = Pico::mMass[0]; // GeV

  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
    // get track pointer

    bool trackeligibleforreplacement = kTRUE;

    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

    // // acceptance and kinematic quality cuts - pt cut is also applied here currently
    // if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) { 
      // get primary track vector
      mTrkMom = trk->pMom(); 
    } else { 
      // get global track vector
      mTrkMom = trk->gMom(mVertex, Bfield); 
    }

    int mctrkid = trk->idTruth();

    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(mctrkid - 1));
    if(!mctrk) {
      trackeligibleforreplacement = kFALSE;
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
    if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
    double eta = mTrkMom.PseudoRapidity();
    double px = mTrkMom.x();
    double py = mTrkMom.y();
    double pz = mTrkMom.z();
    double p = mTrkMom.Mag();
    double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
    short charge = trk->charge();
    int matchedTowerIndex = trk->bemcTowerIndex(); // towerIndex = towerID - 1

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
    int matchedTowerIndex_new = matchedTowerIndex;

    double dca = trk->gDCA(mVertex).Mag();
    int nHitsFit = trk->nHitsFit();
    int nHitsMax = trk->nHitsMax();
    double nHitsRatio = 1.0*nHitsFit/nHitsMax;

    bool goodtrack = (dca < fTrackDCAcut) && (abs(nHitsFit) >= fTracknHitsFit) && (nHitsRatio >= fTracknHitsRatio);

    if (!goodtrack) trackeligibleforreplacement = kFALSE;

    if (fReplaceBadTracks && trackeligibleforreplacement){

      //Track Variables {pt, phi, eta, px, py, pz, p, energy, charge, matchedTowerIndex}

      double relativesmearing = TMath::Sqrt(pow(mctrk->p().Px() - px, 2) + pow(mctrk->p().Py() - py, 2))/(mctrk->p().Pt());

      double fastsimsmearing;

      int pid = mctrk->geantId();

      if(pid ==8 || pid == 9)fastsimsmearing = fPionMomResolution->Eval(mctrk->p().Pt());// Pion
      else if(pid ==11 || pid == 12)fastsimsmearing =  fKaonMomResolution->Eval(mctrk->p().Pt());// Kaon
      else if(pid ==15 || pid == 14)fastsimsmearing =  fProtonMomResolution->Eval(mctrk->p().Pt());// Proton
      else fastsimsmearing = fPionMomResolution->Eval(mctrk->p().Pt());// Catch all: pions

      if (relativesmearing > 3*fastsimsmearing){
        TVector3 fastsimsmearedmom = FastSimMom(mctrk->p(), pid);
        pt_new = fastsimsmearedmom.Perp();
        phi_new = fastsimsmearedmom.Phi();
        if(phi_new < 0.0)    phi_new += 2.0*pi;  // force from 0-2pi
        if(phi_new > 2.0*pi) phi_new -= 2.0*pi;  // force from 0-2pi
        eta_new = fastsimsmearedmom.PseudoRapidity();
        px_new = fastsimsmearedmom.x();
        py_new = fastsimsmearedmom.y();
        pz_new = fastsimsmearedmom.z();
        p_new = fastsimsmearedmom.Mag();
        energy_new = 1.0*TMath::Sqrt(p_new*p_new + pi0mass*pi0mass);
      }

      fTracksSmearedWell.push_back({pt_new, phi_new, eta_new, px_new, py_new, pz_new, p_new, energy_new, charge_new, matchedTowerIndex_new});
    }

    else{
      fTracksSmearedWell.push_back({pt, phi, eta, px, py, pz, p, energy, charge, matchedTowerIndex});
    }

    ////////////////////////// All Kaon Pion Histograms //////////////////////////////////////////
    // if(!IsD0MCVertex(tracksgeantidsfromvertex, iVtx)){
    if (mctrk->geantId() == 8 || mctrk->geantId() == 9){

        hRecoVMCPionPt->Fill(pt, mctrk->pt());
        hRecoVMCPionPhi->Fill(phi, mctrk->p().Phi());
        hRecoVMCPionEta->Fill(eta, mctrk->eta());

        hRecoVMCPionPt_Smeared->Fill(pt_new, mctrk->pt());
        hRecoVMCPionPhi_Smeared->Fill(phi_new, mctrk->p().Phi());
        hRecoVMCPionEta_Smeared->Fill(eta_new, mctrk->eta());
    }

    if (mctrk->geantId() == 11 || mctrk->geantId() == 12){
                
        hRecoVMCKaonPt->Fill(pt, mctrk->pt());
        hRecoVMCKaonPhi->Fill(phi, mctrk->p().Phi());
        hRecoVMCKaonEta->Fill(eta, mctrk->eta());

        hRecoVMCKaonPt_Smeared->Fill(pt_new, mctrk->pt());
        hRecoVMCKaonPhi_Smeared->Fill(phi_new, mctrk->p().Phi());
        hRecoVMCKaonEta_Smeared->Fill(eta_new, mctrk->eta());

    }
    // }
  }
}

////// Helper Functions ////////

Bool_t StTagD0MCEvents::IsD0MCTrack(int trackid){

  bool d0Track = kFALSE;

  for (unsigned short track1 = 0; track1 < mPicoDst->numberOfMcTracks(); track1++){
    StPicoMcTrack *mctrk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track1));

    if (mctrk1->geantId() == 37 || mctrk1->geantId() == 38) {

      double pt = mctrk1->p().Perp();
      double eta = mctrk1->p().PseudoRapidity();

      if (pt < 1. || pt > 10.) continue;
      if (abs(eta) > 1.) continue;

      int idvxstop = mctrk1->idVtxStop() -  1;

      // cout << "ID Vertex Stop : " << idvxstop << endl;

      if (idvxstop < 0) return d0Track; // Need this check. This should never happen for a D0

      StPicoMcVertex *mcvx1 = static_cast<StPicoMcVertex*>(mPicoDst->mcVertex(idvxstop));

      int stopvertexid = mcvx1->id();

      // cout << "ID Vertex Stop : " << stopvertexid << endl;


      StPicoMcTrack *D0Track = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(trackid));

      int idvxstart = D0Track->idVtxStart() -  1;

      // int idvxstop = D0Track->idVtxStop() -  1;

      // if (idvxstop < 0) return d0Track;

      StPicoMcVertex *d0vx = static_cast<StPicoMcVertex*>(mPicoDst->mcVertex(idvxstart));

      int startvertexidfortrack = d0vx->id();

      // cout << "ID Vertex Start : " << startvertexidfortrack << endl;

      if (stopvertexid == startvertexidfortrack) {
        return kTRUE;
      }    
    }

  }

  return d0Track; 

}

/* Returns a matched reco track with our acceptance requirements for an MC track. If there is no matching,
it returns -999. */

int StTagD0MCEvents::GetMatchedRecoTrackFromMCTrack(int mctrkid){
  int recotrackmatch = -999;

  unsigned int ntracks = mPicoDst->numberOfTracks();

  double ratio = -99.;

  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

    // if(!IsAnAcceptableTrack(trk, kFALSE)) continue;

    double dca = trk->gDCA(mVertex).Mag();
    int nHitsFit = trk->nHitsFit();
    int nHitsMax = trk->nHitsMax();
    double nHitsRatio = 1.0*nHitsFit/nHitsMax;

    // additional quality cuts for tracks
    if(dca > fTrackDCAcut)            continue;
    if(nHitsFit < fTracknHitsFit)     continue;
    if(nHitsRatio < fTracknHitsRatio) continue;

    int mctrk = trk->idTruth() - 1;

    if (mctrk == mctrkid) {
      if (nHitsRatio > ratio){
        recotrackmatch = iTracks;
        ratio = nHitsRatio;
      }
    }
  }

  return recotrackmatch;
}

void StTagD0MCEvents::GetAllTracksFromVertex(int vertexid, vector <int> &trackvec){
  if (vertexid < 0) return;
  
  if (fPrintLevel > 1){
    cout << "Called this function for vx " << vertexid << " with ntracks = " << fvertextotrack[vertexid].size() << endl;
    for (int track = 0; track < fvertextotrack[vertexid].size(); track++){
      StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(fvertextotrack[vertexid][track]));
      cout << "Geant ID of tracks = " << mctrk->geantId() << endl;
    }
  }

  for (int track = 0; track < fvertextotrack[vertexid].size(); track++){
    trackvec.push_back(fvertextotrack[vertexid][track]);
    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(fvertextotrack[vertexid][track]));
    // if (mctrk->idVtxStop() != 0) {
      GetAllTracksFromVertex(mctrk->idVtxStop() - 1, trackvec);
    // }
  }

  return;
}

Bool_t StTagD0MCEvents::IsTrackADescendantOfD0Daughters(int trackid, int D0num){ // This is for reco tracks (Needs to be called after getting the list of discarded tracks from MC)
  StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
  int mctrkid = trk->idTruth() - 1;

  if (std::find(fDroppedMCTracks[D0num].begin(), fDroppedMCTracks[D0num].end(), mctrkid) != fDroppedMCTracks[D0num].end()) return kTRUE;

  return kFALSE;
}

Bool_t StTagD0MCEvents::KeepTrack(int trackid) // This is for Reco Tracks
{
  bool keeptrack = kTRUE;

  StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
  
  int particleid = -99;

  if (abs(trk->nSigmaPion()) < 2. && abs(trk->nSigmaKaon()) > 2. && abs(trk->nSigmaProton()) > 2.) particleid = 1*trk->charge();
  else if (abs(trk->nSigmaPion()) > 2. && abs(trk->nSigmaKaon()) < 2. && abs(trk->nSigmaProton()) > 2.) particleid = 2*trk->charge();
  else if (abs(trk->nSigmaPion()) > 2. && abs(trk->nSigmaKaon()) > 2. && abs(trk->nSigmaProton()) < 2.) particleid = 3*trk->charge();

  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, Bfield);
  }

  double pt = mTrkMom.Perp();

  TRandom3 *r = new TRandom3(0);
  double rando = r->Rndm();

  // cout << "Random = " << rando << "\t" << fPionWeight->Eval(pt) << "\t" << fKaonWeight->Eval(pt) << "\t" << fProtonWeight->Eval(pt) << "\t" << fAProtonWeight->Eval(pt) << endl;

  if (particleid == 1 || particleid == -1) keeptrack = (rando > fPionWeight->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == 2 || particleid == -2) keeptrack = (rando > fKaonWeight->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == 3) keeptrack = (rando > fProtonWeight->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == -3) keeptrack = (rando > fAProtonWeight->Eval(pt)) ? kFALSE : kTRUE;

  return keeptrack;
}

TVector3 StTagD0MCEvents::FastSimMom(TVector3 p, int pid)
{
    float pt = p.Perp();
    float pt1 = pt;
    // if(pt1>2) pt1 = 2;//Used for high pt-hat bin smearing test
    if(pt1>10) pt1 = 10;//Used for high pt-hat bin smearing test
    float sPt = -1;

    TRandom3 *r = new TRandom3(0);
    
    if(pid ==8 || pid == 9)sPt = r->Gaus(pt, pt * fPionMomResolution->Eval(pt1));// Pion
    else if(pid ==11 || pid == 12)sPt = r->Gaus(pt, pt * fKaonMomResolution->Eval(pt1));// Kaon
    else if(pid ==15 || pid == 14)sPt = r->Gaus(pt, pt * fProtonMomResolution->Eval(pt1));// Proton
    else sPt = r->Gaus(pt, pt * fPionMomResolution->Eval(pt1));// Catch all: pions

    TVector3 smearedmom(sPt * cos(p.Phi()), sPt * sin(p.Phi()), sPt * sinh(p.PseudoRapidity()));
    
    return smearedmom;
}

