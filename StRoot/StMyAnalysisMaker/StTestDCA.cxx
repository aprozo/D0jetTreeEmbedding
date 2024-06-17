#include "StTestDCA.h"
#include "StMemStat.h"
#include "phys_constants.h"
#include <limits>
#include "math.h"
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
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
// #include "StDcaGeometry.h"
#include "StPhysicalHelix.hh"
#include "StThreeVectorF.hh"


// #include "StarClassLibrary/StRedefinedHelix.hh"


// // KFParticle includes
// #include "StRoot/StEvent/StDcaGeometry.h"
// #include "StRoot/StarRoot/KFParticle.h"

#include "StBTofUtil/tofPathLength.hh"

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

ClassImp(StTestDCA)

//________________________________________________________________________
StTestDCA::StTestDCA(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{ 
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StTestDCA::fRunFlagEnum
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
  // mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0; fMBEventType = 2;
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

  numberofMB5events  = 0;
  numberofMB30events = 0;

}

//
//________________________________________________________________________
StTestDCA::~StTestDCA()
{ 

  
}

//________________________________________________________________________
Int_t StTestDCA::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StTestDCA::Finish() { 

  cout << "StTestDCA::Finish()\n";

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

  cout<<"End of StTestDCA::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StTestDCA::DeclareHistograms() {

  pVtx = new TH3F("pVtx", "pVtx", 100, -2, 2, 100, -2, 2, 300, -6, 6);

  piondcaxy = new TH3F("PionCent_Pt_DCA_XY", "PionCent_Pt_DCA_XY", 9, 0, 9, 200, 0, 20, 200, -0.02, 0.02);
  kaondcaxy = new TH3F("KaonCent_Pt_DCA_XY", "KaonCent_Pt_DCA_XY", 9, 0, 9, 200, 0, 20, 200, -0.02, 0.02);

  piondcaz = new TH3F("PionCent_Pt_DCA_Z", "PionCent_Pt_DCA_Z", 9, 0, 9, 200, 0, 20, 200, -0.02, 0.02);
  kaondcaz = new TH3F("KaonCent_Pt_DCA_Z", "KaonCent_Pt_DCA_Z", 9, 0, 9, 200, 0, 20, 200, -0.02, 0.02);
}
//
// write histograms
//_____________________________________________________________________________
void StTestDCA::WriteHistograms() {

  cout << numberofMB5events  << "\t" << numberofMB30events << endl;

  // for (int i=1; i <= 14; i++){
  //   cuthistogram_event->SetBinContent(i, numberofevents[i-1]);
  //   cuthistogram_event->GetXaxis()->SetBinLabel(i, event_cuts[i-1]);
  //   // cuthistogram_track->SetBinContent(i, numberoftracks[i-1]);
  //   // cuthistogram_track->GetXaxis()->SetBinLabel(i, track_cuts[i-1]);
  // }
  // cuthistogram_event->Write();

  // pVtx->Write();
  // piondcaxy->Write();
  // kaondcaxy->Write();
  // piondcaz->Write();
  // kaondcaz->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StTestDCA::Clear(Option_t *opt) {
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StTestDCA::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;
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

  // cout << "Vertex" << mVertex.X() << "\t" << mVertex.Y() << "\t" << mVertex.Z() << endl;

  // if (mVertex.x() == 0 || mVertex.y() == 0 || mVertex.z() == 0) return kStOK;

  if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5) return kStOK;

  double zVtx_VPD = mPicoEvent->vzVpd();

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
  int cent9 = mCentMaker->GetCent9();
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

    // cout << "Here" << endl;

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

    // check for MB/HT event
    bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
    bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
    bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
    bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);

    bool fRunForMB = kFALSE;
    if(doppAnalysis) fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
    if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

    // ======================== end of Triggers ============================= //

    if (abs(zVtx) > 6.) return kStOk;

    if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > 2.) return kStOK;

    if (abs(zVtx - zVtx_VPD) > 3) return kStOk;

    // int arrMB5_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
    int arrMB5_Run14[] = {450005, 450015, 450025, 450050, 450060};
    int arrMB30_Run14[] = {450008, 450009, 450014, 450018, 450024};
    int arrHT1_Run14[] = {450201, 450211, 460201};
    int arrHT2_Run14[] = {450202, 450212, 460202, 460212};
    int arrHT3_Run14[] = {450203, 450213, 460203};

    if (!doppAnalysis){

      bool matchMB5 = kFALSE;

      for(int i = 0; i < sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14); i++) {
        if(mPicoEvent->isTrigger(arrMB5_Run14[i])) matchMB5 = kTRUE;
        if(matchMB5) break;
      }

      if (matchMB5) numberofMB5events++;

      bool matchMB30 = kFALSE;

      for(int i = 0; i < sizeof(arrMB30_Run14)/sizeof(*arrMB30_Run14); i++) {
        if(mPicoEvent->isTrigger(arrMB30_Run14[i])) matchMB30 = kTRUE;
        if(matchMB30) break;
      }

      if (!matchMB5 && matchMB30) numberofMB30events++;

    }
  }


  // const Int_t ntracks = mPicoDst->numberOfTracks();
  // Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // cout << "Number of tracks found = " << ntracks << endl;

  // cout << Mpion << "\t" << M_PION_PLUS << "\t" << Mkaon << "\t" << M_KAON_PLUS << endl;
  // RunTracks();
  return kStOK;
}


//__________________________________________________________________________________________
  
void StTestDCA::RunTracks(){

  // cout << "Started Run Tracks" << endl;

  // cout << "Using primary tracks : " << std::boolalpha << doUsePrimTracks << endl; 

  pVtx->Fill(mVertex.X(), mVertex.Y(), mVertex.Z());

  const Int_t ntracks = mPicoDst->numberOfTracks();

  for(unsigned short itrk = 0; itrk < ntracks; itrk++){
    // get track pointer

    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
    if(!trk){ continue; }

    TVector3 mTrkMom;
    mTrkMom = trk->gMom(mVertex, Bfield);

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();
    double px = mTrkMom.x();
    double py = mTrkMom.y();
    double pz = mTrkMom.z();
    double p = mTrkMom.Mag();
    short charge = trk->charge();

    if (!trk->isHFTTrack()) continue;
    if (pt < 0.2) continue; 
    if (abs(eta) > 1) continue;
    if (trk->nHitsDedx() < 20) continue;
    if (charge == 0) continue;

    if (trk->nSigmaPion() < 2.0){
      double dcaxy = trk->gDCAxy(mVertex.X(), mVertex.Y());
      double dcaz = trk->gDCAz(mVertex.Z());

      piondcaxy->Fill(mCentMaker->GetCent9(), pt, dcaxy);
      piondcaz->Fill(mCentMaker->GetCent9(), pt, dcaz);
    }

    if (trk->nSigmaKaon() < 2.0){
      double dcaxy = trk->gDCAxy(mVertex.X(), mVertex.Y());
      double dcaz = trk->gDCAz(mVertex.Z());

      kaondcaxy->Fill(mCentMaker->GetCent9(), pt, dcaxy);
      kaondcaz->Fill(mCentMaker->GetCent9(), pt, dcaz);
    }
  }

}