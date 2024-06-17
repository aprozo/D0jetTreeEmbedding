// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StSimJetSaver.h"
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
#include "StRoot/StPicoEvent/StPicoMcVertex.h"
#include "StRoot/StPicoEvent/StPicoMcTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// Bichsel includes
#include "StBichsel/Bichsel.h"


// D0 Includes
#include "StTagD0MCEvents.h"
#include "StMCD0JetMaker.h"
#include "StSimD0EventsJetMaker.h"


// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StSimJetSaver)

//________________________________________________________________________
StSimJetSaver::StSimJetSaver(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* mcjetMakerName = "", const char* recojetMakerName = "", const char *d0taggername = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  mcJets = 0x0 ;
  recoJets = 0x0;
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
  fRunFlag = 0;       // see StSimJetSaver::fRunFlagEnum
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
  fSmearFactor = 2.7;
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
  fMCJetMakerName = mcjetMakerName;
  fRecoJetMakerName = recojetMakerName;
  fD0TaggerName = d0taggername;
  fRhoMakerName = rhoMakerName;

  mcJetsArr.clear();
  recoJetsArr.clear();

  d0MCTrackIndices.clear();
  d0RecoTrackIndices.clear();
  d0Reco4Momenta.clear();

  d0pairs = 0;

  fMCJetTree = {0};
  fRecoJetTree = {0};

  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

}

//
//________________________________________________________________________
StSimJetSaver::~StSimJetSaver()
{ /*  */
  // destructor
  if(hCentrality)  delete hCentrality;
  if(hMultiplicity)delete hMultiplicity;
  if(hJetPt)       delete hJetPt;
  if(hJetCorrPt)   delete hJetCorrPt;

  if (hMCvRecoJetPt) delete hMCvRecoJetPt;
  if (hMCvRecoD0Pt) delete hMCvRecoD0Pt;

  if (hDeltaEtaDeltaPhi) delete hDeltaEtaDeltaPhi;
  if (hDR) delete hDR;

  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StSimJetSaver::Init() {
  StJetFrameworkPicoBase::Init();

  // declare histograms
  DeclareHistograms();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  // fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  mcJets = new TClonesArray("StJet");
  recoJets = new TClonesArray("StJet");
  //fJets->SetName(fJetsName);

  fMCJetTree = {0};
  fRecoJetTree = {0};

  DeclareTree();
  BookTree();

  // BookTree(mcjettree, fMCJetTree);
  // BookTree(recojettree, fRecoJetTree);

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StSimJetSaver::Finish() { 
  cout << "StSimJetSaver::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());

    // fout->mkdir(Form("%s/%s", GetName(), "MC"));
    // fout->mkdir(Form("%s/%s", GetName(), "Reco"));
    // fout->mkdir(Form("%s/%s", GetName(), "Histograms"));

    // fout->cd(Form("%s/%s", GetName(), "MC"));
    fout->cd(GetName());
    WriteTree(mcjettree);

    // fout->cd(Form("%s/%s", GetName(), "Reco"));
    fout->cd(GetName());
    WriteTree(recojettree);

    // fout->cd(Form("%s/%s", GetName(), "Histograms"));
    // fout->cd(GetName());
    // WriteHistograms();

    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StSimJetSaver::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StSimJetSaver::DeclareHistograms() {
  // binning for cent histograms
  int nHistCentBins = 20;

  // binning for mult histograms
  double kHistMultMax = 800.;
  int kHistMultBins = 400;

  // pp specific settings
  if(doppAnalysis) {
    kHistMultMax = 100.;
    kHistMultBins = 100;
  }

  // histograms
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);

  hMCvRecoJetPt = new TH2F("hMCvRecoJetPt", "MC v Reco Jet Pt", 110, -10, 100, 110, -10, 100);
  hRecovMCJetPt = new TH2F("hRecovMCJetPt", "Reco v MC Jet Pt", 110, -10, 100, 110, -10, 100);

  hMCvRecoSmearedJetPt = new TH2F("hMCvRecoSmearedJetPt", "MC v Reco Smeared Jet Pt", 110, -10, 100, 110, -10, 100);
  hRecoSmearedvMCJetPt = new TH2F("hRecoSmearedvMCJetPt", "Reco Smeared v MC Jet Pt", 110, -10, 100, 110, -10, 100);

  hRecoJetPt = new TH1F("hRecoJetPt", "Reco Jet p_{T}", 110, -10, 100);
  hRecoSmearedJetPt = new TH1F("hRecoSmearedJetPt", "Reco Smeared Jet p_{T}", 110, -10, 100);

  hMCJetPt = new TH1F("hMCJetPt", "MC Jet p_{T}", 110, -10, 100);

  hDeltaEtaDeltaPhi = new TH2F("hDeltaEtaDeltaPhi", "Delta Eta v Delta Phi", 100, 0.0, 0.4, 100, 0.0, 0.4);
  hDR = new TH1F("hDR", "Delta R", 100, 0.0, 0.4);

  hMCvRecoD0Pt = new TH2F("hMCvRecoD0Pt", "MC v Reco D0 Pt", 110, -10, 100, 110, -10, 100);

  hMCD0Pt = new TH1F("hMCD0Pt", "MC D0 p_{T}", 100, 0, 100);
  hRecoD0Pt = new TH1F("hRecoD0Pt", "Reco D0 p_{T}", 100, 0, 100);


  SetSumw2();
}
//
// write histograms
//_____________________________________________________________________________
void StSimJetSaver::WriteHistograms() {
  // writing of histograms done here
  // hCentrality->Write();
  // hMultiplicity->Write();
  hJetPt->Write();
  hJetCorrPt->Write();

  hMCvRecoJetPt->Write();
  hRecovMCJetPt->Write();
  hRecoJetPt->Write();

  hMCvRecoSmearedJetPt->Write();
  hRecoSmearedvMCJetPt->Write();
  hRecoSmearedJetPt->Write();

  hMCJetPt->Write();

  hMCvRecoD0Pt->Write();

  hDeltaEtaDeltaPhi->Write();
  hDR->Write();

  hMCD0Pt->Write();
  hRecoD0Pt->Write();
}

void StSimJetSaver::WriteTree(TTree *sometree) {
  sometree->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StSimJetSaver::Clear(Option_t *opt) {
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StSimJetSaver::Make() {
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

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  mD0Tagger = static_cast<StTagD0MCEvents*>(GetMaker(fD0TaggerName));
  if(!mD0Tagger) {
    LOG_WARN << " No D0Tagger! Skip! " << endm;
    return kStWarn;
  }

  d0MCTrackIndices = mD0Tagger->GetD0MCIndices();
  d0RecoTrackIndices = mD0Tagger->GetD0RecoIndices();

  pionReco4Momenta = mD0Tagger->GetPionReco4Momenta();
  kaonReco4Momenta = mD0Tagger->GetKaonReco4Momenta();
  d0Reco4Momenta = mD0Tagger->GetD0Reco4Momenta();

  trackssmearedwithfastsim = mD0Tagger->GetTracksSmearedWell();

  if (d0MCTrackIndices.size() != d0RecoTrackIndices.size()) {
    cout << d0MCTrackIndices.size() << "\t" << d0RecoTrackIndices.size() << endl;
    cout << "Yo wut!!!!!" << endl;
    return kStWarn;
  }

  d0pairs = d0MCTrackIndices.size();

  // if (d0pairs == 0) return kStOK;

  mcJetMaker = static_cast<StMCD0JetMaker*>(GetMaker(fMCJetMakerName));
  const char *mcJetMakerNameCh = fMCJetMakerName;
  if(!mcJetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", mcJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get jet collection associated with JetMaker
  mcJetsArr = mcJetMaker->GetJets();
  

  recoJetMaker = static_cast<StSimD0EventsJetMaker*>(GetMaker(fRecoJetMakerName));
  const char *recoJetMakerNameCh = fRecoJetMakerName;
  if(!recoJetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", recoJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get jet collection associated with JetMaker
  recoJetsArr = recoJetMaker->GetJets();

  // cout << "Jet Array Size = " << mcJetsArr.size() << "\t" << recoJetsArr.size() << endl;

  // RunJets();

  SaveMCJets();
  SaveRecoJets();

  // RunJets();

  return kStOK;
}

void StSimJetSaver::DeclareTree() {
  TString mctreename = "MCJets";
  TString recotreename = "RecoJets";

  mcjettree = new TTree(mctreename.Data(), mctreename.Data());
  recojettree = new TTree(recotreename.Data(), recotreename.Data());
}

void StSimJetSaver::BookTree()
{
  // Branches to save event info
  mcjettree->Branch("RunID", &fMCJetTree.runid, "runid/I");
  mcjettree->Branch("EventId", &fMCJetTree.eventid, "eventid/I");
  mcjettree->Branch("RefMult", &fMCJetTree.refmult, "refmult/F");
  mcjettree->Branch("Centrality", &fMCJetTree.centrality, "centrality/F");
  mcjettree->Branch("Triggers", &fMCJetTree.triggers);
  mcjettree->Branch("PrimaryVertex", &fMCJetTree.primaryvertex);
  mcjettree->Branch("PrimaryVertexErr", &fMCJetTree.primaryvertexerror);
  mcjettree->Branch("JetPt", &fMCJetTree.jetpt, "jetpt/F");
  mcjettree->Branch("JetCorrPt", &fMCJetTree.jetcorrectedpt, "jetcorrectedpt/F");
  mcjettree->Branch("JetEta", &fMCJetTree.jeteta, "jeteta/F");
  mcjettree->Branch("JetPhi", &fMCJetTree.jetphi, "jetphi/F");
  mcjettree->Branch("JetArea", &fMCJetTree.jetarea, "jetarea/F");
  mcjettree->Branch("JetRadius", &fMCJetTree.jetradius, "jetradius/F");
  mcjettree->Branch("JetE", &fMCJetTree.jetenergy, "jetenergy/F");
  mcjettree->Branch("JetNEF", &fMCJetTree.jetnef, "jetnef/F");
  mcjettree->Branch("JetRhoVal", &fMCJetTree.fRhoValforjet, "fRhoValforjet/F");
  mcjettree->Branch("JetHighestTrackPt", &fMCJetTree.jethighesttrackpt, "jethighesttrackpt/F");
  mcjettree->Branch("JetNConst", &fMCJetTree.numberofconstituents, "numberofconstituents/I");
  mcjettree->Branch("D0Mass", &fMCJetTree.d0mass, "d0mass/F");
  mcjettree->Branch("D0Pt", &fMCJetTree.d0pt, "d0pt/F");
  mcjettree->Branch("PionPt", &fMCJetTree.pionpt, "pionpt/F");
  mcjettree->Branch("PionEta", &fMCJetTree.pioneta, "pioneta/F");
  mcjettree->Branch("PionPhi", &fMCJetTree.pionphi, "pionphi/F");
  mcjettree->Branch("PionCharge", &fMCJetTree.pioncharge, "pioncharge/F");
  mcjettree->Branch("KaonPt", &fMCJetTree.kaonpt, "kaonpt/F");
  mcjettree->Branch("KaonEta", &fMCJetTree.kaoneta, "kaoneta/F");
  mcjettree->Branch("KaonPhi", &fMCJetTree.kaonphi, "kaonphi/F");
  mcjettree->Branch("KaonCharge", &fMCJetTree.kaoncharge, "kaoncharge/F");
  mcjettree->Branch("TrackID", &fMCJetTree.mTrackID, "mTrackID[numberofconstituents]/F");
  mcjettree->Branch("TrackPt", &fMCJetTree.mTrackPt, "mTrackPt[numberofconstituents]/F");
  mcjettree->Branch("TrackEta", &fMCJetTree.mTrackEta, "mTrackEta[numberofconstituents]/F");
  mcjettree->Branch("TrackPhi", &fMCJetTree.mTrackPhi, "mTrackPhi[numberofconstituents]/F");
  mcjettree->Branch("TrackPx", &fMCJetTree.mTrackPx, "mTrackPx[numberofconstituents]/F");
  mcjettree->Branch("TrackPy", &fMCJetTree.mTrackPy, "mTrackPy[numberofconstituents]/F");
  mcjettree->Branch("TrackPz", &fMCJetTree.mTrackPz, "mTrackPz[numberofconstituents]/F");
  mcjettree->Branch("TrackCharge", &fMCJetTree.mTrackCharge, "mTrackCharge[numberofconstituents]/F");

  recojettree->Branch("RunID", &fRecoJetTree.runid, "runid/I");
  recojettree->Branch("EventId", &fRecoJetTree.eventid, "eventid/I");
  recojettree->Branch("RefMult", &fRecoJetTree.refmult, "refmult/F");
  recojettree->Branch("Centrality", &fRecoJetTree.centrality, "centrality/F");
  recojettree->Branch("Triggers", &fRecoJetTree.triggers);
  recojettree->Branch("PrimaryVertex", &fRecoJetTree.primaryvertex);
  recojettree->Branch("PrimaryVertexErr", &fRecoJetTree.primaryvertexerror);
  recojettree->Branch("JetPt", &fRecoJetTree.jetpt, "jetpt/F");
  recojettree->Branch("JetCorrPt", &fRecoJetTree.jetcorrectedpt, "jetcorrectedpt/F");
  recojettree->Branch("JetEta", &fRecoJetTree.jeteta, "jeteta/F");
  recojettree->Branch("JetPhi", &fRecoJetTree.jetphi, "jetphi/F");
  recojettree->Branch("JetArea", &fRecoJetTree.jetarea, "jetarea/F");
  recojettree->Branch("JetRadius", &fRecoJetTree.jetradius, "jetradius/F");
  recojettree->Branch("JetE", &fRecoJetTree.jetenergy, "jetenergy/F");
  recojettree->Branch("JetNEF", &fRecoJetTree.jetnef, "jetnef/F");
  recojettree->Branch("JetRhoVal", &fRecoJetTree.fRhoValforjet, "fRhoValforjet/F");
  recojettree->Branch("JetHighestTrackPt", &fRecoJetTree.jethighesttrackpt, "jethighesttrackpt/F");
  recojettree->Branch("JetNConst", &fRecoJetTree.numberofconstituents, "numberofconstituents/I");
  recojettree->Branch("D0Mass", &fRecoJetTree.d0mass, "d0mass/F");
  recojettree->Branch("D0Pt", &fRecoJetTree.d0pt, "d0pt/F");
  recojettree->Branch("PionPt", &fRecoJetTree.pionpt, "pionpt/F");
  recojettree->Branch("PionEta", &fRecoJetTree.pioneta, "pioneta/F");
  recojettree->Branch("PionPhi", &fRecoJetTree.pionphi, "pionphi/F");
  recojettree->Branch("PionCharge", &fRecoJetTree.pioncharge, "pioncharge/F");
  recojettree->Branch("KaonPt", &fRecoJetTree.kaonpt, "kaonpt/F");
  recojettree->Branch("KaonEta", &fRecoJetTree.kaoneta, "kaoneta/F");
  recojettree->Branch("KaonPhi", &fRecoJetTree.kaonphi, "kaonphi/F");
  recojettree->Branch("KaonCharge", &fRecoJetTree.kaoncharge, "kaoncharge/F");
  recojettree->Branch("TrackID", &fRecoJetTree.mTrackID, "mTrackID[numberofconstituents]/F");
  recojettree->Branch("TrackPt", &fRecoJetTree.mTrackPt, "mTrackPt[numberofconstituents]/F");
  recojettree->Branch("TrackEta", &fRecoJetTree.mTrackEta, "mTrackEta[numberofconstituents]/F");
  recojettree->Branch("TrackPhi", &fRecoJetTree.mTrackPhi, "mTrackPhi[numberofconstituents]/F");
  recojettree->Branch("TrackPx", &fRecoJetTree.mTrackPx, "mTrackPx[numberofconstituents]/F");
  recojettree->Branch("TrackPy", &fRecoJetTree.mTrackPy, "mTrackPy[numberofconstituents]/F");
  recojettree->Branch("TrackPz", &fRecoJetTree.mTrackPz, "mTrackPz[numberofconstituents]/F");
  recojettree->Branch("TrackCharge", &fRecoJetTree.mTrackCharge, "mTrackCharge[numberofconstituents]/F");
}


// void StSimJetSaver::BookTree(TTree *jettree, StJetTreeStruct jetstruct)
// {
//   // Branches to save event info
//   jettree->Branch("RunID", &jetstruct.runid, "runid/I");
//   jettree->Branch("EventId", &jetstruct.eventid, "eventid/I");
//   jettree->Branch("RefMult", &jetstruct.refmult, "refmult/F");
//   jettree->Branch("Centrality", &jetstruct.centrality, "centrality/F");
//   jettree->Branch("Triggers", &jetstruct.triggers);
//   jettree->Branch("PrimaryVertex", &jetstruct.primaryvertex);
//   jettree->Branch("PrimaryVertexErr", &jetstruct.primaryvertexerror);
//   jettree->Branch("JetPt", &jetstruct.jetpt, "jetpt/F");
//   jettree->Branch("JetCorrPt", &jetstruct.jetcorrectedpt, "jetcorrectedpt/F");
//   jettree->Branch("JetEta", &jetstruct.jeteta, "jeteta/F");
//   jettree->Branch("JetPhi", &jetstruct.jetphi, "jetphi/F");
//   jettree->Branch("JetArea", &jetstruct.jetarea, "jetarea/F");
//   jettree->Branch("JetRadius", &jetstruct.jetradius, "jetradius/F");
//   jettree->Branch("JetE", &jetstruct.jetenergy, "jetenergy/F");
//   jettree->Branch("JetNEF", &jetstruct.jetnef, "jetnef/F");
//   jettree->Branch("JetRhoVal", &jetstruct.fRhoValforjet, "fRhoValforjet/F");
//   jettree->Branch("JetHighestTrackPt", &jetstruct.jethighesttrackpt, "jethighesttrackpt/F");
//   jettree->Branch("JetNConst", &jetstruct.numberofconstituents, "numberofconstituents/I");
//   jettree->Branch("D0Mass", &jetstruct.d0mass, "d0mass/F");
//   jettree->Branch("PionPt", &jetstruct.pionpt, "pionpt/F");
//   jettree->Branch("PionEta", &jetstruct.pioneta, "pioneta/F");
//   jettree->Branch("PionPhi", &jetstruct.pionphi, "pionphi/F");
//   jettree->Branch("PionCharge", &jetstruct.pioncharge, "pioncharge/F");
//   jettree->Branch("KaonPt", &jetstruct.kaonpt, "kaonpt/F");
//   jettree->Branch("KaonEta", &jetstruct.kaoneta, "kaoneta/F");
//   jettree->Branch("KaonPhi", &jetstruct.kaonphi, "kaonphi/F");
//   jettree->Branch("KaonCharge", &jetstruct.kaoncharge, "kaoncharge/F");
//   jettree->Branch("TrackID", &jetstruct.mTrackID, "mTrackID[numberofconstituents]/F");
//   jettree->Branch("TrackPt", &jetstruct.mTrackPt, "mTrackPt[numberofconstituents]/F");
//   jettree->Branch("TrackEta", &jetstruct.mTrackEta, "mTrackEta[numberofconstituents]/F");
//   jettree->Branch("TrackPhi", &jetstruct.mTrackPhi, "mTrackPhi[numberofconstituents]/F");
//   jettree->Branch("TrackPx", &jetstruct.mTrackPx, "mTrackPx[numberofconstituents]/F");
//   jettree->Branch("TrackPy", &jetstruct.mTrackPy, "mTrackPy[numberofconstituents]/F");
//   jettree->Branch("TrackPz", &jetstruct.mTrackPz, "mTrackPz[numberofconstituents]/F");
//   jettree->Branch("TrackCharge", &jetstruct.mTrackCharge, "mTrackCharge[numberofconstituents]/F");
// }

void StSimJetSaver::SaveMCJets(){

  // cout << " ****************** SAVE MC JETS TEST *******************" << endl;

  for (int jetcollection = 0; jetcollection < d0pairs; jetcollection++){

    if (mcJetsArr.size() == 0) continue;

    mcJets = new TClonesArray("StJet");
    mcJets = mcJetsArr[jetcollection];

    Int_t mcnjets = mcJets->GetEntries();

    if (mcnjets == 0) {
      fMCJetTree = {0};
      fMCJetTree.runid = mPicoEvent->runId();
      fMCJetTree.eventid = mPicoEvent->eventId() + 10000*(jetcollection);
      mcjettree->Fill();
      continue;
    }

    for(int jjet = 0; jjet < mcnjets; jjet++) {
      StJet *mjet = static_cast<StJet*>(mcJets->At(jjet));
      if(!mjet) continue;
      
      // if (mjet->GetNumberOfTracks() <= 1) continue;

      double jetarea = mjet->Area();
      double jetpt = mjet->Pt();
      double corrjetpt = mjet->Pt();
      double jetE = mjet->E();
      double jetEta = mjet->Eta();
      double jetPhi = mjet->Phi();
      double jetNEF = mjet->NEF();

      // get nTracks and maxTrackPt
      double maxtrackpt = mjet->GetMaxTrackPt();
      double NtrackConstit = mjet->GetNumberOfTracks();

      fMCJetTree = {0};

      // Event level information filler

      fMCJetTree.runid = mPicoEvent->runId();
      fMCJetTree.eventid = mPicoEvent->eventId() + 10000*(jetcollection);
      fMCJetTree.refmult = 0.;
      fMCJetTree.centrality = 0.;
      fMCJetTree.triggers = mPicoEvent->triggerIds();

      vector<double> pv;
      pv.clear();
      pv.push_back(mPicoEvent->primaryVertex().X());
      pv.push_back(mPicoEvent->primaryVertex().Y());
      pv.push_back(mPicoEvent->primaryVertex().Z());

      fMCJetTree.primaryvertex = pv;

      vector<double> pverr;
      pverr.clear();
      pverr.push_back(mPicoEvent->primaryVertexError().X());
      pverr.push_back(mPicoEvent->primaryVertexError().Y());
      pverr.push_back(mPicoEvent->primaryVertexError().Z());

      fMCJetTree.primaryvertexerror = pverr;

      // Jet level information filler

      fMCJetTree.jetpt = jetpt;
      fMCJetTree.jetcorrectedpt = corrjetpt;
      fMCJetTree.jeteta = jetEta;
      fMCJetTree.jetphi = jetPhi;
      fMCJetTree.jetarea = jetarea;
      fMCJetTree.jetradius = 0.4;
      fMCJetTree.jetenergy = jetE;
      fMCJetTree.jetnef = jetNEF;
      fMCJetTree.fRhoValforjet = 0.;
      fMCJetTree.jethighesttrackpt = maxtrackpt;
      fMCJetTree.numberofconstituents = NtrackConstit;

      fMCJetTree.d0mass = 1.865;

      for (int jtrk = 0; jtrk < mjet->GetNumberOfTracks(); jtrk++){
        int trackid = mjet->TrackAt(jtrk);

        if (trackid == 99999){

          StPicoMcTrack *trk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(d0MCTrackIndices[jetcollection][0]));
          StPicoMcTrack *trk2 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(d0MCTrackIndices[jetcollection][1]));

          TVector3 mTrk1Mom, mTrk2Mom, mTrkD0Mom;

          mTrk1Mom = trk1->p();
          mTrk2Mom = trk2->p();

          mTrkD0Mom = mTrk1Mom + mTrk2Mom;

          fMCJetTree.d0pt = mTrkD0Mom.Perp(); 
          

          fMCJetTree.pionpt = mTrk1Mom.Perp();
          fMCJetTree.pioneta = mTrk1Mom.PseudoRapidity();
          fMCJetTree.pionphi = standardPhi(mTrk1Mom.Phi());
          fMCJetTree.pioncharge = trk1->charge();

          fMCJetTree.kaonpt = mTrk2Mom.Perp();
          fMCJetTree.kaoneta = mTrk2Mom.PseudoRapidity();
          fMCJetTree.kaonphi = standardPhi(mTrk2Mom.Phi());
          fMCJetTree.kaoncharge = trk2->charge();

          double pt = mTrkD0Mom.Perp();
          double phi = mTrkD0Mom.Phi();
          double eta = mTrkD0Mom.PseudoRapidity();
          double px = mTrkD0Mom.x();
          double py = mTrkD0Mom.y();
          double pz = mTrkD0Mom.z();
          short charge = 0;

          fMCJetTree.mTrackID[jtrk] = 421;
          fMCJetTree.mTrackPt[jtrk] = pt;
          fMCJetTree.mTrackEta[jtrk] = eta;
          fMCJetTree.mTrackPhi[jtrk] = standardPhi(phi);
          fMCJetTree.mTrackPx[jtrk] = px;
          fMCJetTree.mTrackPy[jtrk] = py;
          fMCJetTree.mTrackPz[jtrk] = pz;
          fMCJetTree.mTrackCharge[jtrk] = charge;
        } 

        else{
          StPicoMcTrack *trk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(trackid));

          TVector3 mTrkMom;
          mTrkMom = trk->p();

          double pt = mTrkMom.Perp();
          double phi = mTrkMom.Phi();
          double eta = mTrkMom.PseudoRapidity();
          double px = mTrkMom.x();
          double py = mTrkMom.y();
          double pz = mTrkMom.z();
          short charge = trk->charge();

          fMCJetTree.mTrackID[jtrk] = trk->geantId();
          fMCJetTree.mTrackPt[jtrk] = pt;
          fMCJetTree.mTrackEta[jtrk] = eta;
          fMCJetTree.mTrackPhi[jtrk] = standardPhi(phi);
          fMCJetTree.mTrackPx[jtrk] = px;
          fMCJetTree.mTrackPy[jtrk] = py;
          fMCJetTree.mTrackPz[jtrk] = pz;
          fMCJetTree.mTrackCharge[jtrk] = charge;
        }
      }

      mcjettree->Fill();

      // cout << "Filled " << jjet << " from collection # " << jetcollection << endl;
    }
  }
}

void StSimJetSaver::SaveRecoJets(){

  // cout << " ****************** SAVE RECO JETS TEST *******************" << endl;

  TRandom3 *r3 = new TRandom3(0);

  for (int jetcollection = 0; jetcollection < d0pairs; jetcollection++){
    if (recoJetsArr.size() == 0) continue;

    recoJets = new TClonesArray("StJet");
    recoJets = recoJetsArr[jetcollection];

    Int_t reconjets = recoJets->GetEntries();

    // cout << "Checkpoint 0" << endl;

    if (reconjets == 0) {
      fRecoJetTree = {0};
      fRecoJetTree.runid = mPicoEvent->runId();
      fRecoJetTree.eventid = mPicoEvent->eventId() + 10000*(jetcollection);
      recojettree->Fill();
      continue;
    }

    // cout << "Reco #Jets = " << reconjets << endl;
    // cout << "Checkpoint 1" << endl;


    for(int ijet = 0; ijet < reconjets; ijet++) {  // RECO JET LOOP
      // get jet pointer
      // cout << "Checkpoint 2" << endl;

      StJet *rjet = static_cast<StJet*>(recoJets->At(ijet));
      if(!rjet) continue;

      // if (rjet->GetNumberOfTracks() <= 1) continue;

      // get some jet parameters
      double jetarea = rjet->Area();
      double jetpt = rjet->Pt();
      double corrjetpt = rjet->Pt() + r3->Gaus(0., fSmearFactor);
      double jetE = rjet->E();
      double jetEta = rjet->Eta();
      double jetPhi = rjet->Phi();
      double jetNEF = rjet->NEF();

      // cout << "Jet " << jetpt  << "\t" << jetEta << endl; 
      // cout << "Checkpoint 3" << endl;
      // get nTracks and maxTrackPt
      double maxtrackpt = rjet->GetMaxTrackPt();
      double NtrackConstit = rjet->GetNumberOfTracks();

      fRecoJetTree = {0};

      // Event level information filler

      fRecoJetTree.runid = mPicoEvent->runId();
      fRecoJetTree.eventid = mPicoEvent->eventId() + 10000*(jetcollection);
      fRecoJetTree.refmult = 0.;
      fRecoJetTree.centrality = 0.;
      fRecoJetTree.triggers = mPicoEvent->triggerIds();

      vector<double> pv;
      pv.clear();
      pv.push_back(mPicoEvent->primaryVertex().X());
      pv.push_back(mPicoEvent->primaryVertex().Y());
      pv.push_back(mPicoEvent->primaryVertex().Z());

      fRecoJetTree.primaryvertex = pv;

      vector<double> pverr;
      pverr.clear();
      pverr.push_back(mPicoEvent->primaryVertexError().X());
      pverr.push_back(mPicoEvent->primaryVertexError().Y());
      pverr.push_back(mPicoEvent->primaryVertexError().Z());

      fRecoJetTree.primaryvertexerror = pverr;

      // Jet level information filler

      fRecoJetTree.jetpt = jetpt;
      fRecoJetTree.jetcorrectedpt = corrjetpt;
      fRecoJetTree.jeteta = jetEta;
      fRecoJetTree.jetphi = jetPhi;
      fRecoJetTree.jetarea = jetarea;
      fRecoJetTree.jetradius = 0.4;
      fRecoJetTree.jetenergy = jetE;
      fRecoJetTree.jetnef = jetNEF;
      fRecoJetTree.fRhoValforjet = 0.;
      fRecoJetTree.jethighesttrackpt = maxtrackpt;
      fRecoJetTree.numberofconstituents = NtrackConstit;

      // cout << "Jet " << fRecoJetTree.jetnef << "\t" << fRecoJetTree.jethighesttrackpt << "\t" << fRecoJetTree.numberofconstituents << "\t" << fRecoJetTree.jetpt  << "\t" << fRecoJetTree.jeteta << endl; 
      // cout << "Checkpoint 4" << endl;

      // fRecoJetTree.d0mass = 1.865;

      for (int itrk = 0; itrk < rjet->GetNumberOfTracks(); itrk++){
        int trackid = rjet->TrackAt(itrk);

        if (trackid == 99999){
          if (d0RecoTrackIndices[jetcollection][0] < -999 && d0RecoTrackIndices[jetcollection][1] < -999 ){

            StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(d0RecoTrackIndices[jetcollection][0]));
            StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(d0RecoTrackIndices[jetcollection][1]));

            TVector3 mTrkD0Mom;
            if(doUsePrimTracks) { 
              // get primary track vector
              mTrkD0Mom = trk1->pMom() + trk2->pMom(); 
            } else { 
              // get global track vector
              mTrkD0Mom = trk1->gMom(mVertex, Bfield) + trk2->gMom(mVertex, Bfield); 
            }

            fRecoJetTree.d0mass = 1.865;

            fRecoJetTree.d0pt = mTrkD0Mom.Perp(); 

            fRecoJetTree.pionpt = 0.;
            fRecoJetTree.pioneta = 0.;
            fRecoJetTree.pionphi = 0.;
            fRecoJetTree.pioncharge = 0.;

            fRecoJetTree.kaonpt = 0.;
            fRecoJetTree.kaoneta = 0.;
            fRecoJetTree.kaonphi = 0.;
            fRecoJetTree.kaoncharge = 0.;

            double pt = mTrkD0Mom.Perp();
            double phi = mTrkD0Mom.Phi();
            double eta = mTrkD0Mom.PseudoRapidity();
            double px = mTrkD0Mom.x();
            double py = mTrkD0Mom.y();
            double pz = mTrkD0Mom.z();
            short charge = 0;

            fRecoJetTree.mTrackID[itrk] = 421;
            fRecoJetTree.mTrackPt[itrk] = pt;
            fRecoJetTree.mTrackEta[itrk] = eta;
            fRecoJetTree.mTrackPhi[itrk] = standardPhi(phi);
            fRecoJetTree.mTrackPx[itrk] = px;
            fRecoJetTree.mTrackPy[itrk] = py;
            fRecoJetTree.mTrackPz[itrk] = pz;
            fRecoJetTree.mTrackCharge[itrk] = charge;

          }

          else{

            TVector3 mTrkD0Mom, mPionMom, mKaonMom;
            mPionMom = pionReco4Momenta[jetcollection].Vect();
            mKaonMom = kaonReco4Momenta[jetcollection].Vect();
            mTrkD0Mom = d0Reco4Momenta[jetcollection].Vect();

            fRecoJetTree.d0mass = d0Reco4Momenta[jetcollection].M();
            
            fRecoJetTree.d0pt = mTrkD0Mom.Perp(); 
            // cout << "D0 Mass = " << fRecoJetTree.d0mass << endl;

            fRecoJetTree.pionpt = mPionMom.Perp();
            fRecoJetTree.pioneta = mPionMom.PseudoRapidity();
            fRecoJetTree.pionphi = standardPhi(mPionMom.Phi());
            fRecoJetTree.pioncharge = 1;

            fRecoJetTree.kaonpt = mKaonMom.Perp();
            fRecoJetTree.kaoneta = mKaonMom.PseudoRapidity();
            fRecoJetTree.kaonphi = standardPhi(mKaonMom.Phi());
            fRecoJetTree.kaoncharge = -1;

            double pt = mTrkD0Mom.Perp();
            double phi = mTrkD0Mom.Phi();
            double eta = mTrkD0Mom.PseudoRapidity();
            double px = mTrkD0Mom.x();
            double py = mTrkD0Mom.y();
            double pz = mTrkD0Mom.z();
            short charge = 0;

            fRecoJetTree.mTrackID[itrk] = 421;
            fRecoJetTree.mTrackPt[itrk] = pt;
            fRecoJetTree.mTrackEta[itrk] = eta;
            fRecoJetTree.mTrackPhi[itrk] = standardPhi(phi);
            fRecoJetTree.mTrackPx[itrk] = px;
            fRecoJetTree.mTrackPy[itrk] = py;
            fRecoJetTree.mTrackPz[itrk] = pz;
            fRecoJetTree.mTrackCharge[itrk] = charge;

            // cout << "D0 " << fRecoJetTree.mTrackPt[itrk]  << "\t" << fRecoJetTree.mTrackEta[itrk] << "\t" << fRecoJetTree.mTrackCharge[itrk] << endl; 
            // cout << "Checkpoint 5" << endl;

          }
        }

        else{

          StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
          if(!trk){ continue; }

          double pt = trackssmearedwithfastsim[trackid][0];
          double phi = trackssmearedwithfastsim[trackid][1];
          double eta = trackssmearedwithfastsim[trackid][2];
          double px = trackssmearedwithfastsim[trackid][3];
          double py = trackssmearedwithfastsim[trackid][4];
          double pz = trackssmearedwithfastsim[trackid][5];
          short charge = trackssmearedwithfastsim[trackid][8];

          // Track level information fill

          fRecoJetTree.mTrackID[itrk] = IsWhatParticle(trk);
          fRecoJetTree.mTrackPt[itrk] = pt;
          fRecoJetTree.mTrackEta[itrk] = eta;
          fRecoJetTree.mTrackPhi[itrk] = standardPhi(phi);
          fRecoJetTree.mTrackPx[itrk] = px;
          fRecoJetTree.mTrackPy[itrk] = py;
          fRecoJetTree.mTrackPz[itrk] = pz;
          fRecoJetTree.mTrackCharge[itrk] = charge;
        }

        /*
        else{
          StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
          if(!trk){ continue; }

          // get momentum vector
          TVector3 mTrkMom;
          if(doUsePrimTracks) {
          if(!(trk->isPrimary())) continue; // check if primary
          // get primary track vector
          mTrkMom = trk->pMom();
          } else {
          // get global track vector
          mTrkMom = trk->gMom(mVertex, Bfield);
          }

          // track variables
          double pt = mTrkMom.Perp();
          double phi = mTrkMom.Phi();
          double eta = mTrkMom.PseudoRapidity();
          double px = mTrkMom.x();
          double py = mTrkMom.y();
          double pz = mTrkMom.z();
          short charge = trk->charge();

          // Track level information fill

          fRecoJetTree.mTrackID[itrk] = IsWhatParticle(trk);
          fRecoJetTree.mTrackPt[itrk] = pt;
          fRecoJetTree.mTrackEta[itrk] = eta;
          fRecoJetTree.mTrackPhi[itrk] = standardPhi(phi);
          fRecoJetTree.mTrackPx[itrk] = px;
          fRecoJetTree.mTrackPy[itrk] = py;
          fRecoJetTree.mTrackPz[itrk] = pz;
          fRecoJetTree.mTrackCharge[itrk] = charge;

          // cout << "Tracks " << fRecoJetTree.mTrackID[itrk] << "\t" << fRecoJetTree.mTrackPt[itrk]  << "\t" << fRecoJetTree.mTrackEta[itrk] << "\t" << fRecoJetTree.mTrackCharge[itrk] << endl; 
          // cout << "Checkpoint 6" << endl;

        }
        */
      }

      // cout << "Struct " << fRecoJetTree.jetpt  << "\t" << fRecoJetTree.jeteta << "\t" << fRecoJetTree.mTrackPt[0] << endl;
      // cout << "Checkpoint 7" << endl;
      recojettree->Fill();
    }
  }
}

void StSimJetSaver::IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e){ // NEW PID APPROACH

  pid = 0;
  m = 0.0;
  e = 0.0;

  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, Bfield);
  }

  // track variables
  double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  double px = mTrkMom.x();
  double py = mTrkMom.y();
  double pz = mTrkMom.z();
  double p = mTrkMom.Mag();
  short charge = trk->charge();

  double dedx = trk->dEdx();
  double dedxresolution = trk->dEdxError();

  // bichsel function approximation
  double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
  double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
  double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

  //z - variables
  // double zpi = TMath::Log(dedx/dedxth_pi)/dedxresolution;
  // double zka = TMath::Log(dedx/dedxth_ka)/dedxresolution;
  // double zpr = TMath::Log(dedx/dedxth_pr)/dedxresolution;

  double zpi = trk->nSigmaPion();
  double zka = trk->nSigmaKaon();
  double zpr = trk->nSigmaProton();

  if ((abs(zpi) > 2) && (abs(zka) > 2) && (abs(zpr) > 2)) {pid = 0; m = 0.; e = 0.; return;}

  bool tpc_pion = kFALSE;
  bool tpc_kaon = kFALSE;
  bool tpc_proton = kFALSE;

  if (abs(zpi) <= abs(zka) && abs(zpi) <= abs(zpr)) tpc_pion = kTRUE;
  else if (abs(zka) < abs(zpi) && abs(zka) <= abs(zpr)) tpc_kaon = kTRUE;
  else if (abs(zpr) < abs(zpi) && abs(zpr) < abs(zka)) tpc_proton = kTRUE;

  bool toftrack = kFALSE;

  int tof_loc = trk->bTofPidTraitsIndex();

  if (tof_loc >= 0)
  {
    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    if (tofpointer){toftrack = kTRUE;}
  }

  if(toftrack){

    bool tof_pion = kFALSE;
    bool tof_kaon = kFALSE;
    bool tof_proton = kFALSE;

    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    double invbeta_from_tof = tofpointer->btofBeta();
    invbeta_from_tof = 1/invbeta_from_tof;

    double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
    double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
    double norm_invbeta_pr = TMath::Sqrt(pow(Mproton,2)/pow(p,2) + 1);

    double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.011;
    double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.011;
    double normalisedinvbeta_for_pr = (invbeta_from_tof-norm_invbeta_pr)/0.011;

    if ((abs(normalisedinvbeta_for_pi) > 2) && (abs(normalisedinvbeta_for_ka) > 2) && (abs(normalisedinvbeta_for_pr) > 2)) {pid = 0; m = 0.; e = 0.; return;}

    if (abs(normalisedinvbeta_for_pi) <= abs(normalisedinvbeta_for_ka) && abs(normalisedinvbeta_for_pi) <= abs(normalisedinvbeta_for_pr)) tof_pion = kTRUE;
    else if (abs(normalisedinvbeta_for_ka) < abs(normalisedinvbeta_for_pi) && abs(normalisedinvbeta_for_ka) <= abs(normalisedinvbeta_for_pr)) tof_kaon = kTRUE;
    else if (abs(normalisedinvbeta_for_pr) < abs(normalisedinvbeta_for_pi) && abs(normalisedinvbeta_for_pr) < abs(normalisedinvbeta_for_ka)) tof_proton = kTRUE;



    if (tof_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tof_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else if (tof_proton) {pid = 3*charge; m = Mproton; e = TMath::Sqrt(pow(p,2) + pow(Mproton, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}

  }

  else{
    if (tpc_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tpc_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else if (tpc_proton) {pid = 3*charge; m = Mproton; e = TMath::Sqrt(pow(p,2) + pow(Mproton, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}
  }
}


Int_t StSimJetSaver::IsWhatParticle(StPicoTrack *trk){ // Just to get the PID out
  int pid;
  double m; 
  double e;
  IsWhatParticle(trk, pid, m, e);
  return pid;
}

Double_t StSimJetSaver::standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard > 2*(TMath::Pi())) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}

//
//
//________________________________________________________________________
void StSimJetSaver::RunTracks()
{
  // constants: assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV
  unsigned int ntracks = mPicoDst->numberOfTracks();

  // loop over ALL tracks in PicoDst 
  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // primary track switch: get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = trk->pMom();
    } else {
      // get global track vector
      mTrkMom = trk->gMom(mVertex, Bfield);
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
    short charge = trk->charge();

    // do some things with tracks here
    // ......

    // fill track histograms here
    // .........

  } // track loop

}  // track function

//
//
//________________________________________________________________________
void StSimJetSaver::RunTowers()
{
  // constants: assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // looping over clusters - STAR: matching already done, get # of clusters and set variables
  unsigned int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits();

  // loop over ALL clusters in PicoDst
  for(unsigned short iClus = 0; iClus < nBEmcPidTraits; iClus++){
    StPicoBEmcPidTraits *cluster = static_cast<StPicoBEmcPidTraits*>(mPicoDst->bemcPidTraits(iClus));
    if(!cluster){ cout<<"Cluster pointer does not exist.. iClus = "<<iClus<<endl; continue; }

    // cluster and tower ID
    int clusID = cluster->bemcId();  // index in bemc point array
    int towID = cluster->btowId();   // projected tower Id: 1 - 4800
    int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    if(towID < 0) continue;

    // cluster and tower position - from vertex and ID
    TVector3 towPosition = mEmcPosition->getPosFromVertex(mVertex, towID);
    double towPhi = towPosition.Phi();
    double towEta = towPosition.PseudoRapidity();

    // matched track index
    int trackIndex = cluster->trackIndex();
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackIndex));
    if(!trk) { cout<<"No trk pointer...."<<endl; continue; }
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

  } // BEmc loop

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    // get tower pointer
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // tower ID: tower ID shifted by +1 from array index
    int towerID = itow + 1;
    if(towerID < 0) continue; // double check these aren't still in the event list

    // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    double towerPhi = towerPosition.Phi();
    double towerEta = towerPosition.PseudoRapidity();
    int towerADC = tower->adc();
    double towerE = tower->energy();
    double towerEt = towerE / (1.0*TMath::CosH(towerEta)); // THIS should be USED

    // do stuff with towers and fill histograms here
    // ........


  } // tower loop

}  // run towers function

//
//
// __________________________________________________________________________________
void StSimJetSaver::SetSumw2() {
  hCentrality->Sumw2();
  hMultiplicity->Sumw2();
  hJetPt->Sumw2();
  hJetCorrPt->Sumw2();
}

//
// Function: get relative phi of jet and track (-0.5pi, 1.5pi)
//________________________________________________________________________
Double_t StSimJetSaver::RelativePhi(Double_t mphi,Double_t vphi) const
{ // function to calculate relative PHI
  double dphi = mphi-vphi;

  // set dphi to operate on adjusted scale
  if(dphi < -0.5*pi) dphi += 2.0*TMath::Pi();
  if(dphi >  1.5*pi) dphi -= 2.0*TMath::Pi();

  // test check
  if( dphi < -0.5*pi || dphi > 1.5*pi )
    Form("%s: dPHI not in range [-0.5*Pi, 1.5*Pi]!", GetName());

  return dphi; // dphi in [-0.5Pi, 1.5Pi]                                                                                   
}

//
// Function: calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
//_________________________________________________________________________
Double_t StSimJetSaver::RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{
  Double_t pi = 1.0*TMath::Pi();
  Double_t dphi = 1.0*TMath::Abs(EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi < -pi ){
    dphi = dphi + pi;
  } // this assumes we are doing full jets currently 
 
  if(dphi > 1.5*pi) dphi -= 2.0*pi;
  if((dphi > 1.0*pi) && (dphi < 1.5*pi)) dphi -= 1.0*pi;
  if((dphi > 0.5*pi) && (dphi < 1.0*pi)) dphi -= 1.0*pi;
  dphi = 1.0*TMath::Abs(dphi);

  // test check
  if( dphi < 0.0 || dphi > 0.5*pi ) {
    //Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName());
    cout<<"dPhi not in range [0, 0.5*Pi]!"<<endl;
  }

  return dphi;   // dphi in [0, Pi/2]
}

//
//
//_________________________________________________________________________
void StSimJetSaver::FillEmcTriggers() {
  // number of Emcal Triggers
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // fill for valid triggers
    if(emcTrig->isHT0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(emcTrig->isHT1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(emcTrig->isHT2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(emcTrig->isHT3()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(emcTrig->isJP0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(emcTrig->isJP1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(emcTrig->isJP2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
}
//
// elems: sizeof(myarr)/sizeof(*myarr) prior to passing to function
// upon passing the array collapses to a pointer and can not get size anymore
//________________________________________________________________________
Bool_t StSimJetSaver::DoComparison(int myarr[], int elems) {
  //std::cout << "Length of array = " << (sizeof(myarr)/sizeof(*myarr)) << std::endl;
  bool match = kFALSE;

  // loop over specific physics selection array and compare to specific event trigger
  for(int i=0; i<elems; i++) {
    if(mPicoEvent->isTrigger(myarr[i])) match = kTRUE;
    if(match) break;
  }

  return match;
}
