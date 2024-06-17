// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StSimJetSaver_Test.h"
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

ClassImp(StSimJetSaver_Test)

//________________________________________________________________________
StSimJetSaver_Test::StSimJetSaver_Test(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* mcjetMakerName = "", const char* recojetMakerName = "", const char *d0taggername = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  fRunFlag = 0;       // see StSimJetSaver_Test::fRunFlagEnum
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
StSimJetSaver_Test::~StSimJetSaver_Test()
{ /*  */
  // destructor
  // if(hCentrality)  delete hCentrality;
  // if(hMultiplicity)delete hMultiplicity;
  // if(hJetPt)       delete hJetPt;
  // if(hJetCorrPt)   delete hJetCorrPt;

  // if (hMCvRecoJetPt) delete hMCvRecoJetPt;
  // if (hMCvRecoD0Pt) delete hMCvRecoD0Pt;

  // if (hDeltaEtaDeltaPhi) delete hDeltaEtaDeltaPhi;
  // if (hDR) delete hDR;

  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StSimJetSaver_Test::Init() {
  StJetFrameworkPicoBase::Init();

  // declare histograms
  // DeclareHistograms();

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
Int_t StSimJetSaver_Test::Finish() { 
  cout << "StSimJetSaver_Test::Finish()\n";

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
    WriteTree(jettree);

    // fout->cd(Form("%s/%s", GetName(), "Histograms"));
    // fout->cd(GetName());
    // WriteHistograms();

    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StSimJetSaver_Test::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StSimJetSaver_Test::DeclareHistograms() {
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


  // SetSumw2();
}
//
// write histograms
//_____________________________________________________________________________
void StSimJetSaver_Test::WriteHistograms() {
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

void StSimJetSaver_Test::WriteTree(TTree *sometree) {
  sometree->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StSimJetSaver_Test::Clear(Option_t *opt) {
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StSimJetSaver_Test::Make() {
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

  SaveJets();

  // RunJets();

  return kStOK;
}

void StSimJetSaver_Test::DeclareTree() {
  TString treename = "Jets";

  jettree = new TTree(treename.Data(), treename.Data());
}

void StSimJetSaver_Test::BookTree()
{
  jettree->Branch("MCPrimaryVertex", &fMCJetTree.primaryvertex);
  jettree->Branch("MCJetPt", &fMCJetTree.jetpt, "jetpt/F");
  jettree->Branch("MCJetEta", &fMCJetTree.jeteta, "jeteta/F");
  jettree->Branch("MCJetPhi", &fMCJetTree.jetphi, "jetphi/F");
  jettree->Branch("MCJetArea", &fMCJetTree.jetarea, "jetarea/F");
  jettree->Branch("MCJetE", &fMCJetTree.jetenergy, "jetenergy/F");
  jettree->Branch("MCJetNConst", &fMCJetTree.numberofconstituents, "numberofconstituents/I");
  jettree->Branch("MCD0Pt", &fMCJetTree.d0pt, "d0pt/F");
  jettree->Branch("MCD0Eta", &fMCJetTree.d0eta, "d0eta/F");
  jettree->Branch("MCD0Phi", &fMCJetTree.d0phi, "d0phi/F");
  jettree->Branch("MCPionPt", &fMCJetTree.pionpt, "pionpt/F");
  jettree->Branch("MCPionEta", &fMCJetTree.pioneta, "pioneta/F");
  jettree->Branch("MCPionPhi", &fMCJetTree.pionphi, "pionphi/F");
  jettree->Branch("MCKaonPt", &fMCJetTree.kaonpt, "kaonpt/F");
  jettree->Branch("MCKaonEta", &fMCJetTree.kaoneta, "kaoneta/F");
  jettree->Branch("MCKaonPhi", &fMCJetTree.kaonphi, "kaonphi/F");

  jettree->Branch("RecoJetPt", &fRecoJetTree.jetpt, "jetpt/F");
  jettree->Branch("RecoJetCorrPt", &fRecoJetTree.jetcorrectedpt, "jetcorrectedpt/F");
  jettree->Branch("RecoJetEta", &fRecoJetTree.jeteta, "jeteta/F");
  jettree->Branch("RecoJetPhi", &fRecoJetTree.jetphi, "jetphi/F");
  jettree->Branch("RecoJetArea", &fRecoJetTree.jetarea, "jetarea/F");
  jettree->Branch("RecoJetE", &fRecoJetTree.jetenergy, "jetenergy/F");
  jettree->Branch("RecoJetNConst", &fRecoJetTree.numberofconstituents, "numberofconstituents/I");
  jettree->Branch("RecoD0Pt", &fRecoJetTree.d0pt, "d0pt/F");
  jettree->Branch("RecoD0Eta", &fRecoJetTree.d0eta, "d0eta/F");
  jettree->Branch("RecoD0Phi", &fRecoJetTree.d0phi, "d0phi/F");
  jettree->Branch("RecoPionPt", &fRecoJetTree.pionpt, "pionpt/F");
  jettree->Branch("RecoPionEta", &fRecoJetTree.pioneta, "pioneta/F");
  jettree->Branch("RecoPionPhi", &fRecoJetTree.pionphi, "pionphi/F");
  jettree->Branch("RecoKaonPt", &fRecoJetTree.kaonpt, "kaonpt/F");
  jettree->Branch("RecoKaonEta", &fRecoJetTree.kaoneta, "kaoneta/F");
  jettree->Branch("RecoKaonPhi", &fRecoJetTree.kaonphi, "kaonphi/F");
}

void StSimJetSaver_Test::SaveJets(){

  // cout << " ****************** SAVE MC JETS TEST *******************" << endl;

  for (int jetcollection = 0; jetcollection < d0pairs; jetcollection++){

    mcJets = new TClonesArray("StJet");
    mcJets = mcJetsArr[jetcollection];

    Int_t mcnjets = mcJets->GetEntries();

    recoJets = new TClonesArray("StJet");
    recoJets = recoJetsArr[jetcollection];

    Int_t reconjets = recoJets->GetEntries();

    // cout << "Jets Per D0 = " << mcnjets << "\t" << reconjets << endl;

    // if (mcnjets == 0) continue;
    // if (reconjets == 0)continue;

    fMCJetTree = {0};
    fRecoJetTree = {0};

    StPicoMcTrack *trk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(d0MCTrackIndices[jetcollection][0]));
    StPicoMcTrack *trk2 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(d0MCTrackIndices[jetcollection][1]));

    TVector3 mTrk1Mom, mTrk2Mom, mTrkD0Mom;

    mTrk1Mom = trk1->p();
    mTrk2Mom = trk2->p();

    mTrkD0Mom = mTrk1Mom + mTrk2Mom;

    fMCJetTree.d0pt = mTrkD0Mom.Perp(); 
    fMCJetTree.d0eta = mTrkD0Mom.PseudoRapidity(); 
    fMCJetTree.d0phi = standardPhi(mTrkD0Mom.Phi());
    

    fMCJetTree.pionpt = mTrk1Mom.Perp();
    fMCJetTree.pioneta = mTrk1Mom.PseudoRapidity();
    fMCJetTree.pionphi = standardPhi(mTrk1Mom.Phi());

    fMCJetTree.kaonpt = mTrk2Mom.Perp();
    fMCJetTree.kaoneta = mTrk2Mom.PseudoRapidity();
    fMCJetTree.kaonphi = standardPhi(mTrk2Mom.Phi());

    TVector3 rTrkD0Mom, mPionMom, mKaonMom;
    mPionMom = pionReco4Momenta[jetcollection].Vect();
    mKaonMom = kaonReco4Momenta[jetcollection].Vect();
    rTrkD0Mom = d0Reco4Momenta[jetcollection].Vect();

    fRecoJetTree.d0mass = d0Reco4Momenta[jetcollection].M();
    
    fRecoJetTree.d0pt = rTrkD0Mom.Perp();
    fRecoJetTree.d0eta = rTrkD0Mom.PseudoRapidity(); 
    fRecoJetTree.d0phi = standardPhi(rTrkD0Mom.Phi()); 

    fRecoJetTree.pionpt = mPionMom.Perp();
    fRecoJetTree.pioneta = mPionMom.PseudoRapidity();
    fRecoJetTree.pionphi = standardPhi(mPionMom.Phi());

    fRecoJetTree.kaonpt = mKaonMom.Perp();
    fRecoJetTree.kaoneta = mKaonMom.PseudoRapidity();
    fRecoJetTree.kaonphi = standardPhi(mKaonMom.Phi());
 
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

      // Event level information filler
      vector<double> pv;
      pv.clear();
      pv.push_back(mPicoEvent->primaryVertex().X());
      pv.push_back(mPicoEvent->primaryVertex().Y());
      pv.push_back(mPicoEvent->primaryVertex().Z());

      fMCJetTree.primaryvertex = pv;

      // Jet level information filler

      fMCJetTree.jetpt = jetpt;
      fMCJetTree.jetcorrectedpt = corrjetpt;
      fMCJetTree.jeteta = jetEta;
      fMCJetTree.jetphi = jetPhi;
      fMCJetTree.jetarea = jetarea;
      fMCJetTree.jetenergy = jetE;
      fMCJetTree.numberofconstituents = NtrackConstit;
    }

    TRandom3 *r3 = new TRandom3(0);

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
      double NtrackConstit = rjet->GetNumberOfTracks();


      fRecoJetTree.jetpt = jetpt;
      fRecoJetTree.jetcorrectedpt = corrjetpt;
      fRecoJetTree.jeteta = jetEta;
      fRecoJetTree.jetphi = jetPhi;
      fRecoJetTree.jetarea = jetarea;
      fRecoJetTree.jetenergy = jetE;
      fRecoJetTree.numberofconstituents = NtrackConstit;

    }

    jettree->Fill();
  }
}


Double_t StSimJetSaver_Test::standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard > 2*(TMath::Pi())) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}
