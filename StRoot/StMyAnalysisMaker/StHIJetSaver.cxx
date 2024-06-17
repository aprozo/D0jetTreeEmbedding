// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StHIJetSaver.h"
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
#include "StHIOverlay_Test.h"
#include "StHIMCJets.h"
#include "StHIRecoJets.h"
#include "StHIRho.h"


// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StHIJetSaver)

//________________________________________________________________________
StHIJetSaver::StHIJetSaver(const char* name, StPicoDstMaker *picoMaker, const char* outName1 = "", const char* outName2 = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  mOutName1 = outName1;
  mOutName2 = outName2;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StHIJetSaver::fRunFlagEnum
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
  fout = 0x0;
  fPrintLevel = 0;

  fJetsArrMC.clear();
  fJetsArrReco.clear();
  fJetsArrRecoFromPYTHIAEvent.clear();
  rhovalue.clear();


  fMCJetTree = {0};
  fRecoJetTree = {0};
  fRecoJetTreeFromPYTHIAEvent = {0};

  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

}

//
//________________________________________________________________________
StHIJetSaver::~StHIJetSaver()
{ /*  */

  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StHIJetSaver::Init() {
  StJetFrameworkPicoBase::Init();

  // declare histograms

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  // fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  mcJets = new TClonesArray("StJet");
  recoJets = new TClonesArray("StJet");
  //fJets->SetName(fJetsName);

  DeclareHistograms();

  fMCJetTree = {0};
  fRecoJetTree = {0};
  fRecoJetTreeFromPYTHIAEvent = {0};

  DeclareTree();
  BookTree();

  // BookTree(mcjettree, fMCJetTree);
  // BookTree(recojettree, fRecoJetTree);

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StHIJetSaver::Finish() { 
  cout << "StHIJetSaver::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName1.Data()!="") {
    // TFile *fout = new TFile(mOutName1.Data(), "UPDATE");
    fout->cd();
    // fout->mkdir(GetName());

    // fout->mkdir(Form("%s/%s", GetName(), "MC"));
    // fout->mkdir(Form("%s/%s", GetName(), "Reco"));
    // fout->mkdir(Form("%s/%s", GetName(), "Histograms"));

    // fout->cd(Form("%s/%s", GetName(), "MC"));
    // fout->cd(GetName());
    // WriteTree(mcjettree);

    // // // fout->cd(Form("%s/%s", GetName(), "Reco"));
    // // fout->cd(GetName());
    // WriteTree(recojettree);

    // fout->cd(Form("%s/%s", GetName(), "Histograms"));
    // fout->cd(GetName());
    // WriteHistograms();

    fout->Write();
    fout->Close();
  }

  if (mOutName2.Data() != ""){
    TFile *fout2 = new TFile(mOutName2.Data(), "UPDATE");
    fout2->cd();

    hDiffJetPt->Write();
    fout2->Close();
  }

  cout<<"End of StHIJetSaver::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

void StHIJetSaver::WriteTree(TTree *sometree) {
  sometree->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StHIJetSaver::Clear(Option_t *opt) {
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StHIJetSaver::Make() {
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
  double maxpt = GetMaxTrackPt();
  if(maxpt > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


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


  mHIOverlay = static_cast<StHIOverlay_Test*>(GetMaker("HIOverlay"));
  if (!mHIOverlay){
    LOG_WARN << "No HI Overlay! Skip!" << endm;
    return kStFatal;
  }
  
  fMCTracks = mHIOverlay->GetMCTracks();
  fMCTowers = mHIOverlay->GetMCTowers();

  fRecoTracks = mHIOverlay->GetRecoTracks();
  fRecoTowers = mHIOverlay->GetRecoTowers();

  fOrigin = mHIOverlay->GetMCOrigin();
  fMCD0Information = mHIOverlay->GetMCD0();
  fRecoD0Information = mHIOverlay->GetRecoD0();

  int numberofmcevents = mHIOverlay->GetTheNumberOfEventsToOverLay();

  if (numberofmcevents == 0) return kStOk;

  // // get MC JetMaker
  // mHIMCJets = static_cast<StHIMCJets*>(GetMaker("HIMC"));

  // if(!mHIMCJets) {
  //   LOG_WARN << Form(" No %s! Skip! ", "HIMC") << endm;
  //   return kStWarn;
  // }

  fJetsArrMC = mHIOverlay->GetJets();
  fJetsArrRecoFromPYTHIAEvent = mHIOverlay->GetRecoJets();

  // cout << "Jets Array Size Again = " << fJetsArrMC.size() << "\t" << fJetsArrRecoFromPYTHIAEvent.size() << endl;
  
  // if(fJetsArrMC.size()==0) return kStOk;

  // get Reco JetMaker
  mHIRecoJets = static_cast<StHIRecoJets*>(GetMaker("HIReco"));
  if(!mHIRecoJets) {
    LOG_WARN << Form(" No %s! Skip! ", "HIReco") << endm;
    return kStWarn;
  }

  fJetsArrReco = mHIRecoJets->GetJets();

  // if(fJetsArrReco.size()==0) return kStOk;

  // ============================= RhoMaker ============================== //
  // get RhoMaker pointer from event: old names "StRho_JetsBG"
  mHIRho = static_cast<StHIRho*>(GetMaker("HIRho"));
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!mHIRho) {
    LOG_WARN << Form(" No %s! Skip! ", "HIRho") << endm;
    // return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  if(mHIRho) rhovalue = mHIRho->GetArrayofRhoValue();
  if(rhovalue.size()==0) {
    LOG_WARN << Form("Couldn't get rhovalue object! ") << endm;

    for (int i = 0; i < fJetsArrReco.size(); i++){
      rhovalue.push_back(0);
    }
    // return kStWarn;    
  }

  if (fPrintLevel) 
  cout << "Count again  = "  << numberofmcevents << "\t" << fJetsArrReco.size() << endl;
  // if (fPrintLevel) cout << "Max Track & Tower = " << GetMaxTrackPt() << "\t" << GetMaxTowerEt() << endl;
  if (fPrintLevel) cout << "Max Track & Tower = " << mHIRecoJets->MaxTrackPt() << "\t" << mHIRecoJets->MaxTowerEtBeforeHC() << "\t" << mHIRecoJets->MaxTowerEtAfterHC() << endl;

  for (int i = 0; i < numberofmcevents; i++){
    // SaveMCJets(i);
    // SaveRecoJets(i);
    // cout << "Reco in Jet Saver = " << fJetsArrReco[i]->GetEntries() << endl;
    // StJet *rjet = static_cast<StJet*>(fJetsArrReco[i]->At(0));
    // cout << "pT = " << rjet->Pt() << endl;

    SaveJets(i);
  }
  // RunJets();

  return kStOK;
}

void StHIJetSaver::DeclareTree() {


  // TString mctreename = "MCJets";
  // TString recotreename = "RecoJets";

  TString treename = "Jets";

  TDirectory *writedir;

  if(mOutName1!="") {
    fout = new TFile(mOutName1.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    writedir = (TDirectory *)fout->Get(GetName());
    // fout->Close();
  }

  jettree = new TTree(treename.Data(), treename.Data());
  jettree->SetDirectory(writedir);

  // mcjettree->SetAutoSave(1);
  // recojettree->SetAutoSave(10);
  
}

void StHIJetSaver::DeclareHistograms(){
  hDiffJetPt = new TH1F("hDiffJetPt", "hDiffJetPt", 100, -50, 50);
}

void StHIJetSaver::BookTree()
{
  // Branches to save event info
  jettree->Branch("Centrality", &fRecoJetTree.centrality, "centrality/F");
  jettree->Branch("Weight", &fRecoJetTree.weight, "weight/F");
  jettree->Branch("MCPrimaryVertex", &fMCJetTree.primaryvertex);
  jettree->Branch("RecoPrimaryVertex", &fRecoJetTree.primaryvertex);

  jettree->Branch("RecoMaxTrackPt", &fRecoJetTree.maxtrackpt, "maxtrackpt/F");
  jettree->Branch("RecoMaxTowerEtBeforeHC", &fMCJetTree.maxtowerEt, "maxtowerEt/F");
  jettree->Branch("RecoMaxTowerEtAfterHC", &fRecoJetTree.maxtowerEt, "maxtowerEt/F");

  jettree->Branch("MCJetPt", &fMCJetTree.jetpt, "jetpt/F");
  jettree->Branch("MCJetEta", &fMCJetTree.jeteta, "jeteta/F");
  jettree->Branch("MCJetPhi", &fMCJetTree.jetphi, "jetphi/F");
  jettree->Branch("MCJetArea", &fMCJetTree.jetarea, "jetarea/F");
  jettree->Branch("MCJetE", &fMCJetTree.jetenergy, "jetenergy/F");
  jettree->Branch("MCJetNConst", &fMCJetTree.numberofconstituents, "numberofconstituents/I");
  // jettree->Branch("MCD0Z", &fMCJetTree.d0z, "d0z/F");
  jettree->Branch("MCD0Pt", &fMCJetTree.d0pt, "d0pt/F");
  jettree->Branch("MCD0Eta", &fMCJetTree.d0eta, "d0eta/F");
  jettree->Branch("MCD0Phi", &fMCJetTree.d0phi, "d0phi/F");
  jettree->Branch("MCPionPt", &fMCJetTree.pionpt, "pionpt/F");
  jettree->Branch("MCPionEta", &fMCJetTree.pioneta, "pioneta/F");
  jettree->Branch("MCPionPhi", &fMCJetTree.pionphi, "pionphi/F");
  jettree->Branch("MCKaonPt", &fMCJetTree.kaonpt, "kaonpt/F");
  jettree->Branch("MCKaonEta", &fMCJetTree.kaoneta, "kaoneta/F");
  jettree->Branch("MCKaonPhi", &fMCJetTree.kaonphi, "kaonphi/F");

  jettree->Branch("RecoJetPtFromPYTHIA", &fRecoJetTreeFromPYTHIAEvent.jetpt, "jetpt/F");
  jettree->Branch("RecoJetEtaFromPYTHIA", &fRecoJetTreeFromPYTHIAEvent.jeteta, "jeteta/F");
  jettree->Branch("RecoJetPhiFromPYTHIA", &fRecoJetTreeFromPYTHIAEvent.jetphi, "jetphi/F");
  jettree->Branch("RecoJetEFromPYTHIA", &fRecoJetTreeFromPYTHIAEvent.jetenergy, "jetenergy/F");
  jettree->Branch("RecoJetAreaFromPYTHIA", &fRecoJetTreeFromPYTHIAEvent.jetarea, "jetarea/F");
  jettree->Branch("RecoJetNConstFromPYTHIA", &fRecoJetTreeFromPYTHIAEvent.numberofconstituents, "numberofconstituents/I");

  jettree->Branch("RecoJetPt", &fRecoJetTree.jetpt, "jetpt/F");
  jettree->Branch("RecoJetCorrPt", &fRecoJetTree.jetcorrectedpt, "jetcorrectedpt/F");
  jettree->Branch("RecoJetEta", &fRecoJetTree.jeteta, "jeteta/F");
  jettree->Branch("RecoJetPhi", &fRecoJetTree.jetphi, "jetphi/F");
  jettree->Branch("RecoJetArea", &fRecoJetTree.jetarea, "jetarea/F");
  jettree->Branch("RecoJetE", &fRecoJetTree.jetenergy, "jetenergy/F");
  jettree->Branch("RecoJetRhoVal", &fRecoJetTree.fRhoValforjet, "fRhoValforjet/F");
  jettree->Branch("RecoJetNConst", &fRecoJetTree.numberofconstituents, "numberofconstituents/I");
  // jettree->Branch("RecoD0Z", &fRecoJetTree.d0z, "d0z/F");
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

void StHIJetSaver::SaveJets(int iteration){
  // cout << " ****************** SAVE MC JETS TEST *******************" << endl;

  mcJets = new TClonesArray("StJet");
  mcJets = fJetsArrMC[iteration];

  fMCJetTree = {0}; // Structure to save out jet information
  fRecoJetTree = {0}; // Structure to save out jet information
  fRecoJetTreeFromPYTHIAEvent = {0}; // Structure to save out jet information from reco PYTHIA event

  recoJets = new TClonesArray("StJet");
  recoJets = fJetsArrReco[iteration];
  fRhoVal  = rhovalue[iteration];

  recoJetsFromPYTHIAEvent = new TClonesArray("StJet");
  recoJetsFromPYTHIAEvent = fJetsArrRecoFromPYTHIAEvent[iteration];
  
  Int_t mcnjets = mcJets->GetEntries();
  Int_t reconjets = recoJets->GetEntries();
  Int_t reconjetsfrompythia = recoJetsFromPYTHIAEvent->GetEntries();

  // cout << "Jets = " << "\t" << iteration << "\t" << mcnjets << "\t" << reconjets << "\t" << reconjetsfrompythia << endl;

  // if (mcnjets == 0) return;

  fRecoJetTree.centrality = fCentralityScaled;
  fRecoJetTree.weight = mCentMaker->GetWeight();

  TVector3 mcPion = fMCD0Information[iteration].first;
  TVector3 mcKaon = fMCD0Information[iteration].second;

  TVector3 mcD0 = mcPion+mcKaon;

  fMCJetTree.d0pt = mcD0.Perp();
  fMCJetTree.d0eta = mcD0.PseudoRapidity();
  fMCJetTree.d0phi = standardPhi(mcD0.Phi());

  if (fPrintLevel) cout << "MC D0 Saved = " << fMCJetTree.d0pt << "\t" << fMCJetTree.d0eta << "\t" << fMCJetTree.d0phi << endl;
  
  fMCJetTree.pionpt = mcPion.Perp();
  fMCJetTree.pioneta = mcPion.PseudoRapidity();
  fMCJetTree.pionphi = standardPhi(mcPion.Phi());

  fMCJetTree.kaonpt = mcKaon.Perp();
  fMCJetTree.kaoneta = mcKaon.PseudoRapidity();
  fMCJetTree.kaonphi = standardPhi(mcKaon.Phi());

  TVector3 recoPion = fRecoD0Information[iteration].first;
  TVector3 recoKaon = fRecoD0Information[iteration].second;

  TVector3 recoD0 = recoPion+recoKaon;

  fRecoJetTree.d0pt = recoD0.Perp();
  fRecoJetTree.d0eta = recoD0.PseudoRapidity();
  fRecoJetTree.d0phi = standardPhi(recoD0.Phi());

  // cout << "Reco D0 Saved = " << fRecoJetTree.d0pt << "\t" << fRecoJetTree.d0eta << "\t" << fRecoJetTree.d0phi << endl;
  
  fRecoJetTree.pionpt = recoPion.Perp();
  fRecoJetTree.pioneta = recoPion.PseudoRapidity();
  fRecoJetTree.pionphi = standardPhi(recoPion.Phi());

  fRecoJetTree.kaonpt = recoKaon.Perp();
  fRecoJetTree.kaoneta = recoKaon.PseudoRapidity();
  fRecoJetTree.kaonphi = standardPhi(recoKaon.Phi());

  // fRecoJetTree.maxtrackpt = GetMaxTrackPt();
  // fRecoJetTree.maxtowerEt = GetMaxTowerEt();
  
  fRecoJetTree.maxtrackpt = mHIRecoJets->MaxTrackPt();
  fMCJetTree.maxtowerEt = mHIRecoJets->MaxTowerEtBeforeHC(); // Before Hadronic Correction ET
  fRecoJetTree.maxtowerEt = mHIRecoJets->MaxTowerEtAfterHC(); // Afrer Hadronic Correction ET

  int mcjetid = -99;
  int recojetid = -99;
  int recojetidfrompythia = -99;

  for(int jjet = 0; jjet < mcnjets; jjet++) {
    StJet *mjet = static_cast<StJet*>(mcJets->At(jjet));
    if(!mjet) continue;

    bool d0Jet = kFALSE;

    for (int jtrk = 0; jtrk < mjet->GetNumberOfTracks(); jtrk++){
      int trackid = mjet->TrackAt(jtrk);

      double mass = fMCTracks[iteration][trackid].M();

      if (int(mass)==1) {d0Jet = kTRUE; break;}
    }
    
    if (!d0Jet) continue;

    mcjetid = jjet;
  }

  for(int jjet = 0; jjet < reconjets; jjet++) {
    StJet *rjet = static_cast<StJet*>(recoJets->At(jjet));
    if(!rjet) continue;

    bool d0Jet = kFALSE;

    for (int jtrk = 0; jtrk < rjet->GetNumberOfTracks(); jtrk++){
      int trackid = rjet->TrackAt(jtrk);

      if (trackid >=10000) {
        double mass = fRecoTracks[iteration][trackid-10000].M();

        // if (trackid >= 10000) cout << "Track ID = " << trackid << "\t" << mass << endl;

        if (int(mass)==1) {d0Jet = kTRUE; break;}
      }
    }
    
    if (!d0Jet) continue;

    recojetid = jjet;
  }

  for(int kjet = 0; kjet < reconjetsfrompythia; kjet++) {
    StJet *rpjet = static_cast<StJet*>(recoJetsFromPYTHIAEvent->At(kjet));
    if(!rpjet) continue;

    bool d0Jet = kFALSE;

    for (int jtrk = 0; jtrk < rpjet->GetNumberOfTracks(); jtrk++){
      int trackid = rpjet->TrackAt(jtrk);

      double mass = fRecoTracks[iteration][trackid].M();

      // if (trackid >= 10000) cout << "Track ID = " << trackid << "\t" << mass << endl;

      if (int(mass)==1) {d0Jet = kTRUE; break;}
    }
    
    if (!d0Jet) continue;

    recojetidfrompythia = kjet;
  }

  if (mcjetid != -99){
    StJet *mjet = static_cast<StJet*>(mcJets->At(mcjetid));
    if(!mjet) return;

    vector<double> mcpv;
    mcpv.clear();
    mcpv.push_back(fOrigin[iteration].X());
    mcpv.push_back(fOrigin[iteration].Y());
    mcpv.push_back(fOrigin[iteration].Z());

    if (fPrintLevel) cout << "HIJetSaver Vertex " << mcpv[0] << "\t" << mcpv[1] << "\t" << mcpv[2] << endl;

    fMCJetTree.primaryvertex = mcpv;

    fMCJetTree.jetpt = mjet->Pt();
    fMCJetTree.jeteta = mjet->Eta();
    fMCJetTree.jetphi = mjet->Phi();
    fMCJetTree.jetarea = mjet->Area();
    fMCJetTree.jetenergy = mjet->E();
    fMCJetTree.numberofconstituents = mjet->GetNumberOfTracks();

    if (fPrintLevel) cout << "MC Jet Saved = " << fMCJetTree.jetpt << "\t" << fMCJetTree.jeteta << "\t" << fMCJetTree.jetphi << endl;

    TVector3 mcjet;
    mcjet.SetPtEtaPhi(mjet->Pt(), mjet->Eta(), mjet->Phi());

    double mcz = (mcjet.X()*mcD0.X() + mcjet.Y()*mcD0.Y())/(pow(mjet->Pt(), 2));

    fMCJetTree.d0z = mcz;
  }

  // cout << "Jet info finally saved = " << iteration << "\t" << mjet->Pt() << "\t" << mcD0.Perp() << endl;

  if (recojetid != -99){
    StJet *rjet = static_cast<StJet*>(recoJets->At(recojetid));
    if(!rjet) return;

    vector<double> recopv;
    recopv.clear();
    recopv.push_back(mPicoEvent->primaryVertex().X());
    recopv.push_back(mPicoEvent->primaryVertex().Y());
    recopv.push_back(mPicoEvent->primaryVertex().Z());

    fRecoJetTree.primaryvertex = recopv;

    fRecoJetTree.jetpt = rjet->Pt();
    // cout << "Jet pT = " << rjet->Pt() << endl;
    fRecoJetTree.jetcorrectedpt = rjet->Pt() - fRhoVal*rjet->Area();
    fRecoJetTree.jeteta = rjet->Eta();
    fRecoJetTree.jetphi = rjet->Phi();
    fRecoJetTree.jetarea = rjet->Area();
    fRecoJetTree.jetenergy = rjet->E();
    fRecoJetTree.fRhoValforjet = fRhoVal;
    fRecoJetTree.numberofconstituents = rjet->GetNumberOfTracks();

    // cout << "Reco Jet Saved = " << fRecoJetTree.jetpt << "\t" << fRecoJetTree.jeteta << "\t" << fRecoJetTree.jetphi << endl;

    double RecoJetPx = fRecoJetTree.jetcorrectedpt *TMath::Cos(rjet->Phi());
    double RecoJetPy = fRecoJetTree.jetcorrectedpt *TMath::Sin(rjet->Phi());

    double recoz = (RecoJetPx*recoD0.X() + RecoJetPy*recoD0.Y())/(pow(fRecoJetTree.jetcorrectedpt, 2));

    fRecoJetTree.d0z = recoz;
  }

  // cout << "recojetidfrompythia = " << recojetidfrompythia << endl;

  if (recojetidfrompythia != -99){
    StJet *rjet = static_cast<StJet*>(recoJetsFromPYTHIAEvent->At(recojetidfrompythia));
    if(!rjet) return;

    fRecoJetTreeFromPYTHIAEvent.jetpt = rjet->Pt();
    fRecoJetTreeFromPYTHIAEvent.jeteta = rjet->Eta();
    fRecoJetTreeFromPYTHIAEvent.jetphi = rjet->Phi();
    fRecoJetTreeFromPYTHIAEvent.jetarea = rjet->Area();
    fRecoJetTreeFromPYTHIAEvent.jetenergy = rjet->E();
    fRecoJetTreeFromPYTHIAEvent.numberofconstituents = rjet->GetNumberOfTracks();
  }

  hDiffJetPt->Fill(fRecoJetTree.jetcorrectedpt - fMCJetTree.jetpt);

  jettree->Fill();

}

Double_t StHIJetSaver::standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard > 2*(TMath::Pi())) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}

