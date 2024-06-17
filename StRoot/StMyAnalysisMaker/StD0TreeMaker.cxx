// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StD0TreeMaker.h"
#include "StMemStat.h"


// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <vector>

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
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"
#include "StTagD0Events.h"
#include "StD0EventsJetMaker.h"
#include "StD0Rho.h"


// Bichsel includes
#include "StBichsel/Bichsel.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StD0TreeMaker)

//________________________________________________________________________
StD0TreeMaker::StD0TreeMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  fD0Kind = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StD0TreeMaker::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;

  fMinCentrality = 0;
  fMaxCentrality = 80;

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
  frefCorr2 = -99.;
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
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  d0TrackIndices.clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};
}

//
//________________________________________________________________________
StD0TreeMaker::~StD0TreeMaker()
{ /*  */
  
  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StD0TreeMaker::Init() {
  StJetFrameworkPicoBase::Init();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  //fJets->SetOwner(kTRUE);

  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};

  DeclareTree();
  BookTree();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StD0TreeMaker::Finish() { 
  cout << "StD0TreeMaker::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    // WriteHistograms();

    WriteTree();

    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StD0TreeMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}


void StD0TreeMaker::DeclareTree() {
  TString treename;

  if (fD0Kind == 0){
    treename = "SignalPlusBackground";
  }

  else if (fD0Kind == 1){
    treename = "UnlikeBackground";
  }

  else if (fD0Kind == 2){
    treename = "LikeBackground";
  }

  jettree = new TTree(treename.Data(), treename.Data());
}

void StD0TreeMaker::WriteTree() {
  jettree->Write();
}

//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StD0TreeMaker::Clear(Option_t *opt) {
  fJets->Clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};

}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StD0TreeMaker::Make() {
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

  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
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
  frefCorr2 = refCorr2;
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

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

  mD0Tagger = static_cast<StTagD0Events*>(GetMaker("D0Tagger"));
  if(!mD0Tagger) {
    LOG_WARN << " No D0Tagger! Skip! " << endm;
    return kStWarn;
  }

  if ( fD0Kind==-1  && !mD0Tagger->DoesEventHaveD0()) return kStOK;
  if ( fD0Kind== 1  && !mD0Tagger->DoesEventHaveD0Bg()) return kStOK;
  if ( fD0Kind==-2  && !mD0Tagger->DoesEventHaveD0()) return kStOK;
  if ( fD0Kind== 2  && !mD0Tagger->DoesEventHaveD0Bg()) return kStOK;
  if ( fD0Kind==-3  && !mD0Tagger->DoesEventHaveD0()) return kStOK;
  if ( fD0Kind== 3  && !mD0Tagger->DoesEventHaveD0Bg()) return kStOK;

  if ( fD0Kind==-1 ) d0TrackIndices = mD0Tagger->GetD0Indices();
  if ( fD0Kind== 1 ) d0TrackIndices = mD0Tagger->GetD0BgIndices();
  if ( fD0Kind==-2 ) d0TrackIndices = mD0Tagger->GetTightD0Indices();
  if ( fD0Kind== 2 ) d0TrackIndices = mD0Tagger->GetTightD0BgIndices();
  if ( fD0Kind==-3 ) d0TrackIndices = mD0Tagger->GetLooseD0Indices();
  if ( fD0Kind== 3 ) d0TrackIndices = mD0Tagger->GetLooseD0BgIndices();

  if (d0TrackIndices.size() == 0) return kStOK;

  RunJets();

  // // run Tracks:
  // RunTracks();

  // // run Towers:
  // RunTowers();

  return kStOK;
}

//
//
//_____________________________________________________________________________________________

void StD0TreeMaker::BookTree()
{
  // Branches to save event info
  jettree->Branch("RunID", &fJetTree.runid, "runid/I");
  jettree->Branch("EventId", &fJetTree.eventid, "eventid/I");
  jettree->Branch("RefMult", &fJetTree.refmult, "refmult/F");
  jettree->Branch("Centrality", &fJetTree.centrality, "centrality/F");
  jettree->Branch("Triggers", &fJetTree.triggers);
  jettree->Branch("PrimaryVertex", &fJetTree.primaryvertex);
  jettree->Branch("PrimaryVertexErr", &fJetTree.primaryvertexerror);

  jettree->Branch("D0Mass", &fJetTree.d0mass, "d0mass/F");
  jettree->Branch("D0Pt", &fJetTree.d0pt, "d0pt/F");
  jettree->Branch("D0Eta", &fJetTree.d0eta, "d0eta/F");
  jettree->Branch("D0Phi", &fJetTree.d0phi, "d0phi/F");

  jettree->Branch("PionPt", &fJetTree.pionpt, "pionpt/F");
  jettree->Branch("PionEta", &fJetTree.pioneta, "pioneta/F");
  jettree->Branch("PionPhi", &fJetTree.pionphi, "pionphi/F");
  jettree->Branch("PionCharge", &fJetTree.pioncharge, "pioncharge/F");

  jettree->Branch("KaonPt", &fJetTree.kaonpt, "kaonpt/F");
  jettree->Branch("KaonEta", &fJetTree.kaoneta, "kaoneta/F");
  jettree->Branch("KaonPhi", &fJetTree.kaonphi, "kaonphi/F");
  jettree->Branch("KaonCharge", &fJetTree.kaoncharge, "kaoncharge/F");

}


void StD0TreeMaker::RunJets()
{
  
  for (int jetcollection = 0; jetcollection < d0TrackIndices.size(); jetcollection++){

    fJetTree = {0};

    fJetTree.runid = mPicoEvent->runId();
    fJetTree.eventid = mPicoEvent->eventId();
    fJetTree.refmult = frefCorr2;
    fJetTree.centrality = fCentralityScaled;
    fJetTree.triggers = mPicoEvent->triggerIds();

    vector<double> pv;
    pv.clear();
    pv.push_back(mPicoEvent->primaryVertex().X());
    pv.push_back(mPicoEvent->primaryVertex().Y());
    pv.push_back(mPicoEvent->primaryVertex().Z());

    fJetTree.primaryvertex = pv;

    vector<double> pverr;
    pverr.clear();
    pverr.push_back(mPicoEvent->primaryVertexError().X());
    pverr.push_back(mPicoEvent->primaryVertexError().Y());
    pverr.push_back(mPicoEvent->primaryVertexError().Z());

    fJetTree.primaryvertexerror = pverr;


    iTrack1 = d0TrackIndices[jetcollection][0];
    iTrack2 = d0TrackIndices[jetcollection][1];

    StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(iTrack1));
    if(!trk1){ continue; }

    StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(iTrack2));
    if(!trk2){ continue; }

    TVector3 mTrk1Mom, mTrk2Mom;
    if(doUsePrimTracks) {
      if(!(trk1->isPrimary())) continue; // check if primary
      if(!(trk2->isPrimary())) continue;
      // get primary track vector
      mTrk1Mom = trk1->pMom();
      mTrk1Mom = trk2->pMom();
      } else {
      // get global track vector
      mTrk1Mom = trk1->gMom(mVertex, Bfield);
      mTrk2Mom = trk2->gMom(mVertex, Bfield);
    }

    fJetTree.pionpt = mTrk1Mom.Perp();
    fJetTree.pioneta = mTrk1Mom.PseudoRapidity();
    fJetTree.pionphi = standardPhi(mTrk1Mom.Phi());
    fJetTree.pioncharge = trk1->charge();

    fJetTree.kaonpt = mTrk2Mom.Perp();
    fJetTree.kaoneta = mTrk2Mom.PseudoRapidity();
    fJetTree.kaonphi = standardPhi(mTrk2Mom.Phi());
    fJetTree.kaoncharge = trk2->charge();
    
    TVector3 mTrkD0Mom;
    mTrkD0Mom = mTrk1Mom + mTrk2Mom;


    fJetTree.d0mass = d0TrackIndices[jetcollection][2];
    fJetTree.d0pt = mTrkD0Mom.Perp();
    fJetTree.d0phi = mTrkD0Mom.Phi();
    fJetTree.d0eta = mTrkD0Mom.PseudoRapidity();

    jettree->Fill();

  }

}


Double_t StD0TreeMaker::standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard > 2*(TMath::Pi())) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}


