// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StJetTreeMaker.h"
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

ClassImp(StJetTreeMaker)

//________________________________________________________________________
StJetTreeMaker::StJetTreeMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  fRunFlag = 0;       // see StJetTreeMaker::fRunFlagEnum
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
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  d0TrackIndices.clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};
}

//
//________________________________________________________________________
StJetTreeMaker::~StJetTreeMaker()
{ /*  */
  
  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StJetTreeMaker::Init() {
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
Int_t StJetTreeMaker::Finish() { 
  cout << "StJetTreeMaker::Finish()\n";

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

  cout<<"End of StJetTreeMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}


void StJetTreeMaker::DeclareTree() {
  TString treename = "D0Jets";
  jettree = new TTree(treename.Data(), treename.Data());
}

void StJetTreeMaker::WriteTree() {
  jettree->Write();
}

//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StJetTreeMaker::Clear(Option_t *opt) {
  fJets->Clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};

}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StJetTreeMaker::Make() {
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
  // if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  // if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


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

  // cut on unset centrality, > 80%
  // if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // // cut on centrality for analysis before doing anything
  // if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  // FillEmcTriggers();

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
  // ======================== end of Triggers ============================= //

  // =========================== JetMaker =============================== //
  // get JetMaker pointer
  // JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  // const char *fJetMakerNameCh = fJetMakerName;
  // if(!JetMaker) {
  //   LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
  //   return kStWarn;
  // }

  // // get jet collection associated with JetMaker
  // fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  // if(!fJets) {
  //   LOG_WARN << Form(" No fJets object! Skip! ") << endm;
  //   return kStWarn;
  // }

  // get JetMaker


  mD0Tagger = static_cast<StTagD0Events*>(GetMaker("D0Tagger"));
  if(!mD0Tagger) {
    LOG_WARN << " No D0Tagger! Skip! " << endm;
    return kStWarn;
  }

  if ( fD0Kind==-1  && !mD0Tagger->DoesEventHaveD0()) return kStOK;
  if ( fD0Kind== 1  && !mD0Tagger->DoesEventHaveD0Bg()) return kStOK;
  if ( fD0Kind==-2  && !mD0Tagger->DoesEventHaveTightD0()) return kStOK;
  if ( fD0Kind== 2  && !mD0Tagger->DoesEventHaveTightD0Bg()) return kStOK;
  if ( fD0Kind==-3  && !mD0Tagger->DoesEventHaveLooseD0()) return kStOK;
  if ( fD0Kind== 3  && !mD0Tagger->DoesEventHaveLooseD0Bg()) return kStOK;

  if ( fD0Kind==-1 ) d0TrackIndices = mD0Tagger->GetD0Indices();
  if ( fD0Kind== 1 ) d0TrackIndices = mD0Tagger->GetD0BgIndices();
  if ( fD0Kind==-2 ) d0TrackIndices = mD0Tagger->GetTightD0Indices();
  if ( fD0Kind== 2 ) d0TrackIndices = mD0Tagger->GetTightD0BgIndices();
  if ( fD0Kind==-3 ) d0TrackIndices = mD0Tagger->GetLooseD0Indices();
  if ( fD0Kind== 3 ) d0TrackIndices = mD0Tagger->GetLooseD0BgIndices();

  if (d0TrackIndices.size() == 0) return kStOK;


  D0JetMaker = static_cast<StD0EventsJetMaker*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!D0JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // if we have JetMaker, get jet collection associated with it
  if(D0JetMaker) {
    fJetsArr =  D0JetMaker->GetJets();
    //fJets->SetName("BGJetsRho");  // name is set by Maker who created it
  }
  if(fJetsArr.size()==0) return kStOk;

  // ============================= RhoMaker ============================== //
  // get RhoMaker pointer from event: old names "StRho_JetsBG"
  D0RhoMaker = static_cast<StD0Rho*>(GetMaker(fRhoMakerName));
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!D0RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  rhovalue = D0RhoMaker->GetArrayofRhoValue();
  if(rhovalue.size()==0) {
    LOG_WARN << Form("Couldn't get rhovalue object! ") << endm;

    for (int i = 0; i < fJetsArr.size(); i++){
      rhovalue.push_back(0);
    }

    // return kStWarn;    
  } 
  
  RunJets();

  return kStOK;
}

//
//
//_____________________________________________________________________________________________

void StJetTreeMaker::BookTree()
{
  // Branches to save event info
  jettree->Branch("RunID", &fJetTree.runid, "runid/I");
  jettree->Branch("EventId", &fJetTree.eventid, "eventid/I");
  jettree->Branch("RefMult", &fJetTree.refmult, "refmult/F");
  jettree->Branch("gRefMult", &fJetTree.grefmult, "grefmult/F");
  jettree->Branch("Centrality", &fJetTree.centrality, "centrality/F");
  jettree->Branch("Weight", &fJetTree.weight, "weight/F");
  jettree->Branch("Triggers", &fJetTree.triggers);
  jettree->Branch("PrimaryVertex", &fJetTree.primaryvertex);
  jettree->Branch("PrimaryVertexErr", &fJetTree.primaryvertexerror);

  jettree->Branch("JetPt", &fJetTree.jetpt, "jetpt/F");
  jettree->Branch("JetCorrPt", &fJetTree.jetcorrectedpt, "jetcorrectedpt/F");
  jettree->Branch("JetEta", &fJetTree.jeteta, "jeteta/F");
  jettree->Branch("JetPhi", &fJetTree.jetphi, "jetphi/F");
  jettree->Branch("JetArea", &fJetTree.jetarea, "jetarea/F");
  jettree->Branch("JetRadius", &fJetTree.jetradius, "jetradius/F");
  jettree->Branch("JetE", &fJetTree.jetenergy, "jetenergy/F");
  jettree->Branch("JetNEF", &fJetTree.jetnef, "jetnef/F");
  jettree->Branch("JetRhoVal", &fJetTree.fRhoValforjet, "fRhoValforjet/F");
  jettree->Branch("JetHighestTrackPt", &fJetTree.jethighesttrackpt, "jethighesttrackpt/F");
  jettree->Branch("JetNConst", &fJetTree.numberofconstituents, "numberofconstituents/I");

  jettree->Branch("D0Mass", &fJetTree.d0mass, "d0mass/F");

  jettree->Branch("PionPt", &fJetTree.pionpt, "pionpt/F");
  jettree->Branch("PionEta", &fJetTree.pioneta, "pioneta/F");
  jettree->Branch("PionPhi", &fJetTree.pionphi, "pionphi/F");
  jettree->Branch("PionCharge", &fJetTree.pioncharge, "pioncharge/F");

  jettree->Branch("KaonPt", &fJetTree.kaonpt, "kaonpt/F");
  jettree->Branch("KaonEta", &fJetTree.kaoneta, "kaoneta/F");
  jettree->Branch("KaonPhi", &fJetTree.kaonphi, "kaonphi/F");
  jettree->Branch("KaonCharge", &fJetTree.kaoncharge, "kaoncharge/F");

  // Branches to save jet data info
  jettree->Branch("TrackID", &fJetTree.mTrackID, "mTrackID[numberofconstituents]/F");
  jettree->Branch("TrackPt", &fJetTree.mTrackPt, "mTrackPt[numberofconstituents]/F");
  jettree->Branch("TrackEta", &fJetTree.mTrackEta, "mTrackEta[numberofconstituents]/F");
  jettree->Branch("TrackPhi", &fJetTree.mTrackPhi, "mTrackPhi[numberofconstituents]/F");
  jettree->Branch("TrackPx", &fJetTree.mTrackPx, "mTrackPx[numberofconstituents]/F");
  jettree->Branch("TrackPy", &fJetTree.mTrackPy, "mTrackPy[numberofconstituents]/F");
  jettree->Branch("TrackPz", &fJetTree.mTrackPz, "mTrackPz[numberofconstituents]/F");
  jettree->Branch("TrackCharge", &fJetTree.mTrackCharge, "mTrackCharge[numberofconstituents]/F");
}


void StJetTreeMaker::RunJets()
{
  
  // ====================== Jet loop below ============================
  // loop over Jets in the event: initialize some parameter variables


  for (int jetcollection = 0; jetcollection < fJetsArr.size(); jetcollection++){

    // Jet TClonesArray
    fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it

    fJets = fJetsArr[jetcollection];
    fRhoVal = rhovalue[jetcollection];

    Int_t njets = fJets->GetEntries();

    bool d0Jet = kFALSE;

    // loop over jets
    for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
      // get jet pointer
      StJet *jet = static_cast<StJet*>(fJets->At(ijet));
      if(!jet) continue;

      d0Jet = kFALSE;

      // get some jet parameters
      double jetarea = jet->Area();
      double jetpt = jet->Pt();
      double corrjetpt = jet->Pt() - jetarea*fRhoVal;
      double jetE = jet->E();
      double jetEta = jet->Eta();
      double jetPhi = jet->Phi();
      double jetNEF = jet->NEF();

      // get nTracks and maxTrackPt
      double maxtrackpt = jet->GetMaxTrackPt();
      double NtrackConstit = jet->GetNumberOfTracks();

      for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
        int trackid = jet->TrackAt(itrk);      

        if (trackid == 99999) {
          // cout << "D0 Jet Found." << endl;
          d0Jet = kTRUE;
          break;
        }
      }

      if (d0Jet){

        fJetTree = {0};

        // Event level information filler

        fJetTree.runid = mPicoEvent->runId();
        fJetTree.eventid = mPicoEvent->eventId();
        fJetTree.refmult = mCentMaker->GetrefMult();
        fJetTree.grefmult = mCentMaker->GetgrefMult();
        fJetTree.centrality = mCentMaker->GetCent9();
        fJetTree.weight = mCentMaker->GetWeight();
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

        // Jet level information filler

        fJetTree.jetpt = jetpt;
        fJetTree.jetcorrectedpt = corrjetpt;
        fJetTree.jeteta = jetEta;
        fJetTree.jetphi = jetPhi;
        fJetTree.jetarea = jetarea;
        fJetTree.jetradius = 0.4;
        fJetTree.jetenergy = jetE;
        fJetTree.jetnef = jetNEF;
        fJetTree.fRhoValforjet = fRhoVal;
        fJetTree.jethighesttrackpt = maxtrackpt;
        fJetTree.numberofconstituents = NtrackConstit;

        fJetTree.d0mass = d0TrackIndices[jetcollection][2];

        for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
          int trackid = jet->TrackAt(itrk);

          // fJetTree.mTrackID[itrk] = trackid;

          if (trackid != 99999){
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

            // cout << "=============================================" << endl;

            // cout << "Before filling" << endl;

            // cout << trackid << "\t" << pt << "\t" << phi << "\t" << eta << endl;

            // Track level information fill

            fJetTree.mTrackID[itrk] = 0;
            fJetTree.mTrackPt[itrk] = pt;
            fJetTree.mTrackEta[itrk] = eta;
            fJetTree.mTrackPhi[itrk] = standardPhi(phi);
            fJetTree.mTrackPx[itrk] = px;
            fJetTree.mTrackPy[itrk] = py;
            fJetTree.mTrackPz[itrk] = pz;
            fJetTree.mTrackCharge[itrk] = charge;

          }

          else{
            iTrack1 = d0TrackIndices[jetcollection][0];
            iTrack2 = d0TrackIndices[jetcollection][1];

            StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(iTrack1));
            if(!trk1){ continue; }

            StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(iTrack2));
            if(!trk2){ continue; }

            TVector3 mTrk1Mom, mTrk2Mom;
            mTrk1Mom = trk1->gMom(mVertex, Bfield);
            mTrk2Mom = trk2->gMom(mVertex, Bfield);

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

            double pt = mTrkD0Mom.Perp();
            double phi = mTrkD0Mom.Phi();
            double eta = mTrkD0Mom.PseudoRapidity();
            double px = mTrkD0Mom.x();
            double py = mTrkD0Mom.y();
            double pz = mTrkD0Mom.z();
            short charge = 0;

            // cout << "=============================================" << endl;

            // cout << "Before filling" << endl;

            // cout << trackid << "\t" << pt << "\t" << phi << "\t" << eta << endl;

            fJetTree.mTrackID[itrk] = 421;
            fJetTree.mTrackPt[itrk] = pt;
            fJetTree.mTrackEta[itrk] = eta;
            fJetTree.mTrackPhi[itrk] = standardPhi(phi);
            fJetTree.mTrackPx[itrk] = px;
            fJetTree.mTrackPy[itrk] = py;
            fJetTree.mTrackPz[itrk] = pz;
            fJetTree.mTrackCharge[itrk] = charge;

          }

        }

        jettree->Fill();

      }

    }

  }

}

Double_t StJetTreeMaker::standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard > 2*(TMath::Pi())) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}


