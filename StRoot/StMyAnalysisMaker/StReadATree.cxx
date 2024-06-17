// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StReadATree.h"
#include "StMemStat.h"
#include "phys_constants.h"

// ROOT includes

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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
#include "StRoot/StPicoEvent/StPicoDstReader.h"
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
#include "StRoot/StPicoEvent/StPicoMcTrack.h"
#include "StRoot/StPicoEvent/StPicoMcVertex.h"

#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
// #include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StBTofUtil/tofPathLength.hh"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// extra includes
#include "StJetPicoDefinitions.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StReadATree)

//________________________________________________________________________
StReadATree::StReadATree(const char *name = "", const char* filename = "", const char* outName = "") : StJetFrameworkPicoBase(name)
{
  mFilename = filename;
  mOutName = outName;

  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  Bfield = 0;
  doUsePrimTracks = kTRUE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StTagD0Events::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  doRejectBadRuns = kFALSE;
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

  fMinJetTrackPt = 0.2;
  fMaxJetTrackPt = 30.0; 
  fMinJetClusPt = 0.15;
  fMinJetClusE = 0.2;
  fMinJetTowerE = 0.2;
  fTrackEtaMin = -1.0; 
  fTrackEtaMax = 1.0;
  fTrackPhiMin = 0.0; fTrackPhiMax = 2.0*TMath::Pi();
  fJetTrackEtaMin = -1.0; fJetTrackEtaMax = 1.0;
  fJetTrackPhiMin = 0.0; fJetTrackPhiMax = 2.0*TMath::Pi();
  fJetTrackDCAcut= 3.0;
  fJetTracknHitsFit = 15;
  fJetTracknHitsRatio = 0.52;
  fTrackEfficiency = 1.;
  fJetTowerEMin = 0.2; fJetTowerEMax = 100.0;
  fJetTowerEtaMin = -1.0; fJetTowerEtaMax = 1.0;
  fJetTowerPhiMin = 0.0; fJetTowerPhiMax = 2.0*TMath::Pi();
  mTowerEnergyMin = 0.2;
  mHadronicCorrFrac = 1.;


  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  mEmcPosition = 0x0;

  mBemcGeom = 0x0;
  // mCentMaker = 0x0;
  mBaseMaker = 0x0;

  fAnalysisMakerName = name;

  fEventVectorToFill.clear();
  fTrackVectorToFill.clear();
  fTowerVectorToFill.clear();

}

//
//________________________________________________________________________
StReadATree::~StReadATree()
{ 

}

void StReadATree::DeclareHistograms() {

    hdEta = new TH1D("hdEta", "hdEta", 100, -1, 1);
    hdPhi = new TH1D("hdPhi", "hdPhi", 200, -2*TMath::Pi(), 2*TMath::Pi());
    hdEtadPhi = new TH2D("hdEtadPhi", "hdEtadPhi", 100, -1, 1, 200, -2*TMath::Pi(), 2*TMath::Pi());
    hdTowerIndex = new TH1D("hdTowerIndex", "hdTowerIndex", 9601, -4800.5, 4800.5);
    hPercentTracks = new TH1D("hPercentTracks", "hPercentTracks", 27, 74.5, 100.5);
}

void StReadATree::WriteHistograms() {

    hdEta->Write();
    hdPhi->Write();
    hdEtadPhi->Write();
    hdTowerIndex->Write();
    hPercentTracks->Write();
}
//
//________________________________________________________________________
Int_t StReadATree::Init() {

  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  fEventVectorToFill.clear();
  fTrackVectorToFill.clear();
  fTowerVectorToFill.clear();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  mBemcGeom = StEmcGeom::instance("bemc");

  // if (!AssignDataTree(fVariablesAssignedToTree)) {return kStWarn;}

  cout << "File Name is = " << mFilename << endl;

  if (mFilename == "NULL") return kStOK;

  f = TFile::Open(mFilename);

  if (f->IsOpen()) cout << "File is Open!" << endl;
  fVariablesAssignedToTree = (TTree *)f->Get("PicoDst");

  cout << fVariablesAssignedToTree->GetEntriesFast() << endl;

  fVariablesAssignedToTree->SetMakeClass(1);

  fVariablesAssignedToTree->SetBranchStatus("*",0);

  fVariablesAssignedToTree->SetBranchStatus("Event", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mRunId", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mEventId", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mPrimaryVertexX", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mPrimaryVertexY", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mPrimaryVertexZ", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mPrimaryVertexErrorX", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mPrimaryVertexErrorY", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mPrimaryVertexErrorZ", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mNumberOfGlobalTracks", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mRefMultNeg", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mRefMultPos", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mGRefMult", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mZDCx", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mBBCx", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mId", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mChi2", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mPMomentumX", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mPMomentumY", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mPMomentumZ", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mGMomentumX", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mGMomentumY", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mGMomentumZ", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mOriginX", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mOriginY", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mOriginZ", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mDedx", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mDedxError", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mNHitsFit", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mNHitsMax", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mNHitsDedx", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mNSigmaPion", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mNSigmaKaon", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mNSigmaProton", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mNSigmaElectron", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mTopologyMap[2]", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mBEmcPidTraitsIndex", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mBTofPidTraitsIndex", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mMtdPidTraitsIndex", 1);
  fVariablesAssignedToTree->SetBranchStatus("BTowHit", 1);
  fVariablesAssignedToTree->SetBranchStatus("BTowHit.mAdc", 1);
  fVariablesAssignedToTree->SetBranchStatus("BTowHit.mE", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mTrackIndex", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBemcId", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBemcAdc0", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBemcE0", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBemcE", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBemcZDist", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBemcPhiDist", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBemcSmdNEta", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBemcSmdNPhi", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBtowId", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBtowId23", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBtowE", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBtowE2", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBtowE3", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBtowEtaDist", 1);
  fVariablesAssignedToTree->SetBranchStatus("EmcPidTraits.mBtowPhiDist", 1);

  fVariablesAssignedToTree->SetBranchAddress("Event", &Event_);
  fVariablesAssignedToTree->SetBranchAddress("Event.mRunId", Event_mRunId);
  fVariablesAssignedToTree->SetBranchAddress("Event.mEventId", Event_mEventId);
  fVariablesAssignedToTree->SetBranchAddress("Event.mPrimaryVertexX", Event_mPrimaryVertexX);
  fVariablesAssignedToTree->SetBranchAddress("Event.mPrimaryVertexY", Event_mPrimaryVertexY);
  fVariablesAssignedToTree->SetBranchAddress("Event.mPrimaryVertexZ", Event_mPrimaryVertexZ);
  fVariablesAssignedToTree->SetBranchAddress("Event.mPrimaryVertexErrorX", Event_mPrimaryVertexErrorX);
  fVariablesAssignedToTree->SetBranchAddress("Event.mPrimaryVertexErrorY", Event_mPrimaryVertexErrorY);
  fVariablesAssignedToTree->SetBranchAddress("Event.mPrimaryVertexErrorZ", Event_mPrimaryVertexErrorZ);
  fVariablesAssignedToTree->SetBranchAddress("Event.mNumberOfGlobalTracks", Event_mNumberOfGlobalTracks); 
  fVariablesAssignedToTree->SetBranchAddress("Event.mRefMultNeg", Event_mRefMultNeg);
  fVariablesAssignedToTree->SetBranchAddress("Event.mRefMultPos", Event_mRefMultPos);
  fVariablesAssignedToTree->SetBranchAddress("Event.mGRefMult", Event_mGRefMult); 
  fVariablesAssignedToTree->SetBranchAddress("Event.mZDCx", Event_mZDCx);
  fVariablesAssignedToTree->SetBranchAddress("Event.mBBCx", Event_mBBCx);
  fVariablesAssignedToTree->SetBranchAddress("Track", &Track_);
  fVariablesAssignedToTree->SetBranchAddress("Track.mId", Track_mId);
  fVariablesAssignedToTree->SetBranchAddress("Track.mChi2", Track_mChi2);
  fVariablesAssignedToTree->SetBranchAddress("Track.mPMomentumX", Track_mPMomentumX);
  fVariablesAssignedToTree->SetBranchAddress("Track.mPMomentumY", Track_mPMomentumY);
  fVariablesAssignedToTree->SetBranchAddress("Track.mPMomentumZ", Track_mPMomentumZ);
  fVariablesAssignedToTree->SetBranchAddress("Track.mGMomentumX", Track_mGMomentumX);
  fVariablesAssignedToTree->SetBranchAddress("Track.mGMomentumY", Track_mGMomentumY);
  fVariablesAssignedToTree->SetBranchAddress("Track.mGMomentumZ", Track_mGMomentumZ);
  fVariablesAssignedToTree->SetBranchAddress("Track.mOriginX", Track_mOriginX);
  fVariablesAssignedToTree->SetBranchAddress("Track.mOriginY", Track_mOriginY);
  fVariablesAssignedToTree->SetBranchAddress("Track.mOriginZ", Track_mOriginZ);
  fVariablesAssignedToTree->SetBranchAddress("Track.mDedx", Track_mDedx);
  fVariablesAssignedToTree->SetBranchAddress("Track.mDedxError", Track_mDedxError);
  fVariablesAssignedToTree->SetBranchAddress("Track.mNHitsFit", Track_mNHitsFit);
  fVariablesAssignedToTree->SetBranchAddress("Track.mNHitsMax", Track_mNHitsMax);
  fVariablesAssignedToTree->SetBranchAddress("Track.mNHitsDedx", Track_mNHitsDedx);
  fVariablesAssignedToTree->SetBranchAddress("Track.mNSigmaPion", Track_mNSigmaPion);
  fVariablesAssignedToTree->SetBranchAddress("Track.mNSigmaKaon", Track_mNSigmaKaon);
  fVariablesAssignedToTree->SetBranchAddress("Track.mNSigmaProton", Track_mNSigmaProton);
  fVariablesAssignedToTree->SetBranchAddress("Track.mNSigmaElectron", Track_mNSigmaElectron);
  fVariablesAssignedToTree->SetBranchAddress("Track.mTopologyMap[2]", Track_mTopologyMap);
  fVariablesAssignedToTree->SetBranchAddress("Track.mBEmcPidTraitsIndex", Track_mBEmcPidTraitsIndex);
  fVariablesAssignedToTree->SetBranchAddress("Track.mBTofPidTraitsIndex", Track_mBTofPidTraitsIndex);
  fVariablesAssignedToTree->SetBranchAddress("Track.mMtdPidTraitsIndex", Track_mMtdPidTraitsIndex);
  fVariablesAssignedToTree->SetBranchAddress("BTowHit", &BTowHit_);
  fVariablesAssignedToTree->SetBranchAddress("BTowHit.mAdc", BTowHit_mAdc);
  fVariablesAssignedToTree->SetBranchAddress("BTowHit.mE", BTowHit_mE);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits", &EmcPidTraits_);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mTrackIndex", EmcPidTraits_mTrackIndex);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBemcId", EmcPidTraits_mBemcId);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBemcAdc0", EmcPidTraits_mBemcAdc0);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBemcE0", EmcPidTraits_mBemcE0);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBemcE", EmcPidTraits_mBemcE);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBemcZDist", EmcPidTraits_mBemcZDist);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBemcPhiDist", EmcPidTraits_mBemcPhiDist);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBemcSmdNEta", EmcPidTraits_mBemcSmdNEta);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBemcSmdNPhi", EmcPidTraits_mBemcSmdNPhi);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBtowId", EmcPidTraits_mBtowId);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBtowId23", EmcPidTraits_mBtowId23);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBtowE", EmcPidTraits_mBtowE);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBtowE2", EmcPidTraits_mBtowE2);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBtowE3", EmcPidTraits_mBtowE3);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBtowEtaDist", EmcPidTraits_mBtowEtaDist);
  fVariablesAssignedToTree->SetBranchAddress("EmcPidTraits.mBtowPhiDist", EmcPidTraits_mBtowPhiDist);

  for(int i=0; i<4800; i++) {
    for(int j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

  TString DataFileName(mFilename); 
  
  TString mainfilename = ((TObjString *)DataFileName.Tokenize("/")->Last())->String();
  // mainfilename.ReplaceAll("_3_0_cc_MONASH_TXFile_NoD0Decay_0", "");
  cout << mainfilename.Data() << endl;
  //Here we get the vertex file we need to get the runids

  vertexfilename = mainfilename;
  vertexfilename.ReplaceAll(".picoDst.root", "");
  vertexfilename = TString("/gpfs01/star/pwg/droy1/EMBEDDING_2021_D0Analysis/FullPYTHIA/VERTEXFULL/Vertex_") + vertexfilename + ".txt";
  cout << vertexfilename << endl;

  return kStOK;
}

Int_t StReadATree::Make() {

  if (mFilename == "NULL") return kStOK;
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

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();
  mVertex = mPicoEvent->primaryVertex();

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  GetTOFMatching();
  // FillTheVectorThatReplacesThePicoDst();

  return kStOK;
}


void StReadATree::GetTOFMatching(){
  for(unsigned short iTracks = 0; iTracks < mPicoDst->numberOfTracks(); iTracks++){
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

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

    if (charge == 0) continue;

    if (pt < 0.3) continue; //This is different
    
    if (abs(eta) > 1) continue;

    if (trk->nHitsDedx() < 20) continue; //This is different

    if (float(trk->nHitsDedx())/float(trk->nHitsMax()) < 0.52) { // This is not there in the analysis note. I am a little confused.
      // continue;
    }

    if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 3.) { // This is not there in the analysis note. I am a little confused.
      // continue;
    } 

    if (trk->chi2() > 3) { // This is not there in the analysis note. I am a little confused.
      // continue;
    }

    // HFT Hits only

    if (!trk->isHFTTrack()) continue;
    if (!trk->hasPxl1Hit() && !trk->hasPxl2Hit()) continue;
    if (!trk->hasIstHit()) continue;

    int index2tof = trk->bTofPidTraitsIndex();
    
    float beta = -99.;

    if (index2tof < 0) continue;

    StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(index2tof);
        
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      
      if (beta < 1e-4)
      {
          TVector3 const btofHitPosvec = tofPid->btofHitPos();

          StThreeVectorF const btofHitPos(btofHitPosvec.x(), btofHitPosvec.y(), btofHitPosvec.z()); 
          StThreeVectorF const pVertex(mVertex.x(), mVertex.y(), mVertex.z()); 
          
          StPicoPhysicalHelix helix = trk->helix(Bfield);
          float L = tofPathLength(&pVertex, &btofHitPos, helix.curvature());
          float tof = tofPid->btof();
          if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
          else beta = std::numeric_limits<float>::quiet_NaN();
      }

      // cout << tofPid->btofBeta() - beta << endl;
      if (abs(tofPid->btofBeta() - beta) > 1e-4) cout << "Compare Betas : " << tofPid->btofBeta() << "\t" << beta << endl;
    }
  }
} 

// float StReadATree::getTofBeta(StPicoTrack const* const trk, StPicoDst const* const picoDst, StThreeVectorF const& pVertex) const
// {
//     int index2tof = trk->bTofPidTraitsIndex();
    
//     float beta = std::numeric_limits<float>::quiet_NaN();
    
//     if (index2tof >= 0)
//     {
//         StPicoBTofPidTraits *tofPid = picoDst->btofPidTraits(index2tof);
        
//         if (tofPid)
//         {
//             beta = tofPid->btofBeta();
            
//             if (beta < 1e-4)
//             {
//                 StThreeVectorF const btofHitPos = tofPid->btofHitPos();
                
//                 StPhysicalHelixD helix = trk->helix();
//                 float L = tofPathLength(&pVertex, &btofHitPos, helix.curvature());
//                 float tof = tofPid->btof();
//                 if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
//                 else beta = std::numeric_limits<float>::quiet_NaN();
//             }
//         }
//     }
    
//     return beta;
// }

void StReadATree::FillTheVectorThatReplacesThePicoDst(){


  fEventVectorToFill.clear();
  fTrackVectorToFill.clear();
  fTowerVectorToFill.clear();

  cout << "Number of BEMC Pid Traits = " << mPicoDst->numberOfBEmcPidTraits() << endl;

  int numberOfacceptedtracks = 0;
  int matchedtracks = 0;


  for(unsigned short iTracks = 0; iTracks < mPicoDst->numberOfTracks(); iTracks++){
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

    // if (trk->isPrimary()) continue;

    if (trk->nHitsDedx() < 15) continue;

    if (float(trk->nHitsDedx())/float(trk->nHitsMax()) < 0.52) continue;

    if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 3.) continue;


    double picoeta = trk->gMom(mVertex, Bfield).PseudoRapidity();
    double picopt = trk->gMom(mVertex, Bfield).Perp();

    if (abs(picoeta) > 1) continue;
    if (picopt < 2.0) continue;

    numberOfacceptedtracks++;

    // cout << "Tower ID = " << trk->bemcTowerIndex() <<  "\t" << GetMatchedBtowID(trk) << endl;

    if (trk->bemcTowerIndex() == abs(GetMatchedBtowID(trk)) - 1)matchedtracks++;
    // cout << "New Tower ID = " << GetMatchedBtowID(trk) << endl;

    // double towerE = 0.0;
    // if (trk->bemcTowerIndex() >= 0){
    //   StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(trk->bemcTowerIndex()));
    //   if(!tower) { cout<<"No tower pointer... iTow = "<<trk->bemcTowerIndex()<<endl; continue; }        // corrected energy (hadronically - done below)
    //   towerE = tower->energy();
    // }

    // double emcE = 0.;
    // double emcE2 = 0.;
    // double emcE3 = 0.;
    // double clusmE = 0.;

    // int ID = -99;

    // double towerEnergy = 0;

    // if (trk->bemcPidTraitsIndex() >= 0) {
    //   StPicoBEmcPidTraits * Emc =  mPicoDstMaker->picoDst()->bemcPidTraits(trk->bemcPidTraitsIndex());
    //   emcE = Emc->btowE();
    //   emcE2 = Emc->btowE2();
    //   emcE3 = Emc->btowE3();
    //   clusmE = Emc->bemcE();
    //   ID = Emc->btowId();

    //   StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(trk->bemcTowerIndex()));
    //   if(!tower) { cout<<"No tower pointer... iTow = "<<trk->bemcTowerIndex()<<endl; continue; }        // corrected energy (hadronically - done below)
    //   towerEnergy = tower->energy();
    // }

    // cout << towerE << "\t" << emcE << "\t" << emcE2 << "\t" << emcE3 << "\t" << clusmE << "\t" << towerEnergy << "\t" << ID << endl;
  }

  cout << numberOfacceptedtracks << "\t" << matchedtracks << "\t" << 1.0*matchedtracks/numberOfacceptedtracks*100.0 << endl;

  hPercentTracks->Fill(1.0*matchedtracks/numberOfacceptedtracks*100.0);

  vector<double> vertextocheck = ReadNthLine(vertexfilename.Data(), int(Event_mEventId[0] % 1000)-1);

  // fVariablesAssignedToTree->LoadTree(0);
  // fVariablesAssignedToTree->GetEntryWithIndex(vertextocheck[0],  vertextocheck[1]);

  fVariablesAssignedToTree->GetEntryWithIndex( mPicoEvent->runId(),  mPicoEvent->eventId());
  cout << Event_mRunId[0] << "\t" << Event_mEventId[0] << endl;

  fEventVectorToFill.push_back({Event_mRunId[0], Event_mBBCx[0], Event_mZDCx[0], Event_mGRefMult[0], Event_mRefMultNeg[0]+Event_mRefMultPos[0]});

  TVector3 EventVertex;
  EventVertex.SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

  int ntracks = Event_mNumberOfGlobalTracks[0];

  double pi0mass = Pico::mMass[0]; // GeV

  for(int i=0; i<4800; i++) {
    for(int j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

  int tracks = 0;

  

  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){

    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }


  
    TVector3 TrackVertex;
    TrackVertex.SetXYZ(Track_mOriginX[iTracks], Track_mOriginY[iTracks], Track_mOriginZ[iTracks]);

    TVector3 mTrkMom;
    if(doUsePrimTracks) { 
      // get primary track vector
      mTrkMom.SetXYZ(Track_mPMomentumX[iTracks], Track_mPMomentumY[iTracks], Track_mPMomentumZ[iTracks]); 
    } else { 
      // get global track vector
      mTrkMom.SetXYZ(Track_mGMomentumX[iTracks], Track_mGMomentumY[iTracks], Track_mGMomentumZ[iTracks]);
    }

    double dca = (EventVertex - TrackVertex).Mag();
    int nHitsFit = Track_mNHitsFit[iTracks];
    int nHitsMax = Track_mNHitsMax[iTracks];
    double nHitsRatio = 1.0*nHitsFit/nHitsMax;

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
    if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
    double eta = mTrkMom.PseudoRapidity();

    ////////////////  Accepting the jet track ////////////////

    if(pt < fMinJetTrackPt) continue;
    if(pt > fMaxJetTrackPt) continue; // 20.0 STAR, 100.0 ALICE
    if((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax)) continue;
    if((phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax)) continue;
        
    // additional quality cuts for tracks
    if(dca > fJetTrackDCAcut)            continue;
    if(nHitsFit < fJetTracknHitsFit)     continue;
    if(nHitsRatio < fJetTracknHitsRatio) continue;

    ////////////////  Accepting the jet track ////////////////

    double px = mTrkMom.x();
    double py = mTrkMom.y();
    double pz = mTrkMom.z();
    double p = mTrkMom.Mag();
    double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
    short charge = (Track_mNHitsFit[iTracks] > 0) ? 1 : -1;
    int matchedTowerIndex = Track_mBEmcPidTraitsIndex[iTracks]; // towerIndex = towerID - 1

    double picophi = trk->pMom().Phi();
    if(picophi < 0.0)    picophi += 2.0*pi;  // force from 0-2pi
    if(picophi > 2.0*pi) picophi -= 2.0*pi;  // force from 0-2pi

    if (matchedTowerIndex < 0) {continue;}

    int towerID = EmcPidTraits_mBtowId[matchedTowerIndex];
    double newenergy = BTowHit_mE[towerID-1]/1000.;
    double newenergyfromemc = EmcPidTraits_mBtowE[matchedTowerIndex]/1000.;

    cout << "Energies = " << towerID << "\t" << newenergy << "\t" << newenergyfromemc << "\t" << mPicoDst->bemcPidTraits(trk->bemcPidTraitsIndex())->btowE() << endl;
    cout << "Old Tower ID = " << trk->bemcTowerIndex() << endl;

    if  (trk->isBemcMatchedTrack()){
      // cout << "From Tree    = " << pt << "\t" << phi << "\t" << eta << "\t" << matchedTowerIndex << "\t" << towerID << "\t" << newenergy << endl;
      // cout << "From PicoDst = " << trk->pMom().Perp() << "\t" << picophi << "\t" << trk->pMom().PseudoRapidity() << "\t" << mPicoDst->bemcPidTraits(trk->bemcPidTraitsIndex())->btowE() << endl;


    }

    // cout << "Matched Tower Index = " << matchedTowerIndex << endl;
    // if (matchedTowerIndex >= 0) cout << "Matched Tower Energy " << BTowHit_mE[matchedTowerIndex]/ 1000. << endl;

    //====  matched track index ===
    int trackIndex = iTracks;
    if(trackIndex < 0) { continue; } // can't happen

    mTowerMatchTrkIndex[ matchedTowerIndex ][ mTowerStatusArr[ matchedTowerIndex ]] = trackIndex;
    mTowerStatusArr[ matchedTowerIndex ] = mTowerStatusArr[ matchedTowerIndex ] + 1;  // 1+ means match, 0 for no match

    fTrackVectorToFill.push_back({px, py, pz, energy, tracks + 10000});
    tracks++;
  }

  

  int nTowers = 4800;

  int tower = 0;

  for(int itow = 0; itow < nTowers; itow++) {
    int towerID = itow + 1;
    int towerIndex = towerID - 1;
    if(towerID < 0) continue; // double check these aren't still in the event list

    TVector3 towerPosition = mEmcPosition->getPosFromVertex(EventVertex, towerID);
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

    int towerADC = BTowHit_mAdc[itow];
    double towerEunCorr = BTowHit_mE[itow]/ 1000.;  // uncorrected energy
    double towerE = BTowHit_mE[itow]/ 1000.;        // corrected energy (hadronically - done below)
    double towEtunCorr = towerEunCorr / (1.0*TMath::CosH(towerEta));

    if(towerEunCorr < mTowerEnergyMin) continue;

    // =======================================================================
    // HADRONIC CORRECTION
    double maxEt = 0.;
    double sumEt = 0.;

    // if tower was is matched to a track or multiple, add up the matched track energies - (mult opt.) to then subtract from the corresponding tower
    // August 15: if *have* 1+ matched trk-tow AND uncorrected energy of tower is at least your tower constituent cut, then CONTINUE 
    if(mTowerStatusArr[towerIndex] > 0.5 && towerEunCorr > mTowerEnergyMin) {
      double maxE = 0.0;
      double sumE = 0.0;
      // =======================================================================================================================
      // --- finds max E track matched to tower *AND* the sum of all matched track E and subtract from said tower
      //     USER provides readMacro.C which method to use for their analysis via SetJetHadCorrType(type);
      // loop over ALL matched tracks
      for(int iTracks = 0; iTracks < mTowerStatusArr[towerIndex]; iTracks++) {
        TVector3 mTrkMom;
        if(doUsePrimTracks) { 
          // get primary track vector
          mTrkMom.SetXYZ(Track_mPMomentumX[iTracks], Track_mPMomentumY[iTracks], Track_mPMomentumZ[iTracks]); 
        } else { 
          // get global track vector
          mTrkMom.SetXYZ(Track_mGMomentumX[iTracks], Track_mGMomentumY[iTracks], Track_mGMomentumZ[iTracks]);
        }

        // track variables
        double p = mTrkMom.Mag();
        double E = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
        if(E > maxE) maxE = E;
        sumE = sumE + E;

      } // track loop

      // apply hadronic correction to tower
      maxEt  = (towerEunCorr - (mHadronicCorrFrac * maxE))/(1.0*TMath::CosH(towerEta));
      sumEt  = (towerEunCorr - (mHadronicCorrFrac * sumE))/(1.0*TMath::CosH(towerEta));        

      //=================================================================================================================
    }  // have a track-tower match
    // else - no match so treat towers on their own. Must meet constituent cut

    // Et - correction comparison
    double fMaxEt = (maxEt == 0) ? towerEunCorr / (1.0*TMath::CosH(towerEta)) : maxEt;
    double fSumEt = (sumEt == 0) ? towerEunCorr / (1.0*TMath::CosH(towerEta)) : sumEt;

    // cut on transverse tower energy (more uniform)
    double towerEt = 0.0;
    if(mTowerStatusArr[towerIndex] < 1) { // no matches, use towers uncorrected energy
      towerEt = towerEunCorr / (1.0*TMath::CosH(towerEta)); 
    } else { 
      towerEt = fSumEt;  towerE = fSumEt  * 1.0*TMath::CosH(towerEta);
    }
    if(towerEt == 0) { cout<<"fJetHadCorrType - or - towerE actually 0"<<endl; }  // it was unset, because you provided wrong fJetHadCorrType
    if(towerEt < 0) towerEt = 0.0;
    if(towerEt < mTowerEnergyMin) continue;

    // cout << "Tower E = " << towerID << "\t" << towerE << "\t" << towEtunCorr << endl;

    Double_t energy = BTowHit_mE[itow]/ 1000.;
    Double_t p = TMath::Sqrt(energy*energy - pi0mass*pi0mass);

    double posX = towerPosition.x();
    double posY = towerPosition.y();
    double posZ = towerPosition.z();

    Double_t r = TMath::Sqrt(posX*posX + posY*posY + posZ*posZ);

    TVector3 mTowerMom;
    if(r > 1e-12) {
      mTowerMom.SetX( p*posX/r );
      mTowerMom.SetY( p*posY/r );
      mTowerMom.SetZ( p*posZ/r );
      // energy) ;
    }
    else continue;

    fTowerVectorToFill.push_back({mTowerMom.x(), mTowerMom.y(), mTowerMom.z(), energy, -tower - 10000, towerADC, towerEt, towerE});

    tower++;
  }

}

//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StReadATree::Finish() { 
  
  if (mFilename == "NULL") return kStOK;
  f->Close();

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

  return kStOK;
}

void StReadATree::Clear(Option_t *opt){

}

int StReadATree::GetMatchedBtowID(StPicoTrack *trk){
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

  int trkbemcid = trk->bemcTowerIndex();

  // Check if the track can be projected onto the current radius
  // if not, track can't be matched.
  // By JetCorr request the global track projection to BEMC is used.
  if ( mEmcPosition->projTrack(&bemc_pos, &bemc_mom, trk, mVertex, Bfield, bemc_radius) ) {
    // First, examine track eta. If it falls in two regions:
    // 0 < |eta| < etaMin()
    // etaMax() < |eta| < 1.0
    // then shift the eta for the projection slightly into the neighboring tower
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, trkbemcid + 1);
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

          if (trkbemcid != abs(tow_id) - 1) {
            TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, trkbemcid + 1);

            hdEta->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity());
            hdPhi->Fill(bemc_pos.Phi() - towerPosition.Phi());
            hdEtadPhi->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity(), bemc_pos.Phi() - towerPosition.Phi());
            hdTowerIndex->Fill(trkbemcid - abs(tow_id) + 1);

            cout << "Mismatched Tower = " << trkbemcid << "\t" << abs(tow_id) - 1 << "\t" << bemc_pos.Phi() << "\t" << bemc_pos.PseudoRapidity() << "\t" << towerPosition.Phi() << "\t" << towerPosition.PseudoRapidity() << endl;
          }
          
          return -1*tow_id;
        }

        else{
          if (trkbemcid != abs(tow_id) - 1) {
            TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, trkbemcid + 1);

            hdEta->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity());
            hdPhi->Fill(bemc_pos.Phi() - towerPosition.Phi());
            hdEtadPhi->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity(), bemc_pos.Phi() - towerPosition.Phi());
            hdTowerIndex->Fill(trkbemcid - abs(tow_id) + 1);

            cout << "Mismatched Tower = " << trkbemcid << "\t" << abs(tow_id) - 1 << "\t" << bemc_pos.Phi() << "\t" << bemc_pos.PseudoRapidity() << "\t" << towerPosition.Phi() << "\t" << towerPosition.PseudoRapidity() << endl;          }

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
          if (trkbemcid != abs(tow_id) - 1) {
            TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, trkbemcid + 1);

            hdEta->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity());
            hdPhi->Fill(bemc_pos.Phi() - towerPosition.Phi());
            hdEtadPhi->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity(), bemc_pos.Phi() - towerPosition.Phi());
            hdTowerIndex->Fill(trkbemcid - abs(tow_id) + 1);

            cout << "Mismatched Tower = " << trkbemcid << "\t" << abs(tow_id) - 1 << "\t" << bemc_pos.Phi() << "\t" << bemc_pos.PseudoRapidity() << "\t" << towerPosition.Phi() << "\t" << towerPosition.PseudoRapidity() << endl;          }

          return -1*tow_id;
        }

        else if( mBemcGeom->getBin(bemc_pos_shift_neg.Phi(),bemc_pos_shift_neg.PseudoRapidity(),h_m,h_e,h_s) == 0 && h_s != -1 ){
          mBemcGeom->getId(h_m,h_e,h_s,tow_id);
          if (trkbemcid != abs(tow_id) - 1) {
            TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, trkbemcid + 1);

            hdEta->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity());
            hdPhi->Fill(bemc_pos.Phi() - towerPosition.Phi());
            hdEtadPhi->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity(), bemc_pos.Phi() - towerPosition.Phi());
            hdTowerIndex->Fill(trkbemcid - abs(tow_id) + 1);

            cout << "Mismatched Tower = " << trkbemcid << "\t" << abs(tow_id) - 1 << "\t" << bemc_pos.Phi() << "\t" << bemc_pos.PseudoRapidity() << "\t" << towerPosition.Phi() << "\t" << towerPosition.PseudoRapidity() << endl;          }

          return -1*tow_id;
        } 
      }
    }
  }

  if (trkbemcid != abs(tow_id) - 1) {
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, trkbemcid + 1);

    hdEta->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity());
    hdPhi->Fill(bemc_pos.Phi() - towerPosition.Phi());
    hdEtadPhi->Fill(bemc_pos.PseudoRapidity() - towerPosition.PseudoRapidity(), bemc_pos.Phi() - towerPosition.Phi());
    hdTowerIndex->Fill(trkbemcid - abs(tow_id) + 1);

    cout << "Mismatched Tower = " << trkbemcid << "\t" << abs(tow_id) - 1 << "\t" << bemc_pos.Phi() << "\t" << bemc_pos.PseudoRapidity() << "\t" << towerPosition.Phi() << "\t" << towerPosition.PseudoRapidity() << endl;
  }

  return tow_id;

}


