// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// track and tower input to cluster together using fastjet wrapper and create jets
//      - leading jet tag
//      - switches for defining jet
//      - access to jet constituents
//      - general QA
//      - constitutent subtractor
//
// ################################################################
// $Id$

#include "StHIRecoJets.h"

// ROOT includes
#include "TROOT.h"
#include <TChain.h>
#include <TClonesArray.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include "TFile.h"
#include "TVector3.h"
#include <sstream>
#include <fstream>

// StRoot includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoArrays.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoConstants.h"

// for towers
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcPosition2.h"
class StEmcPosition2;

// jet class and fastjet wrapper and dataset (Run#'s)
#include "StJet.h"
#include "StFJWrapper.h"
#include "SW_UserIndex.h"
#include "StJetFrameworkPicoBase.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistRun14AuAu_P18ih.h"                 // new Run14 AuAu
#include "runlistRun14AuAu_P16id_SL18f_xrootd_MB.h" // Run14 AuAu used by HF group for MB

// centrality includes
#include "StCentMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// D0 Includes
#include "StHIOverlay_Test.h"

// extra includes
#include "StJetPicoDefinitions.h"

// StRoot classes
class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;

// constants - definitions
const Int_t StHIRecoJets::fgkConstIndexShift = 100000;

ClassImp(StHIRecoJets)

    //________________________________________________________________________
    StHIRecoJets::StHIRecoJets() : StMaker(),
                                   fD0Analysis(kFALSE),
                                   doWriteHistos(kFALSE),
                                   doUsePrimTracks(kFALSE),
                                   fDebugLevel(0),
                                   fRunFlag(0), // see StJetFrameworkPicoBase::fRunFlagEnum
                                   doppAnalysis(kFALSE),
                                   fRequireCentSelection(kFALSE),
                                   fMCEventsWithoutCent(kFALSE),
                                   doConstituentSubtr(kFALSE),
                                   fDoEffCorr(kFALSE),
                                   doCorrectTracksforEffBeforeJetReco(kFALSE),
                                   fTrackEfficiencyType(StJetFrameworkPicoBase::kNormalPtEtaBased), // this is default, should not need to change
                                   fEventZVtxMinCut(-40.0),
                                   fEventZVtxMaxCut(40.0),
                                   fCentralitySelectionCut(-99),
                                   fMaxEventTrackPt(30.0),
                                   fMaxEventTowerEt(1000.0), // 30.0
                                   doRejectBadRuns(kFALSE),
                                   Bfield(0.0),
                                   //  mVertex(0x0),
                                   zVtx(0.0),
                                   fRunNumber(0),
                                   fD0Kind(0),
                                   fEmcTriggerEventType(0), // see StJetFrameworkPicoBase::fEmcTriggerFlagEnum
                                   fMBEventType(2),         // kVPDMB5, see StJetFrameworkPicoBase::fMBFlagEnum
                                   fTriggerToUse(0),        // kTriggerAny, see StJetFrameworkPicoBase::fTriggerEventTypeEnum
                                   fCentralityScaled(0.0),
                                   ref16(-99),
                                   ref9(-99),
                                   mOutName(""),
                                   fTracksName(""),
                                   fCaloName(""),
                                   fJetsName(""),
                                   fJetAlgo(1),
                                   fJetType(0),
                                   fRecombScheme(fastjet::BIpt2_scheme), // was BIpt_scheme
                                   fjw("StHIRecoJets", "StHIRecoJets"),
                                   fRadius(0.4),
                                   fMinJetArea(0.001),
                                   fMinJetPt(1.0),
                                   fJetPhiMin(-10.), fJetPhiMax(+10.),
                                   fJetEtaMin(-0.6), fJetEtaMax(0.6),
                                   fGhostArea(0.005),
                                   fMinJetTrackPt(0.2), fMaxJetTrackPt(30.0),
                                   fMinJetClusPt(0.15),
                                   fMinJetClusE(0.2),
                                   fMinJetTowerE(0.2),
                                   fTrackEtaMin(-1.0), fTrackEtaMax(1.0),
                                   fTrackPhiMin(0.0), fTrackPhiMax(2.0 * TMath::Pi()),
                                   fJetTrackEtaMin(-1.0), fJetTrackEtaMax(1.0),
                                   fJetTrackPhiMin(0.0), fJetTrackPhiMax(2.0 * TMath::Pi()),
                                   fJetTrackDCAcut(3.0),
                                   fJetTracknHitsFit(15),
                                   fJetTracknHitsRatio(0.52),
                                   fTrackEfficiency(1.),
                                   fJetTowerEMin(0.2), fJetTowerEMax(100.0),
                                   fJetTowerEtaMin(-1.0), fJetTowerEtaMax(1.0),
                                   fJetTowerPhiMin(0.0), fJetTowerPhiMax(2.0 * TMath::Pi()),
                                   mTowerEnergyMin(0.2),
                                   mHadronicCorrFrac(1.),
                                   fJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks), // default is using all matched Tracks, Aug2019, per Hanseul
                                   fLegacyMode(kFALSE),
                                   fFillGhost(kFALSE),
                                   fJets(0x0),
                                   fJetsBGsub(0x0),
                                   fJetsArr(0x0),
                                   fJetsBGsubArr(0x0),
                                   fFull_Event(0),
                                   fConstituents(0),
                                   mGeom(StEmcGeom::instance("bemc")),
                                   mPicoDstMaker(0x0),
                                   mPicoDst(0x0),
                                   mPicoEvent(0x0),
                                   mHIOverlay(0x0),
                                   mCentMaker(0x0),
                                   mBaseMaker(0x0),
                                   mBemcGeom(0x0),
                                   mEmcPosition(0x0),
                                   grefmultCorr(0x0),
                                   fEfficiencyInputFile(0x0),
                                   fBackground(kFALSE),
                                   fPrintLevel(0)
{
  // Default constructor.
  for (int i = 0; i < 8; i++)
  {
    fEmcTriggerArr[i] = kFALSE;
  }

  for (int i = 0; i < 4800; i++)
  {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;

    for (int j = 0; j < 7; j++)
      mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }
  d0TrackIndices.clear();
  fJetsArr.clear();
  fJetsBGsubArr.clear();
  rhovalue.clear();
  sigmavalue.clear();
}

//________________________________________________________________________
StHIRecoJets::StHIRecoJets(const char *name, double mintrackPt = 0.20, bool doHistos = kFALSE, const char *outName = "") : StMaker(name),
                                                                                                                           fD0Analysis(kFALSE),
                                                                                                                           doWriteHistos(doHistos),
                                                                                                                           doUsePrimTracks(kFALSE),
                                                                                                                           fDebugLevel(0),
                                                                                                                           fRunFlag(0), // see StJetFrameworkPicoBase::fRunFlagEnum
                                                                                                                           doppAnalysis(kFALSE),
                                                                                                                           fRequireCentSelection(kFALSE),
                                                                                                                           fMCEventsWithoutCent(kFALSE),
                                                                                                                           doConstituentSubtr(kFALSE),
                                                                                                                           fDoEffCorr(kFALSE),
                                                                                                                           doCorrectTracksforEffBeforeJetReco(kFALSE),
                                                                                                                           fTrackEfficiencyType(StJetFrameworkPicoBase::kNormalPtEtaBased), // this is default, should not need to change
                                                                                                                           fEventZVtxMinCut(-40.0),
                                                                                                                           fEventZVtxMaxCut(40.0),
                                                                                                                           fCentralitySelectionCut(-99),
                                                                                                                           fMaxEventTrackPt(30.0),
                                                                                                                           fMaxEventTowerEt(1000.0), // 30.0
                                                                                                                           doRejectBadRuns(kFALSE),
                                                                                                                           Bfield(0.0),
                                                                                                                           //  mVertex(0x0),
                                                                                                                           zVtx(0.0),
                                                                                                                           fRunNumber(0),
                                                                                                                           fD0Kind(0),
                                                                                                                           fEmcTriggerEventType(0), // see StJetFrameworkPicoBase::fEMCTriggerFlagEnum
                                                                                                                           fMBEventType(2),         // kVPDMB5, see StJetFrameworkPicoBase::fMBFlagEnum
                                                                                                                           fTriggerToUse(0),        // kTriggerAny, see StJetFrameworkPicoBase::fTriggerEventTypeEnum
                                                                                                                           fCentralityScaled(0.0),
                                                                                                                           ref16(-99),
                                                                                                                           ref9(-99),
                                                                                                                           mOutName(outName),
                                                                                                                           fTracksName("Tracks"),
                                                                                                                           fCaloName("Towers"),
                                                                                                                           fJetsName("Jets"),
                                                                                                                           fJetAlgo(1),
                                                                                                                           fJetType(0),
                                                                                                                           fRecombScheme(fastjet::BIpt2_scheme), // was BIpt2_scheme
                                                                                                                           fjw(name, name),
                                                                                                                           fRadius(0.4),
                                                                                                                           fMinJetArea(0.001),
                                                                                                                           fMinJetPt(1.0),
                                                                                                                           fJetPhiMin(-10), fJetPhiMax(+10),
                                                                                                                           fJetEtaMin(-0.6), fJetEtaMax(0.6),
                                                                                                                           fGhostArea(0.005),
                                                                                                                           fMinJetTrackPt(mintrackPt),
                                                                                                                           fMaxJetTrackPt(30.0),
                                                                                                                           fMinJetClusPt(0.15),
                                                                                                                           fMinJetClusE(0.2),
                                                                                                                           fMinJetTowerE(0.2),
                                                                                                                           fTrackEtaMin(-1.0), fTrackEtaMax(1.0),
                                                                                                                           fTrackPhiMin(0.0), fTrackPhiMax(2.0 * TMath::Pi()),
                                                                                                                           fJetTrackEtaMin(-1.0), fJetTrackEtaMax(1.0),
                                                                                                                           fJetTrackPhiMin(0.0), fJetTrackPhiMax(2.0 * TMath::Pi()),
                                                                                                                           fJetTrackDCAcut(3.0),
                                                                                                                           fJetTracknHitsFit(15),
                                                                                                                           fJetTracknHitsRatio(0.52),
                                                                                                                           fTrackEfficiency(1.),
                                                                                                                           fJetTowerEMin(0.2), fJetTowerEMax(100.0),
                                                                                                                           fJetTowerEtaMin(-1.0), fJetTowerEtaMax(1.0),
                                                                                                                           fJetTowerPhiMin(0.0), fJetTowerPhiMax(2.0 * TMath::Pi()),
                                                                                                                           mTowerEnergyMin(0.2),
                                                                                                                           mHadronicCorrFrac(1.),
                                                                                                                           fJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks), // default is using all matched Tracks, Aug2019, per Hanseul
                                                                                                                           fLegacyMode(kFALSE),
                                                                                                                           fFillGhost(kFALSE),
                                                                                                                           fJets(0x0),
                                                                                                                           fJetsBGsub(0x0),
                                                                                                                           fJetsArr(0x0),
                                                                                                                           fJetsBGsubArr(0x0),
                                                                                                                           fFull_Event(0),
                                                                                                                           fConstituents(0),
                                                                                                                           mGeom(StEmcGeom::instance("bemc")),
                                                                                                                           mPicoDstMaker(0x0),
                                                                                                                           mPicoDst(0x0),
                                                                                                                           mPicoEvent(0x0),
                                                                                                                           mHIOverlay(0x0),
                                                                                                                           mCentMaker(0x0),
                                                                                                                           mBaseMaker(0x0),
                                                                                                                           mBemcGeom(0x0),
                                                                                                                           mEmcPosition(0x0),
                                                                                                                           grefmultCorr(0x0),
                                                                                                                           fEfficiencyInputFile(0x0),
                                                                                                                           fBackground(kFALSE),
                                                                                                                           fPrintLevel(0),
                                                                                                                           fMaxTowerEtBeforeHC(-99),
                                                                                                                           fMaxTowerEtAfterHC(-99)

{
  // Standard constructor.
  for (int i = 0; i < 8; i++)
  {
    fEmcTriggerArr[i] = kFALSE;
  }

  for (int i = 0; i < 4800; i++)
  {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;

    for (int j = 0; j < 7; j++)
      mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

  if (!name)
    return;
  SetName(name);

  fJetsArr.clear();
  fJetsBGsubArr.clear();
  rhovalue.clear();
  sigmavalue.clear();
}

//________________________________________________________________________
StHIRecoJets::~StHIRecoJets()
{

  if (mEmcPosition)
    delete mEmcPosition;

  // track reconstruction efficiency input file
  if (fEfficiencyInputFile)
  {
    fEfficiencyInputFile->Close();
    delete fEfficiencyInputFile;
  }
}
//
//
//________________________________________________________________________
Int_t StHIRecoJets::Init()
{
  // DeclareHistograms();

  // Create user objects.
  fJets = new TClonesArray("StJet");
  fJets->SetName(fJetsName);

  fJetsBGsub = new TClonesArray("StJet");
  fJetsBGsub->SetName(fJetsName + "BGsub");

  fJetsArr.clear();
  fJetsBGsubArr.clear();
  rhovalue.clear();
  sigmavalue.clear();

  // position object for Emc
  mBemcGeom = StEmcGeom::instance("bemc");
  mEmcPosition = new StEmcPosition2();

  // input file - for tracking efficiency: Run14 AuAu and Run12 pp
  const char *input = "";
  // if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) input=Form("./StRoot/StMyAnalysisMaker/Run14_AuAu_200_tracking_efficiency_and_momentum_smearing_dca_3p0_nhit_15_nhitfrac_0p52.root");
  if (fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200)
    input = Form("./StRoot/StMyAnalysisMaker/Run14_efficiencySmaller2D.root");
  if (fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB)
    input = Form("./StRoot/StMyAnalysisMaker/Run14_efficiencySmaller2D.root");
  if (fRunFlag == StJetFrameworkPicoBase::Run12_pp200)
    input = Form("./StRoot/StMyAnalysisMaker/Run12_efficiency_New.root"); // Oct17, 2019 added
  if (fDoEffCorr)
  {
    fEfficiencyInputFile = new TFile(input, "READ");
    if (!fEfficiencyInputFile)
      cout << Form("do not have input file: %s", input);
  }

  // ============================ set jet parameters for fastjet wrapper  =======================
  // recombination schemes:
  // E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, BIpt2_scheme, WTA_pt_scheme, WTA_modp_scheme
  fastjet::RecombinationScheme recombScheme;
  if (fRecombScheme == 0)
    recombScheme = fastjet::E_scheme;
  if (fRecombScheme == 1)
    recombScheme = fastjet::pt_scheme;
  if (fRecombScheme == 2)
    recombScheme = fastjet::pt2_scheme;
  if (fRecombScheme == 3)
    recombScheme = fastjet::Et_scheme;
  if (fRecombScheme == 4)
    recombScheme = fastjet::Et2_scheme;
  if (fRecombScheme == 5)
    recombScheme = fastjet::BIpt_scheme;
  if (fRecombScheme == 6)
    recombScheme = fastjet::BIpt2_scheme;
  if (fRecombScheme == 7)
    recombScheme = fastjet::WTA_pt_scheme;
  if (fRecombScheme == 8)
    recombScheme = fastjet::WTA_modp_scheme;
  if (fRecombScheme == 99)
    recombScheme = fastjet::external_scheme;

  // jet algorithm
  fastjet::JetAlgorithm algorithm;
  if (fJetAlgo == 1)
    algorithm = fastjet::antikt_algorithm;
  if (fJetAlgo == 0)
    algorithm = fastjet::kt_algorithm;
  // extra algorithms
  if (fJetAlgo == 2)
    algorithm = fastjet::cambridge_algorithm;
  if (fJetAlgo == 3)
    algorithm = fastjet::genkt_algorithm;
  if (fJetAlgo == 11)
    algorithm = fastjet::cambridge_for_passive_algorithm;
  if (fJetAlgo == 13)
    algorithm = fastjet::genkt_for_passive_algorithm;
  if (fJetAlgo == 99)
    algorithm = fastjet::plugin_algorithm;
  if (fJetAlgo == 999)
    algorithm = fastjet::undefined_jet_algorithm;
  fastjet::Strategy strategy = fastjet::Best;

  // setup fj wrapper
  fjw.SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw.SetStrategy(strategy);
  fjw.SetGhostArea(fGhostArea);
  fjw.SetR(fRadius);
  fjw.SetAlgorithm(algorithm);       // fJetAlgo);
  fjw.SetRecombScheme(recombScheme); // fRecombScheme);
  fjw.SetMaxRap(1.2);

  // ghost-area specifications
  double ghost_maxrap = 1.2;
  fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);

  // setting legacy mode
  // if(fLegacyMode) { fjw.SetLegacyMode(kTRUE); }

  return kStOK;
}
//
//
//_______________________________________________________________________________________
Int_t StHIRecoJets::Finish()
{

  return kStOK;
}

// Function: clear or delete objects after running
//________________________________________________________________________________
void StHIRecoJets::Clear(Option_t *opt)
{
  fJets->Clear();
  fJetsBGsub->Clear();
}
//
// Function: main loop, called for each event
//________________________________________________________________________________
int StHIRecoJets::Make()
{
  fjw.Clear();
  fFull_Event.clear();

  // cout << "Started with Reco Jets" << endl;
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // ZERO's out the jet array
  fJets->Delete();
  fJetsBGsub->Delete();

  fJetsArr.clear();
  fJetsBGsubArr.clear();
  rhovalue.clear();
  sigmavalue.clear();

  // ZERO these out for double checking they aren't set
  for (int i = 0; i < 4800; i++)
  {
    for (int j = 0; j < 7; j++)
      mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

  // get PicoDstMaker
  mPicoDstMaker = static_cast<StPicoDstMaker *>(GetMaker("picoDst"));
  if (!mPicoDstMaker)
  {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst *>(mPicoDstMaker->picoDst());
  if (!mPicoDst)
  {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent
  mPicoEvent = static_cast<StPicoEvent *>(mPicoDst->event());
  if (!mPicoEvent)
  {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer - this class does not inherit from base class: StJetFrameworkPicoBase, but we want to reduce redundancy
  mBaseMaker = static_cast<StJetFrameworkPicoBase *>(GetMaker("baseClassMaker"));
  if (!mBaseMaker)
  {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();
  deadTowers = mBaseMaker->GetDeadTowers();
  badTowers = mBaseMaker->GetBadTowers();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if (doRejectBadRuns)
  {
    if (!mBaseMaker->IsRunOK(fRunNumber))
      return kStOK;
  }

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  mHIOverlay = static_cast<StHIOverlay_Test *>(GetMaker("HIOverlay"));
  if (!mHIOverlay)
  {
    LOG_WARN << "No HI Overlay! Skip!" << endm;
    return kStFatal;
  }

  fRecoTracks = mHIOverlay->GetRecoTracks();
  fRecoTowers = mHIOverlay->GetRecoTowers();

  int numberofrecoevents = mHIOverlay->GetTheNumberOfEventsToOverLay();

  fMCEventInfo = mHIOverlay->GetMCEventInfo();

  // cout << "==========================================================================" << endl;
  // cout << "The number of events actually overlaid = " << numberofrecoevents << endl;

  if (numberofrecoevents == 0)
    return kStOK;

  for (int iter = 0; iter < numberofrecoevents; iter++)
  {

    if (!fBackground && fPrintLevel)
      cout << "Event Number  = " << iter << endl;

    for (int i = 0; i < 4800; i++)
    {
      fTowerToTriggerTypeHT1[i] = kFALSE;
      fTowerToTriggerTypeHT2[i] = kFALSE;
      fTowerToTriggerTypeHT3[i] = kFALSE;

      for (int j = 0; j < 7; j++)
        mTowerMatchTrkIndex[i][j] = -1;
      mTowerStatusArr[i] = 0;
    }

    // cout << "======================= Event Number After MC = " << fMCEventInfo[iter].first << "\t" << fMCEventInfo[iter].second << "\t" << "=======================" << endl;
    // cout << "Iteration = " << iter << endl;
    FindJets(iter);
    if (!doConstituentSubtr)
      FillJetBranch(iter);
    else
      FillJetBGBranch(iter);
  }

  return kStOK;
}

//
// Main class to FindJets from tracks + towers
//________________________________________________________________________
void StHIRecoJets::FindJets(int iteration)
{
  // clear out existing wrapper object
  fjw.Clear();
  fFull_Event.clear();

  fJets->Delete();
  fJetsBGsub->Delete();

  // assume neutral pion mass
  // additional parameters constructed
  double pi = 1.0 * TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV
  unsigned int ntracks = mPicoDst->numberOfTracks();

  // loop over ALL tracks in PicoDst and add to jet, after acceptance and quality cuts
  if ((fJetType == kFullJet) || (fJetType == kChargedJet))
  {
    fMaxTrackPt = -99;

    for (unsigned short iTracks = 0; iTracks < ntracks; iTracks++)
    {
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack *>(mPicoDst->track(iTracks));
      if (!trk)
      {
        continue;
      }

      // acceptance and kinematic quality cuts - pt cut is also applied here currently
      if (!AcceptJetTrack(trk, Bfield, mVertex))
      {
        continue;
      }

      // get momentum vector of track - global or primary track
      TVector3 mTrkMom;
      if (doUsePrimTracks)
      {
        // get primary track vector
        mTrkMom = trk->pMom();
      }
      else
      {
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield);
      }

      // track variables
      double pt = mTrkMom.Perp();
      double phi = mTrkMom.Phi();
      if (phi < 0.0)
        phi += 2.0 * pi; // force from 0-2pi
      if (phi > 2.0 * pi)
        phi -= 2.0 * pi; // force from 0-2pi
      double eta = mTrkMom.PseudoRapidity();
      double px = mTrkMom.x();
      double py = mTrkMom.y();
      double pz = mTrkMom.z();
      double p = mTrkMom.Mag();
      double energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
      short charge = trk->charge();
      int matchedTowerIndex = abs(GetMatchedBtowID(trk)) - 1; // towerIndex = towerID - 1
      // cout<<"Charge: "<<charge<<"  nHitsFit: "<<trk->nHitsFit()<<endl;

      if (pt > fMaxTrackPt)
        fMaxTrackPt = pt;

      fjw.AddInputVector(px, py, pz, energy, iTracks); // includes E
      fastjet::PseudoJet particle_Track(px, py, pz, energy);
      particle_Track.set_user_index(iTracks);
      fFull_Event.push_back(particle_Track);

      //====  matched track index ===
      int trackIndex = iTracks;
      if (trackIndex < 0)
      {
        continue;
      } // can't happen
      if (matchedTowerIndex < 0)
        continue;
      mTowerMatchTrkIndex[matchedTowerIndex][mTowerStatusArr[matchedTowerIndex]] = trackIndex;
      mTowerStatusArr[matchedTowerIndex] = mTowerStatusArr[matchedTowerIndex] + 1; // 1+ means match, 0 for no match
    }                                                                              // track loop

  } // if full/charged jets

  // full or neutral jets - get towers and apply hadronic correction
  if ((fJetType == kFullJet) || (fJetType == kNeutralJet))
  {

    fMaxTowerEtBeforeHC = -99;
    fMaxTowerEtAfterHC = -99;

    int matchedTowerTrackCounter = 0;

    // print
    // int nTracks = mPicoDst->numberOfTracks();               // number of tracks
    // int nTrigs = mPicoDst->numberOfEmcTriggers();           // number of Emc triggers
    // int nBTowHits = mPicoDst->numberOfBTowHits();           // barrel tower hits, always 4800
    int nTowers = mPicoDst->numberOfBTowHits(); // barrel tower hists, always 4800
    // int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits(); // number of BEMC matched tracks
    // cout<<"nTracks = "<<nTracks<<"  nTrigs = "<<nTrigs<<"  nBTowHits = "<<nBTowHits<<"  nBEmcPidTraits = "<<nBEmcPidTraits<<endl;

    // loop over towers and add input vectors to fastjet
    for (int itow = 0; itow < nTowers; itow++)
    {
      // get tower pointer
      StPicoBTowHit *tower = static_cast<StPicoBTowHit *>(mPicoDst->btowHit(itow));
      if (!tower)
      {
        cout << "No tower pointer... iTow = " << itow << endl;
        continue;
      }

      // tower ID - get from itow shift: which is 1 more than array index
      // itow = matchedTowerIndex
      int towerID = itow + 1;
      int towerIndex = towerID - 1;
      if (towerID < 0)
        continue; // double check these aren't still in the event list

      // tower acceptance cuts - cuts on *bad* towers also
      if (!AcceptJetTower(tower, towerID))
        continue;

      // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
      TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
      double towerPhi = towerPosition.Phi();
      if (towerPhi < 0.0)
        towerPhi += 2.0 * pi; // force from 0-2pi
      if (towerPhi > 2.0 * pi)
        towerPhi -= 2.0 * pi; // force from 0-2pi
      double towerEta = towerPosition.PseudoRapidity();
      int towerADC = tower->adc();
      double towerEunCorr = tower->energy(); // uncorrected energy
      double towerE = tower->energy();       // corrected energy (hadronically - done below)
      double towEtunCorr = towerEunCorr / (1.0 * TMath::CosH(towerEta));

      // cut on min tower energy after filling histos
      if (towerEunCorr < mTowerEnergyMin)
        continue; // if we don't have enough E to start with, why mess around

      if (towerEunCorr > fMaxTowerEtBeforeHC)
        fMaxTowerEtBeforeHC = towerEunCorr;

      // =======================================================================
      // HADRONIC CORRECTION
      double maxEt = 0.;
      double sumEt = 0.;

      // if tower was is matched to a track or multiple, add up the matched track energies - (mult opt.) to then subtract from the corresponding tower
      // August 15: if *have* 1+ matched trk-tow AND uncorrected energy of tower is at least your tower constituent cut, then CONTINUE
      if (mTowerStatusArr[towerIndex] > 0.5 && towerEunCorr > mTowerEnergyMin)
      {
        double maxE = 0.0;
        double sumE = 0.0;
        // =======================================================================================================================
        // --- finds max E track matched to tower *AND* the sum of all matched track E and subtract from said tower
        //     USER provides readMacro.C which method to use for their analysis via SetJetHadCorrType(type);
        // loop over ALL matched tracks
        for (int itrk = 0; itrk < mTowerStatusArr[towerIndex]; itrk++)
        {
          StPicoTrack *trk = static_cast<StPicoTrack *>(mPicoDst->track(mTowerMatchTrkIndex[towerIndex][itrk]));
          if (!trk)
          {
            cout << "No track pointer..." << endl;
            continue;
          }
          if (!AcceptTrack(trk, Bfield, mVertex))
          {
            cout << "track matched back doesn't pass cuts" << endl;
            continue;
          }

          // get track variables to matched tower from 3-vector
          TVector3 mTrkMom;
          if (doUsePrimTracks)
          { // get primary track vector
            mTrkMom = trk->pMom();
          }
          else
          { // get global track vector
            mTrkMom = trk->gMom(mVertex, Bfield);
          }

          // track variables
          double p = mTrkMom.Mag();
          double E = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
          if (E > maxE)
            maxE = E;
          sumE = sumE + E;

          // test statement
          // cout<<itrk<<"  itrkE: "<<E<<"  sumtrkE: "<<sumE<<endl;
        } // track loop

        // apply hadronic correction to tower
        maxEt = (towerEunCorr - (mHadronicCorrFrac * maxE)) / (1.0 * TMath::CosH(towerEta));
        sumEt = (towerEunCorr - (mHadronicCorrFrac * sumE)) / (1.0 * TMath::CosH(towerEta));
        //=================================================================================================================
      } // have a track-tower match
      // else - no match so treat towers on their own. Must meet constituent cut

      // Et - correction comparison
      double fMaxEt = (maxEt == 0) ? towerEunCorr / (1.0 * TMath::CosH(towerEta)) : maxEt;
      double fSumEt = (sumEt == 0) ? towerEunCorr / (1.0 * TMath::CosH(towerEta)) : sumEt;
      // if(mTowerStatusArr[towerIndex] > 0.5) cout<<"towerEunCorr = "<<towerEunCorr<<"  CosH: "<<1.0*TMath::CosH(towerEta)<<"   fMaxEt: "<<fMaxEt<<"   fSumEt: "<<fSumEt<<endl;

      // cut on transverse tower energy (more uniform)
      double towerEt = 0.0;
      if (mTowerStatusArr[towerIndex] < 1)
      { // no matches, use towers uncorrected energy
        towerEt = towerEunCorr / (1.0 * TMath::CosH(towerEta));
      }
      else
      {
        if (fJetHadCorrType == StJetFrameworkPicoBase::kHighestEMatchedTrack)
        {
          towerEt = fMaxEt;
          towerE = fMaxEt * 1.0 * TMath::CosH(towerEta);
        }
        if (fJetHadCorrType == StJetFrameworkPicoBase::kAllMatchedTracks)
        {
          towerEt = fSumEt;
          towerE = fSumEt * 1.0 * TMath::CosH(towerEta);
        }
      }
      if (towerEt == 0)
      {
        cout << "fJetHadCorrType - or - towerE actually 0" << endl;
      } // it was unset, because you provided wrong fJetHadCorrType
      if (towerEt < 0)
        towerEt = 0.0;
      if (towerEt < mTowerEnergyMin)
        continue;

      if (towerEt > fMaxTowerEtAfterHC)
        fMaxTowerEtAfterHC = towerEt;

      // get components from Energy (p - momentum) - the below lines 'converts' the tower energy to momentum:
      TVector3 mom;
      GetMomentum(mom, tower, pi0mass, towerID, towerE);
      double towerPx = mom.x();
      double towerPy = mom.y();
      double towerPz = mom.z();

      // add towers to fastjet - shift tower index (tracks 0+, ghosts = -1, towers < -1)
      int uidTow = -(itow + 2);

      fjw.AddInputVector(towerPx, towerPy, towerPz, towerE, uidTow); // includes E
      fastjet::PseudoJet particle_Track(towerPx, towerPy, towerPz, towerE);
      particle_Track.set_user_index(uidTow);
      fFull_Event.push_back(particle_Track);

    } // tower loop

  } // neutral/full jets

  // Everything upto this will be the same for multiple events. But since we are embedding multiple MC events in one RECO event,
  // from this point on, things will be different for every iteration.

  if (!fBackground)
  {
    // Loop over all saved particles for MC in StHIOverlay_Test vectors and enter them into fastjet wrapper
    for (unsigned short iTracks = 0; iTracks < fRecoTracks[iteration].size(); iTracks++)
    {

      TLorentzVector v;
      v = fRecoTracks[iteration][iTracks];
      // track variables

      double px = v.X();
      double py = v.Y();
      double pz = v.Z();
      double p = v.P();
      double energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
      double mass = v.M();

      if (int(v.M()) == 1 && !fBackground && fPrintLevel)
        cout << "Reco D0 In Jet Loop = " << v.Pt() << "\t" << v.PseudoRapidity() << endl;

      // MC Tracks will start from 10000
      fjw.AddInputVector(px, py, pz, energy, iTracks + 10000); // includes E

      // if (int(v.M()) != 1) cout << Form("Reco Level Track Input = %.2f \t %.2f \t %.2f", px, py, pz) << endl;
      // else cout << Form("Reco Level Track Input D0 = %.2f \t %.2f \t %.2f", px, py, pz) << endl;

    } // track loop

    if (fJetType == kFullJet)
    {
      for (unsigned short iTowers = 0; iTowers < fRecoTowers[iteration].size(); iTowers++)
      {

        TLorentzVector v;
        v = fRecoTowers[iteration][iTowers];
        // tower variables

        double px = v.X();
        double py = v.Y();
        double pz = v.Z();
        double p = v.P();
        double energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
        double mass = v.M();

        // MC Towers will start from
        fjw.AddInputVector(px, py, pz, energy, -1 * iTowers - 10000); // includes E

        // cout << "Reco Level Tower Input = " << px << "\t" << py << "\t" << pz << endl;
        // cout << Form("Reco Level Tower Input = %.2f \t %.2f \t %.2f", px, py, pz) << endl;

      } // tower loop

    } // if full/charged jets

    // cout << "Input Loop Size = " << fRecoTracks[iteration].size() + fRecoTowers[iteration].size() << endl;
    // cout << "Input Loop Size = " << fRecoTracks[iteration].size() + fRecoTowers[iteration].size() << endl;
  }

  // run jet finder
  fjw.Run();
  // fjw.PrintInput();
}

void StHIRecoJets::FillJetBGBranch(int iteration)
{
  fJets->Delete();
  // E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, BIpt2_scheme, WTA_pt_scheme, WTA_modp_scheme
  fastjet::RecombinationScheme recombScheme;
  if (fRecombScheme == 0)
    recombScheme = fastjet::E_scheme;
  if (fRecombScheme == 1)
    recombScheme = fastjet::pt_scheme;
  if (fRecombScheme == 2)
    recombScheme = fastjet::pt2_scheme;
  if (fRecombScheme == 3)
    recombScheme = fastjet::Et_scheme;
  if (fRecombScheme == 4)
    recombScheme = fastjet::Et2_scheme;
  if (fRecombScheme == 5)
    recombScheme = fastjet::BIpt_scheme;
  if (fRecombScheme == 6)
    recombScheme = fastjet::BIpt2_scheme;

  // jet algorithm
  fastjet::JetAlgorithm algorithm;
  if (fJetAlgo == 1)
    algorithm = fastjet::antikt_algorithm;
  if (fJetAlgo == 0)
    algorithm = fastjet::kt_algorithm;
  if (fJetAlgo == 2)
    algorithm = fastjet::cambridge_algorithm;
  if (fJetAlgo == 11)
    algorithm = fastjet::cambridge_for_passive_algorithm;
  fastjet::Strategy strategy = fastjet::Best;

  double jetAbsRapMax = 1.0;

  // create a jet definition for the clustering: We use the anti-kt algorithm with a radius of 0.5
  //----------------------------------------------------------
  fastjet::JetDefinition jet_def(algorithm, fRadius, recombScheme, strategy);

  // create an area definition for the clustering
  // ----------------------------------------------------------
  // ghosts should go up to the acceptance of the detector or (with infinite acceptance) at least 2R beyond the region where you plan to investigate jets.
  double ghost_maxrap = 1.2;
  fastjet::GhostedAreaSpec area_spec(ghost_maxrap, 1, fGhostArea);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);

  // run the jet clustering with the above jet and area definitions for both the hard and full event
  //
  // We retrieve the jets above 7 GeV in both case (note that the 7-GeV cut we be applied again later on after we subtract the jets from the full event)
  // ----------------------------------------------------------
  // fastjet::ClusterSequenceArea clust_seq_full(fFull_Event, jet_def, area_def);

  // minimum jet pt for inclusive jets
  double ptmin = fMinJetPt;
  // Selector sel_jets = SelectorNHardest(2) * SelectorAbsRapMax(3.0);
  // vector<fastjet::PseudoJet> full_jets = sorted_by_pt(clust_seq_full.inclusive_jets(ptmin));
  // vector<fastjet::PseudoJet> full_jets = fjw.GetInclusiveJets();

  // Now turn to the estimation of the background (for the full event)
  //
  // There are different ways to do that. In general, this also requires clustering the particles that will be handled internally in FastJet.
  //
  // The suggested way to proceed is to use a BackgroundEstimator constructed from the following 3 arguments:
  //  - a jet definition used to cluster the particles.
  //    . We strongly recommend using the kt or Cambridge/Aachen algorithm (a warning will be issued otherwise)
  //    . The choice of the radius is a bit more subtle. R=0.4 has been chosen to limit the impact of hard jets; in samples of
  //      dominantly sparse events it may cause the UE/pileup to be underestimated a little, a slightly larger value (0.5 or 0.6) may be better.
  //  - An area definition for which we recommend the use of explicit ghosts (i.e. active_area_explicit_ghosts) As mentionned in the area example (06-area.cc), ghosts should
  //    extend sufficiently far in rapidity to cover the jets used in the computation of the background (see also the comment below)
  //  - A Selector specifying the range over which we will keep the jets entering the estimation of the background (you should
  //    thus make sure the ghosts extend far enough in rapidity to cover the range, a warning will be issued otherwise).
  //    In this particular example, the two hardest jets in the event are removed from the background estimation
  //
  // logical ops:  product*, not !, or ||, and &&
  // ----------------------------------------------------------

  // create what we need for the background estimation
  //----------------------------------------------------------
  fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm, fRadius, recombScheme, strategy);
  // fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetAbsRapMax) * (!fastjet::SelectorNHardest(2));
  fastjet::Selector rho_range_sel = fastjet::SelectorAbsRapMax(3.0); // 3.0
  fastjet::Selector eta_range_sel = fastjet::SelectorAbsEtaMax(0.6); // 0.6
  fastjet::Selector hard_jet_sel = (!fastjet::SelectorNHardest(2));
  fastjet::Selector full_selector = hard_jet_sel * eta_range_sel; // make this the main line as of March 31, 2020 - others should of been aware in past
                                                                  //  fastjet::ClusterSequenceArea clust_seq_rho(fFull_Event, jet_def, area_def); // not used FIXME

  fastjet::JetMedianBackgroundEstimator bge_rho(full_selector, jet_def_for_rho, area_def);
  // TODO next 2 lines commented out to suppress warnings, doesn't affect results - Sept26, 2018
  // fastjet::BackgroundJetScalarPtDensity *scalarPtDensity = new fastjet::BackgroundJetScalarPtDensity();
  // bge_rho.set_jet_density_class(scalarPtDensity); // this changes computation of pt of patches from vector sum to scalar sum. Theor., the scalar sum seems more reasonable.
  //  bge_rho.reset();
  bge_rho.set_particles(fFull_Event);

  //  cout << Form("===== %s =====", GetName()) << endl;
  //  cout << "Event pT Average = " << tmptotal << endl;
  //  cout << bge_rho.description() << endl;

  double rho = bge_rho.rho();
  double sigma = bge_rho.sigma();
  double mean_area = bge_rho.mean_area();

  cout << Form("%s::Rho = %i\t%f", GetName(), fFull_Event.size(), bge_rho.rho()) << endl;

  // tmp_full_event[id0candidate] = fFull_Event;

  // subtractor:
  //----------------------------------------------------------
  fastjet::contrib::ConstituentSubtractor subtractor(&bge_rho);
  subtractor.set_common_bge_for_rho_and_rhom(true); // TODO - need this to omit warning, same results Sept26, 2018

  int d0UserIndex = -999999;
  for (unsigned short iTracks = 0; iTracks < fRecoTracks[iteration].size(); iTracks++)
  {
    TLorentzVector v;
    v = fRecoTracks[iteration][iTracks];

    if (int(v.M()) == 1)
    {
      d0UserIndex = iTracks + 10000; // In HI, because of the way I have assigned D0 user index, I need to figure out which track is the D0 in every iteration.
      break;
    }
  }

  fastjet::Selector notD0 = (!fastjet::SelectorUserIndex(d0UserIndex)); // This is not used here, because only the MB event is used for the background estimation

  subtractor.set_particle_selector(&notD0);

  // this sets the same background estimator to be used for deltaMass density, rho_m, as for pt density, rho:
  // subtractor.use_common_bge_for_rho_and_rhom(true); // for massless input particles it does not make any difference (rho_m is always zero)
  ////cout << subtractor.description() << endl;
  ////cout << "  Giving, for the full event" << "    rho     = " << bge_rho.rho() << "    sigma   = " << bge_rho.sigma() << endl;

  // fill histogram with FastJet calculated rho
  // fHistFJRho->Fill(rho);
  // fHistFJSigma->Fill(sigma);

  // sort jets according to jet pt
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  // ===================================
  // clear constituent array
  fConstituents.clear();

  // loop over FastJet jets
  __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("%d jets found", (Int_t)jets_incl.size()));
  for (UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet)
  {
    Int_t ij = indexes[ijet];
    __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)));

    // PERFORM CUTS ON inclusive JETS before saving
    // cut on min jet pt
    if (jets_incl[ij].perp() < fMinJetPt)
      continue;
    // cut on min jet area
    if (fjw.GetJetArea(ij) < fMinJetArea * TMath::Pi() * fRadius * fRadius)
      continue;
    // cut on eta acceptance
    if ((jets_incl[ij].eta() < fJetEtaMin) || (jets_incl[ij].eta() > fJetEtaMax))
      continue;
    // cut on phi acceptance
    if ((jets_incl[ij].phi() < fJetPhiMin) || (jets_incl[ij].phi() > fJetPhiMax))
      continue;

    // ===========================================================
    // =========  apply subtractor here ===========
    fastjet::PseudoJet subtracted_jet = subtractor(jets_incl[ij]);

    // FIXME - may want to just make this 'fJets'

    // get constituents of jets - the se should have the background (fraction) subtracted
    // fConstituents = fjw.GetJetConstituents(ij);
    // jet->SetJetConstituents(fConstituents);

    // fill jet constituents - these are identified by their index
    vector<fastjet::PseudoJet> constituents = subtracted_jet.constituents();
    //  cout << "Constituent Size = " << constituents.size() << endl;

    //  for (int i = 0; i < constituents.size(); i++){
    //   if (constituents[i].user_index() == 99999){
    //     cout << GetName() << "\t" << "Jet # " << ijet << " Particle with index #" << i << " D0 found with pT " << constituents[i].perp() << endl;
    //   }
    //  }

    bool d0Jet = kFALSE;
    int d0index = -99;
    double deltar = -999;

    for (UInt_t ic = 0; ic < constituents.size(); ++ic)
    {
      // get user defined index
      Int_t uid = constituents[ic].user_index();

      if (uid >= 10000)
      {
        double mass = fRecoTracks[iteration][uid - 10000].M();
        if ((int)mass == 1)
        {
          d0Jet = kTRUE;
          d0index = ic;
          deltar = TMath::Sqrt(pow(subtracted_jet.eta() - constituents[ic].eta(), 2) + pow(dPhi(subtracted_jet.phi(), constituents[ic].phi()), 2));
          break;
        }
      }
    }

    if (!d0Jet)
      continue;

    StJet *jet = new ((*fJets)[jetCount])
        StJet(subtracted_jet.perp(), subtracted_jet.eta(), subtracted_jet.phi(), subtracted_jet.m());

    if (fPrintLevel)
    {
      if (abs(fRecoTracks[iteration][d0UserIndex - 10000].Pt() - constituents[d0index].pt()) > 0.001 || constituents[d0index].pt() > jet->Pt())
      {
        cout << "*************** ERROR ***************" << endl;
      }
      if (doConstituentSubtr)
        cout << GetName() << "\t"
             << "D0 User Index Before And After:" << d0UserIndex << "\t" << constituents[d0index].user_index() << "\t" << fRecoTracks[iteration][d0UserIndex - 10000].Pt() << "\t" << constituents[d0index].pt() << "\t" << jet->Pt() << endl;
    }
    // set label
    jet->SetLabel(ij);

    // area vector and components
    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij)); // FIXME - this MAY not be correct after subtraction, double check
    jet->SetArea(area.perp());                         // same as fjw.GetJetArea(ij)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    jet->SetJetConstituents(constituents); // constituents are the corrected pseudojet objects

    // cout << Form("%s:: Constituent Size = %i", GetName(), jet->GetJetConstituents().size()) << endl;
    FillJetConstituents(jet, constituents, constituents);

    __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("Added jet n. %d, pt = %f, area = %f, constituents = %d", jetCount, jet->Pt(), jet->Area(), jet->GetNumberOfConstituents()));

    jetCount++;

    if (d0Jet && fPrintLevel)
      cout << Form("%s:: Event Info = %i, Rho = %f; Jet pT = %f \t %f;  D0 pT = %f, Delta R = %f", GetName(), fFull_Event.size(), rho, jets_incl[ij].perp(), jet->Pt(), constituents[d0index].perp(), deltar) << endl;
  } // jet loop

  TClonesArray *tempjets = (TClonesArray *)fJets->Clone();
  fJetsArr.push_back(tempjets);

  rhovalue.push_back(rho);
  sigmavalue.push_back(sigma);

  bge_rho.reset();
}

//
/**
 * This method fills the jet output branch (TClonesArray) with the jet found by the FastJet
 * wrapper. Before filling the jet branch, the utilities are prepared. Then the utilities are
 * called for each jet and finally after jet finding the terminate method of all utilities is called.
 */
void StHIRecoJets::FillJetBranch(int iteration)
{

  // cout << "Called FillJetBranch." << endl;
  // get inclusive jets
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();

  fJets->Delete();

  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  // loop over FastJet jets
  __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("%d jets found", (Int_t)jets_incl.size()));

  bool D0Jet = kFALSE;

  for (UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet)
  {
    D0Jet = kFALSE;
    Int_t ij = indexes[ijet];

    __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)));

    if (fPrintLevel)
      cout << "Jets = " << ijet << "\t" << jets_incl[ij].perp() << "\t" << jets_incl[ij].eta() << "\t" << jets_incl[ij].phi() << endl;

    // March 8, 2018 - probably don't need this anymore after coming up with method to pass and extract constituents via index!
    // get constituents of jets
    fConstituents = fjw.GetJetConstituents(ij);

    int D0id = -99;

    for (int jtrk = 0; jtrk < fConstituents.size(); jtrk++)
    {
      Int_t uid = fConstituents[jtrk].user_index();

      if (uid >= 10000)
      {
        double mass = fRecoTracks[iteration][uid - 10000].M();

        if (fPrintLevel)
        {
          if (uid >= 10000)
            cout << "Track ID = " << uid << "\t" << mass << endl;
        }

        if (int(mass) == 1)
        {
          D0id = jtrk;
          D0Jet = kTRUE;
          break;
        }
      }
    }

    // PERFORM CUTS ON inclusive JETS before saving
    // cut on min jet pt
    if (jets_incl[ij].perp() < fMinJetPt)
      continue;
    // cut on min jet area
    if (fjw.GetJetArea(ij) < fMinJetArea * TMath::Pi() * fRadius * fRadius)
      continue;
    // cut on eta acceptance
    if ((jets_incl[ij].eta() < fJetEtaMin) || (jets_incl[ij].eta() > fJetEtaMax))
      continue;
    // cut on phi acceptance
    if ((jets_incl[ij].phi() < fJetPhiMin) || (jets_incl[ij].phi() > fJetPhiMax))
      continue;

    if (!fBackground)
    {
      if (!D0Jet)
        continue;
    }

    // need to figure out how to get m or E from STAR tracks
    StJet *jet = new ((*fJets)[jetCount])
        StJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());

    // if (D0Jet && !fBackground && fPrintLevel) cout << "Found Reco jet with Jet pT = " << jet->Pt() << "\t" << fConstituents.size() << endl;
    // cout << "Found Reco jet with Jet pT = " << jet->Pt() << "\t" << fConstituents.size() << endl;

    // cout << Form("Found Reco jet with Jet pT eta and D0 pT eta = %.2f \t %.2f \t %.2f \t %.2f ", jet->Pt(), jet->Eta(), fConstituents[D0id].perp(), fConstituents[D0id].eta()) << endl;

    // for(UInt_t ic = 0; ic < fConstituents.size(); ++ic) {
    //   if (fConstituents[ic].perp() >= 0.2) cout << Form("Found Constituent pT eta = %.2f \t %.2f ", fConstituents[ic].perp(), fConstituents[ic].eta()) << endl;
    // }

    jet->SetLabel(ij);

    // area vector and components
    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij));
    jet->SetArea(area.perp()); // same as fjw.GetJetArea(ij)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    // fill jet constituents
    vector<fastjet::PseudoJet> constituents = fjw.GetJetConstituents(ij);
    jet->SetJetConstituents(fConstituents);

    FillJetConstituents(jet, constituents, constituents);

    __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("Added jet n. %d, pt = %f, area = %f, constituents = %d", jetCount, jet->Pt(), jet->Area(), jet->GetNumberOfConstituents()));

    jetCount++;
  } // jet loop

  TClonesArray *tempjets = (TClonesArray *)fJets->Clone();

  // cout << "Reco Jet Entries = " << tempjets->GetEntries() << endl;

  fJetsArr.push_back(tempjets);

  // StJet *rjet = static_cast<StJet*>(fJetsArr[iteration]->At(0));

  // cout << "pT = " << rjet->Pt() << endl;

  // cout << "Reco Big Jet Entries = " << fJetsArr.back()->GetEntries() << endl;
}
//
/**
 * This method is called for each jet. It loops over the jet constituents and
 * adds them to the jet object.
 * @param jet Pointer to the AliEmcalJet object where the jet constituents will be added
 * @param constituents List of the jet constituents returned by the FastJet wrapper
 * @param constituents_unsub List of jet constituents before background subtraction
 * @param flag If kTRUE it means that the argument "constituents" is a list of subtracted constituents
 * @param particles_sub Array containing subtracted constituents - not used
 */
void StHIRecoJets::FillJetConstituents(StJet *jet, std::vector<fastjet::PseudoJet> &constituents,
                                       std::vector<fastjet::PseudoJet> &constituents_unsub, Int_t flag, TString particlesSubName)
{
  // initialize some variables/counters
  Double_t neutralE = 0, maxTrack = 0, maxTower = 0;
  Int_t nt = 0; // track counter
  Int_t nc = 0; // tower (cluster) counter
  Int_t ng = 0; // ghost counter
  double pi = 1.0 * TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // initially set track and tower constituent sizes
  jet->SetNumberOfTracks(constituents.size());
  jet->SetNumberOfTowers(constituents.size());

  // loop over constituents for ij'th jet
  for (UInt_t ic = 0; ic < constituents.size(); ++ic)
  {
    // get user defined index
    Int_t uid = constituents[ic].user_index();

    // CHARGED COMPONENT (tracks)
    if (uid >= 0)
    {
      jet->AddTrackAt(uid, nt);
      nt++;
    }

    else if (uid == -1)
    { // Ghosty bois
      ng++;
    }
    else if (uid > -10000 && uid < -1)
    {
      // convert uid to tower index (index of tower - in BTowHit array)
      Int_t towIndex = -(uid + 2); // 1 less than towerID
      jet->AddTowerAt(towIndex, nc);
      nc++;
    }
    else if (uid <= -10000)
    { //
      jet->AddTowerAt(abs(uid), nc);
      nc++;
    }

  } // end of constituent loop

  // set some jet properties
  jet->SetNumberOfTracks(nt);
  jet->SetNumberOfTowers(nc);
}
/**
 * Sorts jets by pT (decreasing)
 * @param[out] indexes This array is used to return the indexes of the jets ordered by pT
 * @param[in] array Vector containing the list of jets obtained by the FastJet wrapper
 * @return kTRUE if at least one jet was found in array; kFALSE otherwise
 */
Bool_t StHIRecoJets::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
{
  static Float_t pt[9999] = {0};
  const Int_t n = (Int_t)array.size();
  if (n < 1)
    return kFALSE;

  for (Int_t i = 0; i < n; i++)
    pt[i] = array[i].perp();

  TMath::Sort(n, pt, indexes);

  return kTRUE;
}
/**
 * An instance of this class can be "locked". Once locked, it cannot be unlocked.
 * If the instance is locked, attempting to change the configuration will throw a
 * fatal and stop the execution of the program. This method checks whether the instance
 * is locked and throw a fatal if it is locked.
 */
Bool_t StHIRecoJets::IsLocked() const
{
  if (fLocked)
  {
    Form("Jet finder task is locked! Changing properties is not allowed.");
    return kStFatal;
  }
  else
  {
    return kStOK;
  }
}
//
// Function: track quality cuts
//________________________________________________________________________
Bool_t StHIRecoJets::AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert)
{
  // constants: assume neutral pion mass
  /// double pi0mass = Pico::mMass[0]; // GeV
  double pi = 1.0 * TMath::Pi();

  // get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  if (doUsePrimTracks)
  {
    if (!(trk->isPrimary()))
      return kFALSE; // check if primary
    // get primary track vector
    mTrkMom = trk->pMom();
  }
  else
  {
    // get global track vector
    mTrkMom = trk->gMom(Vert, B);
  }

  // track variables
  double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  double dca = trk->gDCA(Vert).Mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0 * nHitsFit / nHitsMax;

  // jet track acceptance cuts now - note difference from AcceptJetTrack()
  if ((eta < fTrackEtaMin) || (eta > fTrackEtaMax))
    return kFALSE;
  if (phi < 0.0)
    phi += 2.0 * pi; // force from 0-2pi
  if (phi > 2.0 * pi)
    phi -= 2.0 * pi; // force from 0-2pi
  if ((phi < fTrackPhiMin) || (phi > fTrackPhiMax))
    return kFALSE;

  // additional quality cuts for tracks
  if (dca > fJetTrackDCAcut)
    return kFALSE;
  if (nHitsFit < fJetTracknHitsFit)
    return kFALSE;
  if (nHitsRatio < fJetTracknHitsRatio)
    return kFALSE;

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}
//
// Function: jet track quality cuts
// - this function should only be used for jet constituent tracks
// 	NOT when considering track-tower matches  (Aug 19')
//________________________________________________________________________
Bool_t StHIRecoJets::AcceptJetTrack(StPicoTrack *trk, Float_t B, TVector3 Vert)
{
  // constants: assume neutral pion mass
  /// double pi0mass = Pico::mMass[0]; // GeV
  double pi = 1.0 * TMath::Pi();

  // get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  if (doUsePrimTracks)
  {
    if (!(trk->isPrimary()))
      return kFALSE; // check if primary
    // get primary track vector
    mTrkMom = trk->pMom();
  }
  else
  {
    // get global track vector
    mTrkMom = trk->gMom(Vert, B);
  }

  // track variables
  double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  // double p = mTrkMom.Mag();
  // double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
  // short charge = trk->charge();
  // double dca = (trk->dcaPoint() - mPicoEvent->primaryVertex()).mag();
  double dca = trk->gDCA(Vert).Mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0 * nHitsFit / nHitsMax;

  // jet track acceptance cuts now
  if (pt < fMinJetTrackPt)
    return kFALSE;
  if (pt > fMaxJetTrackPt)
    return kFALSE; // 20.0 STAR, 100.0 ALICE
  if ((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax))
    return kFALSE;
  if (phi < 0.0)
    phi += 2.0 * pi; // force from 0-2pi
  if (phi > 2.0 * pi)
    phi -= 2.0 * pi; // force from 0-2pi
  if ((phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax))
    return kFALSE;

  // additional quality cuts for tracks
  if (dca > fJetTrackDCAcut)
    return kFALSE;
  if (nHitsFit < fJetTracknHitsFit)
    return kFALSE;
  if (nHitsRatio < fJetTracknHitsRatio)
    return kFALSE;

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}
//
// Tower Quality Cuts
//________________________________________________________________________
Bool_t StHIRecoJets::AcceptJetTower(StPicoBTowHit *tower, Int_t towerID)
{
  // constants:
  double pi = 1.0 * TMath::Pi();

  // tower ID - passed into function: make sure some of these aren't still in event array
  if (towerID < 0)
    return kFALSE;

  // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
  TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
  double phi = towerPosition.Phi();
  if (phi < 0.0)
    phi += 2.0 * pi; // force from 0-2pi
  if (phi > 2.0 * pi)
    phi -= 2.0 * pi; // force from 0-2pi
  double eta = towerPosition.PseudoRapidity();

  // check for bad (and dead) towers
  bool TowerOK = mBaseMaker->IsTowerOK(towerID);     // kTRUE means GOOD
  bool TowerDead = mBaseMaker->IsTowerDead(towerID); // kTRUE means BAD
  if (!TowerOK)
  {
    return kFALSE;
  }
  if (TowerDead)
  {
    return kFALSE;
  }

  // jet track acceptance cuts njow
  if ((eta < fJetTowerEtaMin) || (eta > fJetTowerEtaMax))
    return kFALSE;
  if ((phi < fJetTowerPhiMin) || (phi > fJetTowerPhiMax))
    return kFALSE;

  // passed all above cuts - keep tower and fill input vector to fastjet
  return kTRUE;
}

//
// Function: calculate momentum of a tower
//________________________________________________________________________________________________________
Bool_t StHIRecoJets::GetMomentum(TVector3 &mom, const StPicoBTowHit *tower, Double_t mass, Int_t towerID, Double_t CorrectedEnergy) const
{
  // vertex components - only need if below method is used
  // mGeom3->getEtaPhi(towerID,tEta,tPhi);
  // double xVtx = mVertex.x();
  // double yVtx = mVertex.y();
  // double zVtx = mVertex.z();

  // get mass, E, P, ID
  if (mass < 0)
    mass = 0.0;
  // Double_t energy = tower->energy();
  Double_t energy = CorrectedEnergy; // USE THIS, to use the corrected energy to get the momentum components
  Double_t p = 1.0 * TMath::Sqrt(energy * energy - mass * mass);

  // tower ID - passed into function
  // get tower position
  TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
  double posX = towerPosition.x();
  double posY = towerPosition.y();
  double posZ = towerPosition.z();

  // shouldn't need correction with above method
  //  posX-=xVtx;
  //  posY-=yVtx;
  //  posZ-=zVtx;

  // get r, set position components
  Double_t r = TMath::Sqrt(posX * posX + posY * posY + posZ * posZ);
  if (r > 1e-12)
  {
    mom.SetX(p * posX / r);
    mom.SetY(p * posY / r);
    mom.SetZ(p * posZ / r);
  }
  else
  {
    return kFALSE;
  }

  return kTRUE;
}

int StHIRecoJets::GetMatchedBtowID(StPicoTrack *trk)
{
  Double_t bemc_radius = mBemcGeom->Radius();
  // Magnetic field in Tesla
  Double_t mBField_tesla = Bfield / 10.0; // Check this definition. Magnetic fields are minefields of error in STAR

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
  if (mEmcPosition->projTrack(&bemc_pos, &bemc_mom, trk, mVertex, Bfield, bemc_radius))
  {
    // First, examine track eta. If it falls in two regions:
    // 0 < |eta| < etaMin()
    // etaMax() < |eta| < 1.0
    // then shift the eta for the projection slightly into the neighboring tower
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, trkbemcid + 1);
    // cout << bemc_pos.Phi() << "\t" << towerPosition.Phi() << "\t" << bemc_pos.PseudoRapidity() << "\t" << towerPosition.PseudoRapidity() << endl;

    if (fabs(bemc_pos.PseudoRapidity()) < mBemcGeom->EtaMin())
    {
      Double_t unsigned_eta = mBemcGeom->EtaMin() + 0.001;
      // Double_t unsigned_eta = mBemcGeom->EtaMin() + 0.000001;
      Double_t unsigned_theta = 2.0 * atan(exp(-1.0 * unsigned_eta));
      Double_t signed_theta = (bemc_pos.PseudoRapidity() >= 0 ? 1.0 : -1.0) * unsigned_theta;
      bemc_pos.SetTheta(signed_theta);
      close_match = true;
    }
    else if (fabs(bemc_pos.PseudoRapidity()) > mBemcGeom->EtaMax() &&
             fabs(bemc_pos.PseudoRapidity()) < 1.0)
    {
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

    if (mBemcGeom->getBin(bemc_pos.Phi(), bemc_pos.PseudoRapidity(), h_m, h_e, h_s) == 0)
    {
      // If StEmcGeom::getId() == 0: the track was matched successfully. Otherwise,
      // the track was not matched to a tower at this radius, the track was projected
      // into the gap between modules in phi.
      if (h_s != -1)
      {
        mBemcGeom->getId(h_m, h_e, h_s, tow_id);
        if (close_match)
        {
          return -1 * tow_id;
        }

        else
        {
          return tow_id;
        }
      }

      // Track fell in between modules in phi. We will find which module it is closer
      // to by shifting phi slightly.
      else
      {
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

        if (mBemcGeom->getBin(bemc_pos_shift_pos.Phi(), bemc_pos_shift_pos.PseudoRapidity(), h_m, h_e, h_s) == 0 && h_s != -1)
        {
          mBemcGeom->getId(h_m, h_e, h_s, tow_id);
          return -1 * tow_id;
        }

        else if (mBemcGeom->getBin(bemc_pos_shift_neg.Phi(), bemc_pos_shift_neg.PseudoRapidity(), h_m, h_e, h_s) == 0 && h_s != -1)
        {
          mBemcGeom->getId(h_m, h_e, h_s, tow_id);
          return -1 * tow_id;
        }
      }
    }
  }

  return tow_id;
}

Double_t StHIRecoJets::dPhi(Double_t phi1, Double_t phi2)
{
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); // TODO absolute values
  if (deltaPhi > (2 * TMath::Pi()))
    deltaPhi -= 2 * (TMath::Pi());
  if (deltaPhi < (0 * TMath::Pi()))
    deltaPhi += 2 * (TMath::Pi());

  if (deltaPhi > TMath::Pi())
    deltaPhi = 2 * (TMath::Pi()) - deltaPhi;
  return deltaPhi; // dphi in [0, 2Pi]
}