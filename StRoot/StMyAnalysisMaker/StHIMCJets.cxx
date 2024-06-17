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

#include "StHIMCJets.h"

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
#include "StRoot/StPicoEvent/StPicoMcTrack.h"
#include "StRoot/StPicoEvent/StPicoMcVertex.h"
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
#include "StJetFrameworkPicoBase.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistRun14AuAu_P18ih.h" // new Run14 AuAu
#include "runlistRun14AuAu_P16id_SL18f_xrootd_MB.h" // Run14 AuAu used by HF group for MB

// centrality includes
#include "StCentMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// D0 Includes
#include "StHIOverlay.h"

// extra includes
#include "StJetPicoDefinitions.h"

// StRoot classes
class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;

// constants - definitions
const Int_t StHIMCJets::fgkConstIndexShift = 100000;

ClassImp(StHIMCJets)

//________________________________________________________________________
StHIMCJets::StHIMCJets() : 
  StMaker(),
  fD0Analysis(kFALSE),
  doWriteHistos(kFALSE),
  doUsePrimTracks(kFALSE), 
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  doppAnalysis(kFALSE), 
  fRequireCentSelection(kFALSE),
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
  fjw("StHIMCJets", "StHIMCJets"),
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
  fTrackPhiMin(0.0), fTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackEtaMin(-1.0), fJetTrackEtaMax(1.0),
  fJetTrackPhiMin(0.0), fJetTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackDCAcut(3.0),
  fJetTracknHitsFit(15),
  fJetTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fJetTowerEMin(0.2), fJetTowerEMax(100.0),
  fJetTowerEtaMin(-1.0), fJetTowerEtaMax(1.0),
  fJetTowerPhiMin(0.0), fJetTowerPhiMax(2.0*TMath::Pi()),
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
  mCentMaker(0x0),
  mBaseMaker(0x0),
  mHIOverlay(0x0),
  mEmcPosition(0x0),
  grefmultCorr(0x0),
  fEfficiencyInputFile(0x0)
{
  // Default constructor.
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = kFALSE; }

  for(int i=0; i<4800; i++) { 
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE; 

    for(int j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }
  fJetsArr.clear();
  fJetsBGsubArr.clear();
}

//________________________________________________________________________
StHIMCJets::StHIMCJets(const char *name, double mintrackPt = 0.20, bool doHistos = kFALSE, const char* outName = "") : 
  StMaker(name),
  fD0Analysis(kFALSE),
  doWriteHistos(doHistos),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  doppAnalysis(kFALSE),
  fRequireCentSelection(kFALSE),
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
  fTrackPhiMin(0.0), fTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackEtaMin(-1.0), fJetTrackEtaMax(1.0),
  fJetTrackPhiMin(0.0), fJetTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackDCAcut(3.0),
  fJetTracknHitsFit(15),
  fJetTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fJetTowerEMin(0.2), fJetTowerEMax(100.0),
  fJetTowerEtaMin(-1.0), fJetTowerEtaMax(1.0),
  fJetTowerPhiMin(0.0), fJetTowerPhiMax(2.0*TMath::Pi()),
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
  mEmcPosition(0x0),
  grefmultCorr(0x0),
  fEfficiencyInputFile(0x0)
{
  // Standard constructor.
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = kFALSE; }

  for(int i=0; i<4800; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;

    for(int j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

  if (!name) return;
  SetName(name);

  fJetsArr.clear();
  fJetsBGsubArr.clear();
}

//________________________________________________________________________
StHIMCJets::~StHIMCJets()
{
  
  if(mEmcPosition)             delete mEmcPosition;

  // track reconstruction efficiency input file
  if(fEfficiencyInputFile) {
    fEfficiencyInputFile->Close();
    delete fEfficiencyInputFile;
  }
}
//
//
//________________________________________________________________________
Int_t StHIMCJets::Init() {
  // DeclareHistograms();

  // Create user objects.
  fJets = new TClonesArray("StJet");
  fJets->SetName(fJetsName);

  fJetsBGsub = new TClonesArray("StJet");
  fJetsBGsub->SetName(fJetsName+"BGsub");

  fJetsArr.clear();
  fJetsBGsubArr.clear();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // input file - for tracking efficiency: Run14 AuAu and Run12 pp
  const char *input = "";
//if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) input=Form("./StRoot/StMyAnalysisMaker/Run14_AuAu_200_tracking_efficiency_and_momentum_smearing_dca_3p0_nhit_15_nhitfrac_0p52.root");
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) input=Form("./StRoot/StMyAnalysisMaker/Run14_efficiencySmaller2D.root");
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB) input=Form("./StRoot/StMyAnalysisMaker/Run14_efficiencySmaller2D.root");
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)   input=Form("./StRoot/StMyAnalysisMaker/Run12_efficiency_New.root"); // Oct17, 2019 added
  if(fDoEffCorr) {
    fEfficiencyInputFile = new TFile(input, "READ");
    if(!fEfficiencyInputFile) cout<<Form("do not have input file: %s", input);
  }

  // ============================ set jet parameters for fastjet wrapper  =======================
  // recombination schemes:
  // E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, BIpt2_scheme, WTA_pt_scheme, WTA_modp_scheme
  fastjet::RecombinationScheme    recombScheme;
  if (fRecombScheme == 0)     recombScheme = fastjet::E_scheme;
  if (fRecombScheme == 1)     recombScheme = fastjet::pt_scheme; 
  if (fRecombScheme == 2)     recombScheme = fastjet::pt2_scheme;
  if (fRecombScheme == 3)     recombScheme = fastjet::Et_scheme;
  if (fRecombScheme == 4)     recombScheme = fastjet::Et2_scheme;
  if (fRecombScheme == 5)     recombScheme = fastjet::BIpt_scheme;
  if (fRecombScheme == 6)     recombScheme = fastjet::BIpt2_scheme;
  if (fRecombScheme == 7)     recombScheme = fastjet::WTA_pt_scheme;
  if (fRecombScheme == 8)     recombScheme = fastjet::WTA_modp_scheme;
  if (fRecombScheme == 99)    recombScheme = fastjet::external_scheme;

  // jet algorithm
  fastjet::JetAlgorithm          algorithm;
  if (fJetAlgo == 1)      algorithm = fastjet::antikt_algorithm;
  if (fJetAlgo == 0)      algorithm = fastjet::kt_algorithm;
  // extra algorithms
  if (fJetAlgo == 2)      algorithm = fastjet::cambridge_algorithm;
  if (fJetAlgo == 3)      algorithm = fastjet::genkt_algorithm;
  if (fJetAlgo == 11)     algorithm = fastjet::cambridge_for_passive_algorithm;
  if (fJetAlgo == 13)     algorithm = fastjet::genkt_for_passive_algorithm;
  if (fJetAlgo == 99)     algorithm = fastjet::plugin_algorithm;
  if (fJetAlgo == 999)    algorithm = fastjet::undefined_jet_algorithm;
  fastjet::Strategy               strategy = fastjet::Best;

  // cout << "Jet Definition = " << recombScheme << "\t" << algorithm << endl;

  // setup fj wrapper
  fjw.SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw.SetStrategy(strategy);
  fjw.SetGhostArea(fGhostArea);
  fjw.SetR(fRadius);
  fjw.SetAlgorithm(algorithm);        //fJetAlgo);
  fjw.SetRecombScheme(recombScheme);  //fRecombScheme);
  fjw.SetMaxRap(1.2);

  // ghost-area specifications
  double ghost_maxrap = 1.2;
  fastjet::GhostedAreaSpec   area_spec(ghost_maxrap);
  fastjet::AreaDefinition    area_def(fastjet::active_area_explicit_ghosts, area_spec);

  // setting legacy mode
  //if(fLegacyMode) { fjw.SetLegacyMode(kTRUE); }

  // hJetPt = new TH1F("hJetPt", "hJetPt", 40, 0, 40);
  // hJetZ = new TH1F("hJetZ", "hJetZ", 12, 0, 1.2);
  // hD0Pt = new TH1F("hD0Pt", "hD0Pt", 10, 0, 10);

  return kStOK;
}
//
//
//_______________________________________________________________________________________
Int_t StHIMCJets::Finish() {

  // if(mOutName!="") {
  //   TFile *fout = new TFile(mOutName.Data(), "UPDATE");
  //   fout->cd();
  //   fout->mkdir(GetName());
  //   fout->cd(GetName());
  //   hJetPt->Write();
  //   hJetZ->Write();
  //   hD0Pt->Write();
  //   fout->cd();
  //   fout->Write();
  //   fout->Close();
  // }

  return kStOK;
}

//
// Function: clear or delete objects after running
//________________________________________________________________________________
void StHIMCJets::Clear(Option_t *opt) {
  fJets->Clear();
  fJetsBGsub->Clear();
}
//
// Function: main loop, called for each event
//________________________________________________________________________________
int StHIMCJets::Make()
{

  fjw.Clear();
  fFull_Event.clear();

  // cout << "Started with MC Jets" << endl;
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // ZERO's out the jet array
  fJets->Delete();
  fJetsBGsub->Delete();

  fJetsArr.clear();
  fJetsBGsubArr.clear();

  // ZERO these out for double checking they aren't set
  for(int i = 0; i < 4800; i++) {
    for(int j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

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

  // get base class pointer - this class does not inherit from base class: StJetFrameworkPicoBase, but we want to reduce redundancy
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

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  mHIOverlay = static_cast<StHIOverlay*>(GetMaker("HIOverlay"));
  if (!mHIOverlay){
    LOG_WARN << "No HI Overlay! Skip!" << endm;
    return kStFatal;
  }
	
  fMCTracks = mHIOverlay->GetMCTracks();
  fMCTowers = mHIOverlay->GetMCTowers();

  int numberofmcevents = mHIOverlay->GetTheNumberOfEventsToOverLay();

  // cout << "The number of events actually overlaid = " << numberofmcevents << endl;

  if (numberofmcevents == 0) return kStOK;

  for (int i = 0; i < numberofmcevents; i++){
    fJets->Delete();
    fJetsBGsub->Delete();

    // cout << "ITERATION = " << i << endl;
    FindJets(i);
    FillJetBranch(i);
  }

  return kStOK;
}

void StHIMCJets::FindJets(int iteration)
{
  // cout << "This is D0 # " << iteration << "***********************" << endl;

  double pi0mass = Pico::mMass[0]; // GeV
  // clear out existing wrapper object
  fjw.Clear();
  fFull_Event.clear();

  // fjw.PrintJetDescription();

  // Loop over all saved particles for MC in StHIOverlay vectors and enter them into fastjet wrapper 
  for(unsigned short iTracks = 0; iTracks < fMCTracks[iteration].size(); iTracks++){
     
    TLorentzVector v;
    v = fMCTracks[iteration][iTracks];
    // track variables
    
    double px = v.X();
    double py = v.Y();
    double pz = v.Z();
    double p = v.P();
    double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
    double mass = v.M();

    // if (int(mass)==1){  cout << "D0 is in Input loop = " << mass << endl; }

    fjw.AddInputVector(px, py, pz, energy, iTracks ); // includes E
    fastjet::PseudoJet particle_Track(px, py, pz, energy);
    particle_Track.set_user_index(iTracks);
    fFull_Event.push_back(particle_Track);

    double pt = v.Pt();
    double phi = v.Phi();
    if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
    if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
    double eta = v.Eta();

    // cout << "Input = " << pt << "\t" << eta << "\t" << phi << endl;

  } // track loop


  if(fJetType == kFullJet){
    for (unsigned short iTowers = 0; iTowers < fMCTowers[iteration].size(); iTowers++){

      TLorentzVector v;
      v = fMCTowers[iteration][iTowers];
      // tower variables
      
      double px = v.X();
      double py = v.Y();
      double pz = v.Z();
      double p = v.P();
      double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
      double mass = v.M();

      fjw.AddInputVector(px, py, pz, energy, iTowers+10000); // includes E
      fastjet::PseudoJet particle_Track(px, py, pz, energy);
      particle_Track.set_user_index(iTowers+10000);
      fFull_Event.push_back(particle_Track);

      double pt = v.Pt();
      double phi = v.Phi();
      if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
      if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
      double eta = v.Eta();

      // cout << "Input = " << pt << "\t" << eta << "\t" << phi << endl;

    } // tower loop

  }   // if full/charged jets

  // fjw.PrintInput();
  // run jet finder
  fjw.Run();

  // fjw.PrintJetDescription();
}




//
/**
 * This method fills the jet output branch (TClonesArray) with the jet found by the FastJet
 * wrapper. Before filling the jet branch, the utilities are prepared. Then the utilities are
 * called for each jet and finally after jet finding the terminate method of all utilities is called.
 */
void StHIMCJets::FillJetBranch(int iteration)
{

  // cout << "Called FillJetBranch." << endl;
  // get inclusive jets
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();

  // //   // ZERO's out the jet array
  // fJets->Delete();
  // fJetsBGsub->Delete();

  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  // loop over FastJet jets
  __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("%d jets found", (Int_t)jets_incl.size()));
  //for(UInt_t ij = 0, jetCount = 0; ij < jets_incl.size(); ++ij) {

  // cout << jets_incl.size()  << " Jets Found. " << endl;

  // cout << "Cuts = " << fMinJetPt << "\t" << fMinJetArea << "\t" << fJetEtaMin << "\t" << fJetEtaMax << endl;


  for(UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet) {

    Int_t ij = indexes[ijet];
    
    __DEBUG(StJetFrameworkPicoBase::kDebugFillJets,Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)));

    // cout << Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)) << endl;

    // PERFORM CUTS ON inclusive JETS before saving
    // cut on min jet pt
    if(jets_incl[ij].perp() < fMinJetPt) continue;
    // cut on min jet area
    if(fjw.GetJetArea(ij) < fMinJetArea*TMath::Pi()*fRadius*fRadius) continue;
    // cut on eta acceptance
    if((jets_incl[ij].eta() < fJetEtaMin) || (jets_incl[ij].eta() > fJetEtaMax)) continue;
    // cut on phi acceptance 
    if((jets_incl[ij].phi() < fJetPhiMin) || (jets_incl[ij].phi() > fJetPhiMax)) continue;

    // March 8, 2018 - probably don't need this anymore after coming up with method to pass and extract constituents via index!
    // get constituents of jets
    // fConstituents.clear();

    fConstituents = fjw.GetJetConstituents(ij);

    // fill jet constituents
    vector<fastjet::PseudoJet> constituents = fjw.GetJetConstituents(ij);

    bool D0Jet = kFALSE;

    for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
      // get user defined index
      Int_t uid = constituents[ic].user_index();

      if(int(uid) >= 0 && int(uid) < 10000 ){
        TLorentzVector v;
        v = fMCTracks[iteration][uid];

        if (int(v.M())==1) {D0Jet = kTRUE; break;}
      }
    }

    if (!D0Jet) continue;

    // cout << "Found jet with Jet pT = " << jets_incl[ij].perp() << "\t" << jets_incl[ij].eta() << "\t" << jets_incl[ij].phi() << endl;

    // need to figure out how to get m or E from STAR tracks
    StJet *jet = new ((*fJets)[jetCount])
      StJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());

    jet->SetLabel(ij);

    // area vector and components
    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij));
    jet->SetArea(area.perp());  // same as fjw.GetJetArea(ij)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    jet->SetJetConstituents(fConstituents);

    FillJetConstituents(iteration, jet, constituents, constituents);

    __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("Added jet n. %d, pt = %f, area = %f, constituents = %d, tracks = %d, towers = %d", jetCount, jet->Pt(), jet->Area(), jet->GetNumberOfConstituents(), jet->GetNumberOfTracks(), jet->GetNumberOfTowers()));

    jetCount++;
  } // jet loop 

  TClonesArray *tempjets = (TClonesArray *)fJets->Clone();

  fJetsArr.push_back(tempjets);
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

void StHIMCJets::FillJetConstituents(int iteration, StJet *jet, std::vector<fastjet::PseudoJet>& constituents,
    std::vector<fastjet::PseudoJet>& constituents_unsub, Int_t flag, TString particlesSubName)
{
  // initialize some variables/counters
  // Double_t neutralE = 0, maxTrack = 0, maxTower = 0;
  Int_t nt = 0; // track counter
  Int_t nc = 0; // tower (cluster) counter
  Int_t ng = 0; // ghost counter  
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // cout << "Size of constituents = " << constituents.size() << endl;

  // initially set track and tower constituent sizes
  jet->SetNumberOfTracks(constituents.size());
  jet->SetNumberOfTowers(constituents.size());

  // cout << "Jet Number of Tracks = " << jet->GetNumberOfTracks() << endl;
  // cout << "Jet Number of Towers = " << jet->GetNumberOfTowers() << endl;

  // cout << "Jet Pt = " << jet->Pt() << endl;

  for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
    Int_t uid = constituents[ic].user_index();

    if(int(uid) >= 0 && int(uid) < 10000 ){
      // cout << uid  << "\t" << nt << endl;
      jet->AddTrackAt(abs(uid), nt);
      // cout << jet->TrackAt(nt) << endl;
      TLorentzVector v;
      v = fMCTracks[iteration][uid];

      if (int(v.M())==1){  
        // cout << "D0 is in jet = " << v.M() << endl; 
        double mcjetpx = jet->Px();
        double mcjetpy = jet->Py();
        double mcjetpt = jet->Pt();
        double mcD0px = v.Px();
        double mcD0py = v.Py();
        double mcD0pt = v.Pt();

        double z = (mcjetpx*mcD0px + mcjetpy*mcD0py)/(pow(mcjetpt, 2));

        // cout << "Found jet with pT = " << mcjetpt << " D0 pT \t" << mcD0pt << " Z = " << z << endl;
        // hJetPt->Fill(mcjetpt);
        // hJetZ->Fill(z);
        // hD0Pt->Fill(mcD0pt);
      }

      else{
        // cout << "Const = " << v.Pt() << endl;
      }

      nt++;
    }
    else if ( int(uid)>=10000 ){
      // cout << uid  << "\t" << nc << endl;
      jet->AddTowerAt(abs(uid), nc);
      // cout << jet->TowerAt(nc) << endl;
      TLorentzVector v;
      v = fMCTowers[iteration][uid-10000];

      // cout << "Towers = " << v.Pt() << endl;
      nc++;
    }
    else if (uid < 0){
      ng++;
    }
  }

  jet->SetNumberOfTracks(nt);
  jet->SetNumberOfTowers(nc);

}

/**
 * Sorts jets by pT (decreasing)
 * @param[out] indexes This array is used to return the indexes of the jets ordered by pT
 * @param[in] array Vector containing the list of jets obtained by the FastJet wrapper
 * @return kTRUE if at least one jet was found in array; kFALSE otherwise
 */
Bool_t StHIMCJets::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
{
  static Float_t pt[9999] = {0};
  const Int_t n = (Int_t)array.size();
  if(n < 1) return kFALSE;

  for(Int_t i = 0; i < n; i++)
    pt[i] = array[i].perp();

  TMath::Sort(n, pt, indexes);

  return kTRUE;
}

//
// sets errors on histograms up before filling
// set sum weights
//________________________________________________________________________
void StHIMCJets::SetSumw2() {
  TH1::SetDefaultSumw2();
}
