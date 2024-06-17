// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StHIOverlayAngularities.h"
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
#include "StFJWrapper.h"
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

Double_t standardPhi(const Double_t &phi)
{
  Double_t phi_standard = phi;
  if (phi_standard < 0)
    phi_standard += 2 * (TMath::Pi());
  if (phi_standard > 2 * (TMath::Pi()))
    phi_standard -= 2 * (TMath::Pi());
  return phi_standard;
}

ClassImp(StHIOverlayAngularities)

    //________________________________________________________________________
    StHIOverlayAngularities::StHIOverlayAngularities(const char *name, StPicoDstMaker *picoMaker, const char *outName = "", const char *filename = "") : StJetFrameworkPicoBase(name) // StMaker(name),
{
  fLeadingJet = 0x0;
  fSubLeadingJet = 0x0;
  fJets = 0x0;
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
  fRunFlag = 0; // see StHIOverlayAngularities::fRunFlagEnum
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
  fJetAlgo = 1;      // antikt_algorithm
  fJetType = 1;      // kChargedJet
  fRecombScheme = 0; // E_scheme
  fMinJetArea = 0.001;
  fMinJetPt = 0.0;
  fJetPhiMin = 0.;
  fJetPhiMax = 2 * pi;
  fJetEtaMin = -0.6;
  fJetEtaMax = 0.6;
  fGhostArea = 0.01;
  fEventZVtxMinCut = -40.0;
  fEventZVtxMaxCut = 40.0;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fTrackPtMinCut = 0.2;
  fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0;
  fTrackPhiMaxCut = 2.0 * TMath::Pi();
  fTrackEtaMinCut = -1.0;
  fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0;
  fTracknHitsFit = 15;
  fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2;
  fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0;
  fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0;
  fTowerPhiMaxCut = 2.0 * TMath::Pi();
  fCentralityScaled = 0.;
  ref16 = -99;
  ref9 = -99; // FIXME - maybe not make global
  Bfield = 0.0;
  //  mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0;
  fMBEventType = 2;
  fRho = 0x0;
  fRhoVal = 0;
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;

  fMCFileListName = filename;

  for (Int_t i = 0; i < 8; i++)
  {
    fEmcTriggerArr[i] = 0;
  }

  for (Int_t i = 0; i < 4800; i++)
  {
    for (Int_t j = 0; j < 7; j++)
      mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

  mBemcGeom = 0x0;

  fMinJetTrackPt = 0.2;
  fMaxJetTrackPt = 30.0;
  fJetTrackEtaMin = -1.0;
  fJetTrackEtaMax = 1.0;
  fJetTrackPhiMin = 0.0;
  fJetTrackPhiMax = 2.0 * TMath::Pi();
  fJetTrackDCAcut = 3.0;
  fJetTracknHitsFit = 15;
  fJetTracknHitsRatio = 0.52;

  fJetTowerEMin = 0.2;
  fJetTowerEMax = 100.0;
  fJetTowerEtaMin = -1.0;
  fJetTowerEtaMax = 1.0;
  fJetTowerPhiMin = 0.0;
  fJetTowerPhiMax = 2.0 * TMath::Pi();
  mTowerEnergyMin = 0.2;
  mHadronicCorrFrac = 1.;
  fJetHadCorrType = StJetFrameworkPicoBase::kAllMatchedTracks;
  fNumberOfEventsToOverLay = 1;

  for (Int_t i = 0; i < 100; i++)
  {
    fMcEventTracks[i].clear();
    fMcEventTowers[i].clear();
    fRecoMcEventTracks[i].clear();
    fRecoMcEventTowers[i].clear();
    fMcD0Information[i] = {};
    fMcRecoD0Information[i] = {};
    fOrigin[i].SetXYZ(0, 0, 0);
    fMcEventInfo[i] = {};
  }

  fPrintLevel = 0;
  fCentBin = 0;
  fTrackingEfficiency = kFALSE;
  fTrackingEfficiencyPercentage = 1.;
  if (!name)
    return;
  SetName(name);
}

//
//________________________________________________________________________
StHIOverlayAngularities::~StHIOverlayAngularities()
{ /*  */
}

//
//________________________________________________________________________
Int_t StHIOverlayAngularities::Init()
{

  StJetFrameworkPicoBase::Init();

  // position object for Emc
  mBemcGeom = StEmcGeom::instance("bemc");
  mEmcPosition = new StEmcPosition2();

  ifstream filelistforMCEvents(fMCFileListName.Data());

  if (!filelistforMCEvents.is_open())
  {
    LOG_ERROR << "No MC File List! Exiting!" << endm;
    return kStOk;
  }

  string line;

  while (getline(filelistforMCEvents, line))
  {
    TString s(line);
    filenamesforHIOverlay.push_back(s);
  }

  OutputTreeInit();

  // histograms
  //=======================================================================================================//
  hRecoJetPt = new TH1D("hRecoJetPt", "; Jet p_{t}, GeV/c; Counts", 200, 0, 30);
  hMcJetPt = new TH1D("hMcJetPt", "; Jet p_{t}, GeV/c; Counts", 200, 0, 30);
  hMcRecoJetPt = new TH1D("hMcRecoJetPt", "; Jet p_{t}, GeV/c; Counts", 200, 0, 30);
  hMcD0Pt = new TH1D("hMcD0Pt", "; D_{0} p_{t}, GeV/c; Counts", 200, 0, 30);
  hMcRecoD0Pt = new TH1D("hMcRecoD0Pt", "; D_{0} p_{t}, GeV/c; Counts", 200, 0, 30);
  hMcJetZ = new TH1D("hMcJetZ", "; z;  Counts", 200, -1, 2);
  hRecoJetZ = new TH1D("hRecoJetZ", "; z;  Counts", 200, -1, 2);

  hMcD0PtMcRecoD0Pt = new TH2D("hMcD0PtMcRecoD0Pt", ";  Mc D_{0} p_{t}, GeV/c; Reco D_{0}  p_{t}, GeV/c", 200, -5, 30, 200, 0, 30);
  hRecoJetPtMcJetPt = new TH2D("hRecoJetPtMcJetPt", "; HI Reco Jet p_{t}, GeV/c;  Mc Jet p_{t}, GeV/c", 200, -5, 30, 200, 0, 30);
  hRecoJetZMcJetZ = new TH2D("hRecoJetZMcJetZ", "; HI Reco Jet z;  Mc Jet z", 200, -1, 2, 200, -1, 2);

  fractionNeutralToJetPt = new TH1D("fractionNeutralToJetPt", "; Fraction of Neutral to Jet p_{t}; Counts", 200, 0, 1);

  //=======================================================================================================//

  TFile f("/star/u/droy1/Y2019/STAR/Momentum_resolution_SL16d.root");
  fPionMomResolution = (TF1 *)f.Get("fPion")->Clone("fPion");
  fKaonMomResolution = (TF1 *)f.Get("fKaon")->Clone("fKaon");
  fProtonMomResolution = (TF1 *)f.Get("fProton")->Clone("fProton");

  TFile effweight("/star/u/droy1/Y2019/STAR/EffWeightsInCentralityBins.root");

  fPionWeight[0] = (TGraph *)effweight.Get("Pion_0_10");
  fKaonWeight[0] = (TGraph *)effweight.Get("Kaon_0_10");
  fProtonWeight[0] = (TGraph *)effweight.Get("Proton_0_10");
  fAProtonWeight[0] = (TGraph *)effweight.Get("AProton_0_10");

  fPionWeight[1] = (TGraph *)effweight.Get("Pion_10_40");
  fKaonWeight[1] = (TGraph *)effweight.Get("Kaon_10_40");
  fProtonWeight[1] = (TGraph *)effweight.Get("Proton_10_40");
  fAProtonWeight[1] = (TGraph *)effweight.Get("AProton_10_40");

  fPionWeight[2] = (TGraph *)effweight.Get("Pion_40_80");
  fKaonWeight[2] = (TGraph *)effweight.Get("Kaon_40_80");
  fProtonWeight[2] = (TGraph *)effweight.Get("Proton_40_80");
  fAProtonWeight[2] = (TGraph *)effweight.Get("AProton_40_80");

  // ============================ set jet parameters for fastjet wrapper  =======================
  // recombination schemes:
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
  fjw = new StFJWrapper("HIOverlay", "HIOverlay");
  fjw->SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw->SetStrategy(strategy);
  fjw->SetGhostArea(fGhostArea);
  fjw->SetR(fJetRad);
  fjw->SetAlgorithm(algorithm);       // fJetAlgo);
  fjw->SetRecombScheme(recombScheme); // fRecombScheme);
  fjw->SetMaxRap(1.);                // tracks
  fjw->SetMinJetPt(fMinJetPt);

  // // ghost-area specifications
  // Double_t ghost_maxrap = 1.2;
  // fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
  // fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StHIOverlayAngularities::Finish()
{
  cout << "StHIOverlayAngularities::Finish()\n";

  //  Write  to file and close it.
  if (mOutName != "")
  {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    outputTree->Write();
    // write all histograms
    hRecoJetPt->Write();
    fractionNeutralToJetPt->Write();
    hMcJetPt->Write();
    hMcRecoJetPt->Write();
    hMcD0Pt->Write();
    hMcRecoD0Pt->Write();
    hMcJetZ->Write();
    hRecoJetZ->Write();
    hMcD0PtMcRecoD0Pt->Write();
    hRecoJetPtMcJetPt->Write();
    hRecoJetZMcJetZ->Write();

    // fout->Write();
    fout->Close();
  }

  cout << "End of StHIOverlayAngularities::Finish" << endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

// OLD user code says: //  Called every event after Make().
//_____________________________________________________________________________
void StHIOverlayAngularities::Clear(Option_t *opt)
{
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StHIOverlayAngularities::Make()
{
  // ZERO these out for Double_t checking they aren't set
  for (Int_t i = 0; i < 4800; i++)
  {
    for (Int_t j = 0; j < 7; j++)
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
  // get base class pointer
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
  Bfield = mPicoEvent->bField();
  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  Double_t zVtx_VPD = mPicoEvent->vzVpd();
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if ((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut))
    return kStOk;

  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker *>(GetMaker("CentMaker"));
  if (!mCentMaker)
  {
    LOG_WARN << " No CenttMaker! Skip! " << endm;
    return kStWarn;
  }
  // centrality variables
  Int_t grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  Int_t refMult = mCentMaker->GetrefMult();   // see StPicoEvent
  ref9 = mCentMaker->GetRef9();               // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16();             // binning from central -> peripheral
  Int_t cent16 = mCentMaker->GetCent16();     // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  Int_t centbin = mCentMaker->GetRef16();
  Double_t refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  if (fCentralityScaled >= 0 && fCentralityScaled < 10)
    centralitybinforefficiency = 0;
  else if (fCentralityScaled >= 10 && fCentralityScaled < 40)
    centralitybinforefficiency = 1;
  else if (fCentralityScaled >= 40 && fCentralityScaled < 80)
    centralitybinforefficiency = 2;
  // cut on unset centrality, > 80%
  if (cent16 == -1)
    return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them
  // ============================ end of CENTRALITY ============================== //

  if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5)
    return kStOK;
  if ((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut))
    return kStOk;
  Bool_t matchMB = kFALSE;

  set<Int_t> triggers = {450005, 450015, 450025, 450050, 450060, // Run14_AuAu200
                         520001, 520011, 520021, 520031, 520041, 520051, 570002, 570001};

  for (const Int_t &trigger : triggers)
  {
    if (mPicoEvent->isTrigger(trigger))
      matchMB = kTRUE;
    if (matchMB)
      break;
  }

  if (!matchMB)
    return kStOk;
  if (abs(zVtx) > 6.)
    return kStOk;
  if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > 2.)
    return kStOk;
  if (abs(zVtx - zVtx_VPD) > 3)
    return kStOk;
  // ============================ end of VERTEX+TRIGGER ============================== //

  ReadTreeMc();                                // load fMCPico + tree variables
  Int_t nMcEvents = fMCPico->GetEntriesFast(); // get number of events

  Int_t counterEvent = 0;
  Int_t eventsconsidered = 0; // different with counterEvent due to multiple events with D0s in them
  vector<Int_t> randomEventList;

  for (Int_t iD0Event = 0; iD0Event < kMaxNumberOfD0Events; iD0Event++)
  {
    fRecoMcEventTracks[iD0Event].clear();
    fRecoMcEventTowers[iD0Event].clear();
    fMcEventTracks[iD0Event].clear();
    fMcEventTowers[iD0Event].clear();
    fMcD0Information[iD0Event] = {};
    fMcRecoD0Information[iD0Event] = {};
    fOrigin[iD0Event].SetXYZ(0, 0, 0);
    fMcEventInfo[iD0Event] = {};
  }

  while (counterEvent < fNumberOfEventsToOverLay)
  {
    TRandom3 *ran = new TRandom3();
    ran->SetSeed(0);
    Int_t randomMcEvent = ran->Integer(nMcEvents);
    delete ran;
    if (std::find(randomEventList.begin(), randomEventList.end(), randomMcEvent) != randomEventList.end())
      continue;

    randomEventList.push_back(randomMcEvent);
    if (fPrintLevel)
    {
      cout << randomEventList.back() << endl;
      cout << "Size of randomEventList = " << randomEventList.size() << endl;
    }
    if (randomEventList.size() > nMcEvents)
      continue;

    fMCPico->GetEntry(randomMcEvent); // load this MC event into memory
    if (fPrintLevel)
      cout << "======================= Event Number After MC = " << Event_mEventId[0] << "\t" << counterEvent << "=======================" << endl;

    //========================================================================================//
    //================================== Find D0 ids in MC part ==============================//
    //========================================================================================//
    vertexids.clear();
    pionids.clear();
    kaonids.clear();

    matchedpionids.clear();
    matchedkaonids.clear();

    for (Int_t iMcTrack = 0; iMcTrack < McTrack_; iMcTrack++)
    {
      if (McTrack_mGePid[iMcTrack] != 37 && McTrack_mGePid[iMcTrack] != 38)
        continue;
      TLorentzVector particle;
      particle.SetPxPyPzE(McTrack_mPx[iMcTrack], McTrack_mPy[iMcTrack], McTrack_mPz[iMcTrack], McTrack_mE[iMcTrack]);
      // track variables
      Double_t pt = particle.Pt();
      if (pt < 1.0)
        continue; // Only D0s we care about need to have pT > 1 GeV
      Double_t phi = particle.Phi();
      while (phi < 0.0)
        phi += 2.0 * pi; // force from 0-2pi
      while (phi > 2.0 * pi)
        phi -= 2.0 * pi; // force from 0-2pi
      Double_t eta = particle.PseudoRapidity();

      if ((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax))
        continue;

      if (fPrintLevel)
        cout << "MC D0 Found = " << pt << "\t" << eta << "\t" << phi << endl;
      Bool_t isGoodDecay = kTRUE;
      Int_t pionTrack = -99;
      Int_t kaonTrack = -99;

      for (Int_t iMcDaughterTrack = 0; iMcDaughterTrack < McTrack_; iMcDaughterTrack++)
      {
        if (McTrack_mIdVtxStart[iMcDaughterTrack] != McTrack_mIdVtxStop[iMcTrack])
          continue; // We only want the kaon and pion that originated from the D0 we are interested in.

        if (McTrack_mGePid[iMcDaughterTrack] != 8 && McTrack_mGePid[iMcDaughterTrack] != 9 && McTrack_mGePid[iMcDaughterTrack] != 11 && McTrack_mGePid[iMcDaughterTrack] != 12)
        {
          isGoodDecay = kFALSE;
          break; // This should never happen
        }
        // Push the pion and kaon ids into the vector
        if (McTrack_mGePid[iMcDaughterTrack] == 8 || McTrack_mGePid[iMcDaughterTrack] == 9)
        {
          pionTrack = iMcDaughterTrack;
        }
        else if (McTrack_mGePid[iMcDaughterTrack] == 11 || McTrack_mGePid[iMcDaughterTrack] == 12)
        {
          kaonTrack = iMcDaughterTrack;
        }
      }

      if (isGoodDecay)
      {
        if (pionTrack == -99 || kaonTrack == -99)
        {
          cout << "Something wrong with D0. Exiting." << endl;
          continue;
        }

        vertexids.push_back(McTrack_mIdVtxStop[iMcTrack]);
        pionids.push_back(pionTrack);
        kaonids.push_back(kaonTrack);
        matchedpionids.push_back(GetMatchedRecoTrackFromMCTrack(pionTrack));
        matchedkaonids.push_back(GetMatchedRecoTrackFromMCTrack(kaonTrack));
      }
    }
    assert((pionids.size() == kaonids.size()) && "Same number of kaons and pions \n");
    if (fPrintLevel)
      cout << "=============== Number of D0s =============== " << vertexids.size() << "\t" << pionids.size() << endl;

    if (vertexids.size() == 0)
      continue; // While loop continues.

    //========================================================================================//
    //========================================================================================//

    // fill vertices map to loop over further when using GetAllTracksFromVertex()

    fVertexToTracks.clear();

    for (Int_t iMcTrack = 0; iMcTrack < McTrack_; iMcTrack++)
    {
      Int_t idVertexStart = McTrack_mIdVtxStart[iMcTrack] - 1;
      fVertexToTracks[idVertexStart].push_back(iMcTrack);
    }

    //=========================================================================================//
    //================================== Loop over D0 found in event ==========================//
    //=========================================================================================//
    for (UInt_t iD0 = 0; iD0 < vertexids.size(); iD0++)
    {

      if (fPrintLevel)
      {
        cout << "==================================MC part==============================" << endl;
        cout << "=======================================================================" << endl;
      }
      // I only want to use TLorentzVector to propagate the information ahead. All the processing is done within this class.
      // This means I need to find a way to make sure I can identify the D0 track when I save out the information.
      // Since only D0 needs to be identified, I am saving the mass information for the D0 as 1.865.
      // All the other tracks are saved out with pi0mass, because ultimately, jet constituents are assumed to have that mass.

      // How do I propagate charge information though?
      // Do I need it? Mostly, nope! In fact, once I separate track and tower, it should be enough.

      // This loop fills the input vector for the MC Side

      // This function prepares the input list for that event for a particular D0.
      // The event list will of course be a little different for each D0, even for the same event

      fDroppedMCTracks.clear();

      Int_t daughterVertexStop1 = McTrack_mIdVtxStop[pionids[iD0]] - 1;
      Int_t daughterVertexStop2 = McTrack_mIdVtxStop[kaonids[iD0]] - 1;

      GetAllTracksFromVertex(daughterVertexStop1, fDroppedMCTracks);
      GetAllTracksFromVertex(daughterVertexStop2, fDroppedMCTracks);

      if (fPrintLevel)
      {
        cout << "Dropped Track List for D0 # = " << iD0 << endl;
        for (size_t i = 0; i < fDroppedMCTracks.size(); i++)
        {
          TVector3 particle(McTrack_mPx[fDroppedMCTracks[i]], McTrack_mPy[fDroppedMCTracks[i]], McTrack_mPz[fDroppedMCTracks[i]]);
          cout << fDroppedMCTracks[i] << "\t" << McTrack_mGePid[fDroppedMCTracks[i]] << "\t" << particle.Perp() << endl;
        }
        cout << endl;
      }

      //=========================================================================================//
      //================ Loop for preparing MC tracks to fastjet wrapper ========================//
      //=========================================================================================//

      for (Int_t iMcTrack = 0; iMcTrack < McTrack_; iMcTrack++)
      {

        TLorentzVector particle;
        if ((McTrack_mGePid[iMcTrack] == 37 || McTrack_mGePid[iMcTrack] == 38))
          particle.SetXYZM(McTrack_mPx[iMcTrack], McTrack_mPy[iMcTrack], McTrack_mPz[iMcTrack], 1.865); // This is a D0
        else
          particle.SetXYZM(McTrack_mPx[iMcTrack], McTrack_mPy[iMcTrack], McTrack_mPz[iMcTrack], pi0mass);

        // track variables
        Double_t pt = particle.Pt();
        Double_t phi = particle.Phi();
        if (phi < 0.0)
          phi += 2.0 * pi; // force from 0-2pi
        if (phi > 2.0 * pi)
          phi -= 2.0 * pi; // force from 0-2pi
        Double_t eta = particle.PseudoRapidity();

        if (std::find(fDroppedMCTracks.begin(), fDroppedMCTracks.end(), iMcTrack) != fDroppedMCTracks.end())
        {
          if (fPrintLevel)
            cout << "Dropped Track = " << McTrack_mGePid[iMcTrack] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
          continue; // Discard tracks which come from kaon decay after D0 decay
        }

        if (McTrack_mIdVtxStop[iMcTrack] != 0 && (McTrack_mGePid[iMcTrack] != 37 && McTrack_mGePid[iMcTrack] != 38))
          continue; // Only Final State Particles are included in JetMaker. D0s are the only exception

        if (McTrack_mGePid[iMcTrack] == 37 || McTrack_mGePid[iMcTrack] == 38)
        {
          if (McTrack_mIdVtxStop[iMcTrack] != vertexids[iD0])
            continue; // Only consider the current D0
        }
        // Unstable particles which shouldn't make it to the end are discarded by hand. The list provisionally includes:
        /*
          Lambda, Eta, Sigma0, Xi0, Muon, Neutrino, KS0, KL0
        */
        if (McTrack_mGePid[iMcTrack] == 4 || McTrack_mGePid[iMcTrack] == 5 ||
            McTrack_mGePid[iMcTrack] == 6 || McTrack_mGePid[iMcTrack] == 10 ||
            McTrack_mGePid[iMcTrack] == 16 || McTrack_mGePid[iMcTrack] == 17 ||
            McTrack_mGePid[iMcTrack] == 18 || McTrack_mGePid[iMcTrack] == 20 ||
            McTrack_mGePid[iMcTrack] == 22)
          continue;

        // Have to discard Kaons and Pions which come from the Current D0

        if (McTrack_mIdVtxStart[iMcTrack] == vertexids[iD0])
        {
          continue;
        }

        if ((pt < fMinJetTrackPt) || (pt > fMaxJetTrackPt) || (eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax) || (phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax))
          continue;
        if (McTrack_mCharge[iMcTrack] != 0 || McTrack_mGePid[iMcTrack] == 37 || McTrack_mGePid[iMcTrack] == 38)
        {
          fMcEventTracks[counterEvent].push_back(particle); // Neutral particles which are not D0 are saved in towers.
          if (fPrintLevel == 2)
            cout << "Track Constituents = " << McTrack_mGePid[iMcTrack] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
        }
        else
        {
          fMcEventTowers[counterEvent].push_back(particle);
          if (fPrintLevel == 2)
            cout << "Towers = " << McTrack_mGePid[iMcTrack] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
        }
      }
      //=========================================================================================//
      //=========================================================================================//
      fjw->Clear();
      // Loop over all saved particles for MC in vectors and enter them into fastjet wrapper
      for (UInt_t iTrack = 0; iTrack < fMcEventTracks[counterEvent].size(); iTrack++)
      {
        TLorentzVector v = fMcEventTracks[counterEvent][iTrack];
        // track variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
        fjw->AddInputVector(px, py, pz, energy, iTrack + 10000); // includes E
      }                                                          // track loop

      if ((fJetType == kFullJet) || (fJetType == kNeutralJet))
    {
      for (UInt_t iTower = 0; iTower < fMcEventTowers[counterEvent].size(); iTower++)
      {
        TLorentzVector v = fMcEventTowers[counterEvent][iTower];
        // tower variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
        fjw->AddInputVector(px, py, pz, energy, iTower + 20000); // includes E
      }          
      }                                                // tower loop
      // Jet Running
      if (fPrintLevel == 2)
        fjw->PrintInput();

      // This function makes the jets with the MC tracks.
      // If it doesn't find a D0 jet with pt > 5 GeV in the event,
      // we discard the whole event and go to the next case.

      // Make the jets here. If they are not what we want them to be, discard the whole event and keep going.

      StJet *jetMc = DoesItHaveAGoodD0Jet(fMcEventTracks[counterEvent]);
      if (jetMc == NULL)
      {
        fMcEventTracks[counterEvent].clear();
        fMcEventTowers[counterEvent].clear();
        continue;
      }

      fMcJet[counterEvent] = jetMc;

      //=========================================================================================//
      //===================== preparing Reco MC tracks to fastjet wrapper =======================//
      //=========================================================================================//
      if (fPrintLevel)
      {
        cout << "=================================RECO part=============================" << endl;
        cout << "=======================================================================" << endl;
      }
      PrepareSetOfRecoInput(counterEvent, iD0);

      // Loop over all saved particles for MC in StHIOverlayAngularities vectors and enter them into fastjet wrapper
      for (UInt_t iTracks = 0; iTracks < fRecoMcEventTracks[counterEvent].size(); iTracks++)
      {
        TLorentzVector v = fRecoMcEventTracks[counterEvent][iTracks];
        // track variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
        fjw->AddInputVector(px, py, pz, energy, iTracks + 10000); //  10k-20k is Mc tracks range
      }                                                           // track loop

     if ((fJetType == kFullJet) || (fJetType == kNeutralJet))
    {
      for (UInt_t iTowers = 0; iTowers < fRecoMcEventTowers[counterEvent].size(); iTowers++)
      {
        TLorentzVector v = fRecoMcEventTowers[counterEvent][iTowers];
        // tower variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);

        fjw->AddInputVector(px, py, pz, energy, iTowers + 20000); //  //  20k-30k is Mc tower range
      }
      }                                                           // tower loop

      StJet *jetMcReco = DoesItHaveAGoodD0Jet(fRecoMcEventTracks[counterEvent]);
      fMcRecoJet[counterEvent] = jetMcReco;

      if (fPrintLevel == 2)
      {
        cout << "Event Number After Reco = " << Event_mEventId[0] << "\t" << counterEvent << endl;

        for (UInt_t iRecoTrack = 0; iRecoTrack < fRecoMcEventTracks[counterEvent].size(); iRecoTrack++)
        {
          TLorentzVector v = fRecoMcEventTracks[counterEvent][iRecoTrack];
          if ((int)v.M() != 1)
            cout << "Reco Constituents = ";
          else
            cout << "Reco D0 = ";
          cout << v.Pt() << "\t" << v.PseudoRapidity() << "\t" << v.Phi() << "\t" << v.M() << endl;
        }

        for (UInt_t iRecoTower = 0; iRecoTower < fRecoMcEventTowers[counterEvent].size(); iRecoTower++)
        {
          TLorentzVector v = fRecoMcEventTowers[counterEvent][iRecoTower];
          cout << "Reco Tower = " << v.Pt() << "\t" << v.PseudoRapidity() << "\t" << v.Phi() << "\t" << v.M() << endl;
        }
      }
      counterEvent++; // Since one MC event can be invoked multiple times (due to having multiple D0s, this step is necessary.)
    }
    eventsconsidered++; // loop over D0s in event
  }                     // loop over while

  fMCPico->Reset();
  f->Close();
  //========================================================================================//
  //====================================end over MC part====================================//
  //========================================================================================//
  if (fPrintLevel)
  {
    cout << "================================ Real data + embed ====================" << endl;
    cout << "=======================================================================" << endl;
  }
  //=========================================================================================//
  //===================== Loop over Real data and embed PYTHIA tracks =======================//
  //=========================================================================================//

  for (Int_t iMcD0Event = 0; iMcD0Event < counterEvent; iMcD0Event++)
  {
    fFull_Event.clear();

    if (fPrintLevel)
      cout << "Event Number  = " << iMcD0Event << endl;

    fjw->Clear();

    // loop over ALL tracks in PicoDst and add to jet, after acceptance and quality cuts
    if ((fJetType == kFullJet) || (fJetType == kChargedJet))
    {
      for (UInt_t iTrack = 0; iTrack < mPicoDst->numberOfTracks(); iTrack++)
      {
        StPicoTrack *trk = static_cast<StPicoTrack *>(mPicoDst->track(iTrack));
        if (!trk)
        {
          continue;
        }
        // acceptance and kinematic quality cuts - pt cut is also applied here currently
        if (!IsAcceptedTrackAndPt(trk, Bfield, mVertex))
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
        Double_t px = mTrkMom.x();
        Double_t py = mTrkMom.y();
        Double_t pz = mTrkMom.z();
        Double_t p = mTrkMom.Mag();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);

        fjw->AddInputVector(px, py, pz, energy, iTrack); // includes E
        fastjet::PseudoJet particle_Track(px, py, pz, energy);
        particle_Track.set_user_index(iTrack);
        fFull_Event.push_back(particle_Track);
        Int_t matchedTowerIndex = abs(GetMatchedBtowID(trk)) - 1; // towerIndex = towerID - 1
        if (matchedTowerIndex < 0)
          continue;
        mTowerMatchTrkIndex[matchedTowerIndex][mTowerStatusArr[matchedTowerIndex]] = iTrack;
        mTowerStatusArr[matchedTowerIndex] = mTowerStatusArr[matchedTowerIndex] + 1; // 1+ means match, 0 for no match
      }
    }

    // full or neutral jets - get towers and apply hadronic correction
    if ((fJetType == kFullJet) || (fJetType == kNeutralJet))
    {

      fMaxTowerEtBeforeHC = -99;
      fMaxTowerEtAfterHC = -99;

      Int_t matchedTowerTrackCounter = 0;

      // print
      // Int_t nTracks = mPicoDst->numberOfTracks();               // number of tracks
      // Int_t nTrigs = mPicoDst->numberOfEmcTriggers();           // number of Emc triggers
      // Int_t nBTowHits = mPicoDst->numberOfBTowHits();           // barrel tower hits, always 4800
      Int_t nTowers = mPicoDst->numberOfBTowHits(); // barrel tower hists, always 4800
      // Int_t nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits(); // number of BEMC matched tracks
      // cout<<"nTracks = "<<nTracks<<"  nTrigs = "<<nTrigs<<"  nBTowHits = "<<nBTowHits<<"  nBEmcPidTraits = "<<nBEmcPidTraits<<endl;

      // loop over towers and add input vectors to fastjet
      for (Int_t itow = 0; itow < nTowers; itow++)
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
        Int_t towerID = itow + 1;
        Int_t towerIndex = towerID - 1;
        if (towerID < 0)
          continue; // Double_t check these aren't still in the event list

        // tower acceptance cuts - cuts on *bad* towers also
        if (!IsAcceptedTower(tower, towerID))
          continue;

        // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
        TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
        Double_t towerPhi = towerPosition.Phi();
        if (towerPhi < 0.0)
          towerPhi += 2.0 * pi; // force from 0-2pi
        if (towerPhi > 2.0 * pi)
          towerPhi -= 2.0 * pi; // force from 0-2pi
        Double_t towerEta = towerPosition.PseudoRapidity();
        Int_t towerADC = tower->adc();
        Double_t towerEunCorr = tower->energy(); // uncorrected energy
        Double_t towerE = tower->energy();       // corrected energy (hadronically - done below)
        Double_t towEtunCorr = towerEunCorr / (1.0 * TMath::CosH(towerEta));

        // cut on min tower energy after filling histos
        if (towerEunCorr < mTowerEnergyMin)
          continue; // if we don't have enough E to start with, why mess around

        if (towerEunCorr > fMaxTowerEtBeforeHC)
          fMaxTowerEtBeforeHC = towerEunCorr;

        // =======================================================================
        // HADRONIC CORRECTION
        Double_t maxEt = 0.;
        Double_t sumEt = 0.;

        // if tower was is matched to a track or multiple, add up the matched track energies - (mult opt.) to then subtract from the corresponding tower
        // August 15: if *have* 1+ matched trk-tow AND uncorrected energy of tower is at least your tower constituent cut, then CONTINUE
        if (mTowerStatusArr[towerIndex] > 0.5 && towerEunCorr > mTowerEnergyMin)
        {
          Double_t maxE = 0.0;
          Double_t sumE = 0.0;
          // =======================================================================================================================
          // --- finds max E track matched to tower *AND* the sum of all matched track E and subtract from said tower
          //     USER provides readMacro.C which method to use for their analysis via SetJetHadCorrType(type);
          // loop over ALL matched tracks
          for (Int_t itrk = 0; itrk < mTowerStatusArr[towerIndex]; itrk++)
          {
            StPicoTrack *trk = static_cast<StPicoTrack *>(mPicoDst->track(mTowerMatchTrkIndex[towerIndex][itrk]));
            if (!trk)
            {
              cout << "No track pointer..." << endl;
              continue;
            }
            if (!IsAcceptedTrack(trk, Bfield, mVertex))
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
            Double_t p = mTrkMom.Mag();
            Double_t E = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
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
        Double_t fMaxEt = (maxEt == 0) ? towerEunCorr / (1.0 * TMath::CosH(towerEta)) : maxEt;
        Double_t fSumEt = (sumEt == 0) ? towerEunCorr / (1.0 * TMath::CosH(towerEta)) : sumEt;
        // if(mTowerStatusArr[towerIndex] > 0.5) cout<<"towerEunCorr = "<<towerEunCorr<<"  CosH: "<<1.0*TMath::CosH(towerEta)<<"   fMaxEt: "<<fMaxEt<<"   fSumEt: "<<fSumEt<<endl;

        // cut on transverse tower energy (more uniform)
        Double_t towerEt = 0.0;
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
        Double_t towerPx = mom.x();
        Double_t towerPy = mom.y();
        Double_t towerPz = mom.z();

        // add towers to fastjet - shift tower index (tracks 0+, ghosts = -1, towers < -1)
        Int_t uidTow = -(itow + 2);

        fjw->AddInputVector(towerPx, towerPy, towerPz, towerE, uidTow); // includes E
        fastjet::PseudoJet particle_Track(towerPx, towerPy, towerPz, towerE);
        particle_Track.set_user_index(uidTow);
        fFull_Event.push_back(particle_Track);

      } // tower loop

    } // neutral/full jets

    for (UInt_t iMcTrack = 0; iMcTrack < fRecoMcEventTracks[iMcD0Event].size(); iMcTrack++)
    {

      TLorentzVector v = fRecoMcEventTracks[iMcD0Event][iMcTrack];
      // track variables
      Double_t px = v.X();
      Double_t py = v.Y();
      Double_t pz = v.Z();
      Double_t p = v.P();
      Double_t energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);

      if ((int)v.M() == 1 && fPrintLevel)
        cout << "Reco D0 In Jet Loop = " << v.Pt() << "\t" << v.PseudoRapidity() << endl;

      // MC Tracks will start from 10000
      fjw->AddInputVector(px, py, pz, energy, iMcTrack + 10000); // includes E
    }                                                            // track loop

    if (fJetType == kFullJet || (fJetType == kNeutralJet))
    {
      for (UInt_t iTowers = 0; iTowers < fRecoMcEventTowers[iMcD0Event].size(); iTowers++)
      {
        TLorentzVector v;
        v = fRecoMcEventTowers[iMcD0Event][iTowers];
        // tower variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);

        // MC Towers will start from
        fjw->AddInputVector(px, py, pz, energy, -1 * iTowers - 10000); // includes E

      } // tower loop

    } // if full/charged jets

    StJet *jetReco = DoesItHaveAGoodD0Jet(fRecoMcEventTracks[iMcD0Event]); // run fjw and get the jets

    if (jetReco != NULL)
    {
      if (fPrintLevel)
      {
        cout << "==============================background estimate====================" << endl;
        cout << "=======================================================================" << endl;
      }

      fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm, fJetRad);
      fastjet::Selector eta_range_sel = fastjet::SelectorAbsEtaMax(1.0); // 0.6

      Int_t nRejectedJets = 1;

      if (fCentralityScaled >= 0 && fCentralityScaled < 10)
        nRejectedJets = 2;

      fastjet::Selector hard_jet_sel = (!fastjet::SelectorNHardest(nRejectedJets));
      fastjet::Selector pt_min_selector = fastjet::SelectorPtMin(0.01);
      fastjet::Selector full_selector = hard_jet_sel * eta_range_sel * pt_min_selector;

      // Double_t ghost_maxrap = 1.2;
      Double_t ghost_maxrap = 1.0;
      fGhostArea = 0.01;
      // fGhostArea=0.005;
      fastjet::GhostedAreaSpec area_spec(ghost_maxrap, 1, fGhostArea);
      fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);

      fastjet::JetMedianBackgroundEstimator bge_rho(full_selector, jet_def_for_rho, area_def);
      bge_rho.set_particles(fFull_Event);
      Double_t rho = bge_rho.rho();
      fRhoVal = rho;
      Double_t jetcorrptfromarea = jetReco->Pt() - rho * jetReco->Area();

      // cout << "Event Size fjw= " << fjw->GetInputVectors().size() << "\t fFullEvent = " << fFull_Event.size() << endl;
      if (fPrintLevel)
      {
        cout << "Jet Pt = " << jetReco->Pt() << "\t rho = " << rho << "\t area = " << jetReco->Area() << "\t corr pt = " << jetcorrptfromarea << endl;
      }
      jetReco->SetJetPt(jetcorrptfromarea);
    }
    fRecoJet[iMcD0Event] = jetReco;
  } // loop over D0 from MC

  FillTree(counterEvent);
  return kStOK;

} // Make

//
// Function: track quality cuts
//________________________________________________________________________
Bool_t StHIOverlayAngularities::IsAcceptedTrack(StPicoTrack *trk, const Float_t &B, const TVector3 &Vert)
{
  Double_t pi = 1.0 * TMath::Pi();
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
  Double_t pt = mTrkMom.Perp();
  Double_t phi = mTrkMom.Phi();
  Double_t eta = mTrkMom.PseudoRapidity();
  Double_t dca = trk->gDCA(Vert).Mag();
  Int_t nHitsFit = trk->nHitsFit();
  Int_t nHitsMax = trk->nHitsMax();
  Double_t nHitsRatio = 1.0 * nHitsFit / nHitsMax;
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
// Function: jet track quality cuts
// - this function should only be used for jet constituent tracks
// 	NOT when considering track-tower matches  (Aug 19')

Bool_t StHIOverlayAngularities::IsAcceptedTrackAndPt(StPicoTrack *trk, const Float_t &B, const TVector3 &Vert)
{
  if (!IsAcceptedTrack(trk, B, Vert))
    return kFALSE; // first do general track QA cuts
  // get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  if (doUsePrimTracks)
    mTrkMom = trk->pMom(); // get primary track vector
  else
    mTrkMom = trk->gMom(Vert, B); // get global track vector
  Double_t pt = mTrkMom.Perp();
  if (pt < fMinJetTrackPt || pt > fMaxJetTrackPt)
    return kFALSE;
  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}
//
// Tower Quality Cuts
//________________________________________________________________________
Bool_t StHIOverlayAngularities::IsAcceptedTower(StPicoBTowHit *tower, const Int_t &towerID)
{ // constants:
  Double_t pi = 1.0 * TMath::Pi();
  // tower ID - passed into function: make sure some of these aren't still in event array
  if (towerID < 0)
    return kFALSE;
  // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
  TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
  Double_t phi = towerPosition.Phi();
  if (phi < 0.0)
    phi += 2.0 * pi; // force from 0-2pi
  if (phi > 2.0 * pi)
    phi -= 2.0 * pi; // force from 0-2pi
  Double_t eta = towerPosition.PseudoRapidity();
  // check for bad (and dead) towers
  Bool_t TowerOK = mBaseMaker->IsTowerOK(towerID);     // kTRUE means GOOD
  Bool_t TowerDead = mBaseMaker->IsTowerDead(towerID); // kTRUE means BAD
  if (!TowerOK || TowerDead)
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

void StHIOverlayAngularities::PrepareSetOfRecoInput(const Int_t &counterEvent, const Int_t &iD0)
{
  const Int_t numberoftowers = BTowHit_;
  Double_t towerenergy[numberoftowers];
  for (Int_t i = 0; i < numberoftowers; i++)
    towerenergy[i] = 0.;

  vector<Int_t> mTowerToTrack[4800];
  vector<Double_t> mTowerToTrackE[4800];

  TVector3 oVertex;
  oVertex.SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

  if (fPrintLevel == 2)
    cout << "HIOverlay Vertex " << Event_mPrimaryVertexX[0] << "\t" << Event_mPrimaryVertexY[0] << "\t" << Event_mPrimaryVertexZ[0] << endl;

  // This loop fills the input vector for the RECO side
  for (Int_t reco = 0; reco < Track_; reco++)
  {

    TVector3 o;
    o.SetXYZ(Track_mOriginX[reco], Track_mOriginY[reco], Track_mOriginZ[reco]);
    TVector3 g;
    g.SetXYZ(Track_mGMomentumX[reco], Track_mGMomentumY[reco], Track_mGMomentumZ[reco]);

    // track variables
    Double_t pt = g.Perp();
    Double_t phi = g.Phi();
    if (phi < 0.0)
      phi += 2.0 * pi; // force from 0-2pi
    if (phi > 2.0 * pi)
      phi -= 2.0 * pi; // force from 0-2pi
    Double_t eta = g.PseudoRapidity();
    Double_t px = g.x();
    Double_t py = g.y();
    Double_t pz = g.z();
    Double_t p = g.Mag();
    Double_t energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
    short charge = (Track_mNHitsFit[reco] > 0) ? 1 : -1;
    Double_t dca = (oVertex - o).Mag();

    Bool_t goodtrack = (dca < fJetTrackDCAcut) && (abs(Track_mNHitsFit[reco]) >= fTracknHitsFit) && (abs(Double_t(Track_mNHitsFit[reco])) / Double_t(Track_mNHitsMax[reco]) >= fTracknHitsRatio);

    //// Variables For FastSim

    Double_t pt_new = pt;
    Double_t phi_new = phi;
    Double_t eta_new = eta;
    Double_t px_new = px;
    Double_t py_new = py;
    Double_t pz_new = pz;
    Double_t p_new = p;
    Double_t energy_new = energy;
    short charge_new = charge;

    Bool_t mctrackavailable = kTRUE;
    Int_t mcid = Track_mIdTruth[reco] - 1;
    if (mcid < 0)
      mctrackavailable = kFALSE;

    Bool_t isatrackfromD0 = kFALSE;

    // Here, we have two paths to take. If the track needs replacement, we replace it with the fastsim method that is standardised.
    // Else, the pt, eta, phi are sent as is to the final vector.

    if (mctrackavailable && goodtrack)
    {
      if (reco == matchedpionids[iD0] || reco == matchedkaonids[iD0])
        isatrackfromD0 = kTRUE; // Kaons and Pions that come from the current D0 need to be tossed, and replaced by the fast sim version

      TVector3 mg(McTrack_mPx[mcid], McTrack_mPy[mcid], McTrack_mPz[mcid]);
      Double_t relativesmearing = TMath::Sqrt(pow(mg.Px() - px, 2) + pow(mg.Py() - py, 2)) / (mg.Pt());
      Double_t fastsimsmearing;
      Int_t pid = McTrack_mGePid[mcid];

      if (pid == 8 || pid == 9)
        fastsimsmearing = fPionMomResolution->Eval(mg.Pt()); // Pion
      else if (pid == 11 || pid == 12)
        fastsimsmearing = fKaonMomResolution->Eval(mg.Pt()); // Kaon
      else if (pid == 15 || pid == 14)
        fastsimsmearing = fProtonMomResolution->Eval(mg.Pt()); // Proton
      else
        fastsimsmearing = fPionMomResolution->Eval(mg.Pt()); // Catch all: pions

      if (relativesmearing > 3 * fastsimsmearing)
      {
        TVector3 fastsimsmearedmom = FastSimMom(mg, pid);
        pt_new = fastsimsmearedmom.Perp();
        phi_new = fastsimsmearedmom.Phi();
        if (phi_new < 0.0)
          phi_new += 2.0 * pi; // force from 0-2pi
        if (phi_new > 2.0 * pi)
          phi_new -= 2.0 * pi; // force from 0-2pi
        eta_new = fastsimsmearedmom.PseudoRapidity();
        px_new = fastsimsmearedmom.x();
        py_new = fastsimsmearedmom.y();
        pz_new = fastsimsmearedmom.z();
        p_new = fastsimsmearedmom.Mag();
        energy_new = 1.0 * TMath::Sqrt(p_new * p_new + pi0mass * pi0mass);
      }
    }
    if (fPrintLevel == 2)
    {
      if (mctrackavailable)
        cout << Form("Track # %i \t %i \t %i \t %.2f \t %.2f \t %i \t %i \t %.2f \t %i \t %.2f", reco, mcid, McTrack_mGePid[mcid], pt_new, eta_new, isatrackfromD0, abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1, dca, abs(int(Track_mNHitsFit[reco])), abs(Double_t(Track_mNHitsFit[reco])) / Double_t(Track_mNHitsMax[reco])) << endl;
      else
        cout << Form("Track # %i \t %.2f \t %.2f \t %i \t %.2f \t %.2f", reco, pt_new, eta_new, abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1, dca, abs(Double_t(Track_mNHitsFit[reco])) / Double_t(Track_mNHitsMax[reco]), Track_mQATruth[reco]) << endl;
    }

    if (!goodtrack)
      continue;
    // DCA based cuts precede everything else.

    // Track from D0 -> K Pi || D0 is in acceptance range. The KPi do not need to be in acceptance. || The KPi track is projected onto the towers.
    if (isatrackfromD0)
    {
      Int_t matchedTowerIndex = abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1;
      if (matchedTowerIndex >= 0)
      {
        towerenergy[matchedTowerIndex] += energy;
        mTowerToTrack[matchedTowerIndex].push_back(reco);
        mTowerToTrackE[matchedTowerIndex].push_back(energy);
      }
    }

    Int_t particleid = -99;

    Double_t nsigpion = Track_mNSigmaPion[reco] / 1000.;
    Double_t nsigkaon = Track_mNSigmaKaon[reco] / 1000.;
    Double_t nsigproton = Track_mNSigmaProton[reco] / 1000.;

    if (abs(nsigpion) < 2 && abs(nsigkaon) > 2. && abs(nsigproton) > 2.)
      particleid = 1;
    else if (abs(nsigpion) > 2 && abs(nsigkaon) < 2. && abs(nsigproton) > 2.)
      particleid = 2;
    else if (abs(nsigpion) > 2 && abs(nsigkaon) > 2. && abs(nsigproton) < 2.)
      particleid = 3 * charge;

    Bool_t removetrack = kFALSE;
    Bool_t isD0DaugDescendant = kFALSE;
    if (!isatrackfromD0 && pt > fMinJetTrackPt && pt < fMaxJetTrackPt && (eta > fJetTrackEtaMin) && (eta < fJetTrackEtaMax))
    // I am using the old pt for this because the efficiencies were derived with the old pt.
    // I am keeping things consistent with MC
    {
      // KeepTrack returns False if the track is to be discarded // RemoveTrack = True if KeepTrack is False
      removetrack = !KeepTrack(particleid, centralitybinforefficiency, pt);
      // If D0 descendant, record that as well.
      isD0DaugDescendant = kFALSE;
      Int_t mctrkid = Track_mIdTruth[reco] - 1;
      if (std::find(fDroppedMCTracks.begin(), fDroppedMCTracks.end(), mctrkid) != fDroppedMCTracks.end())
        isD0DaugDescendant = kTRUE;
    }

    if (removetrack || isD0DaugDescendant)
    {
      if (fPrintLevel == 2)
      {
        if (mctrackavailable)
          cout << Form("Removed Track Input = %.2f \t %.2f \t %.2f \t %i", px_new, py_new, pz_new, McTrack_mGePid[mcid]) << endl;
        else
          cout << Form("Removed Track Input = %.2f \t %.2f \t %.2f", px_new, py_new, pz_new) << endl;
      }
    }

    // jet track acceptance cuts now
    if (pt_new < fMinJetTrackPt || pt_new > fMaxJetTrackPt)
      continue; // 20.0 STAR, 100.0 ALICE
    if ((eta_new < fJetTrackEtaMin) || (eta_new > fJetTrackEtaMax))
      continue;
    while (phi_new < 0.0)
      phi_new += 2.0 * pi; // force from 0-2pi
    while (phi_new > 2.0 * pi)
      phi_new -= 2.0 * pi; // force from 0-2pi
    if ((phi_new < fJetTrackPhiMin) || (phi_new > fJetTrackPhiMax))
      continue;

    // additional quality cuts for tracks

    // This is questionable. Here, I don't think I should use it. But, in the QM method, I should?
    // I have now included these cuts in both. I think that's the correct thing to do.

    // This place takes care of all energy depositions due to tracks accepted for jet reco.
    // For the reco tracks that we discard, we still need to subtract their contribution from the tower
    // It also includes the kaon and pion from D0.

    Bool_t ignoretrack = removetrack || isatrackfromD0 || isD0DaugDescendant;
    if (ignoretrack)
      continue; // To match the efficiency, we start tossing random tracks.

    Int_t matchedTowerIndex = abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1;
    if (matchedTowerIndex >= 0)
    {
      towerenergy[matchedTowerIndex] += energy;
      mTowerToTrack[matchedTowerIndex].push_back(reco);
      mTowerToTrackE[matchedTowerIndex].push_back(energy);
    }
    TLorentzVector v;
    v.SetXYZM(px_new, py_new, pz_new, pi0mass);
    fRecoMcEventTracks[counterEvent].push_back(v);
    //// Fill the Tower Array with Energy Depositions here
    if (fPrintLevel == 2)
    {
      if (mctrackavailable)
        cout << Form("Reco Level Track Input = %.2f \t %.2f \t %.2f \t %i", px_new, py_new, pz_new, McTrack_mGePid[mcid]) << endl;
      else
        cout << Form("Reco Level Track Input = %.2f \t %.2f \t %.2f", px_new, py_new, pz_new) << endl;
    }
  }

  TVector3 mcKaon(McTrack_mPx[kaonids[iD0]], McTrack_mPy[kaonids[iD0]], McTrack_mPz[kaonids[iD0]]);
  TVector3 mcPion(McTrack_mPx[pionids[iD0]], McTrack_mPy[pionids[iD0]], McTrack_mPz[pionids[iD0]]);

  fMcD0Information[counterEvent] = {mcPion, mcKaon};

  TVector3 recoKaon = FastSimMom(mcKaon, 11);
  TVector3 recoPion = FastSimMom(mcPion, 8);

  fMcRecoD0Information[counterEvent] = {recoPion, recoKaon};

  fOrigin[counterEvent].SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

  TVector3 recoD0;
  recoD0 = recoKaon + recoPion;

  TLorentzVector v;
  v.SetXYZM(recoD0.x(), recoD0.y(), recoD0.z(), 1.865); // This is a D0

  if (fPrintLevel == 2)
    cout << "Reco D0 Entered Jet Loop = " << v.Pt() << "\t" << v.PseudoRapidity() << endl;

  if (fPrintLevel == 2)
    cout << Form("Reco Level Track Input D0 = %.2f \t %.2f \t %.2f", recoD0.x(), recoD0.y(), recoD0.z()) << endl;

  fRecoMcEventTracks[counterEvent].push_back(v);

  // This loop fills the input vector for the towers

  for (Int_t tower = 0; tower < numberoftowers; tower++)
  {
    Int_t towerID = tower + 1;
    if (towerID < 0)
      continue; // Double_t check these aren't still in the event list

    TVector3 towerPosition = mEmcPosition->getPosFromVertex(oVertex, towerID);
    Double_t towerPhi = towerPosition.Phi();
    if (towerPhi < 0.0)
      towerPhi += 2.0 * pi; // force from 0-2pi
    if (towerPhi > 2.0 * pi)
      towerPhi -= 2.0 * pi; // force from 0-2pi
    Double_t towerEta = towerPosition.PseudoRapidity();

    // check for bad (and dead) towers
    Bool_t TowerOK = mBaseMaker->IsTowerOK(towerID);     // kTRUE means GOOD
    Bool_t TowerDead = mBaseMaker->IsTowerDead(towerID); // kTRUE means BAD
    if (!TowerOK)
    {
      continue;
    }
    if (TowerDead)
    {
      continue;
    }

    // jet track acceptance cuts njow
    if ((towerEta < fJetTowerEtaMin) || (towerEta > fJetTowerEtaMax))
      continue;
    if ((towerPhi < fJetTowerPhiMin) || (towerPhi > fJetTowerPhiMax))
      continue;

    Double_t towerEunCorr = Double_t(BTowHit_mE[tower]) / 1000.; // uncorrected energy
    Double_t towerE = Double_t(BTowHit_mE[tower]) / 1000.;       // corrected energy (hadronically - done below)
    Double_t towEtunCorr = towerE / (1.0 * TMath::CosH(towerEta));

    if (fPrintLevel == 2)
      cout << "Tower = " << tower << "\t" << towerE << "\t" << towEtunCorr << endl;

    // cut on min tower energy after filling histos
    if (towerEunCorr < mTowerEnergyMin)
      continue; // if we don't have enough E to start with, why mess around

    // =======================================================================
    // HADRONIC CORRECTION

    Double_t sumEt = (towerEunCorr - towerenergy[tower]) / (1.0 * TMath::CosH(towerEta));
    Double_t towerEt = sumEt;

    if (towerEt < mTowerEnergyMin)
      continue;
    towerE = towerEt * 1.0 * TMath::CosH(towerEta);

    if (fPrintLevel == 2)
      cout << "Tower = " << tower << "\t" << towerE << "\t" << sumEt << endl;

    Double_t p = 1.0 * TMath::Sqrt(towerE * towerE - pi0mass * pi0mass);

    Double_t posX = towerPosition.x();
    Double_t posY = towerPosition.y();
    Double_t posZ = towerPosition.z();

    Double_t r = TMath::Sqrt(posX * posX + posY * posY + posZ * posZ);

    TLorentzVector v;
    v.SetXYZM(p * posX / r, p * posY / r, p * posZ / r, pi0mass);

    fRecoMcEventTowers[counterEvent].push_back(v);

    if (fPrintLevel == 2)
      cout << Form("Reco Level Tower Input = %i \t %.2f \t %.2f \t %.2f \t %.2f", towerID, v.X(), v.Y(), v.Z(), towerE) << endl;
  }
}

//
// Function: calculate momentum of a tower
//________________________________________________________________________________________________________
Bool_t StHIOverlayAngularities::GetMomentum(TVector3 &mom, const StPicoBTowHit *tower, Double_t mass, Int_t towerID, Double_t CorrectedEnergy) const
{
  // get mass, E, P, ID
  if (mass < 0)
    mass = 0.0;
  Double_t energy = CorrectedEnergy; // USE THIS, to use the corrected energy to get the momentum components
  Double_t p = 1.0 * TMath::Sqrt(energy * energy - mass * mass);
  // get tower position
  TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
  Double_t posX = towerPosition.x();
  Double_t posY = towerPosition.y();
  Double_t posZ = towerPosition.z();
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

Int_t StHIOverlayAngularities::GetMatchedBtowID(StPicoTrack *trk)
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
  Int_t trkbemcid = trk->bemcTowerIndex();
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

Bool_t StHIOverlayAngularities::KeepTrack(const Int_t &particleid, const Int_t &centralitybin, const Double_t &pt)
{
  Bool_t keeptrack = kTRUE;
  TRandom3 *r = new TRandom3(0);
  r->SetSeed(0);
  Double_t rando = r->Rndm();
  if (particleid == 1 || particleid == -1)
    keeptrack = (rando > fPionWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE; // Either a pion
  else if (particleid == 2 || particleid == -2)
    keeptrack = (rando > fKaonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == 3)
    keeptrack = (rando > fProtonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == -3)
    keeptrack = (rando > fAProtonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  if (keeptrack)
  {
    if (fTrackingEfficiency)
    {
      keeptrack = (rando > fTrackingEfficiencyPercentage) ? kFALSE : kTRUE;
    }
  }
  delete r;
  return keeptrack;
}

StJet *StHIOverlayAngularities::DoesItHaveAGoodD0Jet(vector<TLorentzVector> &eventTracks)
{
  fjw->Run();
  std::vector<fastjet::PseudoJet> jets_incl = fjw->GetInclusiveJets();
  // sort jets according to jet pt
  Int_t nJets = jets_incl.size();
  Int_t indexes[9999];
  Float_t pt[9999] = {0};
  for (Int_t i = 0; i < nJets; i++)
    pt[i] = jets_incl[i].perp();
  TMath::Sort(nJets, pt, indexes);
  Int_t D0JetIndex = -99;
  for (Int_t iJetUnsorted = 0; iJetUnsorted < nJets; ++iJetUnsorted)
  {

    Int_t iJet = indexes[iJetUnsorted];
    // PERFORM CUTS ON inclusive JETS before saving
    // cut on min jet pt
    if (jets_incl[iJet].perp() < fMinJetPt)
      continue;
    // cut on min jet area
    if (fjw->GetJetArea(iJet) < fMinJetArea * TMath::Pi() * fJetRad * fJetRad)
      continue;
    // cut on eta acceptance
    if ((jets_incl[iJet].eta() < -1. + fJetRad) || (jets_incl[iJet].eta() > 1. - fJetRad))
      continue;
    // // cut on phi acceptance
    if ((jets_incl[iJet].phi() < fJetPhiMin) || (jets_incl[iJet].phi() > fJetPhiMax))
      continue;

    // fill jet constituents
    vector<fastjet::PseudoJet> constituents = fjw->GetJetConstituents(iJet);

    Bool_t D0Jet = kFALSE;

    for (UInt_t iConstituent = 0; iConstituent < constituents.size(); ++iConstituent)
    {
      // get user defined index
      Int_t uid = constituents[iConstituent].user_index();

      if (int(uid) >= 10000 && int(uid) < 20000) // this is a charged MC track id range
      {
        TLorentzVector v = eventTracks[uid - 10000];
        if ((int)v.M() == 1)
        {
          D0Jet = kTRUE;
          D0JetIndex = iJet;
          Double_t D0Pt = v.Pt();
          Double_t D0Eta = v.PseudoRapidity();
          Double_t D0Phi = v.Phi();
          if (fPrintLevel)
            cout << "D0 Found with pT eta phi = " << D0Pt << "\t" << D0Eta << "\t" << D0Phi << endl;
          break;
        }
      }
    }

    if (!D0Jet)
      continue;

    Double_t jet_pt = jets_incl[iJet].perp();
    Double_t neutral_pt = 0.;
    for (UInt_t iConstituent = 0; iConstituent < constituents.size(); ++iConstituent)
    {
      // get user defined index
      Int_t uid = constituents[iConstituent].user_index();
      if (int(uid) < -1 || int(uid) >= 20000) // this is a reco tower
      {
        neutral_pt += constituents[iConstituent].perp();
      }
    }

    fractionNeutralToJetPt->Fill(neutral_pt / jet_pt);

    if ((neutral_pt / jet_pt) > 0.95)
      continue;

    StJet *jet = new StJet(jets_incl[D0JetIndex].perp(), jets_incl[D0JetIndex].eta(), jets_incl[D0JetIndex].phi(), jets_incl[D0JetIndex].m());
    if (fPrintLevel)
      cout << "Jet Found with pT eta phi = " << jet->Pt() << "\t" << jet->Eta() << "\t" << jet->Phi() << endl;

    jet->SetLabel(D0JetIndex);
    // area vector and components
    fastjet::PseudoJet area(fjw->GetJetAreaVector(D0JetIndex));
    jet->SetArea(area.perp()); // same as fjw->GetJetArea(iJet)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());
    constituents = fjw->GetJetConstituents(D0JetIndex);
    jet->SetJetConstituents(constituents);

    fjw->Clear();
    return jet;
  }
  if (fPrintLevel)
  {
    cout << "No D0 Jet Found" << endl;
  }
  fjw->Clear();
  return NULL;
}

void StHIOverlayAngularities::GetAllTracksFromVertex(const Int_t &vertexid, vector<Int_t> &tracksFromVertex)
{
  if (vertexid < 0)
    return;
  if (fPrintLevel == 2)
  {
    cout << "Called this function for vx " << vertexid << " with ntracks = " << fVertexToTracks[vertexid].size() << endl;
    for (UInt_t track = 0; track < fVertexToTracks[vertexid].size(); track++)
    {
      cout << "Geant ID of tracks = " << McTrack_mGePid[fVertexToTracks[vertexid][track]] << endl;
    }
  }
  for (UInt_t track = 0; track < fVertexToTracks[vertexid].size(); track++)
  {
    tracksFromVertex.push_back(fVertexToTracks[vertexid][track]);
    Int_t idvxstop = McTrack_mIdVtxStop[fVertexToTracks[vertexid][track]];
    GetAllTracksFromVertex(idvxstop - 1, tracksFromVertex);
  }
  return;
}

Int_t StHIOverlayAngularities::GetMatchedRecoTrackFromMCTrack(const Int_t &McTrack)
{
  Int_t recotrackmatch = -999;

  TVector3 eventVertex(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);
  Double_t ratio = -99.;
  for (Int_t recoTrack = 0; recoTrack < Track_; recoTrack++)
  {
    // get track pointer
    TVector3 originVertex(Track_mOriginX[recoTrack], Track_mOriginY[recoTrack], Track_mOriginZ[recoTrack]);

    Double_t dca = (eventVertex - originVertex).Mag();
    Int_t nHitsFit = abs(Track_mNHitsFit[recoTrack]);
    Int_t nHitsMax = abs(Track_mNHitsMax[recoTrack]);
    Double_t nHitsRatio = 1.0 * nHitsFit / nHitsMax;

    // additional quality cuts for tracks
    if (dca > fJetTrackDCAcut)
      continue;
    if (nHitsFit < fTracknHitsFit)
      continue;
    if (nHitsRatio < fTracknHitsRatio)
      continue;

    Int_t McTrackFromReco = Track_mIdTruth[recoTrack] - 1;

    if (McTrackFromReco == McTrack)
    {
      if (nHitsRatio > ratio)
      {
        recotrackmatch = recoTrack;
        ratio = nHitsRatio;
      }
    }
  }
  return recotrackmatch;
}

TVector3 StHIOverlayAngularities::FastSimMom(TVector3 p, Int_t pid)
{
  float pt = p.Perp();
  float pt1 = pt;
  if (pt1 > 10)
    pt1 = 10; // Used for high pt-hat bin smearing test
  float sPt = -1;

  TRandom3 *r = new TRandom3();
  r->SetSeed(0);

  if (pid == 8 || pid == 9)
    sPt = r->Gaus(pt, pt * fPionMomResolution->Eval(pt1)); // Pion
  else if (pid == 11 || pid == 12)
    sPt = r->Gaus(pt, pt * fKaonMomResolution->Eval(pt1)); // Kaon
  else if (pid == 15 || pid == 14)
    sPt = r->Gaus(pt, pt * fProtonMomResolution->Eval(pt1)); // Proton
  else
    sPt = r->Gaus(pt, pt * fPionMomResolution->Eval(pt1)); // Catch all: pions

  TVector3 smearedmom(sPt * cos(p.Phi()), sPt * sin(p.Phi()), sPt * sinh(p.PseudoRapidity()));
  delete r;
  return smearedmom;
}

void StHIOverlayAngularities::ReadTreeMc()
{

  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);
  Int_t filenumber = r1->Integer(filenamesforHIOverlay.size());
  delete r1;
  f = new TFile(filenamesforHIOverlay.at(filenumber).Data());
  fMCPico = (TTree *)f->Get("PicoDst");
  if (fPrintLevel)
    cout << filenumber << "\t" << filenamesforHIOverlay.at(filenumber) << "\t" << fMCPico->GetEntriesFast() << "\t" << mCentMaker->GetWeight() << endl;

  fMCPico->SetMakeClass(1);
  fMCPico->SetBranchStatus("*", false);
  fMCPico->SetBranchStatus("Event", true);
  fMCPico->SetBranchStatus("Event.mRunId", true);
  fMCPico->SetBranchStatus("Event.mEventId", true);

  fMCPico->SetBranchStatus("Event.mPrimaryVertexX", true);
  fMCPico->SetBranchStatus("Event.mPrimaryVertexY", true);
  fMCPico->SetBranchStatus("Event.mPrimaryVertexZ", true);

  fMCPico->SetBranchStatus("Track", true);
  fMCPico->SetBranchStatus("Track.mGMomentumX", true);
  fMCPico->SetBranchStatus("Track.mGMomentumY", true);
  fMCPico->SetBranchStatus("Track.mGMomentumZ", true);
  fMCPico->SetBranchStatus("Track.mOriginX", true);
  fMCPico->SetBranchStatus("Track.mOriginY", true);
  fMCPico->SetBranchStatus("Track.mOriginZ", true);
  fMCPico->SetBranchStatus("Track.mNHitsFit", true);
  fMCPico->SetBranchStatus("Track.mNHitsMax", true);
  fMCPico->SetBranchStatus("Track.mNSigmaPion", true);
  fMCPico->SetBranchStatus("Track.mNSigmaKaon", true);
  fMCPico->SetBranchStatus("Track.mNSigmaProton", true);
  fMCPico->SetBranchStatus("Track.mNSigmaElectron", true);
  fMCPico->SetBranchStatus("Track.mBEmcMatchedTowerIndex", true);
  fMCPico->SetBranchStatus("Track.mIdTruth", true);
  fMCPico->SetBranchStatus("Track.mQATruth", true);

  fMCPico->SetBranchStatus("BTowHit", true);
  fMCPico->SetBranchStatus("BTowHit.mE", true);

  fMCPico->SetBranchStatus("McVertex", true);
  fMCPico->SetBranchStatus("McVertex.mId", true);
  fMCPico->SetBranchStatus("McVertex.mNoDaughters", true);
  fMCPico->SetBranchStatus("McVertex.mIdParTrk", true);
  fMCPico->SetBranchStatus("McVertex.mIsInterm", true);
  fMCPico->SetBranchStatus("McVertex.mTime", true);
  fMCPico->SetBranchStatus("McVertex.mVx", true);
  fMCPico->SetBranchStatus("McVertex.mVy", true);
  fMCPico->SetBranchStatus("McVertex.mVz", true);
  fMCPico->SetBranchStatus("McTrack", true);
  fMCPico->SetBranchStatus("McTrack.mId", true);
  fMCPico->SetBranchStatus("McTrack.mGePid", true);
  fMCPico->SetBranchStatus("McTrack.mCharge", true);
  fMCPico->SetBranchStatus("McTrack.mHits[22]", true);
  fMCPico->SetBranchStatus("McTrack.mPx", true);
  fMCPico->SetBranchStatus("McTrack.mPy", true);
  fMCPico->SetBranchStatus("McTrack.mPz", true);
  fMCPico->SetBranchStatus("McTrack.mE", true);
  fMCPico->SetBranchStatus("McTrack.mIsFromShower", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxStart", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxStop", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxItrmd", true);

  fMCPico->SetBranchAddress("Event", &Event_, &b_Event_);
  fMCPico->SetBranchAddress("Event.mRunId", Event_mRunId, &b_Event_mRunId);
  fMCPico->SetBranchAddress("Event.mEventId", Event_mEventId, &b_Event_mEventId);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexX", Event_mPrimaryVertexX, &b_Event_mPrimaryVertexX);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexY", Event_mPrimaryVertexY, &b_Event_mPrimaryVertexY);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexZ", Event_mPrimaryVertexZ, &b_Event_mPrimaryVertexZ);

  fMCPico->SetBranchAddress("Track", &Track_, &b_Track_);
  fMCPico->SetBranchAddress("Track.mGMomentumX", Track_mGMomentumX, &b_Track_mGMomentumX);
  fMCPico->SetBranchAddress("Track.mGMomentumY", Track_mGMomentumY, &b_Track_mGMomentumY);
  fMCPico->SetBranchAddress("Track.mGMomentumZ", Track_mGMomentumZ, &b_Track_mGMomentumZ);
  fMCPico->SetBranchAddress("Track.mOriginX", Track_mOriginX, &b_Track_mOriginX);
  fMCPico->SetBranchAddress("Track.mOriginY", Track_mOriginY, &b_Track_mOriginY);
  fMCPico->SetBranchAddress("Track.mOriginZ", Track_mOriginZ, &b_Track_mOriginZ);
  fMCPico->SetBranchAddress("Track.mNHitsFit", Track_mNHitsFit, &b_Track_mNHitsFit);
  fMCPico->SetBranchAddress("Track.mNHitsMax", Track_mNHitsMax, &b_Track_mNHitsMax);
  fMCPico->SetBranchAddress("Track.mNSigmaPion", Track_mNSigmaPion, &b_Track_mNSigmaPion);
  fMCPico->SetBranchAddress("Track.mNSigmaKaon", Track_mNSigmaKaon, &b_Track_mNSigmaKaon);
  fMCPico->SetBranchAddress("Track.mNSigmaProton", Track_mNSigmaProton, &b_Track_mNSigmaProton);
  fMCPico->SetBranchAddress("Track.mNSigmaElectron", Track_mNSigmaElectron, &b_Track_mNSigmaElectron);
  fMCPico->SetBranchAddress("Track.mBEmcMatchedTowerIndex", Track_mBEmcMatchedTowerIndex, &b_Track_mBEmcMatchedTowerIndex);
  fMCPico->SetBranchAddress("Track.mIdTruth", Track_mIdTruth, &b_Track_mIdTruth);
  fMCPico->SetBranchAddress("Track.mQATruth", Track_mQATruth, &b_Track_mQATruth);

  fMCPico->SetBranchAddress("BTowHit", &BTowHit_, &b_BTowHit_);
  fMCPico->SetBranchAddress("BTowHit.mE", BTowHit_mE, &b_BTowHit_mE);

  fMCPico->SetBranchAddress("McVertex", &McVertex_, &b_McVertex_);
  fMCPico->SetBranchAddress("McVertex.mId", McVertex_mId, &b_McVertex_mId);
  fMCPico->SetBranchAddress("McVertex.mNoDaughters", McVertex_mNoDaughters, &b_McVertex_mNoDaughters);
  fMCPico->SetBranchAddress("McVertex.mIdParTrk", McVertex_mIdParTrk, &b_McVertex_mIdParTrk);
  fMCPico->SetBranchAddress("McVertex.mIsInterm", McVertex_mIsInterm, &b_McVertex_mIsInterm);
  fMCPico->SetBranchAddress("McVertex.mTime", McVertex_mTime, &b_McVertex_mTime);
  fMCPico->SetBranchAddress("McVertex.mVx", McVertex_mVx, &b_McVertex_mVx);
  fMCPico->SetBranchAddress("McVertex.mVy", McVertex_mVy, &b_McVertex_mVy);
  fMCPico->SetBranchAddress("McVertex.mVz", McVertex_mVz, &b_McVertex_mVz);

  fMCPico->SetBranchAddress("McTrack", &McTrack_, &b_McTrack_);
  fMCPico->SetBranchAddress("McTrack.mId", McTrack_mId, &b_McTrack_mId);
  fMCPico->SetBranchAddress("McTrack.mGePid", McTrack_mGePid, &b_McTrack_mGePid);
  fMCPico->SetBranchAddress("McTrack.mCharge", McTrack_mCharge, &b_McTrack_mCharge);
  fMCPico->SetBranchAddress("McTrack.mHits[22]", McTrack_mHits, &b_McTrack_mHits);
  fMCPico->SetBranchAddress("McTrack.mPx", McTrack_mPx, &b_McTrack_mPx);
  fMCPico->SetBranchAddress("McTrack.mPy", McTrack_mPy, &b_McTrack_mPy);
  fMCPico->SetBranchAddress("McTrack.mPz", McTrack_mPz, &b_McTrack_mPz);
  fMCPico->SetBranchAddress("McTrack.mE", McTrack_mE, &b_McTrack_mE);
  fMCPico->SetBranchAddress("McTrack.mIsFromShower", McTrack_mIsFromShower, &b_McTrack_mIsFromShower);
  fMCPico->SetBranchAddress("McTrack.mIdVtxStart", McTrack_mIdVtxStart, &b_McTrack_mIdVtxStart);
  fMCPico->SetBranchAddress("McTrack.mIdVtxStop", McTrack_mIdVtxStop, &b_McTrack_mIdVtxStop);
  fMCPico->SetBranchAddress("McTrack.mIdVtxItrmd", McTrack_mIdVtxItrmd, &b_McTrack_mIdVtxItrmd);
}

void StHIOverlayAngularities::OutputTreeInit()
{
  TString treename = "Jets";

  outputTree = new TTree(treename.Data(), treename.Data());

  // Branches to save event info
  outputTree->Branch("Centrality", &fRecoJetTree.centrality, "centrality/F");
  outputTree->Branch("Weight", &fRecoJetTree.weight, "weight/F");
  outputTree->Branch("RefMult", &fRecoJetTree.refmult, "refmult/F");
  outputTree->Branch("gRefMult", &fRecoJetTree.grefmult, "grefmult/F");
  outputTree->Branch("RefCorr2", &fRecoJetTree.refcorr2, "refcorr2/F");
  outputTree->Branch("McRefMult", &fMcJetTree.refmult, "mcrefmult/F");
  outputTree->Branch("RecoRefMult", &fMcJetTree.grefmult, "recorefmult/F");
  // outputTree->Branch("McPrimaryVertex", &fMcJetTree.primaryvertex);
  // outputTree->Branch("RecoPrimaryVertex", &fRecoJetTree.primaryvertex);

  outputTree->Branch("McD0Pt", &fMcJetTree.d0pt, "d0pt/F");
  outputTree->Branch("McD0Eta", &fMcJetTree.d0eta, "d0eta/F");
  outputTree->Branch("McD0Phi", &fMcJetTree.d0phi, "d0phi/F");
  outputTree->Branch("McPionPt", &fMcJetTree.pionpt, "pionpt/F");
  outputTree->Branch("McPionEta", &fMcJetTree.pioneta, "pioneta/F");
  outputTree->Branch("McPionPhi", &fMcJetTree.pionphi, "pionphi/F");
  outputTree->Branch("McKaonPt", &fMcJetTree.kaonpt, "kaonpt/F");
  outputTree->Branch("McKaonEta", &fMcJetTree.kaoneta, "kaoneta/F");
  outputTree->Branch("McKaonPhi", &fMcJetTree.kaonphi, "kaonphi/F");

  outputTree->Branch("RecoD0Pt", &fRecoJetTree.d0pt, "d0pt/F");
  outputTree->Branch("RecoD0Eta", &fRecoJetTree.d0eta, "d0eta/F");
  outputTree->Branch("RecoD0Phi", &fRecoJetTree.d0phi, "d0phi/F");
  outputTree->Branch("RecoPionPt", &fRecoJetTree.pionpt, "pionpt/F");
  outputTree->Branch("RecoPionEta", &fRecoJetTree.pioneta, "pioneta/F");
  outputTree->Branch("RecoPionPhi", &fRecoJetTree.pionphi, "pionphi/F");
  outputTree->Branch("RecoKaonPt", &fRecoJetTree.kaonpt, "kaonpt/F");
  outputTree->Branch("RecoKaonEta", &fRecoJetTree.kaoneta, "kaoneta/F");
  outputTree->Branch("RecoKaonPhi", &fRecoJetTree.kaonphi, "kaonphi/F");

  outputTree->Branch("McJetPt", &fMcJetTree.jetpt, "jetpt/F");
  outputTree->Branch("McJetEta", &fMcJetTree.jeteta, "jeteta/F");
  outputTree->Branch("McJetPhi", &fMcJetTree.jetphi, "jetphi/F");
  outputTree->Branch("McJetArea", &fMcJetTree.jetarea, "jetarea/F");
  outputTree->Branch("McJetE", &fMcJetTree.jetenergy, "jetenergy/F");
  outputTree->Branch("McJetNConst", &fMcJetTree.numberofconstituents, "numberofconstituents/I");
  outputTree->Branch("McJetLambda_1_half", &fMcJetTree.lambda_1_half, "lambda_1_half/F");
  outputTree->Branch("McJetLambda_1_1", &fMcJetTree.lambda_1_1, "lambda_1_1/F");
  outputTree->Branch("McJetLambda_1_1half", &fMcJetTree.lambda_1_1half, "lambda_1_1half/F");
  outputTree->Branch("McJetLambda_1_2", &fMcJetTree.lambda_1_2, "lambda_1_2/F");
  outputTree->Branch("McJetLambda_1_3", &fMcJetTree.lambda_1_3, "lambda_1_3/F");
  outputTree->Branch("McJetDispersion", &fMcJetTree.dispersion, "dispersion/F");
  outputTree->Branch("McJetD0Z", &fMcJetTree.d0z, "d0z/F");

  outputTree->Branch("McRecoJetPt", &fMcRecoJetTree.jetpt, "jetpt/F");
  outputTree->Branch("McRecoJetEta", &fMcRecoJetTree.jeteta, "jeteta/F");
  outputTree->Branch("McRecoJetPhi", &fMcRecoJetTree.jetphi, "jetphi/F");
  outputTree->Branch("McRecoJetE", &fMcRecoJetTree.jetenergy, "jetenergy/F");
  outputTree->Branch("McRecoJetArea", &fMcRecoJetTree.jetarea, "jetarea/F");
  outputTree->Branch("McRecoJetNConst", &fMcRecoJetTree.numberofconstituents, "numberofconstituents/I");
  outputTree->Branch("McRecoJetLambda_1_half", &fMcRecoJetTree.lambda_1_half, "lambda_1_half/F");
  outputTree->Branch("McRecoJetLambda_1_1", &fMcRecoJetTree.lambda_1_1, "lambda_1_1/F");
  outputTree->Branch("McRecoJetLambda_1_1half", &fMcRecoJetTree.lambda_1_1half, "lambda_1_1half/F");
  outputTree->Branch("McRecoJetLambda_1_2", &fMcRecoJetTree.lambda_1_2, "lambda_1_2/F");
  outputTree->Branch("McRecoJetLambda_1_3", &fMcRecoJetTree.lambda_1_3, "lambda_1_3/F");
  outputTree->Branch("McRecoJetDispersion", &fMcRecoJetTree.dispersion, "dispersion/F");
  outputTree->Branch("McRecoJetD0Z", &fMcRecoJetTree.d0z, "d0z/F");

  outputTree->Branch("RecoJetPt", &fRecoJetTree.jetpt, "jetpt/F");
  outputTree->Branch("RecoJetEta", &fRecoJetTree.jeteta, "jeteta/F");
  outputTree->Branch("RecoJetPhi", &fRecoJetTree.jetphi, "jetphi/F");
  outputTree->Branch("RecoJetArea", &fRecoJetTree.jetarea, "jetarea/F");
  outputTree->Branch("RecoJetE", &fRecoJetTree.jetenergy, "jetenergy/F");
  outputTree->Branch("RecoJetRhoVal", &fRecoJetTree.fRhoValforjet, "fRhoValforjet/F");
  outputTree->Branch("RecoJetNConst", &fRecoJetTree.numberofconstituents, "numberofconstituents/I");
  outputTree->Branch("RecoJetLambda_1_half", &fRecoJetTree.lambda_1_half, "lambda_1_half/F");
  outputTree->Branch("RecoJetLambda_1_1", &fRecoJetTree.lambda_1_1, "lambda_1_1/F");
  outputTree->Branch("RecoJetLambda_1_1half", &fRecoJetTree.lambda_1_1half, "lambda_1_1half/F");
  outputTree->Branch("RecoJetLambda_1_2", &fRecoJetTree.lambda_1_2, "lambda_1_2/F");
  outputTree->Branch("RecoJetLambda_1_3", &fRecoJetTree.lambda_1_3, "lambda_1_3/F");
  outputTree->Branch("RecoJetDispersion", &fRecoJetTree.dispersion, "dispersion/F");
  outputTree->Branch("RecoJetD0Z", &fRecoJetTree.d0z, "d0z/F");
}

void StHIOverlayAngularities::FillTree(const Int_t &numberOfD0Events)
{

  for (Int_t iMcD0Event = 0; iMcD0Event < numberOfD0Events; iMcD0Event++)
  {
    fMcJetTree.Clear();
    fRecoJetTree.Clear();
    fMcRecoJetTree.Clear();

    fRecoJetTree.centrality = fCentralityScaled;
    fRecoJetTree.weight = mCentMaker->GetWeight();

    TVector3 mcPion = fMcD0Information[iMcD0Event].first;
    TVector3 mcKaon = fMcD0Information[iMcD0Event].second;

    TVector3 mcD0 = mcPion + mcKaon;

    fMcJetTree.d0pt = mcD0.Perp();
    fMcJetTree.d0eta = mcD0.PseudoRapidity();
    fMcJetTree.d0phi = standardPhi(mcD0.Phi());

    if (fPrintLevel)
      cout << "MC D0 Saved = " << fMcJetTree.d0pt << "\t" << fMcJetTree.d0eta << "\t" << fMcJetTree.d0phi << endl;

    fMcJetTree.pionpt = mcPion.Perp();
    fMcJetTree.pioneta = mcPion.PseudoRapidity();
    fMcJetTree.pionphi = standardPhi(mcPion.Phi());

    fMcJetTree.kaonpt = mcKaon.Perp();
    fMcJetTree.kaoneta = mcKaon.PseudoRapidity();
    fMcJetTree.kaonphi = standardPhi(mcKaon.Phi());

    fMcJetTree.refmult = fMcEventTracks[iMcD0Event].size();
    fMcJetTree.grefmult = fRecoMcEventTracks[iMcD0Event].size();

    TVector3 recoPion = fMcRecoD0Information[iMcD0Event].first;
    TVector3 recoKaon = fMcRecoD0Information[iMcD0Event].second;

    TVector3 recoD0 = recoPion + recoKaon;

    fRecoJetTree.d0pt = recoD0.Perp();
    fRecoJetTree.d0eta = recoD0.PseudoRapidity();
    fRecoJetTree.d0phi = standardPhi(recoD0.Phi());

    fRecoJetTree.pionpt = recoPion.Perp();
    fRecoJetTree.pioneta = recoPion.PseudoRapidity();
    fRecoJetTree.pionphi = standardPhi(recoPion.Phi());

    fRecoJetTree.kaonpt = recoKaon.Perp();
    fRecoJetTree.kaoneta = recoKaon.PseudoRapidity();
    fRecoJetTree.kaonphi = standardPhi(recoKaon.Phi());

    fRecoJetTree.refmult = mCentMaker->GetrefMult();
    fRecoJetTree.grefmult = mCentMaker->GetgrefMult();
    fRecoJetTree.refcorr2 = mCentMaker->GetRefCorr2();

    StJet *jetMc = fMcJet[iMcD0Event];
    StJet *jetMcReco = fMcRecoJet[iMcD0Event];
    StJet *jetReco = fRecoJet[iMcD0Event];

    FillJet(jetMc, fMcJetTree, mcD0);
    if (fMcRecoJet[iMcD0Event] != NULL)
      FillJet(jetMcReco, fMcRecoJetTree, recoD0);
    if (fRecoJet[iMcD0Event] != NULL)
    {
      fRecoJetTree.fRhoValforjet = fRhoVal;
      FillJet(jetReco, fRecoJetTree, recoD0);
    }
    outputTree->Fill();

    hMcJetPt->Fill(fMcJetTree.jetpt);
    hMcD0Pt->Fill(fMcJetTree.d0pt);
    hMcJetZ->Fill(fMcJetTree.d0z);
    if (fRecoJetTree.numberofconstituents != 0)
    {
      hMcRecoD0Pt->Fill(fRecoJetTree.d0pt);
      hRecoJetPt->Fill(fRecoJetTree.jetpt);
      hMcRecoJetPt->Fill(fMcRecoJetTree.jetpt);
      hRecoJetZ->Fill(fRecoJetTree.d0z);
      hMcD0PtMcRecoD0Pt->Fill(fMcJetTree.d0pt, fRecoJetTree.d0pt);
      hRecoJetPtMcJetPt->Fill(fRecoJetTree.jetpt, fMcJetTree.jetpt);
      hRecoJetZMcJetZ->Fill(fRecoJetTree.d0z, fMcJetTree.d0z);
    }
  }
  if (fPrintLevel)
    cout << "Filling Tree complete" << endl;
}

void StHIOverlayAngularities::FillJet(StJet *jet, StJetTreeStruct &jetTree, const TVector3 &D0)
{
  jetTree.d0z = (jet->Px() * D0.X() + jet->Py() * D0.Y()) / (jet->Pt() * jet->Pt());
  jetTree.jetpt = jet->Pt();
  jetTree.jeteta = jet->Eta();
  jetTree.jetphi = jet->Phi();
  jetTree.jetarea = jet->Area();
  jetTree.jetenergy = jet->E();
  jetTree.numberofconstituents = jet->GetJetConstituents().size();
  jetTree.lambda_1_half = jet->GetAngularityLambda(1., 0.5) / fJetRad;
  jetTree.lambda_1_1 = jet->GetAngularityLambda(1., 1.) / fJetRad;
  jetTree.lambda_1_1half = jet->GetAngularityLambda(1., 1.5) / fJetRad;
  jetTree.lambda_1_2 = jet->GetAngularityLambda(1., 2) / fJetRad;
  jetTree.lambda_1_3 = jet->GetAngularityLambda(1., 3) / fJetRad;
  jetTree.dispersion = sqrt(jet->GetAngularityLambda(2., 0.));
  delete jet;
}