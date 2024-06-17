// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StHIOverlay_Test.h"
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

ClassImp(StHIOverlay_Test)

    //________________________________________________________________________
    StHIOverlay_Test::StHIOverlay_Test(const char *name, StPicoDstMaker *picoMaker, const char *outName = "", const char *filename = "") : StJetFrameworkPicoBase(name) // StMaker(name),
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
  fRunFlag = 0; // see StHIOverlay_Test::fRunFlagEnum
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
  // fjw("HIMC", "HIMC");
  fJets = 0x0;
  // fJetsArr = 0x0;
  fJetsArr.clear();
  fRecoJetsArr.clear();

  fMCFileListName = filename;

  for (int i = 0; i < 8; i++)
  {
    fEmcTriggerArr[i] = 0;
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

  fSetNumberOfEvents = 0;
  fNumberofeventsoverlayed = 0;

  fMCJetPtCut = 5.0;

  for (int i = 0; i < 100; i++)
  {
    fMCEventTracks[i].clear();
    fMCEventTowers[i].clear();
    fRecoEventTracks[i].clear();
    fRecoEventTowers[i].clear();
    fMCD0Information[i] = {};
    fRecoD0Information[i] = {};
    fOrigin[i].SetXYZ(0, 0, 0);
    fMCEventInfo[i] = {};
  }

  randomevents = kTRUE;
  fPrintLevel = 0;

  centraleventsfound = 0;

  fCentBin = 0;

  fTrackingEfficiency = kFALSE;
  fTrackingEfficiencyPercentage = 0.95;

  if (!name)
    return;
  SetName(name);
}

//
//________________________________________________________________________
StHIOverlay_Test::~StHIOverlay_Test()
{ /*  */
}

//
//________________________________________________________________________
Int_t StHIOverlay_Test::Init()
{

  StJetFrameworkPicoBase::Init();

  // position object for Emc
  mBemcGeom = StEmcGeom::instance("bemc");
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  // fJets->SetName(fJetsName);
  fJetsArr.clear();
  fRecoJetsArr.clear();

  ifstream filelistforMCEvents(fMCFileListName.Data());

  if (!filelistforMCEvents.is_open())
  {
    LOG_ERROR << "No MC File List! Exiting!" << endm;
    return kStOk;
  }

  string line;

  // cout << "Files read in: " << endl;

  while (getline(filelistforMCEvents, line))
  {
    TString s(line);

    // cout << s << endl;
    filenamesforHIOverlay.push_back(s);
  }

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

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StHIOverlay_Test::Finish()
{
  cout << "StHIOverlay_Test::Finish()\n";

  //  Write histos to file and close it.
  if (mOutName != "")
  {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout << "End of StHIOverlay_Test::Finish" << endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

// OLD user code says: //  Called every event after Make().
//_____________________________________________________________________________
void StHIOverlay_Test::Clear(Option_t *opt)
{
  // fjw->Clear();
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StHIOverlay_Test::Make()
{
  // fjw->Clear();
  fJets->Delete();

  fJetsArr.clear();
  fRecoJetsArr.clear();
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  for (int i = 0; i < 100; i++)
  {
    fMCEventTracks[i].clear();
    fMCEventTowers[i].clear();
    fRecoEventTracks[i].clear();
    fRecoEventTowers[i].clear();
    fMCD0Information[i] = {};
    fRecoD0Information[i] = {};
    fOrigin[i].SetXYZ(0, 0, 0);
    fMCEventInfo[i] = {};
  }

  fNumberofeventsoverlayed = 0;

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

  // for(int i = 1; i<4801; i++) {  if(!mBaseMaker->IsTowerOK(i))  cout<<"tower: "<<i<<" is not good!!"<<endl;  }
  //  can check for bad towers like this (if needed):
  //  bool isTowOK = mBaseMaker->IsTowerOK(towerID);
  //  ===========================================================================================

  // cut event on max track pt > 30.0 GeV
  if (GetMaxTrackPt() > fMaxEventTrackPt)
    return kStOK;

  // cut event on max tower Et > 30.0 GeV
  // if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  double zVtx_VPD = mPicoEvent->vzVpd();

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
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult = mCentMaker->GetrefMult();   // see StPicoEvent
  ref9 = mCentMaker->GetRef9();             // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16();           // binning from central -> peripheral
  int cent16 = mCentMaker->GetCent16();     // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  // double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  //  for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cout << "Centrality = " << grefMult << "\t" << ref16 << "\t" << cent16 << "\t" << fCentralityScaled << endl;

  if (fCentralityScaled >= 0 && fCentralityScaled < 10)
    centralitybinforefficiency = 0;
  else if (fCentralityScaled >= 10 && fCentralityScaled < 40)
    centralitybinforefficiency = 1;
  else if (fCentralityScaled >= 40 && fCentralityScaled < 80)
    centralitybinforefficiency = 2;

  if (!randomevents)
  {
    cout << "Cent = " << centralitybinforefficiency << "\t" << fCentBin << endl;
    if (centralitybinforefficiency != fCentBin)
      return kStOk;
  } // This is just for test. I want to check it centrality bin-by-bin

  // cut on unset centrality, > 80%
  if (cent16 == -1)
    return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them

  // cout << "1" << endl;

  // cut on centrality for analysis before doing anything
  if (fRequireCentSelection)
  {
    if (!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut))
      return kStOk;
  }

  // cout << "2" << endl;

  // ============================ end of CENTRALITY ============================== //

  if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5)
    return kStOK;

  // cout << "3" << endl;

  if ((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut))
    return kStOk;

  // cout << "4" << endl;

  int arrMB5[] = {450005, 450015, 450025, 450050, 450060,                          // Run14_AuAu200
                  520001, 520011, 520021, 520031, 520041, 520051, 570002, 570001}; // Run16_AuAu200

  bool matchMB = kFALSE;

  for (int i = 0; i < sizeof(arrMB5) / sizeof(*arrMB5); i++)
  {
    if (mPicoEvent->isTrigger(arrMB5[i]))
      matchMB = kTRUE;
    if (matchMB)
      break;
  }

  if (!matchMB)
    return kStOk;

  // cout << "5" << endl;

  if (abs(zVtx) > 6.)
    return kStOk;

  // cout << "6" << endl;

  if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > 2.)
    return kStOk;

  // cout << "7" << endl;

  if (abs(zVtx - zVtx_VPD) > 3)
    return kStOk;

  // cout << "8" << endl;

  // cout << "Centrality = " << fCentralityScaled << endl;
  // pick a random file
  // PickARandomFile();

  SampleMCEvents();
  // cout << "Array Size == " << fNumberofeventsoverlayed << endl;
  // for (int i = 0; i < fNumberofeventsoverlayed; i++){
  //   cout << "Reco Track Tower Size = " << fRecoEventTracks[i].size() << "\t" << fRecoEventTowers[i].size() << endl;
  // }

  // cout << "Jets Array Size (HIOverlay): " << fJetsArr.size() << endl;

  return kStOK;
}

void StHIOverlay_Test::SampleMCEvents()
{

  // cout << "Running sample MC events" << endl;

  fJetsArr.clear();
  fRecoJetsArr.clear();

  TRandom3 *r2 = new TRandom3(0);

  int numberofmceventstoconvulate = fSetNumberOfEvents;

  AssignTreeVariables();

  int nentries = fMCPico->GetEntriesFast();

  vector<int> randomlisttoevents;

  randomlisttoevents.clear();

  int counter = 0;
  int eventsconsidered = 0;

  while (counter < numberofmceventstoconvulate)
  {

    int rand = r2->Integer(nentries);
    if (std::find(randomlisttoevents.begin(), randomlisttoevents.end(), rand) != randomlisttoevents.end())
      continue;
    randomlisttoevents.push_back(rand);
    if (fPrintLevel)
    {
      cout << randomlisttoevents[randomlisttoevents.size() - 1] << endl;
      cout << "Size of randomlisttoevents = " << randomlisttoevents.size() << endl;
    }

    if (randomlisttoevents.size() > nentries)
      continue;

    if (rand >= nentries)
      continue;

    fMCPico->GetEntry(rand);

    if (fPrintLevel)
      cout << "======================= Event Number After MC = " << Event_mEventId[0] << "\t" << counter << "=======================" << endl;

    fvertextotrack = GetMCTrackListForVertex();

    GetD0AndDaughters();

    if (fPrintLevel)
      cout << "=============== Number of D0s =============== " << vertexids.size() << "\t" << pionids.size() << endl;

    if (vertexids.size() == 0)
      continue; // While loop continues.

    int D0Counter = 0;

    // Loop over each D0
    for (int D0 = 0; D0 < vertexids.size(); D0++)
    {

      // cout << "This is D0 # " << D0 << "***********************" << endl;

      // I only want to use TLorentzVector to propagate the information ahead. All the processing is done within this class.
      // This means I need to find a way to make sure I can identify the D0 track when I save out the information.
      // Since only D0 needs to be identified, I am saving the mass information for the D0 as 1.865.
      // All the other tracks are saved out with pi0mass, because ultimately, jet constituents are assumed to have that mass.

      // How do I propagate charge information though?
      // Do I need it? Mostly, nope! In fact, once I separate track and tower, it should be enough.

      // This loop fills the input vector for the MC Side

      // cout << "Number of MC tracks = " << McTrack_ << endl;

      // This function prepares the input list for that event for a particular D0.
      // The event list will of course be a little different for each D0, even for the same event

      PrepareSetOfMCTracks(counter, D0);

      // This function makes the jets with the MC tracks.
      // If it doesn't find a D0 jet with pt > 5 GeV in the event,
      // we discard the whole event and go to the next case.

      if (!DoesItHaveAGoodMCJet(counter))
      {
        fMCEventTracks[counter].clear();
        fMCEventTowers[counter].clear();
        continue;
      }

      fMCEventInfo[counter].first = Event_mRunId[0];
      fMCEventInfo[counter].second = Event_mEventId[0];

      if (fPrintLevel)
        cout << "Event Number After MC = " << Event_mEventId[0] << "\t" << counter << endl;

      PrepareSetOfRecoInput(counter, D0);

      DoesItHaveAGoodRecoJet(counter);

      if (fPrintLevel)
        cout << "Event Number After Reco = " << Event_mEventId[0] << "\t" << counter << endl;

      if (fPrintLevel == -1)
      {
        for (int i = 0; i < fRecoEventTracks[counter].size(); i++)
        {
          TLorentzVector v = fRecoEventTracks[counter][i];
          if (int(v.M()) != 1)
            cout << "Reco Constituents = " << v.Pt() << "\t" << v.PseudoRapidity() << "\t" << v.Phi() << "\t" << v.M() << endl;
          else
            cout << "Reco D0 = " << v.Pt() << "\t" << v.PseudoRapidity() << "\t" << v.Phi() << "\t" << v.M() << endl;
        }

        for (int i = 0; i < fRecoEventTowers[counter].size(); i++)
        {
          TLorentzVector v = fRecoEventTowers[counter][i];
          cout << "Reco Tower = " << v.Pt() << "\t" << v.PseudoRapidity() << "\t" << v.Phi() << "\t" << v.M() << endl;
        }
      }

      counter++; // Since one MC event can be invoked multiple times (due to having multiple D0s, this step is necessary.)
    }
    eventsconsidered++;
    // }
  }

  if (fPrintLevel)
    cout << "Counter = " << counter << endl;
  fNumberofeventsoverlayed = counter;

  fMCPico->Reset();
  f->Close();
}

void StHIOverlay_Test::PrepareSetOfRecoInput(int eventnum, int D0num)
{
  const int numberoftowers = BTowHit_;
  double towerenergy[numberoftowers];
  for (int i = 0; i < numberoftowers; i++)
  {
    towerenergy[i] = 0.;
  }

  vector<int> mTowerToTrack[4800];
  vector<double> mTowerToTrackE[4800];

  for (int i = 0; i < numberoftowers; i++)
  {
    mTowerToTrack[i].clear();
    mTowerToTrackE[i].clear();
  }

  TVector3 oVertex;
  oVertex.SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

  if (fPrintLevel)
    cout << "HIOverlay Vertex " << Event_mPrimaryVertexX[0] << "\t" << Event_mPrimaryVertexY[0] << "\t" << Event_mPrimaryVertexZ[0] << endl;

  // This loop fills the input vector for the RECO side
  for (int reco = 0; reco < Track_; reco++)
  {

    TVector3 o;
    o.SetXYZ(Track_mOriginX[reco], Track_mOriginY[reco], Track_mOriginZ[reco]);
    TVector3 g;
    g.SetXYZ(Track_mGMomentumX[reco], Track_mGMomentumY[reco], Track_mGMomentumZ[reco]);

    // track variables
    double pt = g.Perp();
    double phi = g.Phi();
    if (phi < 0.0)
      phi += 2.0 * pi; // force from 0-2pi
    if (phi > 2.0 * pi)
      phi -= 2.0 * pi; // force from 0-2pi
    double eta = g.PseudoRapidity();
    double px = g.x();
    double py = g.y();
    double pz = g.z();
    double p = g.Mag();
    double energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
    short charge = (Track_mNHitsFit[reco] > 0) ? 1 : -1;

    double dca = (oVertex - o).Mag();

    bool goodtrack = (dca < fJetTrackDCAcut) && (abs(Track_mNHitsFit[reco]) >= fTracknHitsFit) && (abs(double(Track_mNHitsFit[reco])) / double(Track_mNHitsMax[reco]) >= fTracknHitsRatio);

    //// Variables For FastSim

    double pt_new = pt;
    double phi_new = phi;
    double eta_new = eta;
    double px_new = px;
    double py_new = py;
    double pz_new = pz;
    double p_new = p;
    double energy_new = energy;
    short charge_new = charge;

    bool mctrackavailable = kTRUE;

    int mcid = Track_mIdTruth[reco] - 1;

    if (mcid < 0)
      mctrackavailable = kFALSE;

    bool isatrackfromD0 = kFALSE;

    // Here, we have two paths to take. If the track needs replacement, we replace it with the fastsim method that is standardised.
    // Else, the pt, eta, phi are sent as is to the final vector.

    if (mctrackavailable && goodtrack)
    {
      if (reco == matchedpionids[D0num] || reco == matchedkaonids[D0num])
        isatrackfromD0 = kTRUE; // Kaons and Pions that come from the current D0 need to be tossed, and replaced by the fast sim version

      TVector3 mg(McTrack_mPx[mcid], McTrack_mPy[mcid], McTrack_mPz[mcid]);

      double relativesmearing = TMath::Sqrt(pow(mg.Px() - px, 2) + pow(mg.Py() - py, 2)) / (mg.Pt());

      double fastsimsmearing;

      int pid = McTrack_mGePid[mcid];

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
    if (fPrintLevel == -19)
    {
      if (mctrackavailable)
        cout << Form("Track # %i \t %i \t %i \t %.2f \t %.2f \t %i \t %i \t %.2f \t %i \t %.2f", reco, mcid, McTrack_mGePid[mcid], pt_new, eta_new, isatrackfromD0, abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1, dca, abs(int(Track_mNHitsFit[reco])), abs(double(Track_mNHitsFit[reco])) / double(Track_mNHitsMax[reco])) << endl;
      else
        cout << Form("Track # %i \t %.2f \t %.2f \t %i \t %.2f \t %.2f", reco, pt_new, eta_new, abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1, dca, abs(double(Track_mNHitsFit[reco])) / double(Track_mNHitsMax[reco]), Track_mQATruth[reco]) << endl;
    }

    if (!goodtrack)
      continue;
    // DCA based cuts precede everything else.

    // Track from D0 -> K Pi || D0 is in acceptance range. The KPi do not need to be in acceptance. || The KPi track is projected onto the towers.
    if (isatrackfromD0)
    {
      int matchedTowerIndex = abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1;
      if (matchedTowerIndex >= 0)
      {
        towerenergy[matchedTowerIndex] += energy;
        mTowerToTrack[matchedTowerIndex].push_back(reco);
        mTowerToTrackE[matchedTowerIndex].push_back(energy);
      }
    }

    int particleid = -99;

    double nsigpion = Track_mNSigmaPion[reco] / 1000.;
    double nsigkaon = Track_mNSigmaKaon[reco] / 1000.;
    double nsigproton = Track_mNSigmaProton[reco] / 1000.;

    if (abs(nsigpion) < 2 && abs(nsigkaon) > 2. && abs(nsigproton) > 2.)
      particleid = 1;
    else if (abs(nsigpion) > 2 && abs(nsigkaon) < 2. && abs(nsigproton) > 2.)
      particleid = 2;
    else if (abs(nsigpion) > 2 && abs(nsigkaon) > 2. && abs(nsigproton) < 2.)
      particleid = 3 * charge;

    bool removetrack = kFALSE;
    bool isD0DaugDescendant = kFALSE;
    if (!isatrackfromD0 && pt > fMinJetTrackPt && pt < fMaxJetTrackPt && (eta > fJetTrackEtaMin) && (eta < fJetTrackEtaMax))
    // I am using the old pt for this because the efficiencies were derived with the old pt.
    // I am keeping things consistent with MC
    {
      // KeepTrack returns False if the track is to be discarded.
      // RemoveTrack = True if KeepTrack is False
      removetrack = !KeepTrack(particleid, centralitybinforefficiency, pt);
      // If D0 descendant, record that as well.
      isD0DaugDescendant = IsTrackADescendantOfD0Daughters(reco);
    }

    if (removetrack || isD0DaugDescendant)
    {
      if (fPrintLevel == -19)
      {
        if (mctrackavailable)
          cout << Form("Removed Track Input = %.2f \t %.2f \t %.2f \t %i", px_new, py_new, pz_new, McTrack_mGePid[mcid]) << endl;
        else
          cout << Form("Removed Track Input = %.2f \t %.2f \t %.2f", px_new, py_new, pz_new) << endl;
        // if (mctrackavailable) cout << "Removed Track Input = " << px_new << "\t" << py_new << "\t" << pz_new << "\t" << McTrack_mGePid[mcid] << endl;
        // else cout << "Removed Track Input = " << px_new << "\t" << py_new << "\t" << pz_new << endl;
      }
    }

    // if (isD0DaugDescendant) {
    //   if (mctrackavailable) cout << Form("Removed D0 Track Input = %.2f \t %.2f \t %.2f \t %i", px_new, py_new, pz_new, McTrack_mGePid[mcid]) << endl;
    //   else cout << Form("Removed D0 Track Input = %.2f \t %.2f \t %.2f", px_new, py_new, pz_new) << endl;
    //   // if (mctrackavailable) cout << "Removed Track Input = " << px_new << "\t" << py_new << "\t" << pz_new << "\t" << McTrack_mGePid[mcid] << endl;
    //   // else cout << "Removed Track Input = " << px_new << "\t" << py_new << "\t" << pz_new << endl;
    // }

    // if (isatrackfromD0){
    //   if (mctrackavailable) cout << "Removed D0 Reco Level Track = " << px_new << "\t" << py_new << "\t" << pz_new << "\t" << McTrack_mGePid[mcid] << endl;
    //   else cout << "Removed D0 Reco Level Track = " << px_new << "\t" << py_new << "\t" << pz_new << endl;
    // }

    // if (isD0DaugDescendant){
    //   if (mctrackavailable) cout << "Removed D0 Daug Reco Level Track = " << px_new << "\t" << py_new << "\t" << pz_new << "\t" << McTrack_mGePid[mcid] << endl;
    //   else cout << "Removed D0 Daug Reco Level Track = " << px_new << "\t" << py_new << "\t" << pz_new << endl;
    // }

    // jet track acceptance cuts now
    if (pt_new < fMinJetTrackPt)
      continue;
    if (pt_new > fMaxJetTrackPt)
      continue; // 20.0 STAR, 100.0 ALICE
    if ((eta_new < fJetTrackEtaMin) || (eta_new > fJetTrackEtaMax))
      continue;
    if (phi_new < 0.0)
      phi_new += 2.0 * pi; // force from 0-2pi
    if (phi_new > 2.0 * pi)
      phi_new -= 2.0 * pi; // force from 0-2pi
    if ((phi_new < fJetTrackPhiMin) || (phi_new > fJetTrackPhiMax))
      continue;

    // additional quality cuts for tracks

    // This is questionable. Here, I don't think I should use it. But, in the QM method, I should?
    // I have now included these cuts in both. I think that's the correct thing to do.

    // cout << "TRACK = " << px_new << "\t" << py_new << "\t" << pz_new << endl;

    // This place takes care of all energy depositions due to tracks accepted for jet reco.
    // For the reco tracks that we discard, we still need to subtract their contribution from the tower
    // It also includes the kaon and pion from D0.

    // int matchedTowerIndex = abs(GetMatchedBtowID(Track_mBEmcMatchedTowerIndex[reco], g, o, charge_new)) - 1; // towerIndex = towerID - 1

    bool ignoretrack = removetrack || isatrackfromD0 || isD0DaugDescendant;
    // bool ignoretrack = !KeepTrack(particleid, centralitybinforefficiency, pt_new) || isatrackfromD0;

    if (ignoretrack)
      continue; // To match the efficiency, we start tossing random tracks.

    int matchedTowerIndex = abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1;
    if (matchedTowerIndex >= 0)
    {
      towerenergy[matchedTowerIndex] += energy;
      mTowerToTrack[matchedTowerIndex].push_back(reco);
      mTowerToTrackE[matchedTowerIndex].push_back(energy);
    }

    TLorentzVector v;
    v.SetXYZM(px_new, py_new, pz_new, pi0mass);

    fRecoEventTracks[eventnum].push_back(v);
    //// Fill the Tower Array with Energy Depositions here
    if (fPrintLevel == -19)
    {
      if (mctrackavailable)
        cout << Form("Reco Level Track Input = %.2f \t %.2f \t %.2f \t %i", px_new, py_new, pz_new, McTrack_mGePid[mcid]) << endl;
      else
        cout << Form("Reco Level Track Input = %.2f \t %.2f \t %.2f", px_new, py_new, pz_new) << endl;
    }
  }

  TVector3 mcKaon(McTrack_mPx[kaonids[D0num]], McTrack_mPy[kaonids[D0num]], McTrack_mPz[kaonids[D0num]]);
  TVector3 mcPion(McTrack_mPx[pionids[D0num]], McTrack_mPy[pionids[D0num]], McTrack_mPz[pionids[D0num]]);

  fMCD0Information[eventnum] = {mcPion, mcKaon};

  TVector3 recoKaon = FastSimMom(mcKaon, 11);
  TVector3 recoPion = FastSimMom(mcPion, 8);

  fRecoD0Information[eventnum] = {recoPion, recoKaon};

  fOrigin[eventnum].SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

  // fOrigin[counter].push_back(eventorigin);

  TVector3 recoD0;
  recoD0 = recoKaon + recoPion;

  TLorentzVector v;
  v.SetXYZM(recoD0.x(), recoD0.y(), recoD0.z(), 1.865); // This is a D0

  if (fPrintLevel)
    cout << "Reco D0 Entered Jet Loop = " << v.Pt() << "\t" << v.PseudoRapidity() << endl;

  // cout << "Reco Level Track Input D0 = " << recoD0.x() << "\t" << recoD0.y() << "\t" << recoD0.z() << endl;

  if (fPrintLevel == -19)
    cout << Form("Reco Level Track Input D0 = %.2f \t %.2f \t %.2f", recoD0.x(), recoD0.y(), recoD0.z()) << endl;

  fRecoEventTracks[eventnum].push_back(v);

  // This loop fills the input vector for the towers

  for (int tower = 0; tower < numberoftowers; tower++)
  {
    int towerID = tower + 1;
    if (towerID < 0)
      continue; // double check these aren't still in the event list

    TVector3 towerPosition = mEmcPosition->getPosFromVertex(oVertex, towerID);
    double towerPhi = towerPosition.Phi();
    if (towerPhi < 0.0)
      towerPhi += 2.0 * pi; // force from 0-2pi
    if (towerPhi > 2.0 * pi)
      towerPhi -= 2.0 * pi; // force from 0-2pi
    double towerEta = towerPosition.PseudoRapidity();

    // check for bad (and dead) towers
    bool TowerOK = mBaseMaker->IsTowerOK(towerID);     // kTRUE means GOOD
    bool TowerDead = mBaseMaker->IsTowerDead(towerID); // kTRUE means BAD
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

    double towerEunCorr = double(BTowHit_mE[tower]) / 1000.; // uncorrected energy
    double towerE = double(BTowHit_mE[tower]) / 1000.;       // corrected energy (hadronically - done below)
    double towEtunCorr = towerE / (1.0 * TMath::CosH(towerEta));

    if (fPrintLevel > 10)
      cout << "Tower = " << tower << "\t" << towerE << "\t" << towEtunCorr << endl;

    // cut on min tower energy after filling histos
    if (towerEunCorr < mTowerEnergyMin)
      continue; // if we don't have enough E to start with, why mess around

    // cout << Form("Reco Tower = %i \t %.2f \t %.2f", towerID, towEtunCorr, towerE) << endl;

    // =======================================================================
    // HADRONIC CORRECTION

    double sumEt = (towerEunCorr - towerenergy[tower]) / (1.0 * TMath::CosH(towerEta));
    double towerEt = sumEt;
    // cout << Form("FReco Tower = %i \t %.2f \t %.2f", towerID, towerEt, towerenergy[tower]) << endl;
    // for (int i = 0; i < mTowerToTrack[tower].size(); i++){
    //   cout << Form("Track Index matched to Tower = %i \t %i \t %.2f", towerID, mTowerToTrack[tower][i], mTowerToTrackE[tower][i]) << endl;
    // }
    if (towerEt < mTowerEnergyMin)
      continue;
    towerE = towerEt * 1.0 * TMath::CosH(towerEta);

    if (fPrintLevel > 10)
      cout << "Tower = " << tower << "\t" << towerE << "\t" << sumEt << endl;

    Double_t p = 1.0 * TMath::Sqrt(towerE * towerE - pi0mass * pi0mass);

    double posX = towerPosition.x();
    double posY = towerPosition.y();
    double posZ = towerPosition.z();

    Double_t r = TMath::Sqrt(posX * posX + posY * posY + posZ * posZ);

    TLorentzVector v;
    v.SetXYZM(p * posX / r, p * posY / r, p * posZ / r, pi0mass);

    // cout << "TOWER = " << v.X() << "\t" << v.Y() << "\t" << v.Z() << endl;

    fRecoEventTowers[eventnum].push_back(v);

    // cout << "Reco Level Tower Input = " << v.X() << "\t" << v.Y() << "\t" << v.Z() << endl;
    if (fPrintLevel == -19)
      cout << Form("Reco Level Tower Input = %i \t %.2f \t %.2f \t %.2f \t %.2f", towerID, v.X(), v.Y(), v.Z(), towerE) << endl;
  }

  // cout << "Input Loop Size = " << fRecoEventTracks[eventnum].size() + fRecoEventTowers[eventnum].size() << endl;
}

Bool_t StHIOverlay_Test::IsTrackADescendantOfD0Daughters(int trackid)
{ // This is for reco tracks (Needs to be called after getting the list of discarded tracks from MC)
  int mctrkid = Track_mIdTruth[trackid] - 1;

  if (std::find(fDroppedMCTracks.begin(), fDroppedMCTracks.end(), mctrkid) != fDroppedMCTracks.end())
    return kTRUE;

  return kFALSE;
}

Bool_t StHIOverlay_Test::KeepTrack(int particleid, int centralitybin, double pt)
{
  bool keeptrack = kTRUE;

  TRandom3 *r = new TRandom3(0);
  double rando = r->Rndm();

  // cout << "Random = " << rando << "\t" << fPionWeight->Eval(pt) << "\t" << fKaonWeight->Eval(pt) << "\t" << fProtonWeight->Eval(pt) << "\t" << fAProtonWeight->Eval(pt) << endl;

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

  return keeptrack;
}

void StHIOverlay_Test::MCTracksToDiscard(int D0num)
{ // This is for MC tracks

  fDroppedMCTracks.clear();

  int daugstopvx1 = McTrack_mIdVtxStop[pionids[D0num]] - 1;
  int daugstopvx2 = McTrack_mIdVtxStop[kaonids[D0num]] - 1;

  if (fPrintLevel > 1)
    cout << "Pion" << endl;
  GetAllTracksFromVertex(daugstopvx1, fDroppedMCTracks);
  if (fPrintLevel > 1)
    cout << "Kaon" << endl;
  GetAllTracksFromVertex(daugstopvx2, fDroppedMCTracks);

  if (fPrintLevel)
  {
    cout << "Dropped Track List for D0 # = " << D0num << endl;
    for (int i = 0; i < fDroppedMCTracks.size(); i++)
    {
      TVector3 v(McTrack_mPx[fDroppedMCTracks[i]], McTrack_mPy[fDroppedMCTracks[i]], McTrack_mPz[fDroppedMCTracks[i]]);
      cout << fDroppedMCTracks[i] << "\t" << McTrack_mGePid[fDroppedMCTracks[i]] << "\t" << v.Perp() << endl;
    }
    cout << endl;
  }
}

void StHIOverlay_Test::PrepareSetOfMCTracks(int eventnum, int D0num)
{
  /// Here eventnum refers to individual entries. If an MC event has 2 D0s and both are part of jets with pT > 5 GeV, we stage the event as 2 different events.

  MCTracksToDiscard(D0num);

  for (int mc = 0; mc < McTrack_; mc++)
  {

    TLorentzVector v;
    if ((McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38))
      v.SetXYZM(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], 1.865); // This is a D0
    else
      v.SetXYZM(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], pi0mass);

    // track variables
    double pt = v.Pt();
    double phi = v.Phi();
    if (phi < 0.0)
      phi += 2.0 * pi; // force from 0-2pi
    if (phi > 2.0 * pi)
      phi -= 2.0 * pi; // force from 0-2pi
    double eta = v.PseudoRapidity();

    if (std::find(fDroppedMCTracks.begin(), fDroppedMCTracks.end(), mc) != fDroppedMCTracks.end())
    {
      if (fPrintLevel)
        cout << "Dropped Track = " << McTrack_mGePid[mc] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
      continue; // Discard tracks which come from kaon decay after D0 decay
    }

    if (McTrack_mIdVtxStop[mc] != 0 && (McTrack_mGePid[mc] != 37 && McTrack_mGePid[mc] != 38))
      continue; // Only Final State Particles are included in JetMaker. D0s are the only exception

    if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38)
    {
      if (McTrack_mIdVtxStop[mc] != vertexids[D0num])
        continue; // Only consider the current D0
    }
    // Unstable particles which shouldn't make it to the end are discarded by hand. The list provisionally includes:
    /*
      Lambda, Eta, Sigma0, Xi0, Muon, Neutrino, KS0, KL0
    */

    if (McTrack_mGePid[mc] == 4 || McTrack_mGePid[mc] == 5 || McTrack_mGePid[mc] == 6 || McTrack_mGePid[mc] == 10 || McTrack_mGePid[mc] == 16 || McTrack_mGePid[mc] == 17 || McTrack_mGePid[mc] == 18 || McTrack_mGePid[mc] == 20 || McTrack_mGePid[mc] == 22)
      continue;

    // Have to discard Kaons and Pions which come from the Current D0

    if (McTrack_mIdVtxStart[mc] == vertexids[D0num])
    {
      continue;
    }

    if ((pt < fMinJetTrackPt) || (pt > fMaxJetTrackPt) || (eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax) || (phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax))
      continue;
    if (McTrack_mCharge[mc] != 0 || McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38)
    {
      fMCEventTracks[eventnum].push_back(v); // Neutral particles which are not D0 are saved in towers.
      if (fPrintLevel)
        cout << "Constituents = " << McTrack_mGePid[mc] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
    }
    else
    {
      fMCEventTowers[eventnum].push_back(v);
      if (fPrintLevel)
        cout << "Towers = " << McTrack_mGePid[mc] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
    }
  }
}

Bool_t StHIOverlay_Test::DoesItHaveAGoodRecoJet(int eventnum)
{
  // I am saving out the reconstructed jets from the PYTHIA reco events for comparison

  StFJWrapper *fjw = new StFJWrapper("HIRecoNoBG", "HIRecoNoBG"); //! fastjet wrapper
  fjw->Clear();
  fJets->Delete();

  // setup fj wrapper
  fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
  fastjet::JetAlgorithm algorithm = fastjet::antikt_algorithm;
  fastjet::Strategy strategy = fastjet::Best;

  fjw->SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw->SetStrategy(strategy);
  fjw->SetGhostArea(0.005);
  fjw->SetR(0.4);
  fjw->SetAlgorithm(algorithm);       // fJetAlgo);
  fjw->SetRecombScheme(recombScheme); // fRecombScheme);
  fjw->SetMaxRap(1.2);

  // Loop over all saved particles for MC in StHIOverlay_Test vectors and enter them into fastjet wrapper
  for (unsigned short iTracks = 0; iTracks < fRecoEventTracks[eventnum].size(); iTracks++)
  {

    TLorentzVector v;
    v = fRecoEventTracks[eventnum][iTracks];
    // track variables

    double px = v.X();
    double py = v.Y();
    double pz = v.Z();
    double p = v.P();
    double energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
    double mass = v.M();

    fjw->AddInputVector(px, py, pz, energy, iTracks); // includes E

    double pt = v.Pt();
    double phi = v.Phi();
    if (phi < 0.0)
      phi += 2.0 * pi; // force from 0-2pi
    if (phi > 2.0 * pi)
      phi -= 2.0 * pi; // force from 0-2pi
    double eta = v.Eta();

  } // track loop

  for (unsigned short iTowers = 0; iTowers < fRecoEventTowers[eventnum].size(); iTowers++)
  {

    TLorentzVector v;
    v = fRecoEventTowers[eventnum][iTowers];
    // tower variables

    double px = v.X();
    double py = v.Y();
    double pz = v.Z();
    double p = v.P();
    double energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
    double mass = v.M();

    fjw->AddInputVector(px, py, pz, energy, iTowers + 10000); // includes E

    double pt = v.Pt();
    double phi = v.Phi();
    if (phi < 0.0)
      phi += 2.0 * pi; // force from 0-2pi
    if (phi > 2.0 * pi)
      phi -= 2.0 * pi; // force from 0-2pi
    double eta = v.Eta();

  } // tower loop

  // Jet Running

  if (fPrintLevel > 1)
    fjw->PrintInput();

  fjw->Run();

  std::vector<fastjet::PseudoJet> jets_incl = fjw->GetInclusiveJets();

  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  int jetCount = 0; // This is the index of the D0 jet in the jet tree.
  // bool acceptedevent = kFALSE;

  for (UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet)
  {

    Int_t ij = indexes[ijet];

    // PERFORM CUTS ON inclusive JETS before saving
    // cut on min jet pt
    // if(jets_incl[ij].perp() < 5.0) continue;
    // cut on eta acceptance
    if ((jets_incl[ij].eta() < -1.0) || (jets_incl[ij].eta() > 1.0))
      continue;

    // fill jet constituents
    vector<fastjet::PseudoJet> constituents = fjw->GetJetConstituents(ij);

    bool D0Jet = kFALSE;

    double D0Pt, D0Eta, D0Phi;
    for (UInt_t ic = 0; ic < constituents.size(); ++ic)
    {
      // get user defined index
      Int_t uid = constituents[ic].user_index();

      if (int(uid) >= 0 && int(uid) < 10000)
      {
        TLorentzVector v;
        v = fRecoEventTracks[eventnum][uid];

        if (int(v.M()) == 1)
        {
          D0Jet = kTRUE;
          D0Pt = v.Pt();
          D0Eta = v.PseudoRapidity();
          D0Phi = v.Phi();
          break;
        }
      }
    }

    if (!D0Jet)
      continue;

    StJet *jet = new ((*fJets)[jetCount])
        StJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());

    jet->SetLabel(ij);

    // area vector and components
    fastjet::PseudoJet area(fjw->GetJetAreaVector(ij));
    jet->SetArea(area.perp()); // same as fjw->GetJetArea(ij)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    jet->SetJetConstituents(constituents);

    int nt = 0;
    int nc = 0;
    int ng = 0;

    jet->SetNumberOfTracks(constituents.size());
    jet->SetNumberOfTowers(constituents.size());

    for (UInt_t ic = 0; ic < constituents.size(); ++ic)
    {
      Int_t uid = constituents[ic].user_index();

      if (int(uid) >= 0 && int(uid) < 10000)
      {
        jet->AddTrackAt(abs(uid), nt);
        nt++;
      }
      else if (int(uid) >= 10000)
      {
        // cout << uid  << "\t" << nc << endl;
        jet->AddTowerAt(abs(uid), nc);
        nc++;
      }
      else if (uid < 0)
      {
        ng++;
      }
    }
    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfTowers(nc);

    if (fPrintLevel)
      cout << "Jet Found with pT eta phi = " << jet->Pt() << "\t" << jet->Eta() << "\t" << jet->Phi() << endl;
    if (fPrintLevel)
      cout << "D0 Found with pT eta phi = " << D0Pt << "\t" << D0Eta << "\t" << D0Phi << endl;

    break; // Once you encountered a D0 jet to your liking, break.
  }

  TClonesArray *tempjets = (TClonesArray *)fJets->Clone();
  fRecoJetsArr.push_back(tempjets);

  delete fjw;

  return kTRUE; // At this point, I am returning blanket true for this.
}

Bool_t StHIOverlay_Test::DoesItHaveAGoodMCJet(int eventnum)
{
  // Make the jets here. If they are not what we want them to be, discard the whole event and keep going.

  StFJWrapper *fjw = new StFJWrapper("HIMC", "HIMC"); //! fastjet wrapper
  fjw->Clear();
  fJets->Delete();

  // setup fj wrapper
  fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
  fastjet::JetAlgorithm algorithm = fastjet::antikt_algorithm;
  fastjet::Strategy strategy = fastjet::Best;

  fjw->SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw->SetStrategy(strategy);
  fjw->SetGhostArea(0.005);
  fjw->SetR(0.4);
  fjw->SetAlgorithm(algorithm);       // fJetAlgo);
  fjw->SetRecombScheme(recombScheme); // fRecombScheme);
  fjw->SetMaxRap(1.2);

  // Loop over all saved particles for MC in StHIOverlay_Test vectors and enter them into fastjet wrapper
  for (unsigned short iTracks = 0; iTracks < fMCEventTracks[eventnum].size(); iTracks++)
  {

    TLorentzVector v;
    v = fMCEventTracks[eventnum][iTracks];
    // track variables

    double px = v.X();
    double py = v.Y();
    double pz = v.Z();
    double p = v.P();
    double energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
    double mass = v.M();

    fjw->AddInputVector(px, py, pz, energy, iTracks); // includes E

    double pt = v.Pt();
    double phi = v.Phi();
    if (phi < 0.0)
      phi += 2.0 * pi; // force from 0-2pi
    if (phi > 2.0 * pi)
      phi -= 2.0 * pi; // force from 0-2pi
    double eta = v.Eta();

  } // track loop

  for (unsigned short iTowers = 0; iTowers < fMCEventTowers[eventnum].size(); iTowers++)
  {

    TLorentzVector v;
    v = fMCEventTowers[eventnum][iTowers];
    // tower variables

    double px = v.X();
    double py = v.Y();
    double pz = v.Z();
    double p = v.P();
    double energy = 1.0 * TMath::Sqrt(p * p + pi0mass * pi0mass);
    double mass = v.M();

    fjw->AddInputVector(px, py, pz, energy, iTowers + 10000); // includes E

    double pt = v.Pt();
    double phi = v.Phi();
    if (phi < 0.0)
      phi += 2.0 * pi; // force from 0-2pi
    if (phi > 2.0 * pi)
      phi -= 2.0 * pi; // force from 0-2pi
    double eta = v.Eta();

  } // tower loop

  // Jet Running

  if (fPrintLevel > 1)
    fjw->PrintInput();

  fjw->Run();

  std::vector<fastjet::PseudoJet> jets_incl = fjw->GetInclusiveJets();

  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  int jetCount = 0; // This is the index of the D0 jet in the jet tree.
  bool acceptedevent = kFALSE;

  for (UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet)
  {

    Int_t ij = indexes[ijet];

    // PERFORM CUTS ON inclusive JETS before saving
    // cut on min jet pt
    if (jets_incl[ij].perp() < fMCJetPtCut)
      continue;
    // cut on eta acceptance
    if ((jets_incl[ij].eta() < -0.6) || (jets_incl[ij].eta() > 0.6))
      continue;

    // fill jet constituents
    vector<fastjet::PseudoJet> constituents = fjw->GetJetConstituents(ij);

    bool D0Jet = kFALSE;

    double D0Pt, D0Eta, D0Phi;
    for (UInt_t ic = 0; ic < constituents.size(); ++ic)
    {
      // get user defined index
      Int_t uid = constituents[ic].user_index();

      if (int(uid) >= 0 && int(uid) < 10000)
      {
        TLorentzVector v;
        v = fMCEventTracks[eventnum][uid];

        if (int(v.M()) == 1)
        {
          D0Jet = kTRUE;
          D0Pt = v.Pt();
          D0Eta = v.PseudoRapidity();
          D0Phi = v.Phi();
          break;
        }
      }
    }

    if (!D0Jet)
      continue;

    StJet *jet = new ((*fJets)[jetCount])
        StJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());

    jet->SetLabel(ij);

    // area vector and components
    fastjet::PseudoJet area(fjw->GetJetAreaVector(ij));
    jet->SetArea(area.perp()); // same as fjw->GetJetArea(ij)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    jet->SetJetConstituents(constituents);

    Double_t MClambda_1_1 = jet->GetAngularityLambda(1., 1.) / 0.4;
    cout << "MClambda_1_1 = " << MClambda_1_1 << endl;

    int nt = 0;
    int nc = 0;
    int ng = 0;

    jet->SetNumberOfTracks(constituents.size());
    jet->SetNumberOfTowers(constituents.size());

    for (UInt_t ic = 0; ic < constituents.size(); ++ic)
    {
      Int_t uid = constituents[ic].user_index();

      if (int(uid) >= 0 && int(uid) < 10000)
      {
        jet->AddTrackAt(abs(uid), nt);
        nt++;
      }
      else if (int(uid) >= 10000)
      {
        // cout << uid  << "\t" << nc << endl;
        jet->AddTowerAt(abs(uid), nc);
        nc++;
      }
      else if (uid < 0)
      {
        ng++;
      }
    }
    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfTowers(nc);

    TClonesArray *tempjets = (TClonesArray *)fJets->Clone();
    fJetsArr.push_back(tempjets);

    acceptedevent = kTRUE;

    if (fPrintLevel)
      cout << "Jet Found with pT eta phi = " << jet->Pt() << "\t" << jet->Eta() << "\t" << jet->Phi() << endl;
    if (fPrintLevel)
      cout << "D0 Found with pT eta phi = " << D0Pt << "\t" << D0Eta << "\t" << D0Phi << endl;

    break; // Once you encountered a D0 jet to your liking, break.
  }

  delete fjw;

  return acceptedevent;
}

map<int, vector<int>> StHIOverlay_Test::GetMCTrackListForVertex()
{

  map<int, vector<int>> arr;

  for (int mc = 0; mc < McTrack_; mc++)
  {
    int idvxstart = McTrack_mIdVtxStart[mc] - 1;
    arr[idvxstart].push_back(mc);
  }

  return arr;
}

void StHIOverlay_Test::GetAllTracksFromVertex(int vertexid, vector<int> &trackvec)
{
  if (vertexid < 0)
    return;

  if (fPrintLevel > 1)
  {
    cout << "Called this function for vx " << vertexid << " with ntracks = " << fvertextotrack[vertexid].size() << endl;
    for (int track = 0; track < fvertextotrack[vertexid].size(); track++)
    {
      cout << "Geant ID of tracks = " << McTrack_mGePid[fvertextotrack[vertexid][track]] << endl;
    }
  }

  for (int track = 0; track < fvertextotrack[vertexid].size(); track++)
  {
    trackvec.push_back(fvertextotrack[vertexid][track]);

    int idvxstop = McTrack_mIdVtxStop[fvertextotrack[vertexid][track]];
    // if (idvxstop != 0) {
    GetAllTracksFromVertex(idvxstop - 1, trackvec);
    // }
  }

  return;
}

void StHIOverlay_Test::GetD0AndDaughters()
{
  vertexids.clear();
  pionids.clear();
  kaonids.clear();

  matchedpionids.clear();
  matchedkaonids.clear();

  // Saving the D0 vertices
  for (int mc = 0; mc < McTrack_; mc++)
  {
    if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38)
    {
      TLorentzVector v;
      v.SetPxPyPzE(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], McTrack_mE[mc]);
      // track variables
      double pt = v.Pt();
      double phi = v.Phi();
      if (phi < 0.0)
        phi += 2.0 * pi; // force from 0-2pi
      if (phi > 2.0 * pi)
        phi -= 2.0 * pi; // force from 0-2pi
      double eta = v.PseudoRapidity();
      if (pt < 1.0 || pt > 10.0)
        continue; // Only D0s we care about need to have pT > 1 GeV
      if ((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax) || (phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax))
        continue;

      if (fPrintLevel)
        cout << "Number of D0 Daugs = " << McVertex_mNoDaughters[McTrack_mIdVtxStop[mc] - 1] << endl;
      if (fPrintLevel)
        cout << "MC D0 Found = " << pt << "\t" << eta << "\t" << phi << endl;

      bool goodvertex = kTRUE;

      int pionid = -99;
      int kaonid = -99;

      vector<int> tmp;
      tmp.clear();

      for (int daug = 0; daug < McTrack_; daug++)
      {
        if (McTrack_mIdVtxStart[daug] != McTrack_mIdVtxStop[mc])
          continue; // We only want the kaon and pion that originated from the D0 we are interested in.

        TLorentzVector v2;
        v2.SetPxPyPzE(McTrack_mPx[daug], McTrack_mPy[daug], McTrack_mPz[daug], McTrack_mE[daug]);
        double pt2 = v2.Pt();
        double phi2 = v2.Phi();
        if (phi2 < 0.0)
          phi2 += 2.0 * pi; // force from 0-2pi
        if (phi2 > 2.0 * pi)
          phi2 -= 2.0 * pi; // force from 0-2pi
        double eta2 = v2.PseudoRapidity();

        if (McTrack_mGePid[daug] != 8 && McTrack_mGePid[daug] != 9 && McTrack_mGePid[daug] != 11 && McTrack_mGePid[daug] != 12)
        {
          goodvertex = kFALSE;
          break; // This should never happen
        }

        // Push the pion and kaon ids into the vector
        if (McTrack_mGePid[daug] == 8 || McTrack_mGePid[daug] == 9)
        {
          pionid = daug;
        }

        else if (McTrack_mGePid[daug] == 11 || McTrack_mGePid[daug] == 12)
        {
          kaonid = daug;
        }
      }

      if (goodvertex)
      {
        if (pionid == -99 || kaonid == -99)
        {
          cout << "Something wrong with D0. Exiting." << endl;
          continue;
        }
        vertexids.push_back(McTrack_mIdVtxStop[mc]);
        pionids.push_back(pionid);
        kaonids.push_back(kaonid);
        matchedpionids.push_back(GetMatchedRecoTrackFromMCTrack(pionid));
        matchedkaonids.push_back(GetMatchedRecoTrackFromMCTrack(kaonid));
      }

      // cout << "Number of D0 Daugs = " << McVertex_mNoDaughters[McTrack_mIdVtxStop[mc]-1] << endl;
      // cout << "MC D0 Found = " << pt << "\t" << eta << "\t" << phi << endl;
      // cout << "D0 Daugs = " << vertexids.back() << "\t" << pionids.back() << "\t" << kaonids.back() << endl;
    }
  }

  assert((pionids.size() == kaonids.size()) && "Same number of kaons and pions \n");
}

/* Returns a matched reco track with our acceptance requirements for an MC track. If there is no matching,
it returns -999. */

int StHIOverlay_Test::GetMatchedRecoTrackFromMCTrack(int mctrkid)
{
  int recotrackmatch = -999;

  TVector3 oVertex;
  oVertex.SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

  double ratio = -99.;

  for (int reco = 0; reco < Track_; reco++)
  {
    // get track pointer
    TVector3 o;
    o.SetXYZ(Track_mOriginX[reco], Track_mOriginY[reco], Track_mOriginZ[reco]);

    double dca = (oVertex - o).Mag();
    int nHitsFit = abs(Track_mNHitsFit[reco]);
    int nHitsMax = abs(Track_mNHitsMax[reco]);
    double nHitsRatio = 1.0 * nHitsFit / nHitsMax;

    // additional quality cuts for tracks
    if (dca > fJetTrackDCAcut)
      continue;
    if (nHitsFit < fTracknHitsFit)
      continue;
    if (nHitsRatio < fTracknHitsRatio)
      continue;

    int mctrk = Track_mIdTruth[reco] - 1;

    if (mctrk == mctrkid)
    {
      if (nHitsRatio > ratio)
      {
        recotrackmatch = reco;
        ratio = nHitsRatio;
      }
    }
  }

  return recotrackmatch;
}

void StHIOverlay_Test::AssignTreeVariables()
{

  TRandom3 *r1 = new TRandom3(0);
  int filenumber = r1->Integer(filenamesforHIOverlay.size());

  f = new TFile(filenamesforHIOverlay[filenumber].Data());

  fMCPico = (TTree *)f->Get("PicoDst");

  int nentries = fMCPico->GetEntriesFast();

  if (fPrintLevel)
    cout << filenumber << "\t" << filenamesforHIOverlay[filenumber] << "\t" << nentries << "\t" << mCentMaker->GetWeight() << endl;
  // if (fPrintLevel)

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

  // fMCPico->SetBranchStatus("Event", true);

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

// vector <int> StHIOverlay_Test::TracksToToss(int stopvertex, int *)

int StHIOverlay_Test::GetMatchedBtowID(int trkbemcid, TVector3 gMom, TVector3 org, int charge)
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

  // int trkbemcid = trk->bemcTowerIndex();

  // Check if the track can be projected onto the current radius
  // if not, track can't be matched.
  // By JetCorr request the global track projection to BEMC is used.
  if (mEmcPosition->projTrack(&bemc_pos, &bemc_mom, gMom, org, Bfield, charge, bemc_radius))
  {
    // First, examine track eta. If it falls in two regions:
    // 0 < |eta| < etaMin()
    // etaMax() < |eta| < 1.0
    // then shift the eta for the projection slightly into the neighboring tower
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(org, trkbemcid + 1);
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

TVector3 StHIOverlay_Test::FastSimMom(TVector3 p, int pid)
{
  float pt = p.Perp();
  float pt1 = pt;
  // if(pt1>2) pt1 = 2;//Used for high pt-hat bin smearing test
  if (pt1 > 10)
    pt1 = 10; // Used for high pt-hat bin smearing test
  float sPt = -1;

  TRandom3 *r = new TRandom3(0);

  if (pid == 8 || pid == 9)
    sPt = r->Gaus(pt, pt * fPionMomResolution->Eval(pt1)); // Pion
  else if (pid == 11 || pid == 12)
    sPt = r->Gaus(pt, pt * fKaonMomResolution->Eval(pt1)); // Kaon
  else if (pid == 15 || pid == 14)
    sPt = r->Gaus(pt, pt * fProtonMomResolution->Eval(pt1)); // Proton
  else
    sPt = r->Gaus(pt, pt * fPionMomResolution->Eval(pt1)); // Catch all: pions

  TVector3 smearedmom(sPt * cos(p.Phi()), sPt * sin(p.Phi()), sPt * sinh(p.PseudoRapidity()));

  return smearedmom;
}

Bool_t StHIOverlay_Test::IsItASoftEvent()
{
  bool softevent = kTRUE;

  const Int_t ntracks = mPicoDst->numberOfTracks();

  for (unsigned short itrk = 0; itrk < ntracks; itrk++)
  {
    StPicoTrack *trk = static_cast<StPicoTrack *>(mPicoDst->track(itrk));
    if (!trk)
    {
      continue;
    }

    TVector3 mTrkMom;
    mTrkMom = trk->gMom(mVertex, Bfield);

    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    if (abs(eta) > 1.)
      continue;
    if (pt > 5.)
    {
      softevent = kFALSE;
      break;
    }
  }

  return softevent;
}

Bool_t StHIOverlay_Test::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
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
