// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StD0Analysis.h"
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
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"

// KFParticle includes
#include "StRoot/StEvent/StDcaGeometry.h"
#include "StRoot/StarRoot/KFParticle.h"

// Bichsel includes
#include "StBichsel/Bichsel.h"

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

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StD0Analysis)

//________________________________________________________________________
StD0Analysis::StD0Analysis(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kTRUE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StD0Analysis::fRunFlagEnum
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

  for(int i=0; i<10; i++) {numberofevents[i] = 0; numberoftracks[i] = 0;}
  //Mass cut
  fInvMassSignal1 = 1.84;
  fInvMassSignal2 = 1.89;

  fInvMassULBg1 = 1.70;
  fInvMassULBg2 = 2.08;

  fInvMassLSBg1 = 1.70;
  fInvMassLSBg2 = 2.08;

}

//
//________________________________________________________________________
StD0Analysis::~StD0Analysis()
{ 
  if (hCentrality) delete hCentrality;
  if (hMultiplicity) delete hMultiplicity;
  if (cuthistogram_event) delete cuthistogram_event;

  if (decaylengthd0US) delete decaylengthd0US;
  if (distancepikUS) delete distancepikUS;
  if (distanced0PVUS) delete distanced0PVUS;
  if (dcakPVUS) delete dcakPVUS;
  if (dcapiPVUS) delete dcapiPVUS;

  if (decaylengthd0LS) delete decaylengthd0LS;
  if (distancepikLS) delete distancepikLS;
  if (distanced0PVLS) delete distanced0PVLS;
  if (dcakPVLS) delete dcakPVLS;
  if (dcapiPVLS) delete dcapiPVLS;

  if (hnjets) delete hnjets;
  if (hNumberofTracks) delete hNumberofTracks;
  if (hNumberofHFTTracks) delete hNumberofHFTTracks;
  if (hRatioOfAcceptedTracks) delete hRatioOfAcceptedTracks;
  if (hJetPt) delete hJetPt;
  if (hJetPhi) delete hJetPhi;
  if (hJetEta) delete hJetEta;
  if (hConstPt) delete hConstPt;
  if (hConstPhi) delete hConstPhi;
  if (hConstEta) delete hConstEta;

  if (kaonpt) delete kaonpt;
  if (pionpt) delete pionpt;
  if (d0pt) delete d0pt;
  if (kaonpionpt) delete kaonpionpt;

  if (kaonbgpt) delete kaonbgpt;
  if (pionbgpt) delete pionbgpt;
  if (d0bgpt) delete d0bgpt;
  if (kaonpionbgpt) delete kaonpionbgpt;

  for (int topo_bin = 0; topo_bin < 3; topo_bin++){
    for (int jetptbin = 0; jetptbin < 3; jetptbin++){
      for (int ptbin = 0; ptbin < 7; ptbin++){
        if (invmass_ptbin[topo_bin][jetptbin][ptbin]) delete invmass_ptbin[topo_bin][jetptbin][ptbin];
        if (invmassbg_ptbin[topo_bin][jetptbin][ptbin]) delete invmassbg_ptbin[topo_bin][jetptbin][ptbin];
      }
    }
  }

  for (int topo_bin = 0; topo_bin < 3; topo_bin++){
    if (invmass[topo_bin]) delete invmass[topo_bin];
    if (invmassbg[topo_bin]) delete invmassbg[topo_bin];
  }


  if(mEmcPosition) delete mEmcPosition;

}

//________________________________________________________________________
Int_t StD0Analysis::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StD0Analysis::Finish() { 

  cout << "StD0Analysis::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "RECREATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StD0Analysis::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StD0Analysis::DeclareHistograms() {

  //gPtvpPt = new TH2F("gPtvpPt", "gPtvpPt", 1000, 0, 20, 1000, 0, 20);

  // binning for cent histograms
  int nHistCentBins = 20;

  // binning for mult histograms
  double kHistMultMax = 800.;
  int kHistMultBins = 400;

  // pp specific settings
  if(doppAnalysis) {
    kHistMultMax = 100.;
    kHistMultBins = 100.;
  }

  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);

  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  cuthistogram_event = new TH1F("cuthistogram_event", "cuthistogram_event", 10, 0, 10);

  // cuthistogram_track = new TH1F("cuthistogram_track", "cuthistogram_track", 10, 0, 10);

  // z_pi = new TH2F("z_pi", "z_pi", 2000, 0, 20, 2000, -10, 10);

  // z_ka = new TH2F("z_ka", "z_ka", 2000, 0, 20, 2000, -10, 10);

  // normalised_invbetavpT_tof_pi = new TH2F("normalised_invbetavpT_tof_pi", "normalised_invbetavpT_tof_pi", 2000, 0, 20, 2000, -10, 10);

  // normalised_invbetavpT_tof_ka = new TH2F("normalised_invbetavpT_tof_ka", "normalised_invbetavpT_tof_ka", 2000, 0, 20, 2000, -10, 10);

  // invmass = new TH1F("invmass", "invmass", 1000, 0.5, 2.5);
  // invmassbg = new TH1F("invmassbg", "invmassbg", 1000, 0.5, 2.5);
  // kaonpt = new TH1F("kaonpt", "kaonpt", 200, 0, 50);
  // pionpt = new TH1F("pionpt", "pionpt", 200, 0, 50);
  // kaon_pion_pt = new TH2F("kaon_pion_pt", "kaon_pion_pt", 200, 0, 100, 200, 0, 100);
  // d0pt = new TH1F("d0pt", "d0pt", 200, 0, 50);
  // kaonbgpt = new TH1F("kaonbgpt", "kaonbgpt", 200, 0, 50);
  // pionbgpt = new TH1F("pionbgpt", "pionbgpt", 200, 0, 50);
  // kaon_pion_bg_pt = new TH2F("kaon_pion_bg_pt", "kaon_pion_bg_pt", 200, 0, 100, 200, 0, 100);
  // d0bgpt = new TH1F("d0bgpt", "d0bgpt", 200, 0, 50);

  // invmass_wholeevent = new TH1F("invmass_wholeevent", "invmass_wholeevent", 1000, 0.5, 2.5);
  // invmassbg_wholeevent = new TH1F("invmassbg_wholeevent", "invmassbg_wholeevent", 1000, 0.5, 2.5);
  // kaonpt_wholeevent = new TH1F("kaonpt_wholeevent", "kaonpt_wholeevent", 200, 0, 50);
  // pionpt_wholeevent = new TH1F("pionpt_wholeevent", "pionpt_wholeevent", 200, 0, 50);
  // kaon_pion_pt_wholeevent = new TH2F("kaon_pion_pt_wholeevent", "kaon_pion_pt_wholeevent", 200, 0, 100, 200, 0, 100);
  // d0pt_wholeevent = new TH1F("d0pt_wholeevent", "d0pt_wholeevent", 200, 0, 50);
  // kaonbgpt_wholeevent = new TH1F("kaonbgpt_wholeevent", "kaonbgpt_wholeevent", 200, 0, 50);
  // pionbgpt_wholeevent = new TH1F("pionbgpt_wholeevent", "pionbgpt_wholeevent", 200, 0, 50);
  // kaon_pion_bg_pt_wholeevent = new TH2F("kaon_pion_bg_pt_wholeevent", "kaon_pion_bg_pt_wholeevent", 200, 0, 100, 200, 0, 100);
  // d0bgpt_wholeevent = new TH1F("d0bgpt_wholeevent", "d0bgpt_wholeevent", 200, 0, 50);

  // hptdistro = new TH1F("hptdistro", "hptdistro", 2000, 0, 100); // Histogram for Delta Phi

  // hetadistro = new TH1F("hetadistro", "hetadistro", 2000, -5, 5); // Histogram for Delta Phi

  // invmassvpt = new TH2F("invmassvpt", "invmassvpt", 2000, 0, 20, 1000, 0.5, 2.5);


  // for (int jetptbin = 0; jetptbin < 3; jetptbin++){
  //   for (int ptbin = 0; ptbin < numberofptbins; ptbin++){
  //     invmass_ptbin[jetptbin][ptbin] = new TH1F(Form("invmass_ptbin_%i_%i", jetptbin, ptbin), Form("invmass_ptbin_%i_%i", jetptbin, ptbin), 1000, 0.5, 2.5);
  //     invmassbg_ptbin[jetptbin][ptbin] = new TH1F(Form("invmassbg_ptbin_%i_%i", jetptbin, ptbin), Form("invmassbg_ptbin_%i_%i", jetptbin, ptbin), 1000, 0.5, 2.5);
  //   }
  // }

  // Topological histograms



  decaylengthd0US = new TH1F("decaylengthd0US", "decaylengthd0US", 10000, 0, 1);
  distancepikUS = new TH1F("distancepikUS", "distancepikUS", 10000, 0, 1);
  distanced0PVUS = new TH1F("distanced0PVUS", "distanced0PVUS", 10000, 0, 1);
  dcakPVUS = new TH1F("dcakPVUS", "dcakPVUS", 10000, 0, 1);
  dcapiPVUS = new TH1F("dcapiPVUS", "dcapiPVUS", 10000, 0, 1);

  decaylengthd0LS = new TH1F("decaylengthd0LS", "decaylengthd0LS", 10000, 0, 1);
  distancepikLS = new TH1F("distancepikLS", "distancepikLS", 10000, 0, 1);
  distanced0PVLS = new TH1F("distanced0PVLS", "distanced0PVLS", 10000, 0, 1);
  dcakPVLS = new TH1F("dcakPVLS", "dcakPVLS", 10000, 0, 1);
  dcapiPVLS = new TH1F("dcapiPVLS", "dcapiPVLS", 10000, 0, 1);

  // hD0US_pT = new TH1F("hD0US_pT", "hD0US_pT", 1000, 0, 50);
  // hD0LS_pT = new TH1F("hD0LS_pT", "hD0LS_pT", 1000, 0, 50);


  // To check that jet pt, eta, phi distributions are consistent
  hnjets = new TH1F("hnjets", "hnjets", 20, -0.5, 19.5);
  hJetPt = new TH1F("hJetPt", "hJetPt", 1000, 0, 50); // Histogram for Delta Phi
  hJetPtCorr = new TH1F("hJetPtCorr", "hJetPtCorr", 1000, 0, 50); // Histogram for Delta Phi
  hJetEta = new TH1F("hJetEta", "hJetEta", 1000, -2, 2); // Histogram for Delta Phi
  hJetPhi = new TH1F("hJetPhi", "hJetPhi", 1000, -10, 10); // Histogram for Delta Phi
  hNumberofTracks = new TH1F("hNumberofTracks", "hNumberofTracks", 51, -0.5, 49.5);
  hNumberofHFTTracks = new TH1F("hNumberofHFTTracks", "hNumberofHFTTracks", 51, -0.5, 49.5);
  hRatioOfAcceptedTracks = new TH1F("hRatioOfAcceptedTracks", "hRatioOfAcceptedTracks", 102, -0.01, 1.01);
  hConstPt = new TH1F("hConstPt", "hConstPt", 1000, 0, 30); // Histogram for Delta Phi
  hConstPhi = new TH1F("hConstPhi", "hConstPhi", 1000, -10, 10); // Histogram for Delta Phi
  hConstEta = new TH1F("hConstEta", "hConstEta", 1000, -2, 2); // Histogram for Delta Phi


  //Invariant Mass Plots

  kaonpt = new TH1F("kaonpt", "kaonpt", 1000, 0, 30);
  pionpt = new TH1F("pionpt", "pionpt", 1000, 0, 30);
  d0pt = new TH1F("d0pt", "d0pt", 1000, 0, 30);
  kaonpionpt = new TH2F("kaonpionpt", "kaonpionpt", 1000, 0, 30, 1000, 0, 30);
  
  kaonbgpt = new TH1F("kaonbgpt", "kaonbgpt", 1000, 0, 30);
  pionbgpt = new TH1F("pionbgpt", "pionbgpt", 1000, 0, 30);
  d0bgpt = new TH1F("d0bgpt", "d0bgpt", 1000, 0, 30);
  kaonpionbgpt = new TH2F("kaonpionbgpt", "kaonpionbgpt", 1000, 0, 30, 1000, 0, 30);


  for (int topo_bin = 0; topo_bin < 3; topo_bin++){
    for (int jetptbin = 0; jetptbin < 3; jetptbin++){
      for (int ptbin = 0; ptbin < 7; ptbin++){
        invmass_ptbin[topo_bin][jetptbin][ptbin]  = new TH1F(Form("invmass_topo%i_jet%i_d%i", topo_bin, jetptbin, ptbin), Form("invmass_topo%i_jet%i_d%i", topo_bin, jetptbin, ptbin), 1000, 0.5, 2.5);
        invmassbg_ptbin[topo_bin][jetptbin][ptbin] = new TH1F(Form("invmassbg_topo%i_jet%i_d%i", topo_bin, jetptbin, ptbin), Form("invmassbg_topo%i_jet%i_d%i", topo_bin, jetptbin, ptbin), 1000, 0.5, 2.5);
      }
    }
  }

  for (int topo_bin = 0; topo_bin < 3; topo_bin++){
    invmass[topo_bin] = new TH1F(Form("invmass%i", topo_bin), Form("invmass%i", topo_bin), 1000, 0.5, 2.5);
    invmassbg[topo_bin] = new TH1F(Form("invmassbg%i", topo_bin), Form("invmassbg%i", topo_bin), 1000, 0.5, 2.5);
  }
  // Jetshape histograms 

  // for (int cbin = 0; cbin < 4; cbin++){

  //   JetPhi[cbin] = new TH1F(Form("JetPhi %i", cbin),Form("JetPhi %i", cbin), 2000, -10, 10); // Histogram for Delta Pi

  //   ConstPhi[cbin] = new TH1F(Form("ConstPhi %i", cbin),Form("ConstPhi %i", cbin), 200, -10, 10); // Histogram for Delta Pi

  //   JetPt[cbin] = new TH1F(Form("JetPt %i", cbin),Form("JetPt %i", cbin), 2000, 0, 100); // Histogram for Delta Pi

  //   JetPtCorr[cbin] = new TH1F(Form("JetPtCorr %i", cbin),Form("JetPtCorr %i", cbin), 2000, 0, 100); // Histogram for Delta Pi

  //   ConstPt[cbin] = new TH1F(Form("ConstPt %i", cbin),Form("ConstPt %i", cbin), 2000, 0, 100); // Histogram for Delta Pi

  //   DeltaEta[cbin] = new TH1F(Form("DeltaEta %i", cbin),Form("DeltaEta %i", cbin), 200, -1, 1); // Histogram for Delta Eta

  //   DeltaPhi[cbin] = new TH1F(Form("DeltaPhi %i", cbin),Form("DeltaPhi %i", cbin), 1000, -10, 10); // Histogram for Delta Pi\

  //   Rad[cbin] = new TH1F(Form("Rad %i", cbin),Form("Rad %i", cbin), 100, 0, 2); // Histogram for dR


  //   hJetPhid0[cbin] = new TH1F(Form("hJetPhid0 %i", cbin),Form("hJetPhid0 %i", cbin), 2000, -10, 10); // Histogram for Delta Phi

  //   hConstPhid0[cbin] = new TH1F(Form("hConstPhid0 %i", cbin),Form("hConstPhid0 %i", cbin), 200, -10, 10); // Histogram for Delta Phi

  //   hJetPtd0[cbin] = new TH1F(Form("hJetPtd0 %i", cbin),Form("hJetPtd0 %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hJetPtd0Corr[cbin] = new TH1F(Form("hJetPtd0Corr %i", cbin),Form("hJetPtd0Corr %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hConstPtd0[cbin] = new TH1F(Form("hConstPtd0 %i", cbin),Form("hConstPtd0 %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hDeltaEtad0[cbin] = new TH1F(Form("hDeltaEtad0 %i", cbin),Form("hDeltaEtad0 %i", cbin), 200, -1, 1); // Histogram for Delta Eta

  //   hDeltaPhid0[cbin] = new TH1F(Form("hDeltaPhid0 %i", cbin),Form("hDeltaPhid0 %i", cbin), 1000, -10, 10); // Histogram for Delta Phi\

  //   hRadd0[cbin] = new TH1F(Form("hRadd0 %i", cbin),Form("hRadd0 %i", cbin), 100, 0, 2); // Histogram for dR


  //   hJetPhid0BgUL[cbin] = new TH1F(Form("hJetPhid0BgUL %i", cbin),Form("hJetPhid0BgUL %i", cbin), 2000, -10, 10); // Histogram for Delta Phi

  //   hConstPhid0BgUL[cbin] = new TH1F(Form("hConstPhid0BgUL %i", cbin),Form("hConstPhid0BgUL %i", cbin), 200, -10, 10); // Histogram for Delta Phi

  //   hJetPtd0BgUL[cbin] = new TH1F(Form("hJetPtd0BgUL %i", cbin),Form("hJetPtd0BgUL %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hJetPtd0BgULCorr[cbin] = new TH1F(Form("hJetPtd0BgULCorr %i", cbin),Form("hJetPtd0BgULCorr %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hConstPtd0BgUL[cbin] = new TH1F(Form("hConstPtd0BgUL %i", cbin),Form("hConstPtd0BgUL %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hDeltaEtad0BgUL[cbin] = new TH1F(Form("hDeltaEtad0BgUL %i", cbin),Form("hDeltaEtad0BgUL %i", cbin), 200, -1, 1); // Histogram for Delta Eta

  //   hDeltaPhid0BgUL[cbin] = new TH1F(Form("hDeltaPhid0BgUL %i", cbin),Form("hDeltaPhid0BgUL %i", cbin), 1000, -10, 10); // Histogram for Delta Phi\

  //   hRadd0BgUL[cbin] = new TH1F(Form("hRadd0BgUL %i", cbin),Form("hRadd0BgUL %i", cbin), 100, 0, 2); // Histogram for dR


  //   hJetPhid0BgLS[cbin] = new TH1F(Form("hJetPhid0BgLS %i", cbin),Form("hJetPhid0BgLS %i", cbin), 2000, -10, 10); // Histogram for Delta Phi

  //   hConstPhid0BgLS[cbin] = new TH1F(Form("hConstPhid0BgLS %i", cbin),Form("hConstPhid0BgLS %i", cbin), 200, -10, 10); // Histogram for Delta Phi

  //   hJetPtd0BgLS[cbin] = new TH1F(Form("hJetPtd0BgLS %i", cbin),Form("hJetPtd0BgLS %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hJetPtd0BgLSCorr[cbin] = new TH1F(Form("hJetPtd0BgLSCorr %i", cbin),Form("hJetPtd0BgLSCorr %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hConstPtd0BgLS[cbin] = new TH1F(Form("hConstPtd0BgLS %i", cbin),Form("hConstPtd0BgLS %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

  //   hDeltaEtad0BgLS[cbin] = new TH1F(Form("hDeltaEtad0BgLS %i", cbin),Form("hDeltaEtad0BgLS %i", cbin), 200, -1, 1); // Histogram for Delta Eta

  //   hDeltaPhid0BgLS[cbin] = new TH1F(Form("hDeltaPhid0BgLS %i", cbin),Form("hDeltaPhid0BgLS %i", cbin), 1000, -10, 10); // Histogram for Delta Phi\
    
  //   hRadd0BgLS[cbin] = new TH1F(Form("hRadd0BgLS %i", cbin),Form("hRadd0BgLS %i", cbin), 100, 0, 2); // Histogram for dR


  //   hJetShape[cbin] = new TH1F(Form("hJetShape %i", cbin),Form("hJetShape %i", cbin), numberofbins, 0, R);

  //   hJetShaped0[cbin] = new TH1F(Form("hJetShaped0 %i", cbin),Form("hJetShaped0 %i", cbin), numberofbins, 0, R);

  //   hJetShaped0BgUL[cbin] = new TH1F(Form("hJetShaped0BgUL %i", cbin),Form("hJetShaped0BgUL %i", cbin), numberofbins, 0, R);

  //   hJetShaped0BgLS[cbin] = new TH1F(Form("hJetShaped0BgLS %i", cbin),Form("hJetShaped0BgLS %i", cbin), numberofbins, 0, R);

  //   hJetDist[cbin] = new TH1F(Form("hJetDist %i", cbin),Form("hJetDist %i", cbin), numberofbins, 0, R);

  //   hJetDistd0[cbin] = new TH1F(Form("hJetDistd0 %i", cbin),Form("hJetDistd0 %i", cbin), numberofbins, 0, R);

  //   hJetDistd0BgUL[cbin] = new TH1F(Form("hJetDistd0BgUL %i", cbin),Form("hJetDistd0BgUL %i", cbin), numberofbins, 0, R);

  //   hJetDistd0BgLS[cbin] = new TH1F(Form("hJetDistd0BgLS %i", cbin),Form("hJetDistd0BgLS %i", cbin), numberofbins, 0, R);
  // }
}
//
// write histograms
//_____________________________________________________________________________
void StD0Analysis::WriteHistograms() {

  hCentrality->Write();
  hMultiplicity->Write();
  const char *event_cuts[10] = {"Total Events", "MinBias","|V_{z}| < 6 cm", "|V_{z} - V_{z(VPD)}| < 3 cm", "nBEMCMatch > 0", "nBTOFMatch > 0", "", "", "", ""};
  const char *track_cuts[10] = {"Total Tracks", "nHitsDedx > 20", "#frac{nHitsDedx}{nHitsMax} > 0.52", "DCA < 3 cm", "BEMC Hit", "BTOF Hit", "Charged Tracks", "", "", ""};

  for (int i=1; i <= 10; i++){
    cuthistogram_event->SetBinContent(i, numberofevents[i-1]);
    cuthistogram_event->GetXaxis()->SetBinLabel(i, event_cuts[i-1]);
    // cuthistogram_track->SetBinContent(i, numberoftracks[i-1]);
    // cuthistogram_track->GetXaxis()->SetBinLabel(i, track_cuts[i-1]);
  }
  cuthistogram_event->Write();
  // cuthistogram_track->Write();

  // z_pi->Write();
  // z_ka->Write();

  // normalised_invbetavpT_tof_pi->Write();
  // normalised_invbetavpT_tof_ka->Write();

  decaylengthd0US->Write();
  distancepikUS->Write();
  distanced0PVUS->Write();
  dcakPVUS->Write();
  dcapiPVUS->Write();

  decaylengthd0LS->Write();
  distancepikLS->Write();
  distanced0PVLS->Write();
  dcakPVLS->Write();
  dcapiPVLS->Write();

  // hD0US_pT->Write();
  // hD0LS_pT->Write();

  // invmass->Write();
  // invmassbg->Write();

  // kaonpt->Write();
  // pionpt->Write();
  // kaon_pion_pt->Write();
  // d0pt->Write();

  // kaonbgpt->Write();
  // pionbgpt->Write();
  // kaon_pion_bg_pt->Write();
  // d0bgpt->Write();

  // invmass_wholeevent->Write();
  // invmassbg_wholeevent->Write();

  // kaonpt_wholeevent->Write();
  // pionpt_wholeevent->Write();
  // kaon_pion_pt_wholeevent->Write();
  // d0pt_wholeevent->Write();

  // kaonbgpt_wholeevent->Write();
  // pionbgpt_wholeevent->Write();
  // kaon_pion_bg_pt_wholeevent->Write();
  // d0bgpt_wholeevent->Write();

  // for (int jetptbin = 0; jetptbin < 3; jetptbin++){
  //   for (int ptbin = 0; ptbin < numberofptbins; ptbin++){
  //     invmass_ptbin[jetptbin][ptbin]->Write();
  //     invmassbg_ptbin[jetptbin][ptbin]->Write();
  //   }
  // }

  hnjets->Write();
  hNumberofTracks->Write();
  hNumberofHFTTracks->Write();
  hRatioOfAcceptedTracks->Write();
  hJetPt->Write();
  hJetPtCorr->Write();
  hJetPhi->Write();
  hJetEta->Write();
  hConstPt->Write();
  hConstEta->Write();
  hConstPhi->Write();

  kaonpt->Write();
  pionpt->Write();
  d0pt->Write();
  kaonpionpt->Write();

  kaonbgpt->Write();
  pionbgpt->Write();
  d0bgpt->Write();
  kaonpionbgpt->Write();

  for (int topo_bin = 0; topo_bin < 3; topo_bin++){
    for (int jetptbin = 0; jetptbin < 3; jetptbin++){
      for (int ptbin = 0; ptbin < 7; ptbin++){
        invmass_ptbin[topo_bin][jetptbin][ptbin]->Write();
        invmassbg_ptbin[topo_bin][jetptbin][ptbin]->Write();
      }
    }
  }

  for (int topo_bin = 0; topo_bin < 3; topo_bin++){
    invmass[topo_bin]->Write();
    invmassbg[topo_bin]->Write();
  }


}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StD0Analysis::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StD0Analysis::Make() {
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

  // cout << "B = " << Bfield << endl;
  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  double zVtx_VPD = mPicoEvent->vzVpd();

  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


  // HFT Analysis Specific Cuts For Event


  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CentMaker! Skip! " << endm;
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

  //cout << "The centrality bin is " << fCentralityScaled << endl;
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // fill histograms
  hCentrality->Fill(fCentralityScaled);
  hMultiplicity->Fill(refCorr2);

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  //  FillEmcTriggers();

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
  }
  if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<endl;

  numberofevents[0]++;

  // check for MB/HT event
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);

  bool fRunForMB = kFALSE;
  if(doppAnalysis) fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
  if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

  bool match = kFALSE;

  int arrMB5_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};

  for(int i = 0; i < sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14); i++) {
    if(mPicoEvent->isTrigger(arrMB5_Run14[i])) match = kTRUE;
    if(match) break;
  }

  if (!match) return kStOk;

  numberofevents[1]++;
  // ======================== end of Triggers ============================= //

  

  if (abs(zVtx) > 6) return kStOk;
  numberofevents[2]++;

  if (abs(zVtx - zVtx_VPD) > 3) return kStOk;
  numberofevents[3]++;
   
  if (mPicoEvent->nBEMCMatch() == 0) return kStOk;
  numberofevents[4]++;

  if (mPicoEvent->nBTOFMatch() == 0) return kStOk;
  numberofevents[5]++;

  // cout << "MinBias Event Accepted." << endl;

  // =========================== JetMaker =============================== //
  // get JetMaker pointer
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get jet collection associated with JetMaker
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return kStWarn;
  }

  // ============================= RhoMaker ============================== //
  // get RhoMaker pointer from event: old names "StRho_JetsBG"
  RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerName));
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  fRho = static_cast<StRhoParameter*>(RhoMaker->GetRho());
  if(!fRho) {
    LOG_WARN << Form("Couldn't get fRho object! ") << endm;
    return kStWarn;    
  } 
  
  // get rho/area value from rho object     fRho->ls("");
  fRhoVal = fRho->GetVal();
  // =======================================================================

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // cout << "Number of jets found = " << njets << endl;

  
  // TestTracks();

  doHistograms = kTRUE;

  RunJets();
  // ProcessJetForJetShape();

  return kStOK;
}


//__________________________________________________________________________________________
  
// void StD0Analysis::TestTracks(){

//   const Int_t ntracks = mPicoDst->numberOfTracks();

//   for(unsigned short itrk1 = 0; itrk1 < ntracks; itrk1++){
//     // get track pointer
//     StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(itrk1));
//     if(!trk1){ continue; }

//     if (!IsAnAcceptableTrack(trk1)) continue;

//     TVector3 mTrk1Mom;

//     if(doUsePrimTracks) {
//       mTrk1Mom = trk1->pMom();
//     } else {
//       mTrk1Mom = trk1->gMom(mVertex, Bfield);
//     }

//     Double_t pt = mTrk1Mom.Perp();
//     hptdistro->Fill(pt);

//     Double_t eta = mTrk1Mom.PseudoRapidity();
//     hetadistro->Fill(eta);

//     FillPidHistograms(trk1);

//     for(unsigned short itrk2 = itrk1 + 1; itrk2 < ntracks; itrk2++){
//       if (itrk2 >= ntracks) continue;
//       // get track pointer
//       StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(itrk2));
//       if(!trk2){ continue; }

//       if (!IsAnAcceptableTrack(trk2)) continue;

//       if (!TopologicalCuts(trk1, trk2)) continue;

//       TVector3 mTrk2Mom, mResMom;

//       if(doUsePrimTracks) {
//         mTrk2Mom = trk2->pMom();
//       } else {
//         mTrk2Mom = trk2->gMom(mVertex, Bfield);
//       }

//       mResMom = mTrk1Mom + mTrk2Mom;

//       int particle;
//       double mass;
//       TVector3 dummy_vec;
//       InvariantMass(trk1, trk2, particle, mass, dummy_vec);

//       if (particle==-2){
//         kaonpt_wholeevent->Fill(mTrk1Mom.Perp());
//         pionpt_wholeevent->Fill(mTrk2Mom.Perp());
//         kaon_pion_pt_wholeevent->Fill(mTrk1Mom.Perp(), mTrk2Mom.Perp());
//         d0pt_wholeevent->Fill(mResMom.Perp());
//       }

//       if (particle==2){
//         kaonbgpt_wholeevent->Fill(mTrk1Mom.Perp());
//         pionbgpt_wholeevent->Fill(mTrk2Mom.Perp());
//         kaon_pion_bg_pt_wholeevent->Fill(mTrk1Mom.Perp(), mTrk2Mom.Perp());
//         d0bgpt_wholeevent->Fill(mResMom.Perp());
//       }

//       if (particle==-2) invmass_wholeevent->Fill(mass);
//       if (particle==2) invmassbg_wholeevent->Fill(mass);

//       int pid1 = IsWhatParticle(trk1);
//       int pid2 = IsWhatParticle(trk2);

//       if (abs(pid1)==1) dcapiPV->Fill(trk1->gDCA(mPicoEvent->primaryVertex()).Mag());
//       else dcakPV->Fill(trk1->gDCA(mPicoEvent->primaryVertex()).Mag());

//       if (abs(pid2)==2) dcapiPV->Fill(trk2->gDCA(mPicoEvent->primaryVertex()).Mag());
//       else dcakPV->Fill(trk2->gDCA(mPicoEvent->primaryVertex()).Mag());

//     }
//   }
// }

//__________________________________________________________________________________________

void StD0Analysis::RunJets()
{

  double bfield = mPicoEvent->bField();

  Double_t d0ptcuts[7] = {0.0, 0.5, 2.0, 10.0, 20.0, 30.0, 1000.0};


  if(fCorrJetPt) {
    fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
  } else {
    fLeadingJet = GetLeadingJet(fJetMakerName);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
  }

  d0CandidateID = {};
  d0BgCandidateULID = {};
  d0BgCandidateLSID = {};
  
  Int_t njets = fJets->GetEntries();

  hnjets->Fill(njets);

  const Int_t ntracks = mPicoDst->numberOfTracks();

  Int_t d0CandidateJet, d0BgCandidateULJet, d0BgCandidateLSJet;

  for (int ijet = 0; ijet < njets; ijet++){

    int numberofhfttracks = 0;

    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;

    
    // if (abs(jet->Eta()) > 1.) continue;
    double jetpt = jet->Pt();
    double jeteta = jet->Eta();
    double jetphi = jet->Phi();
    double jetarea = jet->Area();
    double corrjetpt = jet->Pt() - jetarea*fRhoVal;
    // if(abs(jeteta) > 1.) continue;
    // if (jet->GetNumberOfTracks() == 1) continue;

    hJetPt->Fill(jetpt);
    hJetPtCorr->Fill(corrjetpt);
    hJetEta->Fill(jeteta);
    hJetPhi->Fill(jetphi);

    int jetptbin = -99;

    if (corrjetpt < 10.0) continue;
    if (corrjetpt < 15.0) jetptbin = 0;
    else if (corrjetpt >= 15.0 && corrjetpt < 25.0) jetptbin = 1;
    else if (corrjetpt >= 25.0) jetptbin = 2;
    else continue;

    hNumberofTracks->Fill(jet->GetNumberOfTracks());

    for(unsigned short itrk1 = 0; itrk1 < jet->GetNumberOfTracks(); itrk1++){
      // get track pointer
      int trackid1 = jet->TrackAt(itrk1);
      StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(trackid1));
      if(!trk1){ continue; }

      TVector3 mTrk1Mom, mTrk2Mom, mResMom;

      if(doUsePrimTracks) {
        mTrk1Mom = trk1->pMom();
      } else {
        mTrk1Mom = trk1->gMom(mVertex, Bfield);
      }

      Double_t phi = standardPhi(mTrk1Mom.Phi());
      Double_t eta = mTrk1Mom.PseudoRapidity();
      Double_t pt = mTrk1Mom.Perp();

      hConstPt->Fill(pt);
      hConstEta->Fill(eta);
      hConstPhi->Fill(phi);

      if (IsAnAcceptableTrack(trk1)) numberofhfttracks++;

      for(unsigned short itrk2 = itrk1 + 1; itrk2 < jet->GetNumberOfTracks(); itrk2++){
      if (itrk2 >= jet->GetNumberOfTracks()) continue;
        // get track pointer
        int trackid2 = jet->TrackAt(itrk2);
        StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(itrk2));
        if(!trk2){ continue; }

        if(doUsePrimTracks) {
          mTrk2Mom = trk2->pMom();
        } else {
          mTrk2Mom = trk2->gMom(mVertex, Bfield);
        }

        int topo_cut = TopologicalCuts(trk1, trk2);

        if (topo_cut == -99) continue;

        mResMom = mTrk1Mom + mTrk2Mom;

        Double_t phi_mom = standardPhi(mResMom.Phi());
        Double_t eta_mom = mResMom.PseudoRapidity();
        Double_t pt_mom = mResMom.Perp();

        // Ensuring that the candidate is within the jet cone

        double deltaetad = dEta(jeteta, eta);
        double deltaphid = dPhi(jetphi, phi);
        double deltard = dR(deltaetad, deltaphid);

        if (deltard > fJetRad) continue;

        int pid1 = IsWhatParticle(trk1);
        int pid2 = IsWhatParticle(trk2);

        double mass = InvariantMass(trk1, trk2);

        if (pid1*pid2 == -2){
          if (abs(pid1)==1){
            pionpt->Fill(mTrk1Mom.Perp());
            kaonpt->Fill(mTrk2Mom.Perp());
            kaonpionpt->Fill(mTrk1Mom.Perp(), mTrk2Mom.Perp());
          }

          else {
            pionpt->Fill(mTrk2Mom.Perp());
            kaonpt->Fill(mTrk1Mom.Perp());
            kaonpionpt->Fill(mTrk2Mom.Perp(), mTrk1Mom.Perp());
          }

          d0pt->Fill(pt_mom);
        }

        if (pid1*pid2 == 2){
          if (abs(pid1)==1){
            pionbgpt->Fill(mTrk1Mom.Perp());
            kaonbgpt->Fill(mTrk2Mom.Perp());
            kaonpionbgpt->Fill(mTrk1Mom.Perp(), mTrk2Mom.Perp());
          }

          else {
            pionbgpt->Fill(mTrk2Mom.Perp());
            kaonbgpt->Fill(mTrk1Mom.Perp());
            kaonpionbgpt->Fill(mTrk2Mom.Perp(), mTrk1Mom.Perp());
          }

          d0bgpt->Fill(pt_mom);
        }

        if (topo_cut == 1){
          for(int ptbin = 0; ptbin < sizeof(d0ptcuts)/sizeof(*d0ptcuts) - 1; ptbin++){
            if ((pt_mom >= d0ptcuts[ptbin]) && (pt_mom < d0ptcuts[ptbin+1])){
              if(pid1*pid2 == -2) invmass_ptbin[0][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == 2) invmassbg_ptbin[0][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == -2) invmass_ptbin[1][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == 2) invmassbg_ptbin[1][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == -2) invmass_ptbin[2][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == 2) invmassbg_ptbin[2][jetptbin][ptbin]->Fill(mass);
            }
          }

          if (pt_mom > 0.5){
            if(pid1*pid2 == -2) invmass[0]->Fill(mass);
            if(pid1*pid2 == 2) invmassbg[0]->Fill(mass);
            if(pid1*pid2 == -2) invmass[1]->Fill(mass);
            if(pid1*pid2 == 2) invmassbg[1]->Fill(mass);
            if(pid1*pid2 == -2) invmass[2]->Fill(mass);
            if(pid1*pid2 == 2) invmassbg[2]->Fill(mass);
          }
        }

        if (topo_cut == 2){
          for(int ptbin = 0; ptbin < sizeof(d0ptcuts)/sizeof(*d0ptcuts) - 1; ptbin++){
            if ((pt_mom >= d0ptcuts[ptbin]) && (pt_mom < d0ptcuts[ptbin+1])){
              if(pid1*pid2 == -2) invmass_ptbin[1][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == 2) invmassbg_ptbin[1][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == -2) invmass_ptbin[2][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == 2) invmassbg_ptbin[2][jetptbin][ptbin]->Fill(mass);
            }
          }

          if (pt_mom > 0.5){
            if(pid1*pid2 == -2) invmass[1]->Fill(mass);
            if(pid1*pid2 == 2) invmassbg[1]->Fill(mass);
            if(pid1*pid2 == -2) invmass[2]->Fill(mass);
            if(pid1*pid2 == 2) invmassbg[2]->Fill(mass);
          }
        }

        if (topo_cut == 3){
          for(int ptbin = 0; ptbin < sizeof(d0ptcuts)/sizeof(*d0ptcuts) - 1; ptbin++){
            if ((pt_mom >= d0ptcuts[ptbin]) && (pt_mom < d0ptcuts[ptbin+1])){
              if(pid1*pid2 == -2) invmass_ptbin[2][jetptbin][ptbin]->Fill(mass);
              if(pid1*pid2 == 2) invmassbg_ptbin[2][jetptbin][ptbin]->Fill(mass);
            }
          }

          if (pt_mom > 0.5){
            if(pid1*pid2 == -2) invmass[2]->Fill(mass);
            if(pid1*pid2 == 2) invmassbg[2]->Fill(mass);
          }
        }
        
      }
    }
    hNumberofHFTTracks->Fill(numberofhfttracks);
    hRatioOfAcceptedTracks->Fill((double)numberofhfttracks/(double)jet->GetNumberOfTracks());
  }
}


Bool_t StD0Analysis::IsAnAcceptableTrack(StPicoTrack *trk){
  if (trk->nHitsDedx() < 15) return kFALSE;
  if (float(trk->nHitsDedx())/float(trk->nHitsMax()) < 0.52) return kFALSE;
  if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 3.) return kFALSE; 
  if (trk->chi2() > 3) return kFALSE;

  // HFT Hits only
  if (!trk->isHFTTrack()) return kFALSE;
  if (!trk->hasPxl1Hit() || !trk->hasPxl2Hit() || !trk->hasIstHit()) return kFALSE;

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

  // bool bemctrack = trk->isBemcMatchedTrack();
  // if ((!bemctrack)) return kFALSE;

  if (pt < 0.5) return kFALSE;
  if (charge == 0) return kFALSE;
  if (abs(eta) > 1) return kFALSE;

  return kTRUE;
}

Bool_t StD0Analysis::IsTrackWithinJet(StJet *jet, StPicoTrack *trk){
  if (!IsAnAcceptableTrack(trk)) return kFALSE;

  double jetpt = jet->Pt();
  double jeteta = jet->Eta();
  double jetphi = jet->Phi();
  double jetarea = jet->Area();

  TVector3 mTrkMom;

  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, Bfield);
  }

  Double_t phi = standardPhi(mTrkMom.Phi());
  Double_t eta = mTrkMom.PseudoRapidity();

  double deltaeta = dEta(jeteta, eta);
  double deltaphi = dPhi(jetphi, phi);
  double deltar = dR(deltaeta, deltaphi);

  // cout << "Delta R" << "\t" << deltar << endl;
  if (deltar < fJetRad) return kTRUE;
  else return kFALSE; 
}

// void StD0Analysis::FillPidHistograms(StPicoTrack *trk){

//   TVector3 mTrkMom;
//   if(doUsePrimTracks) {
//     mTrkMom = trk->pMom();
//   } else {
//     mTrkMom = trk->gMom(mVertex, Bfield);
//   }

//   // track variables
//   double pt = mTrkMom.Perp();
//   double p = mTrkMom.Mag();
//   double eta = mTrkMom.PseudoRapidity();

//   hetadistro->Fill(eta);
//   hptdistro->Fill(pt);

//   bool toftrack = kFALSE;

//   int tof_loc = trk->bTofPidTraitsIndex();

//   if (tof_loc >= 0)
//   {
//     StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
//     if (tofpointer){toftrack = kTRUE;}
//   }

//   if (toftrack){

//     StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
//     double invbeta_from_tof = tofpointer->btofBeta();
//     invbeta_from_tof = 1/invbeta_from_tof;
//     double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
//     double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.011;
//     normalised_invbetavpT_tof_pi->Fill(pt, normalisedinvbeta_for_pi);
//     double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
//     double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.011;
//     normalised_invbetavpT_tof_ka->Fill(pt, normalisedinvbeta_for_ka);
//   }


//   else{
//       double zpi = trk->nSigmaPion();
//       z_pi->Fill(pt, zpi);
//       double zka = trk->nSigmaKaon();
//       z_ka->Fill(pt, zka);
//   }
// }


Int_t StD0Analysis::TopologicalCuts(StPicoTrack *trk1, StPicoTrack *trk2){

  int cutnum = -99;

  int pid1, pid2;
  double m1, m2;
  double e1, e2;

  IsWhatParticle(trk1, pid1, m1, e1);
  IsWhatParticle(trk2, pid2, m2, e2);

  if (abs(pid1*pid2) != 2) return cutnum;

  TVector3 mTrk1Mom;
  if(doUsePrimTracks) {
    mTrk1Mom = trk1->pMom();
  } else {
    mTrk1Mom = trk1->gMom(mVertex, Bfield);
  }

  TVector3 mTrk2Mom;
  if(doUsePrimTracks) {
    mTrk2Mom = trk2->pMom();
  } else {
    mTrk2Mom = trk2->gMom(mVertex, Bfield);
  }

  // to be used for testing with preview II pico production
  StPicoPhysicalHelix kHelix = trk1->helix(Bfield);
  StPicoPhysicalHelix pHelix = trk2->helix(Bfield);

  // move origins of helices to the primary vertex origin
  kHelix.moveOrigin(kHelix.pathLength(mVertex));
  pHelix.moveOrigin(pHelix.pathLength(mVertex));

  // use straight lines approximation to get point of DCA of kaon-pion pair
  TVector3 const kMom = kHelix.momentum(Bfield*(1.e-14));//*(1.e-14);
  TVector3 const pMom = pHelix.momentum(Bfield*(1.e-14));//*(1.e-14);

  // cout << mTrk1Mom.Mag() << "\t" << kMom.Mag() << "\t" << mTrk2Mom.Mag() << "\t" << pMom.Mag() << endl;

  StPicoPhysicalHelix const kStraightLine(kMom, kHelix.origin(), 0, trk1->charge());
  StPicoPhysicalHelix const pStraightLine(pMom, pHelix.origin(), 0, trk2->charge());


  pair<double, double> const ss = kStraightLine.pathLengths(pStraightLine);
  TVector3 const kAtDcaToPion = kStraightLine.at(ss.first);
  TVector3 const pAtDcaToKaon = pStraightLine.at(ss.second);

  // calculate DCA of pion to kaon at their DCA
  Double_t mDcaDaughters = (kAtDcaToPion - pAtDcaToKaon).Mag();

  // calculate Lorentz vector of kaon-pion pair
  TVector3 const kMomAtDca = kHelix.momentumAt(ss.first, Bfield*(1.e-14)); //*(1.e-14);
  TVector3 const pMomAtDca = pHelix.momentumAt(ss.second, Bfield*(1.e-14)); //*(1.e-14);

  TLorentzVector const kFourMom(kMomAtDca, e1);
  TLorentzVector const pFourMom(pMomAtDca, e2);

  TLorentzVector mLorentzVector = kFourMom + pFourMom;

  // calculate cosThetaStar
  //TLorentzVector const kpFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
  //TLorentzVector const kFourMomStar = kFourMom.Boost(kpFourMomReverse);
  //mCosThetaStar = std::cos(kFourMomStar.Vect().Angle(mLorentzVector.Vect()));

  // calculate pointing angle and decay length
  TVector3 const vtxToV0 = (kAtDcaToPion + pAtDcaToKaon) * 0.5 - mVertex;
  Double_t mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());

  //mPointingAngle = vtxToV0.Angle(finalFourVector.vect());

  Double_t mDecayLength = vtxToV0.Mag();

  Double_t perpDcaToVtx = mDecayLength*TMath::Sin(mPointingAngle);

  // calculate DCA of tracks to primary vertex
  Double_t mKaonDca = (kHelix.origin() - mVertex).Mag();
  Double_t mPionDca = (pHelix.origin() - mVertex).Mag();

  if (doHistograms && pid1*pid2==-2) decaylengthd0US->Fill(mDecayLength);
  if (doHistograms && pid1*pid2==-2) distancepikUS->Fill(mDcaDaughters);
  if (doHistograms && pid1*pid2==-2) distanced0PVUS->Fill(perpDcaToVtx);
  if (doHistograms && pid1*pid2==-2) dcakPVUS->Fill(mKaonDca);
  if (doHistograms && pid1*pid2==-2) dcapiPVUS->Fill(mPionDca);

  if (doHistograms && pid1*pid2==2) decaylengthd0LS->Fill(mDecayLength);
  if (doHistograms && pid1*pid2==2) distancepikLS->Fill(mDcaDaughters);
  if (doHistograms && pid1*pid2==2) distanced0PVLS->Fill(perpDcaToVtx);
  if (doHistograms && pid1*pid2==2) dcakPVLS->Fill(mKaonDca);
  if (doHistograms && pid1*pid2==2) dcapiPVLS->Fill(mPionDca);

  // if (mDecayLength < 0.0212) return kFALSE;
  // if (mDcaDaughters > 0.0057) return kFALSE;
  // if (perpDcaToVtx > 0.0038) return kFALSE;
  // if (mKaonDca < 0.0095) return kFALSE;
  // if (mPionDca < 0.0086) return kFALSE;

  if ((mDecayLength > 0.0212) && (mDcaDaughters < 0.0057) && (perpDcaToVtx < 0.0038) && (mKaonDca > 0.0095) && (mPionDca > 0.0086)){
    cutnum = 1;
    return cutnum;
  }

  if ((mDecayLength > 0.5*0.0212) && (mDcaDaughters < 2*0.0057) && (perpDcaToVtx < 2*0.0038) && (mKaonDca > 0.5*0.0095) && (mPionDca > 0.5*0.0086)){
    cutnum = 2;
    return cutnum;
  }

  if ((mDecayLength > 0.25*0.0212) && (mDcaDaughters < 4*0.0057) && (perpDcaToVtx < 4*0.0038) && (mKaonDca > 0.25*0.0095) && (mPionDca > 0.25*0.0086)){
    cutnum = 3;
    return cutnum;
  }

  return cutnum;
  // cout << "Topo Cuts: " << mDecayLength << "\t" << mDcaDaughters << "\t" << perpDcaToVtx << "\t" << mKaonDca << "\t" << mPionDca << endl;
} 


void StD0Analysis::IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e){ // NEW PID APPROACH
  if(!IsAnAcceptableTrack(trk)){pid = 0; m = 0.; e = 0.; return;}

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

  if ((abs(zpi) > 2) && (abs(zka) > 2)) {pid = 0; m = 0.; e = 0.; return;}

  bool tpc_pion = kFALSE;
  bool tpc_kaon = kFALSE;

  tpc_pion = (abs(zpi) <= abs(zka)) ? kTRUE : kFALSE;
  tpc_kaon = (abs(zpi) > abs(zka)) ? kTRUE : kFALSE;

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

    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    double invbeta_from_tof = tofpointer->btofBeta();
    invbeta_from_tof = 1/invbeta_from_tof;

    double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
    double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
    double norm_invbeta_pr = TMath::Sqrt(pow(Mproton,2)/pow(p,2) + 1);

    double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.011;
    double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.011;
    double normalisedinvbeta_for_pr = (invbeta_from_tof-norm_invbeta_pr)/0.011;

    if ((abs(normalisedinvbeta_for_pi) > 2) && (abs(normalisedinvbeta_for_ka) > 2)) {pid = 0; m = 0.; e = 0.; return;}

    tof_pion = (abs(normalisedinvbeta_for_pi) <= abs(normalisedinvbeta_for_ka)) ? kTRUE : kFALSE;
    tof_kaon = (abs(normalisedinvbeta_for_pi) > abs(normalisedinvbeta_for_ka)) ? kTRUE : kFALSE;



    if (tof_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tof_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}

  }

  else{
    if (tpc_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tpc_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}
  }
}


Int_t StD0Analysis::IsWhatParticle(StPicoTrack *trk){ // Just to get the PID out
  int pid;
  double m; 
  double e;
  IsWhatParticle(trk, pid, m, e);
  return pid;
}


void StD0Analysis::IsWhatParticleFlipped(StPicoTrack *trk, int &pidflipped, double &mflipped, double &eflipped){
  int pid;
  double m, e;

  IsWhatParticle(trk, pid, m, e);
  if (pid==1){
    pidflipped==2;
    mflipped = Mkaon;
    eflipped = TMath::Sqrt(pow(e,2) - pow(Mpion, 2) + pow(Mkaon, 2));
  }

  if (pid==-1){
    pidflipped==-2;
    mflipped = Mkaon;
    eflipped = TMath::Sqrt(pow(e,2) - pow(Mpion, 2) + pow(Mkaon, 2));
  }

  if (pid==2){
    pidflipped==1;
    mflipped = Mpion;
    eflipped = TMath::Sqrt(pow(e,2) - pow(Mkaon, 2) + pow(Mpion, 2));
  }

  if (pid==-2){
    pidflipped==-1;
    mflipped = Mpion;
    eflipped = TMath::Sqrt(pow(e,2) - pow(Mkaon, 2) + pow(Mpion, 2));
  }
}


void StD0Analysis::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum){

  int particle1, particle2;
  double mass1, mass2;
  double energy1, energy2;

  // cout << "IsWhatParticle called" << endl;

  IsWhatParticle(trk1, particle1, mass1, energy1);
  IsWhatParticle(trk2, particle2, mass2, energy2);

  TVector3 mTrk1Mom, mTrk2Mom;

  if(doUsePrimTracks) {
    mTrk1Mom = trk1->pMom();
  } else {
    mTrk1Mom = trk1->gMom(mVertex, Bfield);
  }

  if(doUsePrimTracks) {
    mTrk2Mom = trk2->pMom();
  } else {
    mTrk2Mom = trk2->gMom(mVertex, Bfield);
  }

  if ((mass1<0.000001) || (mass2<0.000001) || (energy1<0.000001) || (energy2<0.000001)){
    particle=0;
    invmass=0.;
    momentum = mTrk1Mom + mTrk2Mom;
    return;
  }
  else{
    particle = particle1*particle2;
    invmass = TMath::Sqrt(pow(mass1,2) + pow(mass2,2) + 2*(energy1*energy2 - mTrk1Mom*mTrk2Mom));
    momentum = mTrk1Mom + mTrk2Mom;
  }
}

Double_t StD0Analysis::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2){
  int particle; 
  double invmass; 
  TVector3 momentum;

  InvariantMass(trk1, trk2, particle, invmass, momentum);
  return invmass;
}

void StD0Analysis::InvariantMassFlipped(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum){

  int particle1, particle2;
  double mass1, mass2;
  double energy1, energy2;

  // cout << "IsWhatParticle called" << endl;

  IsWhatParticleFlipped(trk1, particle1, mass1, energy1);
  IsWhatParticleFlipped(trk2, particle2, mass2, energy2);

  TVector3 mTrk1Mom, mTrk2Mom;

  if(doUsePrimTracks) {
    mTrk1Mom = trk1->pMom();
  } else {
    mTrk1Mom = trk1->gMom(mVertex, Bfield);
  }

  if(doUsePrimTracks) {
    mTrk2Mom = trk2->pMom();
  } else {
    mTrk2Mom = trk2->gMom(mVertex, Bfield);
  }

  if ((mass1<0.000001) || (mass2<0.000001) || (energy1<0.000001) || (energy2<0.000001)){
    particle=0;
    invmass=0.;
    momentum = mTrk1Mom + mTrk2Mom;
    return;
  }
  else{
    particle = particle1*particle2;
    invmass = TMath::Sqrt(pow(mass1,2) + pow(mass2,2) + 2*(energy1*energy2 - mTrk1Mom*mTrk2Mom));
    momentum = mTrk1Mom + mTrk2Mom;
  }
}

// void StD0Analysis::IsD0(StPicoTrack *trk1, StPicoTrack *trk2, bool &D0, bool &D0BgUnlike, bool &D0BgLike, double &mass, TVector3 &mom3){

//   int particle;

//   InvariantMass(trk1, trk2, particle, mass, mom3);

//   // cout << "particle = " << particle << "\t" << "Mass = " << mass << endl;

//   D0 = kFALSE; D0BgUnlike = kFALSE; D0BgLike = kFALSE;

//   double pt = mom3.Perp();

//   if(particle == -2){

//     invmass->Fill(mass);

//     invmassvpt->Fill(pt, mass);

//     if ((mass > fInvMassSignal1) && (mass < fInvMassSignal2)){
//       D0 = kTRUE;
//     }
//     if (((mass > fInvMassULBg1) && (mass < fInvMassSignal1)) || ((mass > fInvMassSignal2) && (mass < fInvMassULBg2))){
//       D0BgUnlike = kTRUE;
//     }
//   }

//   if(particle == 2){
//     invmassbg->Fill(mass);
//     invmassbgvpt->Fill(pt, mass);

//     if ((mass > fInvMassLSBg1) && (mass < fInvMassLSBg2)){
//       D0BgLike = kTRUE;
//     }
//   }
// }

// void StD0Analysis::ProcessTrackForKF(StPicoTrack *trk, StPicoTrackCovMatrix *picocov, Double_t params[6], Double_t cov[21]){
//   const Float_t *paramfromcov = picocov->params();
//   const Float_t *sigmafromcov = picocov->sigmas();
//   const Float_t *corrfromcov = picocov->correlations();

//   static StDcaGeometry dcaGeometryPi;

//   Float_t err[15] = {sigmafromcov[0], 
//                      corrfromcov[0], sigmafromcov[1],
//                      corrfromcov[1], corrfromcov[2], sigmafromcov[2],
//                      corrfromcov[3], corrfromcov[4], corrfromcov[5], sigmafromcov[3],
//                      corrfromcov[6], corrfromcov[7], corrfromcov[8], corrfromcov[9], sigmafromcov[4]};

//   dcaGeometryPi.set(paramfromcov, err);
//   dcaGeometryPi.GetXYZ(params, cov);
// }

// Bool_t StD0Analysis::IsD0KFParticle(StPicoTrack *trk1, StPicoTrack *trk2, StPicoTrackCovMatrix *cov1, StPicoTrackCovMatrix *cov2){
//   Double_t params1[6], params2[6];
//   Double_t covar1[21], covar2[21];

//   int charge1 = trk1->charge();
//   int charge2 = trk2->charge();

//   TVector3 mom1;
//   if(doUsePrimTracks) {
//     mom1 = trk1->pMom();
//   } else {
//     mom1 = trk1->gMom(mVertex, 0.0);
//   }

//   TVector3 mom2;
//   if(doUsePrimTracks) {
//     mom2 = trk2->pMom();
//   } else {
//     mom2 = trk2->gMom(mVertex, 0.0);
//   }

//   Int_t pid1, pid2;

//   Double_t m1, m2, e1, e2;

//   IsWhatParticle(trk1, pid1, m1, e1);
//   IsWhatParticle(trk2, pid2, m2, e2);

//   Int_t particle = pid1*pid2;

//   if (particle!=-2){return kFALSE;} 

//   ProcessTrackForKF(trk1, cov1, params1, covar1);
//   ProcessTrackForKF(trk2, cov2, params2, covar2);

//   Int_t pdg1, pdg2;

//   switch(pid1){
//     case 1: pdg1 = 211;     break;
//     case 2: pdg1 = 321;     break;
//     // case 3: pdg1 = 2212;    break;
//     case -1:  pdg1 = -211;   break;
//     case -2:  pdg1 = -321;   break;
//     // case -3:  pdg1 = -2212;  break;
//   }

//   switch(pid2){
//     case 1: pdg2 = 211;     break;
//     case 2: pdg2 = 321;     break;
//     // case 3: pdg2 = 2212;    break;
//     case -1:  pdg2 = -211;   break;
//     case -2:  pdg2 = -321;   break;
//     // case -3:  pdg2 = -2212;  break;
//   }

//   KFParticle d1, d2; 
  

//   d1.Create(params1, covar1, charge1, pdg1);
//   d2.Create(params2, covar2, charge2, pdg2);

//   double dpik = d1.GetDistanceFromParticle(d2);
//   double depik = d1.GetDeviationFromParticle(d2);

//   // cout << "Zero = " << dpik - d2.GetDistanceFromParticle(d1) << endl;

//   distancepik->Fill(dpik);
//   deviationpik->Fill(depik);

//   KFParticle mother = KFParticle(d1, d2);

//   Double_t x_mom = mother.GetX();
//   Double_t y_mom = mother.GetY();
//   Double_t z_mom = mother.GetZ();

//   TVector3 mom_cord(mother.GetX(), mother.GetY(), mother.GetZ());

//   Double_t distanced1fromPV = (trk1->origin() - mPicoEvent->primaryVertex()).Mag();
//   Double_t distanced2fromPV = (trk2->origin() - mPicoEvent->primaryVertex()).Mag();
//   Double_t distanced0fromPV = (mom_cord - mPicoEvent->primaryVertex()).Mag();

//   distanced0PV->Fill(distanced0fromPV);
//   // cout << x_mom << "\t" << y_mom << "\t" << z_mom << "\t" << distancefromPV << endl;

//   double dpid0 = d1.GetDistanceFromParticle(mother);
//   double depid0 = d1.GetDeviationFromParticle(mother);

//   double dkd0 = d2.GetDistanceFromParticle(mother);
//   double dekd0 = d2.GetDeviationFromParticle(mother);

//   distancepid0->Fill(dpid0);
//   deviationpid0->Fill(depid0);

//   distancekd0->Fill(dkd0);
//   deviationkd0->Fill(dekd0);

//   double decaylength, decaylengthsigma;
//   mother.GetDecayLength(decaylength, decaylengthsigma);

//   decaylengthd0->Fill(decaylength);

//   // cout << "Distances " << dpik << "\t" << dpid0 << "\t" << dkd0 << endl;
//   // cout << "Deviations " << depik << "\t" << depid0 << "\t" << dekd0 << endl;
//   // cout << "Decay Length = " << decaylength << endl;

//   double mass1, mass1err;
//   d1.GetMass(mass1, mass1err);
//   double mass2, mass2err;
//   d2.GetMass(mass2, mass2err);

//   double mass_mom, mass_mom_err;
//   mother.GetMass(mass_mom, mass_mom_err);
//   //mother.SetMassConstraint(1.865);

//   kfmass->Fill(mass_mom);

//   if ((distanced0fromPV)<1 && (dpik)<1) return kTRUE;
//   else return kFALSE; 
// }

// void StD0Analysis::InvariantMassWithCuts(StPicoTrack *trk1, StPicoTrack *trk2, StPicoTrackCovMatrix *cov1, StPicoTrackCovMatrix *cov2){

//   double bfield = mPicoEvent->bField();

//   KFParticle::SetField(bfield);

//   Double_t params1[6], params2[6];
//   Double_t covar1[21], covar2[21];

//   int charge1 = trk1->charge();
//   int charge2 = trk2->charge();

//   TVector3 mom1;
//   if(doUsePrimTracks) {
//     mom1 = trk1->pMom();
//   } else {
//     mom1 = trk1->gMom(mVertex, Bfield);
//   }

//   TVector3 mom2;
//   if(doUsePrimTracks) {
//     mom2 = trk2->pMom();
//   } else {
//     mom2 = trk2->gMom(mVertex, Bfield);
//   }

//   TVector3 mom_d0 = mom1 + mom2;
//   double d0pt = mom_d0.Perp();

//   double pt1 = mom1.Perp();
//   double pt2 = mom2.Perp();

//   Int_t pid1, pid2;

//   Double_t m1, m2, e1, e2;

//   IsWhatParticle(trk1, pid1, m1, e1);
//   IsWhatParticle(trk2, pid2, m2, e2);

//   Int_t particle = pid1*pid2;

//   if (abs(particle)!=2){return;} 

//   // ProcessTrackForKF(trk1, cov1, params1, covar1);
//   // ProcessTrackForKF(trk2, cov2, params2, covar2);

//   Int_t pdg1, pdg2;

//   switch(pid1){
//     case 1: pdg1 = 211;     break;
//     case 2: pdg1 = 321;     break;
//     // case 3: pdg1 = 2212;    break;
//     case -1:  pdg1 = -211;   break;
//     case -2:  pdg1 = -321;   break;
//     // case -3:  pdg1 = -2212;  break;
//   }

//   switch(pid2){
//     case 1: pdg2 = 211;     break;
//     case 2: pdg2 = 321;     break;
//     // case 3: pdg2 = 2212;    break;
//     case -1:  pdg2 = -211;   break;
//     case -2:  pdg2 = -321;   break;
//     // case -3:  pdg2 = -2212;  break;
//   }

//   Double_t Cuts[4] = {0.3, 0.4, 0.5, 0.6};
//   Double_t ptcuts[3] = {3.0, 4.0, 5.0};
//   Double_t d0ptcuts[4] = {2.0, 3.0, 5.0, 6.0};

//   //cout << " =========================================================================== " << endl;

//   double dpik = (trk1->origin() - trk2->origin()).Mag();
//   distancepik->Fill(dpik);

//   TVector3 mom_cord = 0.5*(trk1->origin() + trk2->origin());
//   double distanced0fromPV = (mom_cord - mPicoEvent->primaryVertex()).Mag();
//   distanced0PV->Fill(distanced0fromPV);

//   int dummy;
//   double mass, massflipped;
//   TVector3 dummy_vec;
//   InvariantMass(trk1, trk2, dummy, mass, dummy_vec);

//   InvariantMassFlipped(trk1, trk2, dummy, massflipped, dummy_vec);


//   // D0 cuts only
//   if(particle==-2){
//     kpi_ul_pt->Fill(pt1, pt2);

//     hD0pt->Fill(d0pt);

//     invmass->Fill(mass);
//     invmass_flipped->Fill(massflipped);

//     for (int i = 0; i < 4; i++){ // VertexCuts
//       if (dpik > Cuts[i]) continue;
//       for (int j = 0; j < 3; j++){ // PtCuts
//         if (pt1 + pt2 < ptcuts[j]) continue;
//         for (int k = 0; k < 4; k++){ // D0PtCuts
//           if (d0pt < d0ptcuts[k]) continue;
//           invmasscuts[i][j][k]->Fill(mass);
//           invmass_flipped_cuts[i][j][k]->Fill(massflipped);
//         }
//       }
//     }
//   }

//   if(particle==2){
//     kpi_ls_pt->Fill(pt1, pt2);

//     invmassbg->Fill(mass);

//     for (int i = 0; i < 4; i++){ // VertexCuts
//       if (dpik > Cuts[i]) continue;
//       for (int j = 0; j < 3; j++){ // PtCuts
//         if (pt1 + pt2 < ptcuts[j]) continue;
//         for (int k = 0; k < 4; k++){ // D0PtCuts
//           if (d0pt < d0ptcuts[k]) continue;
//           invmassbgcuts[i][j][k]->Fill(mass);
//         }
//       }
//     }
//   }
// }

// void StD0Analysis::ProcessJetForJetShape(){

//   // cout << "Triggered Jet Shape Module" << endl;
//   Int_t centbin = GetFourCentBin(fCentralityScaled);

//   // cout << fCentralityScaled << "\t" << centbin << endl;

//   if (centbin == -99) return;

//   Int_t njets = fJets->GetEntries();
//   unsigned int ntracks = mPicoDst->numberOfTracks();

//   // cout << "Njets = " << njets << endl;

//   float diffdist[numberofbins];
//   float diffd0dist[numberofbins];
//   float diffd0distBgUL[numberofbins];
//   float diffd0distBgLS[numberofbins];

//   float diffshape[numberofbins];
//   float diffd0shape[numberofbins];
//   float diffd0shapeBgUL[numberofbins];
//   float diffd0shapeBgLS[numberofbins];


//   // Differential jet shape for all jets
//   for (int ijet = 0; ijet < njets; ijet++){
//     StJet *jet = static_cast<StJet*>(fJets->At(ijet));
//     if(!jet) continue;

//     double jetpt = jet->Pt();
//     double jeteta = jet->Eta();
//     double jetphi = jet->Phi();
//     double jetarea = jet->Area();
//     double corrjetpt = jet->Pt() - jetarea*fRhoVal;
//     // cout << "Jet Pt from Jet Shape = " << jetpt << "\t" << corrjetpt << endl; 
//     if (corrjetpt < 10.) continue;

//     JetPhi[centbin]->Fill(jetphi);
//     JetPt[centbin]->Fill(jetpt);
//     JetPtCorr[centbin]->Fill(corrjetpt);

//     for (int bin = 0; bin < numberofbins; bin++){
//       diffdist[bin] = 0;
//       diffshape[bin] = 0;
//     }

//     // vector<fastjet::PseudoJet> fConstituents = jet->GetJetConstituents();
//     // for(UInt_t ic = 0; ic < fConstituents.size(); ++ic) {
//     //   // get user defined index
//     //   Int_t uid = fConstituents[ic].user_index();
//     //   double cpt = fConstituents[ic].perp();
//     //   double ceta = fConstituents[ic].eta();
//     //   double cphi = fConstituents[ic].phi();
//     //   cout<<"ic = "<<ic<<", uid = "<<uid<<", cpt = "<<cpt<<", ceta = "<<ceta<<", cphi = "<<cphi<<endl;
//     // }


//     // for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {

//     //   int trackid = jet->TrackAt(itrk);

//     //   StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
//     //   if(!trk){ cout << "No track" << endl; continue; }

//     for(unsigned short itrk = 0; itrk < ntracks; itrk++){
//       // get track pointer
//       StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
//       if(!trk){ continue; }

//       if (!IsTrackWithinJet(jet, trk)) continue;

//       TVector3 mTrkMom;
//       if(doUsePrimTracks) {
//         mTrkMom = trk->pMom();
//       } else {
//         mTrkMom = trk->gMom(mVertex, 0.0);
//       }

//       // track variables
//       double pt = mTrkMom.Perp();
//       double phi = standardPhi(mTrkMom.Phi());
//       double eta = mTrkMom.PseudoRapidity();
//       double px = mTrkMom.x();
//       double py = mTrkMom.y();
//       double pz = mTrkMom.z();
//       double p = mTrkMom.Mag();
//       short charge = trk->charge();

//       ConstPhi[centbin]->Fill(phi);
//       ConstPt[centbin]->Fill(pt);

//       double delphiD = dPhi(jetphi, phi);
//       double deletaD = dEta(jeteta, eta);
//       double delRD = dR(delphiD, deletaD);

//       DeltaPhi[centbin]->Fill(delphiD);
//       DeltaEta[centbin]->Fill(deletaD);
//       Rad[centbin]->Fill(delRD);

//       for (int bin = 0; bin < numberofbins; bin++){
//         if (delRD >= bin*deltar && delRD < (bin+1)*deltar){
//           diffdist[bin]++;
//           diffshape[bin] += pt/jetpt;
//         }
//       } // end of bin loop
//     } // end of track loop

//     for (int bin = 0; bin < numberofbins; bin++){
//       hJetDist[centbin]->Fill((2*bin+1)*deltar/2, diffdist[bin]);
//       hJetShape[centbin]->Fill((2*bin+1)*deltar/2, diffshape[bin]);
//     }
//   }


//   // Jets which have a D0 candidate
//   for (Int_t iJetD0 = 0; iJetD0 < d0CandidateID.size(); iJetD0++){
//     // get jet pointer
//     cout << "iJetD0 = " << d0CandidateID[iJetD0] << endl;
//     StJet *jet = static_cast<StJet*>(fJets->At(d0CandidateID[iJetD0]));
//     if(!jet) continue;

//     double jetpt = jet->Pt();
//     double jeteta = jet->Eta();
//     double jetphi = jet->Phi();
//     double jetarea = jet->Area();
//     double corrjetpt = jet->Pt() - jetarea*fRhoVal;
//     if (corrjetpt < 10.) continue;

//     hJetPhid0[centbin]->Fill(jetphi);
//     hJetPtd0[centbin]->Fill(jetpt);
//     hJetPtd0Corr[centbin]->Fill(corrjetpt);

//     for (int bin = 0; bin < numberofbins; bin++){
//       diffd0dist[bin] = 0;
//       diffd0shape[bin] = 0;
//     }

//     // for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {

//     //   int trackid = jet->TrackAt(itrk);

//     //   StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
//     //   if(!trk){ cout << "No track" << endl; continue; }
//     for(unsigned short itrk = 0; itrk < ntracks; itrk++){
//       // get track pointer
//       StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
//       if(!trk){ continue; }

//       if (!IsTrackWithinJet(jet, trk)) continue;

//       TVector3 mTrkMom;
//       if(doUsePrimTracks) {
//         mTrkMom = trk->pMom();
//       } else {
//         mTrkMom = trk->gMom(mVertex, 0.0);
//       }

//       // track variables
//       double pt = mTrkMom.Perp();
//       double phi = standardPhi(mTrkMom.Phi());
//       double eta = mTrkMom.PseudoRapidity();
//       double px = mTrkMom.x();
//       double py = mTrkMom.y();
//       double pz = mTrkMom.z();
//       double p = mTrkMom.Mag();
//       short charge = trk->charge();

//       hConstPhid0[centbin]->Fill(phi);
//       hConstPtd0[centbin]->Fill(pt);

//       double delphiD = dPhi(jetphi, phi);
//       double deletaD = dEta(jeteta, eta);
//       double delRD = dR(delphiD, deletaD);

//       hDeltaPhid0[centbin]->Fill(delphiD);
//       hDeltaEtad0[centbin]->Fill(deletaD);
//       hRadd0[centbin]->Fill(delRD);

//       for (int bin = 0; bin < numberofbins; bin++){
//         if (delRD >= bin*deltar && delRD < (bin+1)*deltar){
//           diffd0dist[bin]++;
//           diffd0shape[bin] += pt/jetpt;
//         }
//       } // end of bin loop
//     } // end of track loop

//     for (int bin = 0; bin < numberofbins; bin++){
//       hJetDistd0[centbin]->Fill((2*bin+1)*deltar/2, diffd0dist[bin]);
//       hJetShaped0[centbin]->Fill((2*bin+1)*deltar/2, diffd0shape[bin]);
//     }
//   } //end of jet loop

//   // Jets which have a D0 Bg from UL 
//   for (Int_t iJetD0ul = 0; iJetD0ul < d0BgCandidateULID.size(); iJetD0ul++){
//     // get jet pointer
//     cout << "iJetD0 UL = " << d0BgCandidateULID[iJetD0ul] << endl;

//     StJet *jet = static_cast<StJet*>(fJets->At(d0BgCandidateULID[iJetD0ul]));
//     if(!jet) continue;
//     numberoftracks[0]+=jet->GetNumberOfTracks();

//     double jetpt = jet->Pt();
//     double jeteta = jet->Eta();
//     double jetphi = jet->Phi();
//     double jetarea = jet->Area();
//     double corrjetpt = jet->Pt() - jetarea*fRhoVal;
//     if (corrjetpt < 10.) continue;


//     hJetPhid0BgUL[centbin]->Fill(jetphi);
//     hJetPtd0BgUL[centbin]->Fill(jetpt);
//     hJetPtd0BgULCorr[centbin]->Fill(corrjetpt);

//     for (int bin = 0; bin < numberofbins; bin++){
//       diffd0distBgUL[bin] = 0;
//       diffd0shapeBgUL[bin] = 0;
//     }

//     // for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {

//     //   int trackid = jet->TrackAt(itrk);

//     //   StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
//     //   if(!trk){ cout << "No track" << endl; continue; }

//     for(unsigned short itrk = 0; itrk < ntracks; itrk++){
//       // get track pointer
//       StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
//       if(!trk){ continue; }

//       if (!IsTrackWithinJet(jet, trk)) continue;

//       TVector3 mTrkMom;
//       if(doUsePrimTracks) {
//         mTrkMom = trk->pMom();
//       } else {
//         mTrkMom = trk->gMom(mVertex, 0.0);
//       }

//       // track variables
//       double pt = mTrkMom.Perp();
//       double phi = mTrkMom.Phi();
//       double eta = mTrkMom.PseudoRapidity();
//       double px = mTrkMom.x();
//       double py = mTrkMom.y();
//       double pz = mTrkMom.z();
//       double p = mTrkMom.Mag();
//       short charge = trk->charge();

//       hConstPhid0BgUL[centbin]->Fill(phi);
//       hConstPtd0BgUL[centbin]->Fill(pt);

//       double delphiD = dPhi(jetphi, phi);
//       double deletaD = dEta(jeteta, eta);
//       double delRD = dR(delphiD, deletaD);

//       hDeltaPhid0BgUL[centbin]->Fill(delphiD);
//       hDeltaEtad0BgUL[centbin]->Fill(deletaD);
//       hRadd0BgUL[centbin]->Fill(delRD);

//       for (int bin = 0; bin < numberofbins; bin++){
//         if (delRD >= bin*deltar && delRD < (bin+1)*deltar){
//           diffd0distBgUL[bin]++;
//           diffd0shapeBgUL[bin]+= pt/jetpt;
//         }
//       } // end of bin loop
//     } // end of track loop

//     for (int bin = 0; bin < numberofbins; bin++){
//       hJetDistd0BgUL[centbin]->Fill((2*bin+1)*deltar/2, diffd0distBgUL[bin]);
//       hJetShaped0BgUL[centbin]->Fill((2*bin+1)*deltar/2, diffd0shapeBgUL[bin]);
//     }
//   } //end of jet loop

//   // Jets which have a D0 Bg from LS
//   for (Int_t iJetD0ls = 0; iJetD0ls < d0BgCandidateLSID.size(); iJetD0ls++){
//     // get jet pointer
//     cout << "iJetD0 LS = " << d0BgCandidateLSID[iJetD0ls] << endl;

//     StJet *jet = static_cast<StJet*>(fJets->At(d0BgCandidateLSID[iJetD0ls]));
//     if(!jet) continue;
//     numberoftracks[0]+=jet->GetNumberOfTracks();

//     double jetpt = jet->Pt();
//     double jeteta = jet->Eta();
//     double jetphi = jet->Phi();
//     double jetarea = jet->Area();
//     double corrjetpt = jet->Pt() - jetarea*fRhoVal;
//     if (corrjetpt < 10.) continue;

//     hJetPhid0BgLS[centbin]->Fill(jetphi);
//     hJetPtd0BgLS[centbin]->Fill(jetpt);
//     hJetPtd0BgLSCorr[centbin]->Fill(corrjetpt);

//     for (int bin = 0; bin < numberofbins; bin++){
//       diffd0distBgLS[bin] = 0;
//       diffd0shapeBgLS[bin] = 0;
//     }

//     // for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {

//     //   int trackid = jet->TrackAt(itrk);

//     //   StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
//     //   if(!trk){ cout << "No track" << endl; continue; }

//     for(unsigned short itrk = 0; itrk < ntracks; itrk++){
//       // get track pointer
//       StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
//       if(!trk){ continue; }

//       if (!IsTrackWithinJet(jet, trk)) continue;

//       TVector3 mTrkMom;
//       if(doUsePrimTracks) {
//         mTrkMom = trk->pMom();
//       } else {
//         mTrkMom = trk->gMom(mVertex, 0.0);
//       }

//       // track variables
//       double pt = mTrkMom.Perp();
//       double phi = mTrkMom.Phi();
//       double eta = mTrkMom.PseudoRapidity();
//       double px = mTrkMom.x();
//       double py = mTrkMom.y();
//       double pz = mTrkMom.z();
//       double p = mTrkMom.Mag();
//       short charge = trk->charge();

//       hConstPhid0BgLS[centbin]->Fill(phi);
//       hConstPtd0BgLS[centbin]->Fill(pt);

//       double delphiD = dPhi(jetphi, phi);
//       double deletaD = dEta(jeteta, eta);
//       double delRD = dR(delphiD, deletaD);

//       hDeltaPhid0BgLS[centbin]->Fill(delphiD);
//       hDeltaEtad0BgLS[centbin]->Fill(deletaD);
//       hRadd0BgLS[centbin]->Fill(delRD);

//       for (int bin = 0; bin < numberofbins; bin++){
//         if (delRD >= bin*deltar && delRD < (bin+1)*deltar){
//           diffd0distBgLS[bin]++;
//           diffd0shapeBgLS[bin] += pt/jetpt;
//         }
//       } // end of bin loop
//     } // end of track loop

//     for (int bin = 0; bin < numberofbins; bin++){
//       hJetDistd0BgLS[centbin]->Fill((2*bin+1)*deltar/2, diffd0distBgLS[bin]);
//       hJetShaped0BgLS[centbin]->Fill((2*bin+1)*deltar/2, diffd0shapeBgLS[bin]);
//     }
//   } //end of jet loop
// } //end of function

// Centrality bin getter
Int_t StD0Analysis::GetFourCentBin(Double_t scaledCent) const{
  int centbin = -99;
  // get centrality bin number
  if(scaledCent >= 0 && scaledCent <  10.0)  { centbin = 0; }
  else if(scaledCent >= 10.0 && scaledCent <  20.0)                { centbin = 1; }
  else if(scaledCent >= 20.0 && scaledCent <  50.0)                { centbin = 2; }
  else if(scaledCent >= 50.0 && scaledCent <= 80.0)                { centbin = 3; }

  return centbin;
}


// jet function

// void StD0Analysis::WriteInfo(){

//   std::ofstream normalisedinvbetaVp_pi;
//   normalisedinvbetaVp_pi.open("normalisedinvbetaVp_pi.csv", ios::out | ios::app | ios::binary);
//   normalisedinvbetaVp_pi << track_p << "\t" << normalisedinvbeta_for_pi << endl;
//   normalisedinvbetaVp_pi.close();

//   std::ofstream normalisedinvbetaVp_ka;
//   normalisedinvbetaVp_ka.open("normalisedinvbetaVp_ka.csv", ios::out | ios::app | ios::binary);
//   normalisedinvbetaVp_ka << track_p << "\t" << normalisedinvbeta_for_ka << endl;
//   normalisedinvbetaVp_ka.close();

// }
