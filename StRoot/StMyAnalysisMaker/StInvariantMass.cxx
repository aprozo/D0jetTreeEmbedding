// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StInvariantMass.h"
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

ClassImp(StInvariantMass)

//________________________________________________________________________
StInvariantMass::StInvariantMass(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  fRunFlag = 0;       // see StInvariantMass::fRunFlagEnum
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
StInvariantMass::~StInvariantMass()
{ 
  //if (gPtvpPt) delete gPtvpPt;
  if (dEdXvpT) delete dEdXvpT;
  if (invbetavpT) delete invbetavpT;
  // if (event_cuts) delete event_cuts;
  // if (track_cuts) delete track_cuts;

  if (cuthistogram_event) delete cuthistogram_event;
  if (cuthistogram_track) delete cuthistogram_track;
  //gPtvpPt->Write();
  // if (dEdXvpT) delete dEdXvpT;
  if (dEdXvp) delete dEdXvp;

  if (dEdXvpT_pion) delete dEdXvpT_pion;
  if (dEdXvpT_kaon) delete dEdXvpT_kaon;
  if (dEdXvpT_proton) delete dEdXvpT_proton;
  if (dEdXvpT_electron) delete dEdXvpT_electron;

  if (dEdXvp_pion) delete dEdXvp_pion;
  if (dEdXvp_kaon) delete dEdXvp_kaon;
  if (dEdXvp_proton) delete dEdXvp_proton;
  if (dEdXvp_electron) delete dEdXvp_electron;

  if (dEdXthvp_pi) delete dEdXthvp_pi;
  if (dEdXthvp_ka) delete dEdXthvp_ka;
  if (dEdXthvp_pr) delete dEdXthvp_pr;

  if (z_pi) delete z_pi;
  if (z_ka) delete z_ka;
  if (z_pr) delete z_pr;

  // if (invbetavpT) delete invbetavpT;
  if (invbetavpT_tof) delete invbetavpT_tof;
  if (normalised_invbetavpT_tof_pi) delete normalised_invbetavpT_tof_pi;
  if (normalised_invbetavpT_tof_ka) delete normalised_invbetavpT_tof_ka;
  if (normalised_invbetavpT_tof_pr) delete normalised_invbetavpT_tof_pr;

  if (mvpT) delete mvpT;
  if (EvP) delete EvP;

  for (int cbin = 0; cbin < 4; cbin++){
    if (hJetShape[cbin]) delete hJetShape[cbin];
    if (hJetShaped0[cbin]) delete hJetShaped0[cbin];
    if (hJetShaped0BgUL[cbin]) delete hJetShaped0BgUL[cbin];
    if (hJetShaped0BgLS[cbin]) delete hJetShaped0BgLS[cbin];

    if (hJetDist[cbin]) delete hJetDist[cbin];
    if (hJetDistd0[cbin]) delete hJetDistd0[cbin];
    if (hJetDistd0BgUL[cbin]) delete hJetDistd0BgUL[cbin];
    if (hJetDistd0BgLS[cbin]) delete hJetDistd0BgLS[cbin];    
  }



}

//________________________________________________________________________
Int_t StInvariantMass::Init() {
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
Int_t StInvariantMass::Finish() { 

  cout << "StInvariantMass::Finish()\n";

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

  cout<<"End of StInvariantMass::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StInvariantMass::DeclareHistograms() {

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

  cuthistogram_track = new TH1F("cuthistogram_track", "cuthistogram_track", 10, 0, 10);

  dEdXvpT = new TH2F("dEdXvpT", "dEdXvpT", 2000, 0, 20, 1000, 0, 40);

  dEdXvp = new TH2F("dEdXvp", "dEdXvp", 2000, 0, 20, 1000, 0, 40);

  dEdXvpT_pion = new TH2F("dEdXvpT_pion", "dEdXvpT_pion", 2000, 0, 20, 1000, 0, 40);

  dEdXvpT_kaon = new TH2F("dEdXvpT_kaon", "dEdXvpT_kaon", 2000, 0, 20, 1000, 0, 40);

  dEdXvpT_proton = new TH2F("dEdXvpT_proton", "dEdXvpT_proton", 2000, 0, 20, 1000, 0, 40);

  dEdXvpT_electron = new TH2F("dEdXvpT_electron", "dEdXvpT_electron", 2000, 0, 20, 1000, 0, 40);

  dEdXvp_pion = new TH2F("dEdXvp_pion", "dEdXvp_pion", 2000, 0, 20, 1000, 0, 40);

  dEdXvp_kaon = new TH2F("dEdXvp_kaon", "dEdXvp_kaon", 2000, 0, 20, 1000, 0, 40);

  dEdXvp_proton = new TH2F("dEdXvp_proton", "dEdXvp_proton", 2000, 0, 20, 1000, 0, 40);

  dEdXvp_electron = new TH2F("dEdXvp_electron", "dEdXvp_electron", 2000, 0, 20, 1000, 0, 40);

  dEdXthvp_pi = new TH2F("dEdXthvp_pi", "dEdXthvp_pi", 2000, 0, 20, 1000, 0, 40);

  dEdXthvp_ka = new TH2F("dEdXthvp_ka", "dEdXthvp_ka", 2000, 0, 20, 1000, 0, 40);

  dEdXthvp_pr = new TH2F("dEdXthvp_pr", "dEdXthvp_pr", 2000, 0, 20, 1000, 0, 40);

  z_pi = new TH2F("z_pi", "z_pi", 2000, 0, 20, 2000, -10, 10);

  z_ka = new TH2F("z_ka", "z_ka", 2000, 0, 20, 2000, -10, 10);

  z_pr = new TH2F("z_pr", "z_pr", 2000, 0, 20, 2000, -10, 10);

  z_pi_jetconst = new TH2F("z_pi_jetconst", "z_pi_jetconst", 2000, 0, 20, 2000, -10, 10);

  z_ka_jetconst = new TH2F("z_ka_jetconst", "z_ka_jetconst", 2000, 0, 20, 2000, -10, 10);

  invbetavpT = new TH2F("invbetavpT", "invbetavpT", 2000, 0, 20, 1000, 0, 5);

  invbetavpT_tof = new TH2F("invbetavpT_tof", "invbetavpT_tof", 2000, 0, 20, 1000, 0, 10);

  normalised_invbetavpT_tof_pi = new TH2F("normalised_invbetavpT_tof_pi", "normalised_invbetavpT_tof_pi", 2000, 0, 20, 2000, -10, 10);

  normalised_invbetavpT_tof_ka = new TH2F("normalised_invbetavpT_tof_ka", "normalised_invbetavpT_tof_ka", 2000, 0, 20, 2000, -10, 10);

  normalised_invbetavpT_tof_pr = new TH2F("normalised_invbetavpT_tof_pr", "normalised_invbetavpT_tof_pr", 2000, 0, 20, 2000, -10, 10);

  normalised_invbetavpT_tof_pi_jetconst = new TH2F("normalised_invbetavpT_tof_pi_jetconst", "normalised_invbetavpT_tof_pi_jetconst", 2000, 0, 20, 2000, -10, 10);

  normalised_invbetavpT_tof_ka_jetconst = new TH2F("normalised_invbetavpT_tof_ka_jetconst", "normalised_invbetavpT_tof_ka_jetconst", 2000, 0, 20, 2000, -10, 10);

  chosen_pi_from_z_vs_pt = new TH2F("chosen_pi_from_z_vs_pt", "chosen_pi_from_z_vs_pt", 2000, 0, 20, 2000, -10, 10);

  chosen_ka_from_z_vs_pt = new TH2F("chosen_ka_from_z_vs_pt", "chosen_ka_from_z_vs_pt", 2000, 0, 20, 2000, -10, 10);

  chosen_pi_from_invbeta_vs_pt = new TH2F("chosen_pi_from_invbeta_vs_pt", "chosen_pi_from_invbeta_vs_pt", 2000, 0, 20, 2000, -10, 10);

  chosen_ka_from_invbeta_vs_pt = new TH2F("chosen_ka_from_invbeta_vs_pt", "chosen_ka_from_invbeta_vs_pt", 2000, 0, 20, 2000, -10, 10);

  mvpT = new TH2F("mvp", "mvp", 2000, -5, 5, 1000, -1.5, 1.5);

  EvP = new TH2F("EvP", "EvP", 1000, 0, 10, 1000, 0, 10);

  invmass = new TH1F("invmass", "invmass", 1000, 0.5, 2.5);

  invmass_flipped = new TH1F("invmass_flipped", "invmass_flipped", 1000, 0.5, 2.5); 

  kpi_ul_pt = new TH2F("kpi_ul_pt", "kpi_ul_pt", 1000, 0, 50, 1000, 0, 50);

  kpi_ls_pt = new TH2F("kpi_ls_pt", "kpi_ls_pt", 1000, 0, 50, 1000, 0, 50);

  hD0pt = new TH1F("hD0pt", "hD0pt", 1000, 0, 50);

  trackchi2byndf = new TH1F("trackchi2byndf", "trackchi2byndf", 100, 0, 10);
  // Different cuts to try out for now:
  // D0 from PV: Try values: {0.5, 1.0, 1.5, 2.0, 2.5, 3}
  // K from Pi: Try values: {0.5, 1.0, 1.5, 2.0}
  // Combined cuts: Try values: {{0.5, 0.5}, {1,1}, {1.2, 1.2}, {1.5, 1.5}, {2,2}}

  for (int i = 0; i < 6; i++){
    nskavnspi[i] = new TH2F(Form("nskavnspi_pt%i", i),Form("nskavnspi_pt%i", i), 100, -0.5, 3.5, 100, -0.5, 3.5);
    dedxkavdedxpi[i] = new TH2F(Form("dedxkavdedxpi_pt%i", i),Form("dedxkavdedxpi_pt%i", i), 100, -0.5, 3.5, 100, -0.5, 3.5);
    invbetakavinvbetapi[i] = new TH2F(Form("invbetakavinvbetapi_pt%i", i),Form("invbetakavinvbetapi_pt%i", i), 100, -0.5, 3.5, 100, -0.5, 3.5);
  }

  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 4; j++){
      for (int k = 0; k < 4; k++){
        invmasscuts[i][j][k] = new TH1F(Form("invmasscuts_%i_%i_%i", i, j, k),Form("invmasscuts_%i_%i_%i", i, j, k), 1000, 0.5, 2.5); // Histogram for Delta Pi
        invmass_flipped_cuts[i][j][k] = new TH1F(Form("invmass_flipped_cuts%i_%i_%i", i, j, k),Form("invmass_flipped_cuts%i_%i_%i", i, j, k), 1000, 0.5, 2.5); // Histogram for Delta Pi
        invmassbgcuts[i][j][k] = new TH1F(Form("invmassbgcuts%i_%i_%i", i, j, k),Form("invmassbgcuts%i_%i_%i", i, j, k), 1000, 0.5, 2.5); // Histogram for Delta Pi
      }
    }
  }
  // for (int cutnum = 0; cutnum < 27; cutnum++){

  //   invmasscuts[cutnum] = new TH1F(Form("invmasscuts %i", cutnum),Form("invmasscuts %i", cutnum), 1000, 0.5, 2.5); // Histogram for Delta Pi
  //   invmass_flipped_cuts[cutnum] = new TH1F(Form("invmass_flipped_cuts %i", cutnum),Form("invmass_flipped_cuts %i", cutnum), 1000, 0.5, 2.5); // Histogram for Delta Pi
  //   invmassbgcuts[cutnum] = new TH1F(Form("invmassbgcuts %i", cutnum),Form("invmassbgcuts %i", cutnum), 1000, 0.5, 2.5); // Histogram for Delta Pi

  // } 

  invmassvpt = new TH2F("invmassvpt", "invmassvpt", 2000, 0, 20, 1000, 0.5, 2.5);

  invmassbg = new TH1F("invmassbg", "invmassbg", 1000, 0.5, 2.5);

  invmassbgvpt = new TH2F("invmassbgvpt", "invmassbgvpt", 2000, 0, 20, 1000, 0.5, 2.5);

  kfmass = new TH1F("kfmass", "kfmass", 1000, 0.5, 2.5);

  distancepik = new TH1F("distancepik", "distancepik", 10000, -0.05, 0.5);

  deviationpik = new TH1F("deviationpik", "deviationpik", 10000, -0.5, 1000.5);

  distancepid0 = new TH1F("distancepid0", "distancepid0", 10000, -0.5, 1000.5);

  deviationpid0 = new TH1F("deviationpid0", "deviationpid0", 10000, -0.5, 1000.5);

  distancekd0 = new TH1F("distancekd0", "distancekd0", 10000, -0.5, 1000.5);

  deviationkd0 = new TH1F("deviationkd0", "deviationkd0", 10000, -0.5, 1000.5);

  decaylengthd0 = new TH1F("decaylengthd0", "decaylengthd0", 10000, -0.05, 0.5);

  distanced0PV = new TH1F("distanced0PV", "distanced0PV", 10000, -0.05, 0.5);


  // Jetshape histograms 

  for (int cbin = 0; cbin < 4; cbin++){

    JetPhi[cbin] = new TH1F(Form("JetPhi %i", cbin),Form("JetPhi %i", cbin), 2000, -10, 10); // Histogram for Delta Pi

    ConstPhi[cbin] = new TH1F(Form("ConstPhi %i", cbin),Form("ConstPhi %i", cbin), 200, -10, 10); // Histogram for Delta Pi

    JetPt[cbin] = new TH1F(Form("JetPt %i", cbin),Form("JetPt %i", cbin), 2000, 0, 100); // Histogram for Delta Pi

    JetPtCorr[cbin] = new TH1F(Form("JetPtCorr %i", cbin),Form("JetPtCorr %i", cbin), 2000, 0, 100); // Histogram for Delta Pi

    ConstPt[cbin] = new TH1F(Form("ConstPt %i", cbin),Form("ConstPt %i", cbin), 2000, 0, 100); // Histogram for Delta Pi

    DeltaEta[cbin] = new TH1F(Form("DeltaEta %i", cbin),Form("DeltaEta %i", cbin), 200, -1, 1); // Histogram for Delta Eta

    DeltaPhi[cbin] = new TH1F(Form("DeltaPhi %i", cbin),Form("DeltaPhi %i", cbin), 1000, -10, 10); // Histogram for Delta Pi\

    Rad[cbin] = new TH1F(Form("Rad %i", cbin),Form("Rad %i", cbin), 100, 0, 2); // Histogram for dR


    hJetPhid0[cbin] = new TH1F(Form("hJetPhid0 %i", cbin),Form("hJetPhid0 %i", cbin), 2000, -10, 10); // Histogram for Delta Phi

    hConstPhid0[cbin] = new TH1F(Form("hConstPhid0 %i", cbin),Form("hConstPhid0 %i", cbin), 200, -10, 10); // Histogram for Delta Phi

    hJetPtd0[cbin] = new TH1F(Form("hJetPtd0 %i", cbin),Form("hJetPtd0 %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hJetPtd0Corr[cbin] = new TH1F(Form("hJetPtd0Corr %i", cbin),Form("hJetPtd0Corr %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hConstPtd0[cbin] = new TH1F(Form("hConstPtd0 %i", cbin),Form("hConstPtd0 %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hDeltaEtad0[cbin] = new TH1F(Form("hDeltaEtad0 %i", cbin),Form("hDeltaEtad0 %i", cbin), 200, -1, 1); // Histogram for Delta Eta

    hDeltaPhid0[cbin] = new TH1F(Form("hDeltaPhid0 %i", cbin),Form("hDeltaPhid0 %i", cbin), 1000, -10, 10); // Histogram for Delta Phi\

    hRadd0[cbin] = new TH1F(Form("hRadd0 %i", cbin),Form("hRadd0 %i", cbin), 100, 0, 2); // Histogram for dR


    hJetPhid0BgUL[cbin] = new TH1F(Form("hJetPhid0BgUL %i", cbin),Form("hJetPhid0BgUL %i", cbin), 2000, -10, 10); // Histogram for Delta Phi

    hConstPhid0BgUL[cbin] = new TH1F(Form("hConstPhid0BgUL %i", cbin),Form("hConstPhid0BgUL %i", cbin), 200, -10, 10); // Histogram for Delta Phi

    hJetPtd0BgUL[cbin] = new TH1F(Form("hJetPtd0BgUL %i", cbin),Form("hJetPtd0BgUL %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hJetPtd0BgULCorr[cbin] = new TH1F(Form("hJetPtd0BgULCorr %i", cbin),Form("hJetPtd0BgULCorr %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hConstPtd0BgUL[cbin] = new TH1F(Form("hConstPtd0BgUL %i", cbin),Form("hConstPtd0BgUL %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hDeltaEtad0BgUL[cbin] = new TH1F(Form("hDeltaEtad0BgUL %i", cbin),Form("hDeltaEtad0BgUL %i", cbin), 200, -1, 1); // Histogram for Delta Eta

    hDeltaPhid0BgUL[cbin] = new TH1F(Form("hDeltaPhid0BgUL %i", cbin),Form("hDeltaPhid0BgUL %i", cbin), 1000, -10, 10); // Histogram for Delta Phi\

    hRadd0BgUL[cbin] = new TH1F(Form("hRadd0BgUL %i", cbin),Form("hRadd0BgUL %i", cbin), 100, 0, 2); // Histogram for dR


    hJetPhid0BgLS[cbin] = new TH1F(Form("hJetPhid0BgLS %i", cbin),Form("hJetPhid0BgLS %i", cbin), 2000, -10, 10); // Histogram for Delta Phi

    hConstPhid0BgLS[cbin] = new TH1F(Form("hConstPhid0BgLS %i", cbin),Form("hConstPhid0BgLS %i", cbin), 200, -10, 10); // Histogram for Delta Phi

    hJetPtd0BgLS[cbin] = new TH1F(Form("hJetPtd0BgLS %i", cbin),Form("hJetPtd0BgLS %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hJetPtd0BgLSCorr[cbin] = new TH1F(Form("hJetPtd0BgLSCorr %i", cbin),Form("hJetPtd0BgLSCorr %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hConstPtd0BgLS[cbin] = new TH1F(Form("hConstPtd0BgLS %i", cbin),Form("hConstPtd0BgLS %i", cbin), 2000, 0, 100); // Histogram for Delta Phi

    hDeltaEtad0BgLS[cbin] = new TH1F(Form("hDeltaEtad0BgLS %i", cbin),Form("hDeltaEtad0BgLS %i", cbin), 200, -1, 1); // Histogram for Delta Eta

    hDeltaPhid0BgLS[cbin] = new TH1F(Form("hDeltaPhid0BgLS %i", cbin),Form("hDeltaPhid0BgLS %i", cbin), 1000, -10, 10); // Histogram for Delta Phi\
    
    hRadd0BgLS[cbin] = new TH1F(Form("hRadd0BgLS %i", cbin),Form("hRadd0BgLS %i", cbin), 100, 0, 2); // Histogram for dR


    hJetShape[cbin] = new TH1F(Form("hJetShape %i", cbin),Form("hJetShape %i", cbin), numberofbins, 0, R);

    hJetShaped0[cbin] = new TH1F(Form("hJetShaped0 %i", cbin),Form("hJetShaped0 %i", cbin), numberofbins, 0, R);

    hJetShaped0BgUL[cbin] = new TH1F(Form("hJetShaped0BgUL %i", cbin),Form("hJetShaped0BgUL %i", cbin), numberofbins, 0, R);

    hJetShaped0BgLS[cbin] = new TH1F(Form("hJetShaped0BgLS %i", cbin),Form("hJetShaped0BgLS %i", cbin), numberofbins, 0, R);

    hJetDist[cbin] = new TH1F(Form("hJetDist %i", cbin),Form("hJetDist %i", cbin), numberofbins, 0, R);

    hJetDistd0[cbin] = new TH1F(Form("hJetDistd0 %i", cbin),Form("hJetDistd0 %i", cbin), numberofbins, 0, R);

    hJetDistd0BgUL[cbin] = new TH1F(Form("hJetDistd0BgUL %i", cbin),Form("hJetDistd0BgUL %i", cbin), numberofbins, 0, R);

    hJetDistd0BgLS[cbin] = new TH1F(Form("hJetDistd0BgLS %i", cbin),Form("hJetDistd0BgLS %i", cbin), numberofbins, 0, R);
  }
}
//
// write histograms
//_____________________________________________________________________________
void StInvariantMass::WriteHistograms() {

  hCentrality->Write();
  hMultiplicity->Write();

  const char *event_cuts[10] = {"Total Events", "|V_{z} - V_{z(VPD)}| < 6 cm", "nBEMCMatch > 0", "nBTOFMatch > 0", "|V_{z}| < 60 cm", "", "", "", "", ""};
  const char *track_cuts[10] = {"Total Tracks", "nHitsDedx > 20", "#frac{nHitsDedx}{nHitsMax} > 0.52", "DCA < 3 cm", "BEMC Hit", "BTOF Hit", "Charged Tracks", "", "", ""};

  for (int i=1; i <= 10; i++){
    cuthistogram_event->SetBinContent(i, numberofevents[i-1]);
    cuthistogram_event->GetXaxis()->SetBinLabel(i, event_cuts[i-1]);
    cuthistogram_track->SetBinContent(i, numberoftracks[i-1]);
    cuthistogram_track->GetXaxis()->SetBinLabel(i, track_cuts[i-1]);
  }
  // cuthistogram_event->GetXaxis()->LabelsOption(“v”);
  // cuthistogram_track->GetXaxis()->LabelsOption(“v”);

  cuthistogram_event->Write();
  cuthistogram_track->Write();

  trackchi2byndf->Write();
  for (int i = 0; i < 6; i++){
    nskavnspi[i]->Write();
    dedxkavdedxpi[i]->Write();
    invbetakavinvbetapi[i]->Write();
  }
  //gPtvpPt->Write();
  dEdXvpT->Write();
  dEdXvp->Write();

  dEdXvpT_pion->Write();
  dEdXvpT_kaon->Write();
  dEdXvpT_proton->Write();
  dEdXvpT_electron->Write();

  dEdXvp_pion->Write();
  dEdXvp_kaon->Write();
  dEdXvp_proton->Write();
  dEdXvp_electron->Write();

  dEdXthvp_pi->Write();
  dEdXthvp_ka->Write();
  dEdXthvp_pr->Write();

  z_pi->Write();
  z_ka->Write();
  // z_pr->Write();

  z_pi_jetconst->Write();
  z_ka_jetconst->Write();

  invbetavpT->Write();
  invbetavpT_tof->Write();
  normalised_invbetavpT_tof_pi->Write();
  normalised_invbetavpT_tof_ka->Write();
  // normalised_invbetavpT_tof_pr->Write();

  normalised_invbetavpT_tof_pi_jetconst->Write();
  normalised_invbetavpT_tof_ka_jetconst->Write();

  mvpT->Write();
  EvP->Write();

  chosen_pi_from_z_vs_pt->Write();
  chosen_ka_from_z_vs_pt->Write();
  chosen_pi_from_invbeta_vs_pt->Write();
  chosen_ka_from_invbeta_vs_pt->Write();

  invmass->Write();
  invmass_flipped->Write();
  invmassvpt->Write();
  invmassbg->Write();
  invmassbgvpt->Write();
  kfmass->Write();

  distancepik->Write();
  deviationpik->Write();
  distancepid0->Write();
  deviationpid0->Write();
  distancekd0->Write();
  deviationkd0->Write();
  decaylengthd0->Write();
  distanced0PV->Write();


  kpi_ul_pt->Write();
  kpi_ls_pt->Write();

  hD0pt->Write();

  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 3; j++){
      for (int k = 0; k < 4; k++){
        invmasscuts[i][j][k]->Write();
        invmass_flipped_cuts[i][j][k]->Write();
        invmassbgcuts[i][j][k]->Write();
      }
    }
  } 

  for (int cbin = 0; cbin < 4; cbin++){

    JetPhi[cbin]->Write(); ConstPhi[cbin]->Write(); JetPt[cbin]->Write(); JetPtCorr[cbin]->Write(); ConstPt[cbin]->Write(); DeltaEta[cbin]->Write(); DeltaPhi[cbin]->Write(); Rad[cbin]->Write();


    hJetPhid0[cbin]->Write(); hConstPhid0[cbin]->Write(); hJetPtd0[cbin]->Write(); hJetPtd0Corr[cbin]->Write(); hConstPtd0[cbin]->Write(); hDeltaEtad0[cbin]->Write(); hDeltaPhid0[cbin]->Write(); hRadd0[cbin]->Write();


    hJetPhid0BgUL[cbin]->Write(); hConstPhid0BgUL[cbin]->Write(); hJetPtd0BgUL[cbin]->Write(); hJetPtd0BgULCorr[cbin]->Write(); hConstPtd0BgUL[cbin]->Write(); hDeltaEtad0BgUL[cbin]->Write(); hDeltaPhid0BgUL[cbin]->Write(); hRadd0BgUL[cbin]->Write();


    hJetPhid0BgLS[cbin]->Write(); hConstPhid0BgLS[cbin]->Write(); hJetPtd0BgLS[cbin]->Write(); hJetPtd0BgLSCorr[cbin]->Write(); hConstPtd0BgLS[cbin]->Write(); hDeltaEtad0BgLS[cbin]->Write(); hDeltaPhid0BgLS[cbin]->Write(); hRadd0BgLS[cbin]->Write();


    hJetShape[cbin]->Write(); hJetShaped0[cbin]->Write(); hJetShaped0BgUL[cbin]->Write(); hJetShaped0BgLS[cbin]->Write(); hJetDist[cbin]->Write(); hJetDistd0[cbin]->Write(); hJetDistd0BgUL[cbin]->Write(); hJetDistd0BgLS[cbin]->Write();
  }

}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StInvariantMass::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StInvariantMass::Make() {
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

  // check for MB/HT event
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);

  // // Only look for D0s in MB events
  // bool fRunForMB = kFALSE;  // used to differentiate pp and AuAu   
  // if(doppAnalysis)  fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;   
  // if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE; 
  // if(!fRunForMB) return kStOK; 
  
  // ======================== end of Triggers ============================= //

  numberofevents[0]++;

  if (abs(mPicoEvent->primaryVertex().z()- mPicoEvent->vzVpd() > 6)) return kStOk;
  numberofevents[1]++;
   
  if (mPicoEvent->nBEMCMatch()==0) return kStOk;
  numberofevents[2]++;

  if (mPicoEvent->nBTOFMatch() == 0) return kStOk;
  numberofevents[3]++;

  if (abs(mPicoEvent->primaryVertex().z() > 60)) return kStOk;
  numberofevents[4]++;

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

  // TestTracks();

  RunJets();
  // ProcessJetForJetShape();

  return kStOK;
}

// //________________________________________________________________________
// void StInvariantMass::TestTracks(){

//   if(fCorrJetPt) {
//     fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
//     fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
//   } else {
//     fLeadingJet = GetLeadingJet(fJetMakerName);
//     fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
//   }

//   const Int_t ntracks = mPicoDst->numberOfTracks();

//   for(unsigned short itrk = 0; itrk < ntracks; itrk++){
//     // get track pointer
//     StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
//     if(!trk){ continue; }
//     if(!IsAnAcceptableTrack(trk)){ continue; }

//     FillPidHistograms(trk);
//   }

//   for(unsigned short itrk1 = 0; itrk1 < ntracks; itrk1++){
//     // get track pointer
//     StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(itrk1));
//     if(!trk1){ continue; }

//     for(unsigned short itrk2 = itrk1 + 1; itrk2 < ntracks; itrk2++){
//       if (itrk2 >= ntracks) continue;
//       // get track pointer
//       StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(itrk2));
//       if(!trk2){ continue; }

//       StPicoTrackCovMatrix *cov1 = static_cast<StPicoTrackCovMatrix*>(mPicoDst->trackCovMatrix(itrk1));
//       StPicoTrackCovMatrix *cov2 = static_cast<StPicoTrackCovMatrix*>(mPicoDst->trackCovMatrix(itrk2));

//       bool D0, D0BgUnlike, D0BgLike;
//       double mass;
//       TVector3 mom3;

//       IsD0(trk1, trk2, D0, D0BgUnlike, D0BgLike, mass, mom3);
//       IsD0KFParticle(trk1, trk2, cov1, cov2);
//     }
//   }
// }

//________________________________________________________________________
void StInvariantMass::RunJets()
{

  double bfield = mPicoEvent->bField();
 //cout << "Magnetic Field = " << bfield << endl;

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
  
  // const Int_t ntracks = mPicoDst->numberOfTracks();

  // for(unsigned short itrk = 0; itrk < ntracks; itrk++){
  //   // get track pointer
  //   StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
  //   if(!trk){ continue; }
  //   FillPidHistogramsAllTracks(trk);
  // }

  Int_t njets = fJets->GetEntries();
  //if (njets!=0) cout << njets << endl;
  //cout << njets << endl;

  Int_t d0CandidateJet, d0BgCandidateULJet, d0BgCandidateLSJet;

  for (int ijet = 0; ijet < njets; ijet++){

    Bool_t D0Candidate = kFALSE;
    Bool_t D0BgCandidateUL = kFALSE;
    Bool_t D0BgCandidateLS = kFALSE;

    //cout << "Jet #" << ijet << endl;
    // get jet pointer
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;

    //if (abs(jet->Eta()) > 1.) continue;
    double jetpt = jet->Pt();
    double jeteta = jet->Eta();
    double jetphi = jet->Phi();
    double jetarea = jet->Area();
    double corrjetpt = jet->Pt() - jetarea*fRhoVal;
    //cout << "Jet #" << ijet << "\t" << "Jet Pt = " << jetpt << "\t" << corrjetpt << endl; 
    if (!testrun)
    {
      if (corrjetpt < 15.) continue;
      if (jet->GetNumberOfTracks() == 1) continue;
    }

    numberoftracks[0]+=jet->GetNumberOfTracks();

    const Int_t ntracks = mPicoDst->numberOfTracks();
    Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

    // for(unsigned short itrk = 0; itrk < ntracks; itrk++){
    for(unsigned short itrk = 0; itrk < nglobaltracks; itrk++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
      if(!trk){ continue; }

      if (!IsTrackWithinJet(jet, trk)) continue;
      // FillPidHistograms(trk);
      trackchi2byndf->Fill(trk->chi2());

      if (trk->chi2() > 3) continue;

      FillPidHistograms(trk);

      
    }

    // for(int jettrk = 0; jettrk < jet->GetNumberOfTracks(); jettrk++) {
    //   int trackid = jet->TrackAt(jettrk);      

    //   // get jet track pointer
    //   StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
    //   if(!trk){ continue; }
    //   FillPidHistogramsForJetConst(trk);
    // } // track constit loop

    // loop over ALL tracks in PicoDst 
    // for(unsigned short itrk1 = 0; itrk1 < ntracks; itrk1++){
    for(unsigned short itrk1 = 0; itrk1 < nglobaltracks; itrk1++){
      // get track pointer
      StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(itrk1));
      if(!trk1){ continue; }

      if (!IsTrackWithinJet(jet, trk1)) continue;
      if (trk1->chi2() > 3) continue;

      // for(unsigned short itrk2 = itrk1 + 1; itrk2 < ntracks; itrk2++){
      for(unsigned short itrk2 = itrk1 + 1; itrk2 < nglobaltracks; itrk2++){
        if (itrk2 >= ntracks) continue;
        // get track pointer
        StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(itrk2));
        if(!trk2){ continue; }

        if (!IsTrackWithinJet(jet, trk2)) continue;
        if (trk2->chi2() > 3) continue;

        StPicoTrackCovMatrix *cov1 = static_cast<StPicoTrackCovMatrix*>(mPicoDst->trackCovMatrix(itrk1));
        StPicoTrackCovMatrix *cov2 = static_cast<StPicoTrackCovMatrix*>(mPicoDst->trackCovMatrix(itrk2));

        TVector3 mTrk1Mom, mTrk2Mom, mResMom;

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

        mResMom = mTrk1Mom + mTrk2Mom;


        if (mResMom.Perp() < 0.5) continue;
        Double_t phi = standardPhi(mResMom.Phi());
        Double_t eta = mResMom.PseudoRapidity();

        double deltaetad = dEta(jeteta, eta);
        double deltaphid = dPhi(jetphi, phi);
        double deltard = dR(deltaetad, deltaphid);

        if (deltard > fJetRad) continue;

        InvariantMassWithCuts(trk1, trk2, cov1, cov2);
        
        // if(!IsD0KFParticle(trk1, trk2, cov1, cov2)) continue;

        // bool D0, D0BgUnlike, D0BgLike;
        // double mass;
        // TVector3 mom3;

        // IsD0(trk1, trk2, D0, D0BgUnlike, D0BgLike, mass, mom3);
        

        // if ((!D0) || (!D0BgUnlike) || (!D0BgLike)) continue;

        // if (D0) D0Candidate = kTRUE;
        // if (D0BgUnlike) D0BgCandidateUL = kTRUE;
        // if (D0BgLike) D0BgCandidateLS = kTRUE;

        // if ((D0Candidate) && (D0BgCandidateUL) && (D0BgCandidateLS)) break;
      }
      // if ((D0Candidate) && (D0BgCandidateUL) && (D0BgCandidateLS)) break;
    }

    // if (D0Candidate){ 
    //   d0CandidateID.push_back(ijet); 
    //   cout << "ijet D0= " << ijet << endl; 
    //   cout << "Found a D0 Jet" << endl;
    // }
    // if (D0BgCandidateUL){ 
    //   d0BgCandidateULID.push_back(ijet); 
    //   cout << "ijet D0 UL= " << ijet << endl; 
    //   cout << "Found a UL Jet" << endl;
    // }
    // if (D0BgCandidateLS){ 
    //   d0BgCandidateLSID.push_back(ijet); 
    //   cout << "ijet D0 LS= " << ijet << endl; 
    //   cout << "Found a LS Jet" << endl;
    // }
  }
}


Bool_t StInvariantMass::IsAnAcceptableTrack(StPicoTrack *trk){
  if (trk->nHitsDedx() < 15) return kFALSE;
  if (float(trk->nHitsDedx())/float(trk->nHitsMax()) < 0.52) return kFALSE;
  if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 3.) return kFALSE; 
  // if (trk->chi2() > 3) return kFALSE;

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

Bool_t StInvariantMass::IsTrackWithinJet(StJet *jet, StPicoTrack *trk){
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

  // cout << "Eta = " << jeteta << "\t" << eta << endl;
  // cout << "Phi = " << jetphi << "\t" << phi << endl;

  double deltaeta = dEta(jeteta, eta);
  double deltaphi = dPhi(jetphi, phi);
  double deltar = dR(deltaeta, deltaphi);

  // cout << "Delta R" << "\t" << deltar << endl;
  if (deltar < fJetRad) return kTRUE;
  else return kFALSE; 
}

void StInvariantMass::FillPidHistograms(StPicoTrack *trk){
  if (!IsAnAcceptableTrack(trk)) return;

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



  // dedx info
  double dedx = trk->dEdx();
  double dedxresolution = trk->dEdxError();

  // bichsel function approximation
  double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
  double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
  double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

  //z - variables
  double zpi = TMath::Log(dedx/dedxth_pi)/dedxresolution;
  double zka = TMath::Log(dedx/dedxth_ka)/dedxresolution;
  double zpr = TMath::Log(dedx/dedxth_pr)/dedxresolution;

  dEdXvpT->Fill(pt, dedx);
  z_pi->Fill(p, zpi);
  z_ka->Fill(p, zka);
  z_pr->Fill(p, zpr);

  bool toftrack = kFALSE;

  int tof_loc = trk->bTofPidTraitsIndex();

  if (tof_loc >= 0)
  {
    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    if (tofpointer){toftrack = kTRUE;}
  }

  if(toftrack){
    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    double invbeta_from_tof = tofpointer->btofBeta();
    
    // cout << "Beta " << invbeta_from_tof << "\n"; 
    invbeta_from_tof = 1/invbeta_from_tof;

    invbetavpT_tof->Fill(pt, invbeta_from_tof);

    double energy = invbeta_from_tof*p;
    double masssq = (pow(invbeta_from_tof, 2) - 1)*pow(p,2);

    mvpT->Fill(p*charge, masssq);

    double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
    double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
    double norm_invbeta_pr = TMath::Sqrt(pow(Mproton,2)/pow(p,2) + 1);

    double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.011;
    double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.011;
    double normalisedinvbeta_for_pr = (invbeta_from_tof-norm_invbeta_pr);

    normalised_invbetavpT_tof_pi->Fill(p, normalisedinvbeta_for_pi);
    normalised_invbetavpT_tof_ka->Fill(p, normalisedinvbeta_for_ka);
    normalised_invbetavpT_tof_pr->Fill(p, normalisedinvbeta_for_pr);

    if (abs(zpi) < 2 && abs(normalisedinvbeta_for_pi) < 2){chosen_pi_from_invbeta_vs_pt->Fill(p, normalisedinvbeta_for_pi); chosen_pi_from_z_vs_pt->Fill(p, zpi); return;}
    if (abs(zka) < 2 && abs(normalisedinvbeta_for_ka) < 2){chosen_ka_from_invbeta_vs_pt->Fill(p, normalisedinvbeta_for_ka); chosen_ka_from_z_vs_pt->Fill(p, zka); return;}

    for (int i = 0; i < 6; i++){

      double pt_low = i*5;
      double pt_high = (i+1)*5;

      if ((pt > pt_low) && (pt < pt_high)){ 
        invbetakavinvbetapi[i]->Fill(abs(normalisedinvbeta_for_pi), abs(normalisedinvbeta_for_ka));
      }
    }

  }

  else{
    if (abs(zpi) < 2){chosen_pi_from_z_vs_pt->Fill(p, zpi); return;}
    if (abs(zka) < 2){chosen_ka_from_z_vs_pt->Fill(p, zka); return;}
  }

  for (int i = 0; i < 6; i++){

    double pt_low = i*5;
    double pt_high = (i+1)*5;

    if ((pt > pt_low) && (pt < pt_high)){ 
      nskavnspi[i]->Fill(abs(trk->nSigmaPion()), abs(trk->nSigmaKaon()));
      dedxkavdedxpi[i]->Fill(abs(zpi), abs(zka));
    }
  }


}


// void StInvariantMass::FillPidHistogramsForJetConst(StPicoTrack *trk){
//   if (!IsAnAcceptableTrack(trk)) return;

//   TVector3 mTrkMom;
//   if(doUsePrimTracks) {
//     mTrkMom = trk->pMom();
//   } else {
//     mTrkMom = trk->gMom(mVertex, 0.0);
//   }

//   // track variables
//   double pt = mTrkMom.Perp();
//   double phi = mTrkMom.Phi();
//   double eta = mTrkMom.PseudoRapidity();
//   double px = mTrkMom.x();
//   double py = mTrkMom.y();
//   double pz = mTrkMom.z();
//   double p = mTrkMom.Mag();
//   short charge = trk->charge();

//   // dedx info
//   double dedx = trk->dEdx();
//   double dedxresolution = trk->dEdxError();

//   // bichsel function approximation
//   double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
//   double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
//   double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

//   //z - variables
//   double zpi = TMath::Log(dedx/dedxth_pi)/dedxresolution;
//   double zka = TMath::Log(dedx/dedxth_ka)/dedxresolution;
//   double zpr = TMath::Log(dedx/dedxth_pr)/dedxresolution;

//   z_pi_jetconst->Fill(pt, zpi);
//   z_ka_jetconst->Fill(pt, zka);

//   // int tower_loc = trk->bemcTowerIndex();
//   int tof_loc = trk->bTofPidTraitsIndex();

//   if (tof_loc >= 0)
//   {
//     StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));

//     if (tofpointer){
//       // cout << "TOF track found" << endl;
//       double invbeta_from_tof = tofpointer->btofBeta();
      
//       // cout << "Beta " << invbeta_from_tof << "\n"; 
//       invbeta_from_tof = 1/invbeta_from_tof;

//       invbetavpT_tof->Fill(pt*charge, invbeta_from_tof);

//       double energy = invbeta_from_tof*p;
//       double masssq = (pow(invbeta_from_tof, 2) - 1)*pow(p,2);

//       mvpT->Fill(p*charge, masssq);

//       double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
//       double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
//       double norm_invbeta_pr = TMath::Sqrt(pow(Mproton,2)/pow(p,2) + 1);

//       double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.011;
//       double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.012;
//       double normalisedinvbeta_for_pr = (invbeta_from_tof-norm_invbeta_pr);
      
//       normalised_invbetavpT_tof_pi_jetconst->Fill(pt, normalisedinvbeta_for_pi);
//       normalised_invbetavpT_tof_ka_jetconst->Fill(pt, normalisedinvbeta_for_ka);
      
//     }
//   }
// }


void StInvariantMass::IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e){ // NEW PID APPROACH
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

    tof_pion = (abs(normalisedinvbeta_for_pi) <= abs(normalisedinvbeta_for_ka)) ? kTRUE : kFALSE;
    tof_kaon = (abs(normalisedinvbeta_for_pi) > abs(normalisedinvbeta_for_ka)) ? kTRUE : kFALSE;

    if ((tpc_pion != tof_pion) || (tpc_kaon != tof_kaon) || ((abs(normalisedinvbeta_for_pi) > 2) && (abs(normalisedinvbeta_for_ka) > 2))){
      if (tof_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
      else if (tof_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
      else {pid = 0; m = 0.; e = 0.; return;}
    }

    else if (tpc_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tpc_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}
  }

  else{
    if (tpc_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tpc_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}
  }
}

// void StInvariantMass::IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e){
//   if(!IsAnAcceptableTrack(trk)){pid = 0; m = 0.; e = 0.; return;}

//   pid = 0;
//   m = 0.0;
//   e = 0.0;

//   TVector3 mTrkMom;
//   if(doUsePrimTracks) {
//     mTrkMom = trk->pMom();
//   } else {
//     mTrkMom = trk->gMom(mVertex, 0.0);
//   }

//   // track variables
//   double pt = mTrkMom.Perp();
//   double phi = mTrkMom.Phi();
//   double eta = mTrkMom.PseudoRapidity();
//   double px = mTrkMom.x();
//   double py = mTrkMom.y();
//   double pz = mTrkMom.z();
//   double p = mTrkMom.Mag();
//   short charge = trk->charge();

//   double dedx = trk->dEdx();
//   double dedxresolution = trk->dEdxError();

//   // bichsel function approximation
//   double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
//   double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
//   double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

//   //z - variables
//   double zpi = TMath::Log(dedx/dedxth_pi)/dedxresolution;
//   double zka = TMath::Log(dedx/dedxth_ka)/dedxresolution;
//   double zpr = TMath::Log(dedx/dedxth_pr)/dedxresolution;

//   bool toftrack = kFALSE;

//   int tof_loc = trk->bTofPidTraitsIndex();

//   if (tof_loc >= 0)
//   {
//     StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
//     if (tofpointer){toftrack = kTRUE;}
//   }

//   if(toftrack){
//     StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
//     double invbeta_from_tof = tofpointer->btofBeta();
//     invbeta_from_tof = 1/invbeta_from_tof;

//     double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
//     double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
//     double norm_invbeta_pr = TMath::Sqrt(pow(Mproton,2)/pow(p,2) + 1);

//     double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.011;
//     double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.011;
//     double normalisedinvbeta_for_pr = (invbeta_from_tof-norm_invbeta_pr);

//     if (abs(zpi) < 2 && abs(normalisedinvbeta_for_pi) < 2 && charge==1){pid = 1; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
//     if (abs(zka) < 2 && abs(normalisedinvbeta_for_ka) < 2 && charge==1){pid = 2; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
//     if (abs(zpi) < 2 && abs(normalisedinvbeta_for_pi) < 2 && charge==-1){pid = -1; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
//     if (abs(zka) < 2 && abs(normalisedinvbeta_for_ka) < 2 && charge==-1){pid = -2; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
//   }

//   else{
//     if (abs(zpi) < 2 && charge==1){pid = 1; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
//     if (abs(zka) < 2 && charge==1){pid = 2; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
//     if (abs(zpi) < 2 && charge==-1){pid = -1; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
//     if (abs(zka) < 2 && charge==-1){pid = -2; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
//   }
// }


void StInvariantMass::IsWhatParticleFlipped(StPicoTrack *trk, int &pidflipped, double &mflipped, double &eflipped){
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


void StInvariantMass::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum){

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

void StInvariantMass::InvariantMassFlipped(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum){

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

// void StInvariantMass::IsD0(StPicoTrack *trk1, StPicoTrack *trk2, bool &D0, bool &D0BgUnlike, bool &D0BgLike, double &mass, TVector3 &mom3){

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

// void StInvariantMass::ProcessTrackForKF(StPicoTrack *trk, StPicoTrackCovMatrix *picocov, Double_t params[6], Double_t cov[21]){
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

// Bool_t StInvariantMass::IsD0KFParticle(StPicoTrack *trk1, StPicoTrack *trk2, StPicoTrackCovMatrix *cov1, StPicoTrackCovMatrix *cov2){
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

void StInvariantMass::InvariantMassWithCuts(StPicoTrack *trk1, StPicoTrack *trk2, StPicoTrackCovMatrix *cov1, StPicoTrackCovMatrix *cov2){

  double bfield = mPicoEvent->bField();

  KFParticle::SetField(bfield);

  Double_t params1[6], params2[6];
  Double_t covar1[21], covar2[21];

  int charge1 = trk1->charge();
  int charge2 = trk2->charge();

  TVector3 mom1;
  if(doUsePrimTracks) {
    mom1 = trk1->pMom();
  } else {
    mom1 = trk1->gMom(mVertex, Bfield);
  }

  TVector3 mom2;
  if(doUsePrimTracks) {
    mom2 = trk2->pMom();
  } else {
    mom2 = trk2->gMom(mVertex, Bfield);
  }

  TVector3 mom_d0 = mom1 + mom2;
  double d0pt = mom_d0.Perp();

  double pt1 = mom1.Perp();
  double pt2 = mom2.Perp();

  Int_t pid1, pid2;

  Double_t m1, m2, e1, e2;

  IsWhatParticle(trk1, pid1, m1, e1);
  IsWhatParticle(trk2, pid2, m2, e2);

  Int_t particle = pid1*pid2;

  if (abs(particle)!=2){return;} 

  // ProcessTrackForKF(trk1, cov1, params1, covar1);
  // ProcessTrackForKF(trk2, cov2, params2, covar2);

  Int_t pdg1, pdg2;

  switch(pid1){
    case 1: pdg1 = 211;     break;
    case 2: pdg1 = 321;     break;
    // case 3: pdg1 = 2212;    break;
    case -1:  pdg1 = -211;   break;
    case -2:  pdg1 = -321;   break;
    // case -3:  pdg1 = -2212;  break;
  }

  switch(pid2){
    case 1: pdg2 = 211;     break;
    case 2: pdg2 = 321;     break;
    // case 3: pdg2 = 2212;    break;
    case -1:  pdg2 = -211;   break;
    case -2:  pdg2 = -321;   break;
    // case -3:  pdg2 = -2212;  break;
  }

  Double_t Cuts[4] = {0.3, 0.4, 0.5, 0.6};
  Double_t ptcuts[3] = {3.0, 4.0, 5.0};
  Double_t d0ptcuts[4] = {2.0, 3.0, 5.0, 6.0};

  //cout << " =========================================================================== " << endl;

  double dpik = (trk1->origin() - trk2->origin()).Mag();
  distancepik->Fill(dpik);

  TVector3 mom_cord = 0.5*(trk1->origin() + trk2->origin());
  double distanced0fromPV = (mom_cord - mPicoEvent->primaryVertex()).Mag();
  distanced0PV->Fill(distanced0fromPV);

  int dummy;
  double mass, massflipped;
  TVector3 dummy_vec;
  InvariantMass(trk1, trk2, dummy, mass, dummy_vec);

  InvariantMassFlipped(trk1, trk2, dummy, massflipped, dummy_vec);


  // D0 cuts only
  if(particle==-2){
    kpi_ul_pt->Fill(pt1, pt2);

    hD0pt->Fill(d0pt);

    invmass->Fill(mass);
    invmass_flipped->Fill(massflipped);

    for (int i = 0; i < 4; i++){ // VertexCuts
      if (dpik > Cuts[i]) continue;
      for (int j = 0; j < 3; j++){ // PtCuts
        if (pt1 + pt2 < ptcuts[j]) continue;
        for (int k = 0; k < 4; k++){ // D0PtCuts
          if (d0pt < d0ptcuts[k]) continue;
          invmasscuts[i][j][k]->Fill(mass);
          invmass_flipped_cuts[i][j][k]->Fill(massflipped);
        }
      }
    }
  }

  if(particle==2){
    kpi_ls_pt->Fill(pt1, pt2);

    invmassbg->Fill(mass);

    for (int i = 0; i < 4; i++){ // VertexCuts
      if (dpik > Cuts[i]) continue;
      for (int j = 0; j < 3; j++){ // PtCuts
        if (pt1 + pt2 < ptcuts[j]) continue;
        for (int k = 0; k < 4; k++){ // D0PtCuts
          if (d0pt < d0ptcuts[k]) continue;
          invmassbgcuts[i][j][k]->Fill(mass);
        }
      }
    }
  }
}

// void StInvariantMass::ProcessJetForJetShape(){

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
Int_t StInvariantMass::GetFourCentBin(Double_t scaledCent) const{
  int centbin = -99;
  // get centrality bin number
  if(scaledCent >= 0 && scaledCent <  10.0)  { centbin = 0; }
  else if(scaledCent >= 10.0 && scaledCent <  20.0)                { centbin = 1; }
  else if(scaledCent >= 20.0 && scaledCent <  50.0)                { centbin = 2; }
  else if(scaledCent >= 50.0 && scaledCent <= 80.0)                { centbin = 3; }

  return centbin;
}


// jet function

// void StInvariantMass::WriteInfo(){

//   std::ofstream normalisedinvbetaVp_pi;
//   normalisedinvbetaVp_pi.open("normalisedinvbetaVp_pi.csv", ios::out | ios::app | ios::binary);
//   normalisedinvbetaVp_pi << track_p << "\t" << normalisedinvbeta_for_pi << endl;
//   normalisedinvbetaVp_pi.close();

//   std::ofstream normalisedinvbetaVp_ka;
//   normalisedinvbetaVp_ka.open("normalisedinvbetaVp_ka.csv", ios::out | ios::app | ios::binary);
//   normalisedinvbetaVp_ka << track_p << "\t" << normalisedinvbeta_for_ka << endl;
//   normalisedinvbetaVp_ka.close();

// }
