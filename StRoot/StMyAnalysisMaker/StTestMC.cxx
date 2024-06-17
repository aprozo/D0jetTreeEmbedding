// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StTestMC.h"
#include "StMemStat.h"
#include "phys_constants.h"
#include <limits>
#include "math.h"

// ROOT includes
#include "TF1.h"
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

#include "StPhysicalHelix.hh"
#include "StThreeVectorF.hh"

#include "StBTofUtil/tofPathLength.hh"

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

// #include "StBichsel/Bichsel.h"


// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StTestMC)

//________________________________________________________________________
StTestMC::StTestMC(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StTestMC::fRunFlagEnum
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
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  fKaonMomResolution = NULL;
  fPionMomResolution = NULL;
  fProtonMomResolution = NULL;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

}

//
//________________________________________________________________________
StTestMC::~StTestMC()
{ /*  */
  // destructor
  if(hCentrality)  delete hCentrality;
  if(hMultiplicity)delete hMultiplicity;
  if(hJetPt)       delete hJetPt;
  if(hJetCorrPt)   delete hJetCorrPt;

  if(mEmcPosition) delete mEmcPosition;

//______________________Test Code #FIXME ________________________________
  if(hTrackPt) delete hTrackPt;
  if(hTrackEt) delete hTrackEt;

}

//
//________________________________________________________________________
Int_t StTestMC::Init() {
  StJetFrameworkPicoBase::Init();

  // declare histograms
  DeclareHistograms();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  //fJets->SetOwner(kTRUE);

  // TString vertexOutName = mOutName;
  // vertexOutName.ReplaceAll(".root", "");

  // vertexfile.open(Form("Vertex_%s.txt", vertexOutName.Data()));

  TFile f("/star/u/droy1/Y2019/STAR/Momentum_resolution_SL16d.root");
  fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPion");
  fKaonMomResolution = (TF1*)f.Get("fKaon")->Clone("fKaon");
  fProtonMomResolution = (TF1*)f.Get("fProton")->Clone("fProton");
    
  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StTestMC::Finish() { 
  cout << "StTestMC::Finish()\n";
  // vertexfile.close();

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
//    fout->mkdir(GetName());
//    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }
  
  cout<<"End of StTestMC::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StTestMC::DeclareHistograms() {
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

  // histograms
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);

  hD0Daughters = new TH1F("hD0Daughters", "Number of D0 Daughters", 11, -0.5, 10.5);
  hD0DecayLength = new TH1F("hD0DecayLength", "D0 Decay Length", 10000, 0.0, 0.1);
  

  hD0MCMass = new TH1F ("hD0MCMass","hD0MCMass",400,1.0, 5.0);

  //___________________TEST Code #FIXME____________________________________

  
  hTrackPt = new TH2F("hTrackPt", "No. of tracks", 200, -2.0, 2.0, 200, -5.0, 5.0);
  hTrackEt = new TH2F("hTrackEt", "No. of tracks", 200, -2.0, 2.0, 200, -5.0, 5.0);

  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 200, -100., 100.);

  hDCA = new TH1F("hDCA", "hDCA", 200, 0, 2);

  hTPCPt = new TH1F("hTPCPt", "TPC hits vs Pt", 32, -0.5, 15.5);
  hHFTPt = new TH1F("hHFTPt", "HFT hits vs Pt", 32, -0.5, 15.5);
  hMCPt = new TH1F("hMCPt", "MC hits vs Pt", 32, -0.5, 15.5);
  hGenPt = new TH1F("hGenPt", "Gen vs Pt", 32, -0.5, 15.5);

  hTPCEta = new TH1F("hTPCEta", "TPC hits vs Eta", 40, -1.0, 1.0);
  hHFTEta = new TH1F("hHFTEta", "HFT hits vs Eta", 40, -1.0, 1.0);
  hMCEta = new TH1F("hMCEta", "MC hits vs Eta", 40, -1.0, 1.0);
  hGenEta = new TH1F("hGenEta", "Gen vs Eta", 40, -1.0, 1.0);

  hTPCPhi = new TH1F("hTPCPhi", "TPC hits vs Phi", 80, -10, 10);
  hHFTPhi = new TH1F("hHFTPhi", "HFT hits vs Phi", 80, -10, 10);
  hMCPhi = new TH1F("hMCPhi", "MC hits vs Phi", 80, -10, 10);
  hGenPhi = new TH1F("hGenPhi", "Gen vs Phi", 80, -10, 10);

  hMatchedRecoTrack = new TH1F("hMatchedRecoTrack", "hMatchedRecoTrack", 100, 0, 100);
  hRecoTrack = new TH1F("hRecoTrack", "hRecoTrack", 100, 0, 100);

  ///////////////////////// D0 Kaon Pion Histograms ///////////////////////////////////////////////

  hMCD0KaonPt = new TH1F("hMCD0KaonPt", "MC D0 Kaon Pt", 32, -0.5, 15.5);
  hMCD0PionPt = new TH1F("hMCD0PionPt", "MC D0 Pion Pt", 32, -0.5, 15.5);

  hRecoD0KaonPt = new TH1F("hRecoD0KaonPt", "Reco D0 Kaon Pt", 32, -0.5, 15.5);
  hRecoD0PionPt = new TH1F("hRecoD0PionPt", "Reco D0 Pion Pt", 32, -0.5, 15.5);

  hRecoVMCD0KaonPt = new TH2F("hRecoVMCD0KaonPt", "RecovMC D0 Kaon Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);
  hRecoVMCD0PionPt = new TH2F("hRecoVMCD0PionPt", "RecovMC D0 Pion Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);

  hMCD0KaonPhi = new TH1F("hMCD0KaonPhi", "MC D0 Kaon Phi", 80, -10, 10);
  hMCD0PionPhi = new TH1F("hMCD0PionPhi", "MC D0 Pion Phi", 80, -10, 10);

  hRecoD0KaonPhi  = new TH1F("hRecoD0KaonPhi", "Reco D0 Kaon Phi", 80, -10, 10);
  hRecoD0PionPhi  = new TH1F("hRecoD0PionPhi", "Reco D0 Pion Phi", 80, -10, 10);

  hRecoVMCD0KaonPhi = new TH2F("hRecoVMCD0KaonPhi", "RecovMC D0 Kaon Phi", 80, -10, 10, 80, -10, 10);
  hRecoVMCD0PionPhi = new TH2F("hRecoVMCD0PionPhi", "RecovMC D0 Pion Phi", 80, -10, 10, 80, -10, 10);

  hMCD0KaonEta = new TH1F("hMCD0KaonEta", "MC D0 Kaon Eta", 40, -1.0, 1.0);
  hMCD0PionEta = new TH1F("hMCD0PionEta", "MC D0 Pion Eta", 40, -1.0, 1.0);

  hRecoD0KaonEta = new TH1F("hRecoD0KaonEta", "Reco D0 Kaon Eta", 40, -1.0, 1.0);
  hRecoD0PionEta = new TH1F("hRecoD0PionEta", "Reco D0 Pion Eta", 40, -1.0, 1.0);

  hRecoVMCD0KaonEta = new TH2F("hRecoVMCD0KaonEta", "RecovMC D0 Kaon Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);
  hRecoVMCD0PionEta = new TH2F("hRecoVMCD0PionEta", "RecovMC D0 Pion Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);


  ////////////////////////// All Kaon Pion Histograms /////////////////////////////////////////////

  hMCD0Eta = new TH1F("hMCD0Eta", "MC D0 Eta", 200, -2, 2);

  hMCKaonPt = new TH1F("hMCKaonPt", "MC Kaon Pt", 32, -0.5, 15.5);
  hMCPionPt = new TH1F("hMCPionPt", "MC Pion Pt", 32, -0.5, 15.5);

  hRecoKaonPt = new TH1F("hRecoKaonPt", "Reco Kaon Pt", 32, -0.5, 15.5);
  hRecoPionPt = new TH1F("hRecoPionPt", "Reco Pion Pt", 32, -0.5, 15.5);

  hRecoVMCKaonPt = new TH2F("hRecoVMCKaonPt", "RecovMC Kaon Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);
  hRecoVMCPionPt = new TH2F("hRecoVMCPionPt", "RecovMC Pion Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);

  hMCKaonPhi = new TH1F("hMCKaonPhi", "MC Kaon Phi", 80, -10, 10);
  hMCPionPhi = new TH1F("hMCPionPhi", "MC Pion Phi", 80, -10, 10);

  hRecoKaonPhi  = new TH1F("hRecoKaonPhi", "Reco Kaon Phi", 80, -10, 10);
  hRecoPionPhi  = new TH1F("hRecoPionPhi", "Reco Pion Phi", 80, -10, 10);

  hRecoVMCKaonPhi = new TH2F("hRecoVMCKaonPhi", "RecovMC Kaon Phi", 80, -10, 10, 80, -10, 10);
  hRecoVMCPionPhi = new TH2F("hRecoVMCPionPhi", "RecovMC Pion Phi", 80, -10, 10, 80, -10, 10);

  hMCKaonEta = new TH1F("hMCKaonEta", "MC Kaon Eta", 40, -1.0, 1.0);
  hMCPionEta = new TH1F("hMCPionEta", "MC Pion Eta", 40, -1.0, 1.0);

  hRecoKaonEta = new TH1F("hRecoKaonEta", "Reco Kaon Eta", 40, -1.0, 1.0);
  hRecoPionEta = new TH1F("hRecoPionEta", "Reco Pion Eta", 40, -1.0, 1.0);

  hRecoVMCKaonEta = new TH2F("hRecoVMCKaonEta", "RecovMC Kaon Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);
  hRecoVMCPionEta = new TH2F("hRecoVMCPionEta", "RecovMC Pion Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);

  ////////////////////////// All Kaon Pion (TPC) PID Histograms //////////////////////////////////////////

  hRecoTPCKaonPt = new TH1F("hRecoTPCKaonPt", "Reco TPC Kaon Pt", 32, -0.5, 15.5);
  hRecoTPCPionPt = new TH1F("hRecoTPCPionPt", "Reco TPC Pion Pt", 32, -0.5, 15.5);

  hRecoVMCTPCKaonPt = new TH2F("hRecoVMCTPCKaonPt", "RecovMC TPC Kaon Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);
  hRecoVMCTPCPionPt = new TH2F("hRecoVMCTPCPionPt", "RecovMC TPC Pion Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);

  hRecoTPCKaonPhi  = new TH1F("hRecoTPCKaonPhi", "Reco TPC Kaon Phi", 80, -10, 10);
  hRecoTPCPionPhi  = new TH1F("hRecoTPCPionPhi", "Reco TPC Pion Phi", 80, -10, 10);

  hRecoVMCTPCKaonPhi = new TH2F("hRecoVMCTPCKaonPhi", "RecovMC TPC Kaon Phi", 80, -10, 10, 80, -10, 10);
  hRecoVMCTPCPionPhi = new TH2F("hRecoVMCTPCPionPhi", "RecovMC TPC Pion Phi", 80, -10, 10, 80, -10, 10);

  hRecoTPCKaonEta = new TH1F("hRecoTPCKaonEta", "Reco TPC Kaon Eta", 40, -1.0, 1.0);
  hRecoTPCPionEta = new TH1F("hRecoTPCPionEta", "Reco TPC Pion Eta", 40, -1.0, 1.0);

  hRecoVMCTPCKaonEta = new TH2F("hRecoVMCTPCKaonEta", "RecovMC TPC Kaon Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);
  hRecoVMCTPCPionEta = new TH2F("hRecoVMCTPCPionEta", "RecovMC TPC Pion Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);

  ////////////////////////// All Kaon Pion (TPC + TOF) PID Histograms //////////////////////////////////////////

  hRecoTOFKaonPt = new TH1F("hRecoTOFKaonPt", "Reco TOF Kaon Pt", 32, -0.5, 15.5);
  hRecoTOFPionPt = new TH1F("hRecoTOFPionPt", "Reco TOF Pion Pt", 32, -0.5, 15.5);

  hRecoVMCTOFKaonPt = new TH2F("hRecoVMCTOFKaonPt", "RecovMC TOF Kaon Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);
  hRecoVMCTOFPionPt = new TH2F("hRecoVMCTOFPionPt", "RecovMC TOF Pion Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);

  hRecoTOFKaonPhi  = new TH1F("hRecoTOFKaonPhi", "Reco TOF Kaon Phi", 80, -10, 10);
  hRecoTOFPionPhi  = new TH1F("hRecoTOFPionPhi", "Reco TOF Pion Phi", 80, -10, 10);

  hRecoVMCTOFKaonPhi = new TH2F("hRecoVMCTOFKaonPhi", "RecovMC TOF Kaon Phi", 80, -10, 10, 80, -10, 10);
  hRecoVMCTOFPionPhi = new TH2F("hRecoVMCTOFPionPhi", "RecovMC TOF Pion Phi", 80, -10, 10, 80, -10, 10);

  hRecoTOFKaonEta = new TH1F("hRecoTOFKaonEta", "Reco TOF Kaon Eta", 40, -1.0, 1.0);
  hRecoTOFPionEta = new TH1F("hRecoTOFPionEta", "Reco TOF Pion Eta", 40, -1.0, 1.0);

  hRecoVMCTOFKaonEta = new TH2F("hRecoVMCTOFKaonEta", "RecovMC TOF Kaon Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);
  hRecoVMCTOFPionEta = new TH2F("hRecoVMCTOFPionEta", "RecovMC TOF Pion Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);

  ////////////////////////// All Kaon Pion (TPC + HFT) PID Histograms //////////////////////////////////////////

  hRecoHFTKaonPt = new TH1F("hRecoHFTKaonPt", "Reco HFT Kaon Pt", 32, -0.5, 15.5);
  hRecoHFTPionPt = new TH1F("hRecoHFTPionPt", "Reco HFT Pion Pt", 32, -0.5, 15.5);

  hRecoVMCHFTKaonPt = new TH2F("hRecoVMCHFTKaonPt", "RecovMC HFT Kaon Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);
  hRecoVMCHFTPionPt = new TH2F("hRecoVMCHFTPionPt", "RecovMC HFT Pion Pt", 64, -0.5, 15.5, 64, -0.5, 15.5);

  hRecoHFTKaonPhi  = new TH1F("hRecoHFTKaonPhi", "Reco HFT Kaon Phi", 80, -10, 10);
  hRecoHFTPionPhi  = new TH1F("hRecoHFTPionPhi", "Reco HFT Pion Phi", 80, -10, 10);

  hRecoVMCHFTKaonPhi = new TH2F("hRecoVMCHFTKaonPhi", "RecovMC HFT Kaon Phi", 80, -10, 10, 80, -10, 10);
  hRecoVMCHFTPionPhi = new TH2F("hRecoVMCHFTPionPhi", "RecovMC HFT Pion Phi", 80, -10, 10, 80, -10, 10);

  hRecoHFTKaonEta = new TH1F("hRecoHFTKaonEta", "Reco HFT Kaon Eta", 40, -1.0, 1.0);
  hRecoHFTPionEta = new TH1F("hRecoHFTPionEta", "Reco HFT Pion Eta", 40, -1.0, 1.0);

  hRecoVMCHFTKaonEta = new TH2F("hRecoVMCHFTKaonEta", "RecovMC HFT Kaon Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);
  hRecoVMCHFTPionEta = new TH2F("hRecoVMCHFTPionEta", "RecovMC HFT Pion Eta", 40, -1.0, 1.0, 40, -1.0, 1.0);

  hD0Mass= new TH1F ("hD0Mass","hD0Mass",400,1.7,2.1);

}
//
// write histograms
//_____________________________________________________________________________
void StTestMC::WriteHistograms() {
  // writing of histograms done here
  // hCentrality->Write();
  // hMultiplicity->Write();
  // hJetPt->Write();
  // hJetCorrPt->Write();
  // hTrackPt->SetOption("LEGO2");
  // hTrackEt->SetOption("LEGO2");
  // hTrackPt->Write();
  // hTrackEt->Write();

  hEventZVertex->Write();
  hDCA->Write();

  hD0Daughters->Write();
  hD0DecayLength->Write();

  hMCD0Eta->Write();

  hD0MCMass->Write();

  // hTPCPt->Write();
  // hTPCEta->Write();
  // hTPCPhi->Write();

  // hHFTPt->Write();
  // hHFTEta->Write();
  // hHFTPhi->Write();

  // hMCPt->Write();
  // hMCEta->Write();
  // hMCPhi->Write();

  // hGenPt->Write();
  // hGenEta->Write();
  // hGenPhi->Write();
  hD0Mass->Write();
  hMatchedRecoTrack->Write();
  hRecoTrack->Write();
///////////////////////// D0 Kaon Pion Histograms ///////////////////////////////////////////////

    hMCD0KaonPt->Write();
    hMCD0PionPt->Write();

    hRecoD0KaonPt->Write();
    hRecoD0PionPt->Write();

    hRecoVMCD0KaonPt->Write();
    hRecoVMCD0PionPt->Write();

    hMCD0KaonPhi->Write();
    hMCD0PionPhi->Write();

    hRecoD0KaonPhi->Write();
    hRecoD0PionPhi->Write();

    hRecoVMCD0KaonPhi->Write();
    hRecoVMCD0PionPhi->Write();

    hMCD0KaonEta->Write();
    hMCD0PionEta->Write();

    hRecoD0KaonEta->Write();
    hRecoD0PionEta->Write();

    hRecoVMCD0KaonEta->Write();
    hRecoVMCD0PionEta->Write();

    ////////////////////////// All Kaon Pion Histograms //////////////////////////////////////////

    hMCKaonPt->Write();
    hMCPionPt->Write();

    hRecoKaonPt->Write();
    hRecoPionPt->Write();

    hRecoVMCKaonPt->Write();
    hRecoVMCPionPt->Write();

    hMCKaonPhi->Write();
    hMCPionPhi->Write();

    hRecoKaonPhi->Write();
    hRecoPionPhi->Write();

    hRecoVMCKaonPhi->Write();
    hRecoVMCPionPhi->Write();

    hMCKaonEta->Write();
    hMCPionEta->Write();

    hRecoKaonEta->Write();
    hRecoPionEta->Write();

    hRecoVMCKaonEta->Write();
    hRecoVMCPionEta->Write();

    ////////////////////////// All Kaon Pion (TPC) PID Histograms //////////////////////////////////////////

    hRecoTPCKaonPt->Write();
    hRecoTPCPionPt->Write();

    hRecoVMCTPCKaonPt->Write();
    hRecoVMCTPCPionPt->Write();

    hRecoTPCKaonPhi->Write();
    hRecoTPCPionPhi->Write();

    hRecoVMCTPCKaonPhi->Write();
    hRecoVMCTPCPionPhi->Write();

    hRecoTPCKaonEta->Write();
    hRecoTPCPionEta->Write();

    hRecoVMCTPCKaonEta->Write();
    hRecoVMCTPCPionEta->Write();

    ////////////////////////// All Kaon Pion (TPC + TOF) PID Histograms //////////////////////////////////////////

    hRecoTOFKaonPt->Write();
    hRecoTOFPionPt->Write();

    hRecoVMCTOFKaonPt->Write();
    hRecoVMCTOFPionPt->Write();

    hRecoTOFKaonPhi->Write();
    hRecoTOFPionPhi->Write();

    hRecoVMCTOFKaonPhi->Write();
    hRecoVMCTOFPionPhi->Write();

    hRecoTOFKaonEta->Write();
    hRecoTOFPionEta->Write();

    hRecoVMCTOFKaonEta->Write();
    hRecoVMCTOFPionEta->Write();

    ////////////////////////// All Kaon Pion (TPC + HFT) PID Histograms //////////////////////////////////////////

    hRecoHFTKaonPt->Write();
    hRecoHFTPionPt->Write();

    hRecoVMCHFTKaonPt->Write();
    hRecoVMCHFTPionPt->Write();

    hRecoHFTKaonPhi->Write();
    hRecoHFTPionPhi->Write();

    hRecoVMCHFTKaonPhi->Write();
    hRecoVMCHFTPionPhi->Write();

    hRecoHFTKaonEta->Write();
    hRecoHFTPionEta->Write();

    hRecoVMCHFTKaonEta->Write();
    hRecoVMCHFTPionEta->Write();

}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StTestMC::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StTestMC::Make() {
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

  double zVtx_VPD = mPicoEvent->vzVpd();

  // cout << "VPD Vtx: " << zVtx_VPD << endl;

  hEventZVertex->Fill(zVtx);

  // -0.224548 -0.0349843  15.9442

  // if (mVertex.x() != -0.224548 || mVertex.y() != -0.0349843 || mVertex.z() != 15.9442) return kStOk;

  // if (mPicoEvent->eventId() != 1983595) return kStOk;
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


  if (abs(zVtx) > 6) return kStOk;
  // numberofevents[2]++;

  // if (abs(zVtx - zVtx_VPD) > 3) return kStOk;
  // numberofevents[3]++;

  /*

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
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 
  */
  // cout << mPicoEvent->runId() << "\t" << mPicoEvent->eventId()<< endl;

  // vertexfile << mPicoEvent->runId() << "\t" << mPicoEvent->eventId() << "\t" << mVertex.x() << "\t" << mVertex.y() << "\t" << mVertex.z() << endl;
  // cout << "Primary Vertex : " << mVertex.x() << "\t" << mVertex.y() << "\t" << mVertex.z() << endl;

  // cout << "Centrality : " << cent16 << endl;

  // vertexfile << mPicoEvent->runId() << "\t" << mPicoEvent->eventId() << "\t" << mVertex.x() << "\t" << mVertex.y() << "\t" << mVertex.z() << endl;

  // cout << " ====================================================================================================== " << endl;

  /*
  // fill histograms
  hCentrality->Fill(fCentralityScaled);
  hMultiplicity->Fill(refCorr2);

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // ============================ end of CENTRALITY ============================== //

  */

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  // FillEmcTriggers();

  // // get trigger IDs from PicoEvent class and loop over them
  // vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  // if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"EventTriggers: ";
  // for(unsigned int i=0; i<mytriggers.size(); i++) {
  //   if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
  // }
  // if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<endl;

  // // check for MB/HT event
  // bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  // bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  // bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  // bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);
  // ======================== end of Triggers ============================= //

  // // =========================== JetMaker =============================== //
  // // get JetMaker pointer
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

  // // ============================= RhoMaker ============================== //
  // // get RhoMaker pointer from event: old names "StRho_JetsBG"
  // RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerName));
  // const char *fRhoMakerNameCh = fRhoMakerName;
  // if(!RhoMaker) {
  //   LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
  //   return kStWarn;
  // }

  // // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  // fRho = static_cast<StRhoParameter*>(RhoMaker->GetRho());
  // if(!fRho) {
  //   LOG_WARN << Form("Couldn't get fRho object! ") << endm;
  //   return kStWarn;    
  // } 
  
  // // get rho/area value from rho object     fRho->ls("");
  // fRhoVal = fRho->GetVal();
  // // =======================================================================

  // // get number of jets, tracks, and global tracks in events
  // Int_t njets = fJets->GetEntries();
  // const Int_t ntracks = mPicoDst->numberOfTracks();
  // Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // // run Jets:
  // RunJets();

  // run Tracks:
  RunTracks();

  // run Towers:
  // RunTowers();

  return kStOK;
}

//
//
//________________________________________________________________________
void StTestMC::RunTracks()
{
  // constants: assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV
  unsigned int ntracks = mPicoDst->numberOfTracks();
  unsigned int nmctracks = mPicoDst->numberOfMcTracks();
  const int nmcvertices = mPicoDst->numberOfMcVertices();
  vector<int> tracksgeantidsfromvertex[nmcvertices];

  // cout << ntracks << "\t" << nmctracks << "\t" << nmcvertices << endl;

  // cout << "PRIMARY VERTEX : " << mPicoEvent->primaryVertex().X() << "\t" << mPicoEvent->primaryVertex().Y() << "\t" << mPicoEvent->primaryVertex().Z() << endl;

  for(unsigned short iTracks = 0; iTracks < nmctracks; iTracks++){
    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(iTracks));

    // cout << "iTracks = " << iTracks << "\t" << mctrk->id() << "\t" << mctrk->p().X() << "\t" << mctrk->p().Y() << endl;

    int idvx = mctrk->idVtxStart();
    // if (idvx <= 1) continue;
    if (idvx >= nmcvertices) continue;

    // if (mctrk->isFromShower()) continue;

    tracksgeantidsfromvertex[idvx].push_back(mctrk->geantId());
  }

  int arrSize = sizeof(tracksgeantidsfromvertex)/sizeof(tracksgeantidsfromvertex[0]);


   for (unsigned short track1 = 0; track1 < mPicoDst->numberOfMcTracks(); track1++){
      StPicoMcTrack *mctrk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track1));

      if (mctrk1->geantId() == 37 || mctrk1->geantId() == 38) {

        cout << mctrk1->idVtxStart() << "\t" << mctrk1->idVtxStop() << "\t" << tracksgeantidsfromvertex[mctrk1->idVtxStop()].size() << endl;

        hMCD0Eta->Fill(mctrk1->eta());

        hD0Daughters->Fill(tracksgeantidsfromvertex[mctrk1->idVtxStop()].size()); 

        StPicoMcVertex *d0vxstart = static_cast<StPicoMcVertex*>(mPicoDst->mcVertex(mctrk1->idVtxStart() - 1));
        StPicoMcVertex *d0vxstop  = static_cast<StPicoMcVertex*>(mPicoDst->mcVertex(mctrk1->idVtxStop() - 1));

        hD0DecayLength->Fill((d0vxstart->position() - d0vxstop->position()).Mag());

      }
    }

  for (unsigned short track1 = 0; track1 < nmctracks; track1++){
    StPicoMcTrack *mctrk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track1));

    if (mctrk1->geantId() == 8 || mctrk1->geantId() == 9 || mctrk1->geantId() == 11 || mctrk1->geantId() == 12) {

      if (IsD0MCTrack(track1)){
        // cout << "TRACK : "  << mctrk1->geantId() << "\t" << mctrk1->pt() << "\t" << mctrk1->eta() <<  endl;

        int idvxstart = mctrk1->idVtxStart() -  1;

        StPicoMcVertex *d0vx = static_cast<StPicoMcVertex*>(mPicoDst->mcVertex(idvxstart));

        // cout << "VERTEX : " << d0vx->position().X() << "\t" << d0vx->position().Y() << "\t" << d0vx->position().Z() << endl;

        // cout << "DCA : " << (mPicoEvent->primaryVertex() - d0vx->position()).Mag() << endl;

        hDCA->Fill((mPicoEvent->primaryVertex() - d0vx->position()).Mag());

      } 

    }

    // cout << "Number of D0s = " << D0Counter << endl;

    if (mctrk1->geantId() == 8){

      for (unsigned short track2 = 0; track2 < nmctracks; track2++){
        StPicoMcTrack *mctrk2 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track2));

        if (mctrk2->geantId() == 12){

          double mass1 = Mpion;
          double mass2 = Mkaon;

          TLorentzVector mTrk1Mom = mctrk1->fourMomentum();
          TLorentzVector mTrk2Mom = mctrk2->fourMomentum();

          double invmass = TMath::Sqrt(pow(mass1,2) + pow(mass2,2) + 2*(mTrk1Mom*mTrk2Mom));

          hD0MCMass->Fill(invmass);

          // cout << "InvMass = " << invmass << endl;
        }

      }

    }

    if (mctrk1->geantId() == 9){

      for (unsigned short track2 = 0; track2 < nmctracks; track2++){
        StPicoMcTrack *mctrk2 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track2));

        if (mctrk2->geantId() == 11){

          double mass1 = Mpion;
          double mass2 = Mkaon;

          TLorentzVector mTrk1Mom = mctrk1->fourMomentum();
          TLorentzVector mTrk2Mom = mctrk2->fourMomentum();

          double invmass = TMath::Sqrt(pow(mass1,2) + pow(mass2,2) + 2*(mTrk1Mom*mTrk2Mom));

          hD0MCMass->Fill(invmass);

          // cout << "InvMass = " << invmass << endl;
        }

      }

    }

  }

  double px = 0.0;
  double py = 0.0;

  for(unsigned short iTracks = 0; iTracks < nmctracks; iTracks++){
    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(iTracks));

    int idvx = mctrk->idVtxStart();

    if (idvx != 1) continue;

    px += mctrk->p().X();
    py += mctrk->p().Y();

    // cout << px << "\t" << py << endl;

  }

  // cout << "MC pT = " << TMath::Sqrt(pow(px, 2) + pow(py, 2)) << endl;

  double pxReco = 0.0;
  double pyReco = 0.0;

  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));

    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = trk->pMom();
    } else {
      // get global track vector
      mTrkMom = trk->gMom(mVertex, Bfield);
    }

    pxReco += mTrkMom.x();
    pyReco += mTrkMom.y();

  }

  // cout << "Reco pT = " << TMath::Sqrt(pow(pxReco, 2) + pow(pyReco, 2)) << endl;

    
  for(unsigned short iTracks = 0; iTracks < nmctracks; iTracks++){
    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(iTracks));

    int idvx = mctrk->idVtxStart();
    // if (idvx <= 1) continue;
    // if (idvx >= nmcvertices) continue;
    // if (mctrk->isFromShower()) continue;

    int mctrkid = mctrk->id();

    // if (mctrkid > 200) cout << "MC Track ID = " << mctrkid << endl;

    // Track acceptance cut                                                                                                                                                                                                                                    
    if (mctrk->pt() < 0.2 || abs(mctrk->eta()) > 1) continue;


    // ============== This is the D0 Fast Sim test ====                                                                                                                                                           
    if (IsD0MCVertex(tracksgeantidsfromvertex, idvx)) {
    	for(unsigned short jTracks = iTracks+1; jTracks < nmctracks; jTracks++){
    	    StPicoMcTrack *mctrk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(jTracks));
    	    int idvx1 = mctrk1->idVtxStart();
    	    if (mctrk1->pt() < 0.2 || abs(mctrk1->eta()) > 1) continue;
    	    if (IsD0MCVertex(tracksgeantidsfromvertex, idvx1)) {
    		if(idvx!=idvx1)continue;//Same D0 decay
    		// Note this calc assumes D0 is forced into Kpi only and that PID is perfect
    		StPicoTrack *d1 = SmearMom(mctrk->p(),mctrk->geantId());
    		StPicoTrack *d2 = SmearMom(mctrk1->p(),mctrk1->geantId());
    		TVector3 gmom1 = d1->gMom();
    		TVector3 gmom2 = d2->gMom();
    		TLorentzVector *pp1 = new TLorentzVector();
    		TLorentzVector *pp2 = new TLorentzVector();
    		if(mctrk->geantId() == 8 || mctrk->geantId() == 9){
    		    pp1->SetXYZM(gmom1.Px(),gmom1.Py(),gmom1.Pz(),0.1396);
    		    pp2->SetXYZM(gmom2.Px(),gmom2.Py(),gmom2.Pz(),0.4937);
    		}
    		else{
    		    pp1->SetXYZM(gmom1.Px(),gmom1.Py(),gmom1.Pz(),0.4937);
            pp2->SetXYZM(gmom2.Px(),gmom2.Py(),gmom2.Pz(),0.1396);
    		}
    		TLorentzVector *pD0 = new TLorentzVector(*pp1+*pp2);
    		hD0Mass->Fill(pD0->M());
    		delete pp1; delete pp2;delete pD0;
    	    }
    	}
    }
   // ============== End D0 Fast Sim test ============     

    ////////////////////////// All Kaon Pion Histograms //////////////////////////////////////////
  // if (!IsD0MCVertex(tracksgeantidsfromvertex, idvx)) {
  	if (mctrk->geantId() == 8 || mctrk->geantId() == 9){   
	    hMCPionPt->Fill(mctrk->pt());
	    hMCPionEta->Fill(mctrk->eta());
	    hMCPionPhi->Fill(mctrk->p().Phi());

      if (IsD0MCTrack(iTracks)){
        // if (!IsD0MCVertex(tracksgeantidsfromvertex, idvx)) continue;
        hMCD0PionPt->Fill(mctrk->pt());
        hMCD0PionEta->Fill(mctrk->eta());
        hMCD0PionPhi->Fill(mctrk->p().Phi());
      }
  	}

  	if (mctrk->geantId() == 11 || mctrk->geantId() == 12){
	    hMCKaonPt->Fill(mctrk->pt());
	    hMCKaonEta->Fill(mctrk->eta());
	    hMCKaonPhi->Fill(mctrk->p().Phi());

      if (IsD0MCTrack(iTracks)){
        // if (!IsD0MCVertex(tracksgeantidsfromvertex, idvx)) continue;
        hMCD0KaonPt->Fill(mctrk->pt());
        hMCD0KaonEta->Fill(mctrk->eta());
        hMCD0KaonPhi->Fill(mctrk->p().Phi());
      }
    }
  
  // }
    ///////////////////////// D0 Kaon Pion Histograms ///////////////////////////////////////////////

    // if (IsD0MCVertex(tracksgeantidsfromvertex, idvx)) {
    //   if (mctrk->geantId() == 8 || mctrk->geantId() == 9){
    //     hMCD0PionPt->Fill(mctrk->pt());
    //     hMCD0PionEta->Fill(mctrk->eta());
    //     hMCD0PionPhi->Fill(mctrk->p().Phi());
    //   }

    //   if (mctrk->geantId() == 11 || mctrk->geantId() == 12){
    //     hMCD0KaonPt->Fill(mctrk->pt());
    //     hMCD0KaonEta->Fill(mctrk->eta());
    //     hMCD0KaonPhi->Fill(mctrk->p().Phi());
    //   }
    // }

  }

  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // primary track switch: get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
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

    if (pt < 0.2 || abs(eta) > 1.0) continue;
    
    int mctrkid = trk->idTruth();

    if (mctrkid > mPicoDst->numberOfMcTracks()) continue;

    cout << "MC Track ID = " << mctrkid << endl;

    hRecoTrack->Fill(pt);

    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(mctrkid - 1));
    if(!mctrk) continue;

    hMatchedRecoTrack->Fill(pt);

    //if (mctrk->isFromShower()) continue;
    
    cout << "MC Track ID = " << mctrkid << "\t" << mctrk->id() << endl;

    int iVtx = mctrk->idVtxStart();

    //if (iVtx >= nmcvertices) continue;


    ////////////////////////// All Kaon Pion Histograms //////////////////////////////////////////
    // if(!IsD0MCVertex(tracksgeantidsfromvertex, iVtx)){
    	if (mctrk->geantId() == 8 || mctrk->geantId() == 9){
    	    hRecoPionPt->Fill(pt);//pt);
    	    hRecoPionEta->Fill(eta);
    	    hRecoPionPhi->Fill(phi);  
    	    
    	    hRecoVMCPionPt->Fill(pt, mctrk->pt());
    	    hRecoVMCPionPhi->Fill(phi, mctrk->p().Phi());
    	    hRecoVMCPionEta->Fill(eta, mctrk->eta());
    	}

    	if (mctrk->geantId() == 11 || mctrk->geantId() == 12){
    	    
    	    hRecoKaonPt->Fill(pt);
    	    hRecoKaonEta->Fill(eta);
    	    hRecoKaonPhi->Fill(phi);
    	    
    	    hRecoVMCKaonPt->Fill(pt, mctrk->pt());
    	    hRecoVMCKaonPhi->Fill(phi, mctrk->p().Phi());
    	    hRecoVMCKaonEta->Fill(eta, mctrk->eta());
    	}
    // }

    ///////////////////////// D0 Kaon Pion Histograms ///////////////////////////////////////////////

    if (IsD0MCVertex(tracksgeantidsfromvertex, iVtx)) {
      if (mctrk->geantId() == 8 || mctrk->geantId() == 9){
        hRecoD0PionPt->Fill(pt);
        hRecoD0PionEta->Fill(eta);
        hRecoD0PionPhi->Fill(phi);

        hRecoVMCD0PionPt->Fill(pt, mctrk->pt());
        hRecoVMCD0PionEta->Fill(eta, mctrk->eta());
        hRecoVMCD0PionPhi->Fill(phi, mctrk->p().Phi());
      }

      if (mctrk->geantId() == 11 || mctrk->geantId() == 12){
        hRecoD0KaonPt->Fill(pt);
        hRecoD0KaonEta->Fill(eta);
        hRecoD0KaonPhi->Fill(phi);

        hRecoVMCD0KaonPt->Fill(pt, mctrk->pt());
        hRecoVMCD0KaonEta->Fill(eta, mctrk->eta());
        hRecoVMCD0KaonPhi->Fill(phi, mctrk->p().Phi());
      }
    }


    ////////////////////////// All Kaon Pion (TPC) PID Histograms //////////////////////////////////////////

    double zpi = trk->nSigmaPion();
    double zka = trk->nSigmaKaon();
    double zpr = trk->nSigmaProton();

    bool k_tpc_pion = IsTpcPion(trk);
    bool k_tpc_kaon = IsTpcKaon(trk);

    bool k_pion = IsPion(trk);
    bool k_kaon = IsKaon(trk);

    // k_tpc_pion = (abs(zpi) <= abs(zka)) ? kTRUE : kFALSE;
    // k_tpc_kaon = (abs(zpi) > abs(zka)) ? kTRUE : kFALSE;

    // if ((abs(zpi) > 2) && (abs(zka) > 2)){
    //   k_tpc_pion = kFALSE;
    //   k_tpc_kaon = kFALSE;
    // }

    // if ((abs(zpi) <= 2)) k_tpc_pion = kTRUE;
    // if ((abs(zka) <= 2)) k_tpc_kaon = kTRUE;

    // if (k_tpc_pion) k_tpc_kaon = kFALSE;
    // if (k_tpc_kaon) k_tpc_pion = kFALSE;

    if (k_pion){
      hRecoTPCPionPt->Fill(pt);
      hRecoTPCPionEta->Fill(eta);
      hRecoTPCPionPhi->Fill(phi);

      hRecoVMCTPCPionPt->Fill(pt, mctrk->pt());
      hRecoVMCTPCPionEta->Fill(eta, mctrk->eta());
      hRecoVMCTPCPionPhi->Fill(phi, mctrk->p().Phi());
    }

    if (k_kaon){
      hRecoTPCKaonPt->Fill(pt);
      hRecoTPCKaonEta->Fill(eta);
      hRecoTPCKaonPhi->Fill(phi);

      hRecoVMCTPCKaonPt->Fill(pt, mctrk->pt());
      hRecoVMCTPCKaonEta->Fill(eta, mctrk->eta());
      hRecoVMCTPCKaonPhi->Fill(phi, mctrk->p().Phi());
    }

    ////////////////////////// All Kaon Pion (TPC + TOF) PID Histograms //////////////////////////////////////////

    // int pid = abs(IsWhatParticle(trk));

    if (k_pion){
      hRecoTOFPionPt->Fill(pt);
      hRecoTOFPionEta->Fill(eta);
      hRecoTOFPionPhi->Fill(phi);

      hRecoVMCTOFPionPt->Fill(pt, mctrk->pt());
      hRecoVMCTOFPionEta->Fill(eta, mctrk->eta());
      hRecoVMCTOFPionPhi->Fill(phi, mctrk->p().Phi());
    }

    if (k_kaon){
      hRecoTOFKaonPt->Fill(pt);
      hRecoTOFKaonEta->Fill(eta);
      hRecoTOFKaonPhi->Fill(phi);

      hRecoVMCTOFKaonPt->Fill(pt, mctrk->pt());
      hRecoVMCTOFKaonEta->Fill(eta, mctrk->eta());
      hRecoVMCTOFKaonPhi->Fill(phi, mctrk->p().Phi());
    }

    ////////////////////////// All Kaon Pion (TPC + HFT) PID Histograms //////////////////////////////////////////

    if (trk->hasPxl1Hit() || trk->hasPxl2Hit() || trk->hasIstHit()){
      if (k_pion){
        hRecoHFTPionPt->Fill(pt);
        hRecoHFTPionEta->Fill(eta);
        hRecoHFTPionPhi->Fill(phi);

        hRecoVMCHFTPionPt->Fill(pt, mctrk->pt());
        hRecoVMCHFTPionEta->Fill(eta, mctrk->eta());
        hRecoVMCHFTPionPhi->Fill(phi, mctrk->p().Phi());
      }

      if (k_kaon){
        hRecoHFTKaonPt->Fill(pt);
        hRecoHFTKaonEta->Fill(eta);
        hRecoHFTKaonPhi->Fill(phi);

        hRecoVMCHFTKaonPt->Fill(pt, mctrk->pt());
        hRecoVMCHFTKaonEta->Fill(eta, mctrk->eta());
        hRecoVMCHFTKaonPhi->Fill(phi, mctrk->p().Phi());
      }
    }
  }

}  // track function


// Bool_t StTestMC::IsD0KaonPion()
// {
//   if ()
// }


Bool_t StTestMC::IsD0MCTrack(int trackid){

  bool d0Track = kFALSE;

  for (unsigned short track1 = 0; track1 < mPicoDst->numberOfMcTracks(); track1++){
    StPicoMcTrack *mctrk1 = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(track1));

    if (mctrk1->geantId() == 37 || mctrk1->geantId() == 38) {

      int idvxstop = mctrk1->idVtxStop() -  1;

      // cout << "ID Vertex Stop : " << idvxstop << endl;

      if (idvxstop < 0) return d0Track; // Need this check

      StPicoMcVertex *mcvx1 = static_cast<StPicoMcVertex*>(mPicoDst->mcVertex(idvxstop));

      int stopvertexid = mcvx1->id();

      // cout << "ID Vertex Stop : " << stopvertexid << endl;


      StPicoMcTrack *D0Track = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(trackid));

      int idvxstart = D0Track->idVtxStart() -  1;

      StPicoMcVertex *d0vx = static_cast<StPicoMcVertex*>(mPicoDst->mcVertex(idvxstart));

      int startvertexidfortrack = d0vx->id();

      // cout << "ID Vertex Start : " << startvertexidfortrack << endl;

      if (stopvertexid == startvertexidfortrack) d0Track = kTRUE;
    
    }

  }

  return d0Track; 

}

Bool_t StTestMC::IsD0MCVertex(vector<int> arr[], int idvx)
{

  bool ismcvertex = kFALSE;

  // int arrSize = sizeof(arr)/sizeof(arr[0]);

  // cout << "Array Size = " << sizeof(arr)/sizeof(arr[0]) << endl;
  // if (idvx >= arrSize) return ismcvertex;

  // if (arr[idvx].size() == 2 && ((std::find(arr[idvx].begin(), arr[idvx].end(), 8) && std::find(arr[idvx].begin(), arr[idvx].end(), 12)) || (std::find(arr[idvx].begin(), arr[idvx].end(), 9) && std::find(arr[idvx].begin(), arr[idvx].end(), 11)))) ismcvertex = kTRUE;
  // else ismcvertex = kFALSE;

  if (arr[idvx].size() == 2){
    
    if (arr[idvx][0] == 8 && arr[idvx][1] == 12) ismcvertex = kTRUE;
    else if (arr[idvx][0] == 12 && arr[idvx][1] == 8) ismcvertex = kTRUE;
    else if (arr[idvx][0] == 9 && arr[idvx][1] == 11) ismcvertex = kTRUE;
    else if (arr[idvx][0] == 11 && arr[idvx][1] == 9) ismcvertex = kTRUE;
    else ismcvertex = kFALSE;
    // if (ismcvertex) cout << idvx << "\t" << arr[idvx][0] << "\t" << arr[idvx][1] << endl;
  }

  return ismcvertex;
}


bool StTestMC::IsTpcPion(StPicoTrack *trk){
  double zpi = trk->nSigmaPion();
  if (abs(zpi) < 3.) return true;
  return false;
}

bool StTestMC::IsTpcKaon(StPicoTrack *trk){
  double zka = trk->nSigmaKaon();
  if (abs(zka) < 2.) return true;
  return false;
}

bool StTestMC::IsPion(StPicoTrack *trk){
  if (!IsTpcPion(trk)) return false;
  float beta = GetTofBeta(trk);
  if (isnan(beta) || beta < 0) return true;
  double p = trk->gMom(mVertex, Bfield).Mag();
  float oneOverBetaExpected = sqrt(M_PION_PLUS*M_PION_PLUS / p / p + 1);
  double nsigma = abs(1./beta - oneOverBetaExpected);
  if (nsigma < 0.03) return true;
  return false;
}

bool StTestMC::IsKaon(StPicoTrack *trk){
  if (!IsTpcKaon(trk)) return false;
  float beta = GetTofBeta(trk);
  if (isnan(beta) || beta < 0) return true;
  double p = trk->gMom(mVertex, Bfield).Mag();
  float oneOverBetaExpected = sqrt(M_KAON_PLUS*M_KAON_PLUS / p / p + 1);
  double nsigma = abs(1./beta - oneOverBetaExpected);
  if (nsigma < 0.03) return true;
  return false;
}

float StTestMC::GetTofBeta(StPicoTrack *trk){
  int index2tof = trk->bTofPidTraitsIndex();
    
  float beta = std::numeric_limits<double>::quiet_NaN();

  if (index2tof < 0) return beta;

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
    // if (abs(tofPid->btofBeta() - beta) > 1e-4) cout << "Compare Betas : " << tofPid->btofBeta() << "\t" << beta << endl;
  }

  return beta;

}

void StTestMC::IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e){ // NEW PID APPROACH
  if(!AcceptTrack(trk, Bfield, mVertex)){pid = 0; m = 0.; e = 0.; return;}

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
  // double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
  // double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
  // double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

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


Int_t StTestMC::IsWhatParticle(StPicoTrack *trk){ // Just to get the PID out
  int pid;
  double m; 
  double e;
  IsWhatParticle(trk, pid, m, e);
  return pid;
}


//
//
//________________________________________________________________________
void StTestMC::RunTowers()
{
  // constants: assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // looping over clusters - STAR: matching already done, get # of clusters and set variables
  unsigned int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits();

  // loop over ALL clusters in PicoDst and add to jet //TODO
  for(unsigned short iClus = 0; iClus < nBEmcPidTraits; iClus++){
    StPicoBEmcPidTraits *cluster = static_cast<StPicoBEmcPidTraits*>(mPicoDst->bemcPidTraits(iClus));
    if(!cluster){ cout<<"Cluster pointer does not exist.. iClus = "<<iClus<<endl; continue; }

    // cluster and tower ID
    int clusID = cluster->bemcId();  // index in bemc point array
    int towID = cluster->btowId();   // projected tower Id: 1 - 4800
    int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    if(towID < 0) continue;

    // cluster and tower position - from vertex and ID
    TVector3 towPosition = mEmcPosition->getPosFromVertex(mVertex, towID);
    double towPhi = towPosition.Phi();
    double towEta = towPosition.PseudoRapidity();

    // matched track index
    int trackIndex = cluster->trackIndex();
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackIndex));
    if(!trk) { cout<<"No trk pointer...."<<endl; continue; }
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

  } // BEmc loop

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    // get tower pointer
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // tower ID: tower ID shifted by +1 from array index
    int towerID = itow + 1;
    if(towerID < 0) continue; // double check these aren't still in the event list

    // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    double towerPhi = towerPosition.Phi();
    double towerEta = towerPosition.PseudoRapidity();
    int towerADC = tower->adc();
    double towerE = tower->energy();
    double towerEt = towerE / (1.0*TMath::CosH(towerEta)); // THIS should be USED

    // do stuff with towers and fill histograms here
    // ........
    hTrackEt->Fill(towerEta, towerPhi, towerEt);


  } // tower loop

}  // run towers function

//
//
// __________________________________________________________________________________
void StTestMC::SetSumw2() {
  hCentrality->Sumw2();
  hMultiplicity->Sumw2();
  hJetPt->Sumw2();
  hJetCorrPt->Sumw2();
}

//
// Function: get relative phi of jet and track (-0.5pi, 1.5pi)
//________________________________________________________________________
Double_t StTestMC::RelativePhi(Double_t mphi,Double_t vphi) const
{ // function to calculate relative PHI
  double dphi = mphi-vphi;

  // set dphi to operate on adjusted scale
  if(dphi < -0.5*pi) dphi += 2.0*TMath::Pi();
  if(dphi >  1.5*pi) dphi -= 2.0*TMath::Pi();

  // test check
  if( dphi < -0.5*pi || dphi > 1.5*pi )
    Form("%s: dPHI not in range [-0.5*Pi, 1.5*Pi]!", GetName());

  return dphi; // dphi in [-0.5Pi, 1.5Pi]                                                                                   
}

//
// Function: calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
//_________________________________________________________________________
Double_t StTestMC::RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{
  Double_t pi = 1.0*TMath::Pi();
  Double_t dphi = 1.0*TMath::Abs(EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi < -pi ){
    dphi = dphi + pi;
  } // this assumes we are doing full jets currently 
 
  if(dphi > 1.5*pi) dphi -= 2.0*pi;
  if((dphi > 1.0*pi) && (dphi < 1.5*pi)) dphi -= 1.0*pi;
  if((dphi > 0.5*pi) && (dphi < 1.0*pi)) dphi -= 1.0*pi;
  dphi = 1.0*TMath::Abs(dphi);

  // test check
  if( dphi < 0.0 || dphi > 0.5*pi ) {
    //Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName());
    cout<<"dPhi not in range [0, 0.5*Pi]!"<<endl;
  }

  return dphi;   // dphi in [0, Pi/2]
}

//
//
//_________________________________________________________________________
void StTestMC::FillEmcTriggers() {
  // number of Emcal Triggers
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // fill for valid triggers
    if(emcTrig->isHT0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(emcTrig->isHT1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(emcTrig->isHT2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(emcTrig->isHT3()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(emcTrig->isJP0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(emcTrig->isJP1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(emcTrig->isJP2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
}
//
// elems: sizeof(myarr)/sizeof(*myarr) prior to passing to function
// upon passing the array collapses to a pointer and can not get size anymore
//________________________________________________________________________
Bool_t StTestMC::DoComparison(int myarr[], int elems) {
  //std::cout << "Length of array = " << (sizeof(myarr)/sizeof(*myarr)) << std::endl;
  bool match = kFALSE;

  // loop over specific physics selection array and compare to specific event trigger
  for(int i=0; i<elems; i++) {
    if(mPicoEvent->isTrigger(myarr[i])) match = kTRUE;
    if(match) break;
  }

  return match;
}
StPicoTrack * StTestMC::SmearMom(TVector3 p, int pid)
{
    float pt = p.Perp();
    float pt1 = pt;
    if(pt1>10) pt1 = 10;//Used for high pt-hat bin smearing test
    float sPt = -1;
    
    if(pid ==8 || pid == 9)sPt = gRandom->Gaus(pt, pt * fPionMomResolution->Eval(pt1));// Pion
    else if(pid ==11 || pid == 12)sPt = gRandom->Gaus(pt, pt * fKaonMomResolution->Eval(pt1));// Kaon
    else if(pid ==15 || pid == 14)sPt = gRandom->Gaus(pt, pt * fProtonMomResolution->Eval(pt1));// Proton
    else sPt = gRandom->Gaus(pt, pt * fPionMomResolution->Eval(pt1));// Catch all: pions
    
    StPicoTrack *sMom = new StPicoTrack();
    sMom->setGlobalMomentum(sPt * cos(p.Phi()), sPt * sin(p.Phi()), sPt * sinh(p.PseudoRapidity()));
    return sMom;
}
