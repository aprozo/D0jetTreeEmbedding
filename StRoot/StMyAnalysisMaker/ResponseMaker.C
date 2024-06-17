#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TChain.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>
#ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
// #include "StJetTreeStruct.h"
#include <vector>

#pragma link C++ class vector<int> +;

using namespace std;

// #pragma link C++ class StJetTreeStruct+;

// #pragma link C++ class vector<float> +;
// #pragma link C++ class vector<vector<float> >+;
// #pragma link C++ class vector<int> +;
// #pragma link C++ class vector<vector<int> >+;
#endif



const double R = 0.4;
const double deltar = 0.05;
const int numberofbins = R/deltar;


// Important functions that are used repeatedly in the calculations

// function to calculate relative phi between 2 objects and shift between 0 and 2pi 
//___________________________________________________________________________________________
Double_t dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); //TODO absolute values
  if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

  if (deltaPhi > TMath::Pi()) deltaPhi= 2*(TMath::Pi()) - deltaPhi;
  return deltaPhi;   // dphi in [0, 2Pi]
}

Double_t standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dEta(Double_t eta1, Double_t eta2) {
  Double_t deltaEta;
  deltaEta = eta1 - eta2;

  return deltaEta;
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dR(Double_t delphi, Double_t deleta) {
  Double_t dRad;
  dRad = TMath::Sqrt(pow(delphi,2) + pow(deleta,2));

  return dRad;
}

// function to calculate invariant mass of a pair of objects
//___________________________________________________________________________________________
Double_t Mass(Double_t m1, Double_t m2, Double_t E1, Double_t E2, Double_t px1, Double_t px2, Double_t py1, Double_t py2, Double_t pz1, Double_t pz2) {
  Double_t m;
  m = TMath::Sqrt(pow(m1,2) + pow(m2,2) + 2*(E1*E2 - px1*px2 - py1*py2 - pz1*pz2));

  return m;
}

// function to calculate invariant mass of a pair of objects
//___________________________________________________________________________________________
Double_t pTforD0(Double_t px1, Double_t px2, Double_t py1, Double_t py2) {
  Double_t m;
  m = TMath::Sqrt(pow(px1 + px2, 2) + pow(py1 + py2, 2));

  return m;
}

Double_t pX(Double_t pT, Double_t phi){
  Double_t px;
  px = pT*TMath::Cos(phi);

  return px;
}

Double_t pY(Double_t pT, Double_t phi){
  Double_t py;
  py = pT*TMath::Sin(phi);

  return py;
}

Double_t p(Double_t E, Double_t m){
   
  Double_t pmag;
  pmag = TMath::Sqrt(pow(E,2) - pow(m,2));

  return pmag;
}

Double_t p(Double_t px, Double_t py, Double_t pz){
   
  Double_t pmag;
  pmag = TMath::Sqrt(pow(px,2) + pow(py,2) + pow(pz, 2));

  return pmag;
}

Double_t pZ(Double_t E, Double_t m, Double_t eta){
   
  Double_t pz;
  pz = p(E,m)*(TMath::Exp(2*eta)-1)/(TMath::Exp(2*eta)+1);

  return pz;
}

int GetD0PtBin(Double_t pt){
  if (pt > 0 && pt < 1.0) return 0;
  if (pt >= 1.0 && pt < 3.0) return 1;
  if (pt >= 3.0 && pt < 5.0) return 2;
  if (pt >= 5.0 && pt < 7.0) return 3;
  if (pt >= 7.0 && pt < 10.0) return 4;
  if (pt >= 10.0 && pt < 15.0) return 5;
  if (pt >= 15.0 && pt < 20.0) return 6;
  if (pt >= 20.0 && pt < 30.0) return 7;
  if (pt >= 30.0) return 8;
}

int GetJetPtBin(Double_t pt){
  int jetbin = -99;

  if (pt < 0.0) jetbin = 0;
  else if (pt >= 0.0 && pt < 5.0) jetbin = 1;
  else if (pt >= 5.0 && pt < 10.0) jetbin = 2;
  else if (pt >= 10.0 && pt < 15.0) jetbin = 3;
  else if (pt >= 15.0 && pt < 20.0) jetbin = 4;
  else if (pt >= 20.0 && pt < 30.0) jetbin = 5;
  else if (pt >= 30.0 && pt < 50.0) jetbin = 6;
  else jetbin = 7;

  return jetbin;
}

// void JetAnalysis(){}

void ResponseMaker(TFile *outfile = 0x0){

  TStopwatch timer;

  timer.Start();

  // Declaring Global variables

  const int maximum = 5000;

  Int_t           mcRunID;
  Int_t           mcEventId;
  Float_t         mcRefMult;
  Float_t         mcCentrality;
  vector<int> *mcTriggers = new vector<int>;
  vector<double>  *mcPrimaryVertex = new vector<double>;
  vector<double>  *mcPrimaryVertexErr = new vector<double>;
  Float_t         mcJetPt;
  Float_t         mcJetCorrPt;
  Float_t         mcJetEta;
  Float_t         mcJetPhi;
  Float_t         mcJetArea;
  Float_t         mcJetRadius;
  Float_t         mcJetE;
  Float_t         mcJetNEF;
  Float_t         mcJetRhoVal;
  Float_t         mcJetHighestTrackPt;
  Int_t           mcJetNConst;
  Float_t         mcD0Mass;
  Float_t         mcPionPt;
  Float_t         mcPionEta;
  Float_t         mcPionPhi;
  Float_t         mcPionCharge;
  Float_t         mcKaonPt;
  Float_t         mcKaonEta;
  Float_t         mcKaonPhi;
  Float_t         mcKaonCharge;
  Float_t         mcTrackID[maximum];   //[numberofconstituents]
  Float_t         mcTrackPt[maximum];   //[numberofconstituents]
  Float_t         mcTrackEta[maximum];   //[numberofconstituents]
  Float_t         mcTrackPhi[maximum];   //[numberofconstituents]
  Float_t         mcTrackPx[maximum];   //[numberofconstituents]
  Float_t         mcTrackPy[maximum];   //[numberofconstituents]
  Float_t         mcTrackPz[maximum];   //[numberofconstituents]
  Float_t         mcTrackCharge[maximum];   //[numberofconstituents]

  Int_t           recoRunID;
  Int_t           recoEventId;
  Float_t         recoRefMult;
  Float_t         recoCentrality;
  vector<int> *recoTriggers = new vector<int>;
  vector<double>  *recoPrimaryVertex = new vector<double>;
  vector<double>  *recoPrimaryVertexErr = new vector<double>;
  Float_t         recoJetPt;
  Float_t         recoJetCorrPt;
  Float_t         recoJetEta;
  Float_t         recoJetPhi;
  Float_t         recoJetArea;
  Float_t         recoJetRadius;
  Float_t         recoJetE;
  Float_t         recoJetNEF;
  Float_t         recoJetRhoVal;
  Float_t         recoJetHighestTrackPt;
  Int_t           recoJetNConst;
  Float_t         recoD0Mass;
  Float_t         recoPionPt;
  Float_t         recoPionEta;
  Float_t         recoPionPhi;
  Float_t         recoPionCharge;
  Float_t         recoKaonPt;
  Float_t         recoKaonEta;
  Float_t         recoKaonPhi;
  Float_t         recoKaonCharge;
  Float_t         recoTrackID[maximum];   //[numberofconstituents]
  Float_t         recoTrackPt[maximum];   //[numberofconstituents]
  Float_t         recoTrackEta[maximum];   //[numberofconstituents]
  Float_t         recoTrackPhi[maximum];   //[numberofconstituents]
  Float_t         recoTrackPx[maximum];   //[numberofconstituents]
  Float_t         recoTrackPy[maximum];   //[numberofconstituents]
  Float_t         recoTrackPz[maximum];   //[numberofconstituents]
  Float_t         recoTrackCharge[maximum];   //[numberofconstituents]


  TFile *filein;

  filein = TFile::Open("/Users/diptanilroy/Desktop/Y2020/STAR/STARAnalysis/2021/EMBEDDING/RAGHAVSFRAMEWORK/ResponsePlots/SimJetTree.root")
  
  TString mcinputstring = "SimJetSaver/MCJets";
  TString recoinputstring = "SimJetSaver/RecoJets";

  // TFile *f = TFile::Open("test.root");

  TTree *mcjettree = (TTree*)filein->Get(mcinputstring.Data());
  TTree *recojettree = (TTree*)filein->Get(recoinputstring.Data());

  // // jettree->SetBranchStatus("Triggers", 0);
  // jettree->SetBranchStatus("PrimaryVertex", 0);
  // jettree->SetBranchStatus("PrimaryVertexErr", 0);

  mcjettree->SetBranchAddress("RunID", &mcRunID);
  mcjettree->SetBranchAddress("EventId", &mcEventId);
  mcjettree->SetBranchAddress("RefMult", &mcRefMult);
  mcjettree->SetBranchAddress("Centrality", &mcCentrality);
  mcjettree->SetBranchAddress("Triggers", &mcTriggers);
  mcjettree->SetBranchAddress("PrimaryVertex", &mcPrimaryVertex);
  mcjettree->SetBranchAddress("PrimaryVertexErr", &mcPrimaryVertexErr);
  mcjettree->SetBranchAddress("JetPt", &mcJetPt);
  mcjettree->SetBranchAddress("JetCorrPt", &mcJetCorrPt);
  mcjettree->SetBranchAddress("JetEta", &mcJetEta);
  mcjettree->SetBranchAddress("JetPhi", &mcJetPhi);
  mcjettree->SetBranchAddress("JetArea", &mcJetArea);
  mcjettree->SetBranchAddress("JetRadius", &mcJetRadius);
  mcjettree->SetBranchAddress("JetE", &mcJetE);
  mcjettree->SetBranchAddress("JetNEF", &mcJetNEF);
  mcjettree->SetBranchAddress("JetRhoVal", &mcJetRhoVal);
  mcjettree->SetBranchAddress("JetHighestTrackPt", &mcJetHighestTrackPt);
  mcjettree->SetBranchAddress("JetNConst", &mcJetNConst);
  mcjettree->SetBranchAddress("D0Mass", &mcD0Mass);
  mcjettree->SetBranchAddress("PionPt", &mcPionPt);
  mcjettree->SetBranchAddress("PionEta", &mcPionEta);
  mcjettree->SetBranchAddress("PionPhi", &mcPionPhi);
  mcjettree->SetBranchAddress("PionCharge", &mcPionCharge);
  mcjettree->SetBranchAddress("KaonPt", &mcKaonPt);
  mcjettree->SetBranchAddress("KaonEta", &mcKaonEta);
  mcjettree->SetBranchAddress("KaonPhi", &mcKaonPhi);
  mcjettree->SetBranchAddress("KaonCharge", &mcKaonCharge);
  mcjettree->SetBranchAddress("TrackID", mcTrackID);
  mcjettree->SetBranchAddress("TrackPt", mcTrackPt);
  mcjettree->SetBranchAddress("TrackEta", mcTrackEta);
  mcjettree->SetBranchAddress("TrackPhi", mcTrackPhi);
  mcjettree->SetBranchAddress("TrackPx", mcTrackPx);
  mcjettree->SetBranchAddress("TrackPy", mcTrackPy);
  mcjettree->SetBranchAddress("TrackPz", mcTrackPz);
  mcjettree->SetBranchAddress("TrackCharge", mcTrackCharge);

  recojettree->SetBranchAddress("RunID", &recoRunID);
  recojettree->SetBranchAddress("EventId", &recoEventId);
  recojettree->SetBranchAddress("RefMult", &recoRefMult);
  recojettree->SetBranchAddress("Centrality", &recoCentrality);
  recojettree->SetBranchAddress("Triggers", &recoTriggers);
  recojettree->SetBranchAddress("PrimaryVertex", &recoPrimaryVertex);
  recojettree->SetBranchAddress("PrimaryVertexErr", &recoPrimaryVertexErr);
  recojettree->SetBranchAddress("JetPt", &recoJetPt);
  recojettree->SetBranchAddress("JetCorrPt", &recoJetCorrPt);
  recojettree->SetBranchAddress("JetEta", &recoJetEta);
  recojettree->SetBranchAddress("JetPhi", &recoJetPhi);
  recojettree->SetBranchAddress("JetArea", &recoJetArea);
  recojettree->SetBranchAddress("JetRadius", &recoJetRadius);
  recojettree->SetBranchAddress("JetE", &recoJetE);
  recojettree->SetBranchAddress("JetNEF", &recoJetNEF);
  recojettree->SetBranchAddress("JetRhoVal", &recoJetRhoVal);
  recojettree->SetBranchAddress("JetHighestTrackPt", &recoJetHighestTrackPt);
  recojettree->SetBranchAddress("JetNConst", &recoJetNConst);
  recojettree->SetBranchAddress("D0Mass", &recoD0Mass);
  recojettree->SetBranchAddress("PionPt", &recoPionPt);
  recojettree->SetBranchAddress("PionEta", &recoPionEta);
  recojettree->SetBranchAddress("PionPhi", &recoPionPhi);
  recojettree->SetBranchAddress("PionCharge", &recoPionCharge);
  recojettree->SetBranchAddress("KaonPt", &recoKaonPt);
  recojettree->SetBranchAddress("KaonEta", &recoKaonEta);
  recojettree->SetBranchAddress("KaonPhi", &recoKaonPhi);
  recojettree->SetBranchAddress("KaonCharge", &recoKaonCharge);
  recojettree->SetBranchAddress("TrackID", recoTrackID);
  recojettree->SetBranchAddress("TrackPt", recoTrackPt);
  recojettree->SetBranchAddress("TrackEta", recoTrackEta);
  recojettree->SetBranchAddress("TrackPhi", recoTrackPhi);
  recojettree->SetBranchAddress("TrackPx", recoTrackPx);
  recojettree->SetBranchAddress("TrackPy", recoTrackPy);
  recojettree->SetBranchAddress("TrackPz", recoTrackPz);
  recojettree->SetBranchAddress("TrackCharge", recoTrackCharge);


  ////////////////////////////////// Make Jet Shape Histograms /////////////////////////////////

  int mcentries = mcjettree->GetEntries();
  int recoentries = recojettree->GetEntries();

  cout << mcentries << "\t" << recoentries << endl;


  for (int event = 0; event < mcentries;  event++){
    jettree->GetEntry(event);

    if (mode != 1){ if (D0Mass < 1.80 || D0Mass > 1.92) continue; }
    
    int jetbin = -99;
    jetbin = GetJetPtBin(JetCorrPt);

    int d0bin = -99;

    int d0index = -99;

    for (int itrk = 0; itrk < JetNConst; itrk++){
      if (TrackID[itrk] == 99999) {

        d0bin = GetD0PtBin(TrackPt[itrk]);
        d0index = itrk;
        break;
      }
    }

    if (Centrality == 0) Centrality += 0.5;

    if (Centrality < fCentralityMin || Centrality >= fCentralityMax) continue;

    double d0efficiencyvalue = 1.0;

    if (Centrality == 0.5) aCentrality[d0bin][jetbin]->Fill(Centrality-0.5);
    else aCentrality[d0bin][jetbin]->Fill(Centrality);

    aMultiplicity[d0bin][jetbin]->Fill(RefMult);

    aSelectedD0JetPt[d0bin][jetbin]->Fill(JetPt, 1./d0efficiencyvalue);
    aSelectedD0CorrJetPt[d0bin][jetbin]->Fill(JetCorrPt, 1./d0efficiencyvalue);
    aSelectedD0JetPhi[d0bin][jetbin]->Fill(JetPhi, 1./d0efficiencyvalue);
    aSelectedD0JetEta[d0bin][jetbin]->Fill(JetEta, 1./d0efficiencyvalue);

    aSelectedPionPt[d0bin][jetbin]->Fill(PionPt, 1./d0efficiencyvalue);
    aSelectedPionEta[d0bin][jetbin]->Fill(PionEta, 1./d0efficiencyvalue);
    aSelectedPionPhi[d0bin][jetbin]->Fill(PionPhi, 1./d0efficiencyvalue);

    aSelectedKaonPt[d0bin][jetbin]->Fill(KaonPt, 1./d0efficiencyvalue);
    aSelectedKaonEta[d0bin][jetbin]->Fill(KaonEta, 1./d0efficiencyvalue);
    aSelectedKaonPhi[d0bin][jetbin]->Fill(KaonPhi, 1./d0efficiencyvalue);

    aSelectedD0Pt[d0bin][jetbin]->Fill(TrackPt[d0index], 1./d0efficiencyvalue);
    aSelectedD0Eta[d0bin][jetbin]->Fill(TrackEta[d0index], 1./d0efficiencyvalue);
    aSelectedD0Phi[d0bin][jetbin]->Fill(standardPhi(TrackPhi[d0index]), 1./d0efficiencyvalue);

    aSelectedD0Invmass[d0bin][jetbin]->Fill(D0Mass, 1./d0efficiencyvalue);
    // if (Centrality == 0.5) cout << Centrality << endl;

    aSelectedD0PtvSelectedD0JetPt[d0bin][jetbin]->Fill(TrackPt[d0index], JetCorrPt, 1./d0efficiencyvalue);

    double deltaD0phi = dPhi(JetPhi, standardPhi(TrackPhi[d0index]));
    double deltaD0eta = dEta(JetEta, TrackEta[d0index]);
    double deltaD0R = dR(deltaD0phi, deltaD0eta);

    double JetPx = JetCorrPt*TMath::Cos(JetPhi);
    double JetPy = JetCorrPt*TMath::Sin(JetPhi);

    double z = (JetPx*TrackPx[d0index] + JetPy*TrackPy[d0index])/pow(JetCorrPt, 2);

    // if (TrackPt[d0index] > 3) cout << z << endl;

    adphiD0Jet[d0bin][jetbin]->Fill(deltaD0phi, 1./d0efficiencyvalue);
    adetaD0Jet[d0bin][jetbin]->Fill(deltaD0eta, 1./d0efficiencyvalue);
    adphivdetaD0Jet[d0bin][jetbin]->Fill(deltaD0phi, deltaD0eta, 1./d0efficiencyvalue);
    adRD0Jet[d0bin][jetbin]->Fill(deltaD0R, 1./d0efficiencyvalue);
    aZD0Jet[d0bin][jetbin]->Fill(z, 1./d0efficiencyvalue);
    adRD0Jet_Zweighted[d0bin][jetbin]->Fill(deltaD0R, z/d0efficiencyvalue);

    aZvdR[d0bin][jetbin]->Fill(deltaD0R, z, 1./d0efficiencyvalue);

    for(int itrk = 0; itrk < JetNConst; itrk++) {

      if (TrackID[itrk] != 99999){
      
        double deltaphi = dPhi(JetPhi, standardPhi(TrackPhi[itrk]));
        double deltaeta = dEta(JetEta, TrackEta[itrk]);
        double deltar = dR(deltaphi, deltaeta);

        aJetDistD0[d0bin][jetbin]->Fill(deltar);
        aJetShapeD0[d0bin][jetbin]->Fill(deltar, TrackPt[itrk]/JetPt);
      }
    }

    aJetDistD0[d0bin][jetbin]->Fill(deltaD0R, 1./d0efficiencyvalue);
    aJetShapeD0[d0bin][jetbin]->Fill(deltaD0R, TrackPt[d0index]/(JetPt*d0efficiencyvalue));
  }


  TString outputstring;

  outputstring = "D0JetShape";

  if (mode == 0) outputstring += "SigBg_";
  if (mode == 1) outputstring += "Unlike_";
  if (mode == 2) outputstring += "Like_";

  outputstring += TString(Form("%i",fCentralityMin));
  outputstring += "-";
  outputstring += TString(Form("%i",fCentralityMax));

  if (!outfile) return; 

  outfile->mkdir(outputstring.Data());
  outfile->cd(outputstring.Data());


  for (int d0ptbin = 0; d0ptbin < 11; d0ptbin++){
    for (int jetptbin = 0; jetptbin < 8; jetptbin++){
      // The next histograms are before dR < 0.4 cuts
      aCentrality[d0ptbin][jetptbin]->Write();
      aMultiplicity[d0ptbin][jetptbin]->Write();
      adphiD0Jet[d0ptbin][jetptbin]->Write();
      adetaD0Jet[d0ptbin][jetptbin]->Write();
      adphivdetaD0Jet[d0ptbin][jetptbin]->Write();
      adRD0Jet[d0ptbin][jetptbin]->Write();
      aZD0Jet[d0ptbin][jetptbin]->Write();
      adRD0Jet_Zweighted[d0ptbin][jetptbin]->Write();
      aZvdR[d0ptbin][jetptbin]->Write();

      // The next histograms are after dR < 0.4 cuts
      aSelectedPionPt[d0ptbin][jetptbin]->Write();
      aSelectedKaonPt[d0ptbin][jetptbin]->Write();
      aSelectedD0Pt[d0ptbin][jetptbin]->Write();
      aSelectedPionPhi[d0ptbin][jetptbin]->Write();
      aSelectedKaonPhi[d0ptbin][jetptbin]->Write();
      aSelectedD0Phi[d0ptbin][jetptbin]->Write();
      aSelectedPionEta[d0ptbin][jetptbin]->Write();
      aSelectedKaonEta[d0ptbin][jetptbin]->Write();
      aSelectedD0Eta[d0ptbin][jetptbin]->Write();
      aSelectedD0Invmass[d0ptbin][jetptbin]->Write();

      aSelectedD0PtvSelectedD0JetPt[d0ptbin][jetptbin]->Write();

      aSelectedD0JetPt[d0ptbin][jetptbin]->Write();
      aSelectedD0CorrJetPt[d0ptbin][jetptbin]->Write();
      aSelectedD0JetPhi[d0ptbin][jetptbin]->Write();
      aSelectedD0JetEta[d0ptbin][jetptbin]->Write();
      aJetDistD0[d0ptbin][jetptbin]->Write();
      aJetShapeD0[d0ptbin][jetptbin]->Write();
    }
  }
  
  for (int d0ptbin = 0; d0ptbin < 11; d0ptbin++){
    for (int jetptbin = 0; jetptbin < 8; jetptbin++){
      if (aCentrality[d0ptbin][jetptbin]) delete aCentrality[d0ptbin][jetptbin];
      if (aMultiplicity[d0ptbin][jetptbin]) delete aMultiplicity[d0ptbin][jetptbin];
      if (adphiD0Jet[d0ptbin][jetptbin]) delete adphiD0Jet[d0ptbin][jetptbin];
      if (adetaD0Jet[d0ptbin][jetptbin]) delete adetaD0Jet[d0ptbin][jetptbin];
      if (adphivdetaD0Jet[d0ptbin][jetptbin]) delete adphivdetaD0Jet[d0ptbin][jetptbin];
      if (adRD0Jet[d0ptbin][jetptbin]) delete adRD0Jet[d0ptbin][jetptbin];
      if (aZD0Jet[d0ptbin][jetptbin]) delete aZD0Jet[d0ptbin][jetptbin];
      if (adRD0Jet_Zweighted[d0ptbin][jetptbin]) delete adRD0Jet_Zweighted[d0ptbin][jetptbin];
      if (aZvdR[d0ptbin][jetptbin]) delete aZvdR[d0ptbin][jetptbin];

      if (aSelectedPionPt[d0ptbin][jetptbin]) delete aSelectedPionPt[d0ptbin][jetptbin];
      if (aSelectedKaonPt[d0ptbin][jetptbin]) delete aSelectedKaonPt[d0ptbin][jetptbin];
      if (aSelectedD0Pt[d0ptbin][jetptbin]) delete aSelectedD0Pt[d0ptbin][jetptbin];
      if (aSelectedPionPhi[d0ptbin][jetptbin]) delete aSelectedPionPhi[d0ptbin][jetptbin];
      if (aSelectedKaonPhi[d0ptbin][jetptbin]) delete aSelectedKaonPhi[d0ptbin][jetptbin];
      if (aSelectedD0Phi[d0ptbin][jetptbin]) delete aSelectedD0Phi[d0ptbin][jetptbin];
      if (aSelectedPionEta[d0ptbin][jetptbin]) delete aSelectedPionEta[d0ptbin][jetptbin];
      if (aSelectedKaonEta[d0ptbin][jetptbin]) delete aSelectedKaonEta[d0ptbin][jetptbin];
      if (aSelectedD0Eta[d0ptbin][jetptbin]) delete aSelectedD0Eta[d0ptbin][jetptbin];
      if (aSelectedD0Invmass[d0ptbin][jetptbin]) delete aSelectedD0Invmass[d0ptbin][jetptbin];

      if (aSelectedD0PtvSelectedD0JetPt[d0ptbin][jetptbin]) delete aSelectedD0PtvSelectedD0JetPt[d0ptbin][jetptbin];

      if (aSelectedD0JetPt[d0ptbin][jetptbin]) delete aSelectedD0JetPt[d0ptbin][jetptbin];
      if (aSelectedD0CorrJetPt[d0ptbin][jetptbin]) delete aSelectedD0CorrJetPt[d0ptbin][jetptbin];
      if (aSelectedD0JetPhi[d0ptbin][jetptbin]) delete aSelectedD0JetPhi[d0ptbin][jetptbin];
      if (aSelectedD0JetEta[d0ptbin][jetptbin]) delete aSelectedD0JetEta[d0ptbin][jetptbin];
      if (aJetDistD0[d0ptbin][jetptbin]) delete aJetDistD0[d0ptbin][jetptbin];
      if (aJetShapeD0[d0ptbin][jetptbin]) delete aJetShapeD0[d0ptbin][jetptbin];

    }
  }

  jettree->Reset();
  filein->Close();

  timer.Stop();
  cout << "Real Time Used: " << timer.RealTime()/60 << "m" << endl;
}

void D0Jets_AuAu_2021(bool vertexdefault = kTRUE){
  const char *outfile;
  if (vertexdefault) outfile = "VertexDefault/AuAu_JetShape_Aug13.root";
  else outfile = "VertexVpdOrDefault/AuAu_JetShape_Aug13.root";

  TFile *histos = new TFile(outfile, "RECREATE");

  // Mode = 0 (Signal plus background)
  // Mode = 1 (Unlike background)
  // Mode = 2 (Like background)

  JetAnalysis(0, 00, 10, histos, vertexdefault);
  JetAnalysis(0, 10, 20, histos, vertexdefault);
  JetAnalysis(0, 20, 40, histos, vertexdefault);
  JetAnalysis(0, 40, 60, histos, vertexdefault);
  JetAnalysis(0, 60, 80, histos, vertexdefault);
  JetAnalysis(0, 00, 80, histos, vertexdefault);


  JetAnalysis(1, 00, 10, histos, vertexdefault);
  JetAnalysis(1, 10, 20, histos, vertexdefault);
  JetAnalysis(1, 20, 40, histos, vertexdefault);
  JetAnalysis(1, 40, 60, histos, vertexdefault);
  JetAnalysis(1, 60, 80, histos, vertexdefault);
  JetAnalysis(1, 00, 80, histos, vertexdefault);

  JetAnalysis(2, 00, 10, histos, vertexdefault);
  JetAnalysis(2, 10, 20, histos, vertexdefault);
  JetAnalysis(2, 20, 40, histos, vertexdefault);
  JetAnalysis(2, 40, 60, histos, vertexdefault);
  JetAnalysis(2, 60, 80, histos, vertexdefault);
  JetAnalysis(2, 00, 80, histos, vertexdefault);

  histos->Close();


}
