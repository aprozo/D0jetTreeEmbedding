#ifndef StHIJetSaverCombined_h
#define StHIJetSaverCombined_h

#include "StJetFrameworkPicoBase.h"

#include "StJetTreeStruct.h"
class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3;
class THnSparse;
class TProfile;
class TString;
class TVector3;
class TLorentzVector;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoMcTrack;
class StPicoMcVertex;
class StRefMultCorr;

// jet-framework classes
class StCentMaker;
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;

class StHIOverlay_Test;
class StHIRecoJetsCombined;

#ifndef __ObjectsToPass__H__
#define __ObjectsToPass__H__
struct ObjectsToPass
{
  ObjectsToPass()
  {
    lambda_1_1 = -1;
    lambda_1_1half = -1;
    lambda_1_2 = -1;
    lambda_1_3 = -1;
  };
  Double_t lambda_1_1;
  Double_t lambda_1_1half;
  Double_t lambda_1_2;
  Double_t lambda_1_3;
};
#endif //!__ObjectsToPass__H__

class StHIJetSaverCombined : public StJetFrameworkPicoBase
{
public:
  StHIJetSaverCombined(const char *name, StPicoDstMaker *picoMaker, const char *outName1, const char *outName2);
  virtual ~StHIJetSaverCombined();

  // class required functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *opt = "");
  virtual Int_t Finish();

  void DeclareHistograms();
  void DeclareTree();
  void BookTree();
  void WriteTree(TTree *sometree);

  virtual void SetPrintLevel(int i) { fPrintLevel = i; }

protected:
  void SaveJets(int iteration);
  void SaveMCJets(int iteration);
  void SaveRecoJets(int iteration);

  Double_t standardPhi(Double_t phi);

  const char *fD0TaggerName;

  // switches
  Bool_t doPrintEventCounter; // print event # switch
  Bool_t fDoEffCorr;          // efficiency correction to tracks

  // event selection types
  UInt_t fEmcTriggerEventType; // Physics selection of event used for signal
  UInt_t fMBEventType;         // Physics selection of event used for MB
  Int_t fEmcTriggerArr[8];     // EMCal triggers array: used to select signal and do QA

private:
  Int_t fPrintLevel;
  // variables
  Int_t fRunNumber;

  const double Mpion = 0.139570;
  const double Mkaon = 0.493677;
  const double Mproton = 0.938272;

  TString mOutName1;
  TString mOutName2;

  double fRhoval;
  // Rho objects
  StRhoParameter *GetRhoFromEvent(const char *name);

  // bad and dead tower list
  std::set<Int_t> badTowers;
  std::set<Int_t> deadTowers;

  // bad run list
  std::set<Int_t> badRuns;

  // base class pointer object
  StJetFrameworkPicoBase *mBaseMaker;

  vector<TLorentzVector> *fMCTracks;
  vector<TLorentzVector> *fMCTowers;

  vector<TLorentzVector> *fRecoTracks;
  vector<TLorentzVector> *fRecoTowers;

  TVector3 *fOrigin;
  pair<TVector3, TVector3> *fMCD0Information;
  pair<TVector3, TVector3> *fRecoD0Information;

  std::vector<double> fRhoValue;
  std::vector<double> rhovalueCS;
  std::vector<double> sigvalueCS;

  double fRhoValCS;
  double fSigValCS;

  TTree *mcjettree;
  TTree *recojettree;

  TTree *jettree;

  StJetTreeStruct fMCJetTree;
  StJetTreeStruct fRecoJetTree;
  StJetTreeStruct fRecoJetTreeFromPYTHIAEvent;
  StJetTreeStruct fRecoJetTreeCS;
  StJetTreeStruct fRecoJetTreeCS2;

  // maker names
  TString fAnalysisMakerName;
  TString fEventMixerMakerName;

  TString fMCJetMakerName;
  TString fRecoJetMakerName;

  StHIOverlay_Test *mHIOverlay;
  StHIRecoJetsCombined *mHIRecoJetsCombined;

  TClonesArray *mcJets;
  TClonesArray *recoJets;
  TClonesArray *recoJetsFromPYTHIAEvent;
  TClonesArray *recoJetsCS;
  TClonesArray *recoJetsCS2;

  TFile *fout;

  std::vector<TClonesArray *> fJetsArrMC;
  std::vector<TClonesArray *> fJetsArrRecoFromPYTHIAEvent; // These are just RECO jets from PYTHIA events without any background at all. If an event has been saved out, there will be an entry on this.
  std::vector<TClonesArray *> fJetsAreaArr;
  std::vector<TClonesArray *> fJetsCS1Arr;
  std::vector<TClonesArray *> fJetsCS2Arr;

  std::vector<double> fD0CorrectedPt;
  std::vector<ObjectsToPass> fObjectsToPass;

  TH1F *hDiffJetPt;

  ClassDef(StHIJetSaverCombined, 1)
};
#endif
