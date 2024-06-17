#ifndef StHIOverlayAngularities_h
#define StHIOverlayAngularities_h

#include "StJetFrameworkPicoBase.h"

#include "StEmcUtil/geometry/StEmcGeom.h"

#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

#include "TGraph.h"
#include "StJetTreeStruct.h"
#include "StFJWrapper.h"
#include "FJ_includes.h"
#include "StJet.h"

// old file kept
#include "StPicoConstants.h"

class StFJWrapper;
class StJetFrameworkPicoBase;

// ROOT classes
class StJet;

class TString;
class TVector3;
class TLorentzVector;
class TString;
class TTree;
class TBranch;
class TGraph;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

// jet-framework classes
class StCentMaker;
class StJetMakerTask;

class StJet;
class StRho;
class StRhoParameter;

class StHIOverlayAngularities : public StJetFrameworkPicoBase
{
public:
  StHIOverlayAngularities(const char *name, StPicoDstMaker *picoMaker, const char *outName, const char *filename);
  virtual ~StHIOverlayAngularities();

  // class required functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *opt = "");
  virtual Int_t Finish();

  void SetNumberOfEventsToOverLay(Int_t a) { fNumberOfEventsToOverLay = a; }

  virtual void SetACentralityBinForTest(Int_t i)
  {
    if (i >= 0 && i <= 2)
      fCentBin = i;
  }
  virtual void SetPrintLevel(Int_t i) { fPrintLevel = i; }

  void SetJetTrackDCAcut(Double_t d) { fJetTrackDCAcut = d; }

  void DoTrackingEfficiency(const Double_t &percentage)
  {
    fTrackingEfficiency = kTRUE;
    fTrackingEfficiencyPercentage = percentage;
  } // percentage of tracks to keep

  void SetJetAlgo(Int_t a) { fJetAlgo = a; }
  void SetJetType(Int_t t) { fJetType = t; }

  void SetMinJetTrackPt(Double_t min) { fMinJetTrackPt = min; }
  void SetMaxJetTrackPt(Double_t max) { fMaxJetTrackPt = max; }
  void SetJetTrackEtaRange(Double_t etaMin, Double_t etaMax)
  {
    fJetTrackEtaMin = etaMin;
    fJetTrackEtaMax = etaMax;
  }
  void SetJetTrackPhiRange(Double_t phiMin, Double_t phiMax)
  {
    fJetTrackPhiMax = phiMin;
    fJetTrackPhiMax = phiMax;
  }

  void SetRecombScheme(Int_t scheme) { fRecombScheme = scheme; }
  void SetMinJetArea(Double_t a) { fMinJetArea = a; }
  void SetMinJetPt(Double_t j) { fMinJetPt = j; }
  void SetGhostArea(Double_t gharea) { fGhostArea = gharea; }
  void SetJetEtaRange(Double_t etaMin, Double_t etaMax)
  {
    fJetEtaMin = etaMin;
    fJetEtaMax = etaMax;
  }
  void SetJetPhiRange(Double_t phiMin, Double_t phiMax)
  {
    fJetPhiMin = phiMin;
    fJetPhiMax = phiMax;
  }

  void SetMinJetTowerE(Double_t min) { mTowerEnergyMin = min; }
  void SetJetTowerERange(Double_t energyMin, Double_t energyMax)
  {
    fJetTowerEMin = energyMin;
    fJetTowerEMax = energyMax;
  }
  void SetJetTowerEtaRange(Double_t etaMin, Double_t etaMax)
  {
    fJetTowerEtaMin = etaMin;
    fJetTowerEtaMax = etaMax;
  }
  void SetJetTowerPhiRange(Double_t phiMin, Double_t phiMax)
  {
    fJetTowerPhiMin = phiMin;
    fJetTowerPhiMax = phiMax;
  }

  // set hadronic correction fraction and type for matched tracks to towers
  void SetHadronicCorrFrac(float frac) { mHadronicCorrFrac = frac; }
  void SetJetHadCorrType(Int_t hct) { fJetHadCorrType = hct; }

  virtual void SetEventZVtxRange(Double_t zmi, Double_t zma)
  {
    fEventZVtxMinCut = zmi;
    fEventZVtxMaxCut = zma;
  }
  virtual void SetEmcTriggerEventType(UInt_t te) { fEmcTriggerEventType = te; }
  virtual void SetMBEventType(UInt_t mbe) { fMBEventType = mbe; }
  virtual void SetTriggerToUse(UInt_t ttu) { fTriggerToUse = ttu; }
  virtual void SetMaxEventTrackPt(Double_t mxpt) { fMaxEventTrackPt = mxpt; }
  virtual void SetMaxEventTowerEt(Double_t mxEt) { fMaxEventTowerEt = mxEt; }
  virtual void SetRejectBadRuns(Bool_t rj) { doRejectBadRuns = rj; }

  //// Tree Variables

  // Fixed size dimensions of array or collections stored in the TTree if any.
  static const Int_t kMaxEvent = 1;
  static const Int_t kMaxTrack = 1000;
  static const Int_t kMaxEmcTrigger = 89;
  static const Int_t kMaxMtdTrigger = 1;
  static const Int_t kMaxBTowHit = 4800;

  static const Int_t kMaxMcVertex = 6000;
  static const Int_t kMaxMcTrack = 6000;


  TTree *fMCPico;
  // Declaration of leaf types
  Int_t Event_;
  Int_t Event_mRunId[kMaxEvent];   //[Event_]
  Int_t Event_mEventId[kMaxEvent]; //[Event_]
  Float_t Event_mPrimaryVertexX[kMaxEvent];
  Float_t Event_mPrimaryVertexY[kMaxEvent];
  Float_t Event_mPrimaryVertexZ[kMaxEvent];

  Int_t Track_;
  Float_t Track_mGMomentumX[kMaxTrack]; //[Track_]
  Float_t Track_mGMomentumY[kMaxTrack]; //[Track_]
  Float_t Track_mGMomentumZ[kMaxTrack]; //[Track_]
  Float_t Track_mOriginX[kMaxTrack];    //[Track_]
  Float_t Track_mOriginY[kMaxTrack];    //[Track_]
  Float_t Track_mOriginZ[kMaxTrack];    //[Track_]

  Char_t Track_mNHitsFit[kMaxTrack];  //[Track_]
  UChar_t Track_mNHitsMax[kMaxTrack]; //[Track_]

  Short_t Track_mNSigmaPion[kMaxTrack];            //[Track_]
  Short_t Track_mNSigmaKaon[kMaxTrack];            //[Track_]
  Short_t Track_mNSigmaProton[kMaxTrack];          //[Track_]
  Short_t Track_mNSigmaElectron[kMaxTrack];        //[Track_]
  Short_t Track_mBEmcMatchedTowerIndex[kMaxTrack]; //[Track_]
  UShort_t Track_mIdTruth[kMaxTrack];              //[Track_]
  UShort_t Track_mQATruth[kMaxTrack];              //[Track_]

  Int_t BTowHit_;
  Short_t BTowHit_mE[kMaxBTowHit]; //[BTowHit_]

  Int_t McVertex_;
  Int_t McVertex_mId[kMaxMcVertex];             //[McVertex_]
  UShort_t McVertex_mNoDaughters[kMaxMcVertex]; //[McVertex_]
  Int_t McVertex_mIdParTrk[kMaxMcVertex];       //[McVertex_]
  Int_t McVertex_mIsInterm[kMaxMcVertex];       //[McVertex_]
  Float_t McVertex_mTime[kMaxMcVertex];         //[McVertex_]
  Float_t McVertex_mVx[kMaxMcVertex];           //[McVertex_]
  Float_t McVertex_mVy[kMaxMcVertex];           //[McVertex_]
  Float_t McVertex_mVz[kMaxMcVertex];           //[McVertex_]
  Int_t McTrack_;
  UShort_t McTrack_mId[kMaxMcTrack];         //[McTrack_]
  Int_t McTrack_mGePid[kMaxMcTrack];         //[McTrack_]
  Char_t McTrack_mCharge[kMaxMcTrack];       //[McTrack_]
  UChar_t McTrack_mHits[kMaxMcTrack][22];    //[McTrack_]
  Float_t McTrack_mPx[kMaxMcTrack];          //[McTrack_]
  Float_t McTrack_mPy[kMaxMcTrack];          //[McTrack_]
  Float_t McTrack_mPz[kMaxMcTrack];          //[McTrack_]
  Float_t McTrack_mE[kMaxMcTrack];           //[McTrack_]
  Bool_t McTrack_mIsFromShower[kMaxMcTrack]; //[McTrack_]
  Short_t McTrack_mIdVtxStart[kMaxMcTrack];  //[McTrack_]
  Short_t McTrack_mIdVtxStop[kMaxMcTrack];   //[McTrack_]
  Short_t McTrack_mIdVtxItrmd[kMaxMcTrack];  //[McTrack_]

  // List of branches
  TBranch *b_Event_;                //!
  TBranch *b_Event_mRunId;          //!
  TBranch *b_Event_mEventId;        //!
  TBranch *b_Event_mPrimaryVertexX; //!
  TBranch *b_Event_mPrimaryVertexY; //!
  TBranch *b_Event_mPrimaryVertexZ; //!

  TBranch *b_Track_;            //!
  TBranch *b_Track_mGMomentumX; //!
  TBranch *b_Track_mGMomentumY; //!
  TBranch *b_Track_mGMomentumZ; //!
  TBranch *b_Track_mOriginX;    //!
  TBranch *b_Track_mOriginY;    //!
  TBranch *b_Track_mOriginZ;    //!

  TBranch *b_Track_mNHitsFit;              //!
  TBranch *b_Track_mNHitsMax;              //!
  TBranch *b_Track_mNSigmaPion;            //!
  TBranch *b_Track_mNSigmaKaon;            //!
  TBranch *b_Track_mNSigmaProton;          //!
  TBranch *b_Track_mNSigmaElectron;        //!
  TBranch *b_Track_mTopologyMap;           //!
  TBranch *b_Track_mBEmcMatchedTowerIndex; //!
  TBranch *b_Track_mIdTruth;               //!
  TBranch *b_Track_mQATruth;               //!

  TBranch *b_BTowHit_;   //!
  TBranch *b_BTowHit_mE; //!

  TBranch *b_McVertex_;             //!
  TBranch *b_McVertex_mId;          //!
  TBranch *b_McVertex_mNoDaughters; //!
  TBranch *b_McVertex_mIdParTrk;    //!
  TBranch *b_McVertex_mIsInterm;    //!
  TBranch *b_McVertex_mTime;        //!
  TBranch *b_McVertex_mVx;          //!
  TBranch *b_McVertex_mVy;          //!
  TBranch *b_McVertex_mVz;          //!
  TBranch *b_McTrack_;              //!
  TBranch *b_McTrack_mId;           //!
  TBranch *b_McTrack_mGePid;        //!
  TBranch *b_McTrack_mCharge;       //!
  TBranch *b_McTrack_mHits;         //!
  TBranch *b_McTrack_mPx;           //!
  TBranch *b_McTrack_mPy;           //!
  TBranch *b_McTrack_mPz;           //!
  TBranch *b_McTrack_mE;            //!
  TBranch *b_McTrack_mIsFromShower; //!
  TBranch *b_McTrack_mIdVtxStart;   //!
  TBranch *b_McTrack_mIdVtxStop;    //!
  TBranch *b_McTrack_mIdVtxItrmd;   //!

  void ReadTreeMc();
  void FillTree(const Int_t &numberOfD0Events);
  void FillJet(StJet *jet, StJetTreeStruct &jetTree, const TVector3 & D0);
  void OutputTreeInit();

  void GetAllTracksFromVertex(const Int_t &vertexid, vector<Int_t> &trackvec);
  StJet *DoesItHaveAGoodD0Jet(vector<TLorentzVector> &eventTracks);

  void PrepareSetOfRecoInput(const Int_t & counterEvent, const Int_t & iD0);
  Int_t GetMatchedRecoTrackFromMCTrack(const Int_t & mctrkid);
  Bool_t IsAcceptedTrack(StPicoTrack *trk, const Float_t &B, const TVector3 &Vert);
  Bool_t IsAcceptedTrackAndPt(StPicoTrack *trk, const Float_t &B, const TVector3 &Vert);
  Bool_t IsAcceptedTower(StPicoBTowHit *tower, const Int_t &towerID);

  TVector3 FastSimMom(TVector3 p, Int_t pid);
  Bool_t GetMomentum(TVector3 &mom, const StPicoBTowHit *tower, Double_t mass, Int_t towerID, Double_t CorrectedEnergy) const;

  Int_t GetMatchedBtowID(StPicoTrack *trk);
  Bool_t KeepTrack(const Int_t & particleid, const Int_t & centralitybin, const Double_t &  pt);
  Bool_t GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;
  Int_t fPrintLevel;

  // switches
  Bool_t doPrintEventCounter; // prInt_t event # switch
  Bool_t fDoEffCorr;          // efficiency correction to tracks

  // event selection types
  UInt_t fTriggerToUse;        // trigger to use for analysis
  UInt_t fEmcTriggerEventType; // Physics selection of event used for signal
  UInt_t fMBEventType;         // Physics selection of event used for MB
  Int_t fEmcTriggerArr[8];     // EMCal triggers array: used to select signal and do QA
  Int_t fCentBin;

  Double_t pi0mass = Pico::mMass[0]; // GeV

  // variables
  Int_t fRunNumber;
  Int_t centralitybinforefficiency;
  Double_t fRhoVal;
  vector<Int_t> vertexids;
  vector<Int_t> pionids;
  vector<Int_t> kaonids;
  vector<Int_t> matchedpionids;
  vector<Int_t> matchedkaonids;
  map<int, vector<Int_t>> fVertexToTracks;
  // Array of vectors which saves the MC track IDs matched to each vertex (START) // This can be private as I am not calling this function outside this class
  //  The dropped MC tracks will be different for each D0 because we only want to drop final state tracks which are from the KPi from D0 decaying.
  vector<Int_t> fDroppedMCTracks; // Dropped MC Tracks (This will be tracks which came from the KPi from D0 decaying. We don't want them.)

  Int_t fNumberOfEventsToOverLay;
  Double_t fMinJetArea; // min area to keep jet in output
  Double_t fMinJetPt;   // min jet pt to keep jet in output
  Double_t fJetPhiMin;  // minimum phi to keep jet in output
  Double_t fJetPhiMax;  // maximum phi to keep jet in output
  Double_t fJetEtaMin;  // minimum eta to keep jet in output
  Double_t fJetEtaMax;  // maximum eta to keep jet in output
  Double_t fGhostArea;  // ghost area

  Double_t fMinJetTrackPt; // min jet track transverse momentum cut
  Double_t fMaxJetTrackPt; // max jet track transverse momentum cut

  Double_t fJetTrackEtaMin;     // min jet track eta cut
  Double_t fJetTrackEtaMax;     // max jet track eta cut
  Double_t fJetTrackPhiMin;     // min jet track phi cut
  Double_t fJetTrackPhiMax;     // max jet track phi cut
  Double_t fJetTrackDCAcut;     // max jet track dca cut
  Int_t fJetTracknHitsFit;      // requirement for track hits
  Double_t fJetTracknHitsRatio; // requirement for nHitsFit / nHitsMax

  Double_t fJetTowerEMin;   // min jet tower energy cut
  Double_t fJetTowerEMax;   // max jet tower energy cut
  Double_t fJetTowerEtaMin; // min jet tower eta cut
  Double_t fJetTowerEtaMax; // max jet tower eta cut
  Double_t fJetTowerPhiMin; // min jet tower phi cut
  Double_t fJetTowerPhiMax; // max jet tower phi cut
  Double_t mTowerEnergyMin; // min jet tower energy cut

  bool fTrackingEfficiency;               // track reconstruction efficiency
  Double_t fTrackingEfficiencyPercentage; // percentage of tracks to keep

  Double_t fMaxTowerEtBeforeHC;
  Double_t fMaxTowerEtAfterHC;
  Float_t mHadronicCorrFrac; // hadronic correction fraction from 0.0 to 1.0
  Int_t fJetHadCorrType;     // hadronic correction type to be used
  Double_t mTowerMatchTrkIndex[4800][7];
  Int_t mTowerStatusArr[4800];
  // bad and dead tower list

  // histos
  TH1F *hCentrality;   //!
  TH1F *hMultiplicity; //!
  TH1F *hJetPt;        //!
  TH1F *hJetCorrPt;    //!

  TFile *f;
  TFile *fout;

  // This is where I will store MC event by event information. All the processing will be done in this class, and once finished, we will have no access to the picodsts we got the files from.
  // All the track selection cuts for jets are done here.
  static const Int_t kMaxNumberOfD0Events = 100;
  vector<TLorentzVector> fMcEventTracks[kMaxNumberOfD0Events];         // For charged particles
  vector<TLorentzVector> fMcEventTowers[kMaxNumberOfD0Events];         // For neutral particles
  vector<TLorentzVector> fRecoMcEventTracks[kMaxNumberOfD0Events];     // For charged tracks
  vector<TLorentzVector> fRecoMcEventTowers[kMaxNumberOfD0Events];     // For neutral towers
  pair<TVector3, TVector3> fMcD0Information[kMaxNumberOfD0Events];     // Pion momenta, Kaon momenta (3 components each)
  pair<TVector3, TVector3> fMcRecoD0Information[kMaxNumberOfD0Events]; // Pion momenta, Kaon momenta (3 components each)
  TVector3 fOrigin[kMaxNumberOfD0Events];                              // MC Event Origin Information
  pair<int, int> fMcEventInfo[kMaxNumberOfD0Events];                   // RunID, EventID
  vector<TLorentzVector> fRecoTracks[kMaxNumberOfD0Events];
  vector<TLorentzVector> fRecoTowers[kMaxNumberOfD0Events];
  StJet *fMcJet[kMaxNumberOfD0Events];
  StJet *fMcRecoJet[kMaxNumberOfD0Events];
  StJet *fRecoJet[kMaxNumberOfD0Events];

  vector<fastjet::PseudoJet> fFull_Event; //! jet input vectors

  TTree *outputTree;

  StJetTreeStruct fMcJetTree;
  StJetTreeStruct fRecoJetTree;
  StJetTreeStruct fMcRecoJetTree;

TH1D *hRecoJetPt ;   
TH1D *hMcJetPt;
TH1D *hMcRecoJetPt;
TH1D *hMcD0Pt ;
TH1D *hMcRecoD0Pt;
TH1D *hMcJetZ;
TH1D *hRecoJetZ;
TH2D *hMcD0PtMcRecoD0Pt;
TH2D *hRecoJetPtMcJetPt;
TH2D *hRecoJetZMcJetZ;

TH1D *fractionNeutralToJetPt;

  std::set<Int_t> badTowers;
  std::set<Int_t> deadTowers;
  // bad run list
  std::set<Int_t> badRuns;

  // base class pointer object
  StJetFrameworkPicoBase *mBaseMaker;

  // position objection
  StEmcGeom *mBemcGeom;

  // maker names
  TString fAnalysisMakerName;
  TString fEventMixerMakerName;

  TString fMCFileListName;

  vector<TString> filenamesforHIOverlay;

  TF1 *fKaonMomResolution;
  TF1 *fPionMomResolution;
  TF1 *fProtonMomResolution;

  TGraph *fPionWeight[3];
  TGraph *fKaonWeight[3];
  TGraph *fProtonWeight[3];
  TGraph *fAProtonWeight[3];

  Int_t fJetAlgo;      // jet algorithm (kt, akt, etc)
  Int_t fJetType;      // jet type (full, charged, neutral)
  Int_t fRecombScheme; // recombination scheme used by fastjet

  StFJWrapper *fjw;     //! fastjet wrapper
                        // jet attributes
 
  ClassDef(StHIOverlayAngularities, 1)
};
#endif
