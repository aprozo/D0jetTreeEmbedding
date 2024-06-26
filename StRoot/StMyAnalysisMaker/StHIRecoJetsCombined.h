#ifndef StHIRecoJetsCombined_H
#define StHIRecoJetsCombined_H

// $Id$

#include <set>

// for clusters
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
class StEmcGeom;
class StEmcPosition2;

// ROOT classes
class TClonesArray;
class TObjArray;
class TList;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3;
class TProfile;
class TString;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoBTowHit;
class StEmcGeom;

// Jet classes
class StFJWrapper;
class StJetUtility;

// Centrality class
class StCentMaker;
class StRefMultCorr;
// D0 Includes
class StHIOverlay_Test;

// STAR includes
#include "StFJWrapper.h"
#include "FJ_includes.h"
#include "StJet.h"
#include "StJetFrameworkPicoBase.h"

namespace fastjet
{
  class PseudoJet;
}

/**
 * @brief General jet finder task implementing a wrapper for FastJet
 *
 * This class implements a wrapper for the FastJet jet finder. It allows to set a jet definition (jet algorithm, recombination scheme) and the
 * list of jet constituents. The jet finding is delegated to the class StFJWrapper which implements an interface to FastJet.
 *
 * The below is not functional yet:
 * The FastJet contrib utilities are available via the StJetUtility base class
 * and its derived classes. Utilities can be added via the AddUtility(StJetUtility*) method.
 * All the utilities added in the list will be executed. Users can implement new utilities
 * deriving a new class from StJetUtility to interface functionalities of the FastJet contribs.
 */
#ifndef __ObjectsToPass__H__
#define __ObjectsToPass__H__
struct ObjectsToPass
{ ObjectsToPass(){
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
#endif  //!__ObjectsToPass__H__


class StHIRecoJetsCombined : public StMaker
{
public:
  // jet type enumerator
  enum EJetType_t
  {
    kFullJet,
    kChargedJet,
    kNeutralJet
  };

  StHIRecoJetsCombined();
  StHIRecoJetsCombined(const char *name, double mintrackPt, bool dohistos, const char *outName);
  virtual ~StHIRecoJetsCombined();

  // needed class functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *opt = "");
  virtual Int_t Finish();

  int GetMatchedBtowID(StPicoTrack *trk);

  // booking of histograms (optional)
  void DeclareHistograms();
  void WriteHistograms();
  void WriteHadCorrQAHistograms();

  void FindJets(int iteration);

  // switches
  virtual void SetUsePrimaryTracks(Bool_t P) { doUsePrimTracks = P; }
  virtual void SetDebugLevel(Int_t l) { fDebugLevel = l; }
  virtual void SetRunFlag(Int_t f) { fRunFlag = f; }
  virtual void SetdoppAnalysis(Bool_t pp) { doppAnalysis = pp; }
  virtual void SetTurnOnCentSelection(Bool_t o) { fRequireCentSelection = o; }
  virtual void SetCentralityBinCut(Int_t c) { fCentralitySelectionCut = c; }
  virtual void SetdoConstituentSubtr(Bool_t c) { doConstituentSubtr = c; }

  virtual void SetMCAnalysisWithoutCent() { fMCEventsWithoutCent = kTRUE; }

  void SetForBackground() { fBackground = kTRUE; }
  virtual void SetPrintLevel(int i) { fPrintLevel = i; }
  // event setters
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

  // track setters
  void SetTrackEtaRange(Double_t etmi, Double_t etma)
  {
    fTrackEtaMin = etmi;
    fTrackEtaMax = etma;
  }
  void SetTrackPhiRange(Double_t ptmi, Double_t ptma)
  {
    fTrackPhiMax = ptmi;
    fTrackPhiMax = ptma;
  }
  virtual void SetDoEffCorr(Bool_t effcorr) { fDoEffCorr = effcorr; }
  virtual void SetDoCorrectTracksforEffBeforeJetReco(Bool_t ec) { doCorrectTracksforEffBeforeJetReco = ec; }
  virtual void SetTrackEfficiencyType(Int_t t) { fTrackEfficiencyType = t; }

  // common setters
  void SetClusName(const char *n) { fCaloName = n; }
  void SetTracksName(const char *n) { fTracksName = n; }
  void SetJetsName(const char *n) { fJetsName = n; }
  void SetJetAlgo(Int_t a) { fJetAlgo = a; }
  void SetJetType(Int_t t) { fJetType = t; }
  void SetRecombScheme(Int_t scheme) { fRecombScheme = scheme; }
  void SetMinJetArea(Double_t a) { fMinJetArea = a; }
  void SetMinJetPt(Double_t j) { fMinJetPt = j; }
  void SetRadius(Double_t r) { fRadius = r; }
  void SetGhostArea(Double_t gharea) { fGhostArea = gharea; }
  void SetJetEtaRange(Double_t emi, Double_t ema)
  {
    fJetEtaMin = emi;
    fJetEtaMax = ema;
  }
  void SetJetPhiRange(Double_t pmi, Double_t pma)
  {
    fJetPhiMin = pmi;
    fJetPhiMax = pma;
  }

  // set Jets for D0s, D0s Bg US, D0s Bg LS
  void SetD0Kind(Int_t kind) { fD0Kind = kind; }

  void SetMinJetTrackPt(Double_t min) { fMinJetTrackPt = min; }
  void SetMaxJetTrackPt(Double_t max) { fMaxJetTrackPt = max; }
  void SetJetTrackEtaRange(Double_t etmi, Double_t etma)
  {
    fJetTrackEtaMin = etmi;
    fJetTrackEtaMax = etma;
  }
  void SetJetTrackPhiRange(Double_t ptmi, Double_t ptma)
  {
    fJetTrackPhiMax = ptmi;
    fJetTrackPhiMax = ptma;
  }
  void SetJetTrackDCAcut(Double_t d) { fJetTrackDCAcut = d; }
  void SetJetTracknHitsFit(Double_t h) { fJetTracknHitsFit = h; }
  void SetJetTracknHitsRatio(Double_t r) { fJetTracknHitsRatio = r; }
  void SetMinJetTowerE(Double_t min) { mTowerEnergyMin = min; }
  void SetJetTowerERange(Double_t enmi, Double_t enmx)
  {
    fJetTowerEMin = enmi;
    fJetTowerEMax = enmx;
  }
  void SetJetTowerEtaRange(Double_t temi, Double_t temx)
  {
    fJetTowerEtaMin = temi;
    fJetTowerEtaMax = temx;
  }
  void SetJetTowerPhiRange(Double_t tpmi, Double_t tpmx)
  {
    fJetTowerPhiMin = tpmi;
    fJetTowerPhiMax = tpmx;
  }
  void SetMinJetClusPt(Double_t min) { fMinJetClusPt = min; }
  void SetMinJetClusE(Double_t min) { fMinJetClusE = min; }

  void SetLocked() { fLocked = kTRUE; }
  void SetTrackEfficiency(Double_t t) { fTrackEfficiency = t; }
  void SetLegacyMode(Bool_t mode) { fLegacyMode = mode; }
  void SetFillGhost(Bool_t b = kTRUE) { fFillGhost = b; }

  // for jet substructure routines
  StJetUtility *AddUtility(StJetUtility *utility);
  TObjArray *GetUtilities() { return fUtilities; }

  // jets
  // TClonesArray*          GetJets()                        { return fJets; }
  // TClonesArray*          GetJetsBGsub()                   { return fJetsBGsub; }

  std::vector<TClonesArray *> GetAreaSubtractedJets() { return fJetsAreaArr; }
  std::vector<TClonesArray *> GetCS1Jets() { return fJetsCS1Arr; }
  std::vector<TClonesArray *> GetCS2Jets() { return fJetsCS2Arr; }
  std::vector<double> GetD0CorrectedPt() { return fD0CorrectedPt; }
  std::vector<ObjectsToPass> GetObjectsToPass() { return fObjectsToPass; }
  std::vector<double> GetArrayOfRhoValue() { return fRhoValue; }

  // getters
  Double_t GetGhostArea() { return fGhostArea; }
  const char *GetJetsName() { return fJetsName.Data(); }
  const char *GetJetsTag() { return fJetsTag.Data(); }
  Double_t GetJetEtaMin() { return fJetEtaMin; }
  Double_t GetJetEtaMax() { return fJetEtaMax; }
  Double_t GetJetPhiMin() { return fJetPhiMin; }
  Double_t GetJetPhiMax() { return fJetPhiMax; }
  UInt_t GetJetType() { return fJetType; }
  UInt_t GetJetAlgo() { return fJetAlgo; }
  Int_t GetRecombScheme() { return fRecombScheme; }
  Bool_t GetLegacyMode() { return fLegacyMode; }
  Double_t GetMinJetArea() { return fMinJetArea; }
  Double_t GetMinJetPt() { return fMinJetPt; }
  Double_t GetRadius() { return fRadius; }
  Double_t GetJetTrackEtaMin() { return fJetTrackEtaMin; }
  Double_t GetJetTrackEtaMax() { return fJetTrackEtaMax; }
  Double_t GetJetTrackPhiMin() { return fJetTrackPhiMin; }
  Double_t GetJetTrackPhiMax() { return fJetTrackPhiMax; }
  Double_t GetTrackEfficiency() { return fTrackEfficiency; }
  Double_t GetJetTowerEtaMin() { return fJetTowerEtaMin; }
  Double_t GetJetTowerEtaMax() { return fJetTowerEtaMax; }
  Double_t GetJetTowerPhiMin() { return fJetTowerPhiMin; }
  Double_t GetJetTowerPhiMax() { return fJetTowerPhiMax; }

  Bool_t IsLocked() const;

  std::set<Int_t> GetBadTowers() { return badTowers; }
  std::set<Int_t> GetDeadTowers() { return deadTowers; }

  double MaxTrackPt() { return fMaxTrackPt; }
  double MaxTowerEtBeforeHC() { return fMaxTowerEtBeforeHC; }
  double MaxTowerEtAfterHC() { return fMaxTowerEtAfterHC; }

  // set hadronic correction fraction and type for matched tracks to towers
  void SetHadronicCorrFrac(float frac) { mHadronicCorrFrac = frac; }
  void SetJetHadCorrType(Int_t hct) { fJetHadCorrType = hct; }

protected:
  // this 1st version is deprecated as the parameters are global for the class and already set
  void FindJets(TObjArray *tracks, TObjArray *towers, Int_t algo, Double_t radius);
  void FindJets();
  // Int_t FindJets(); // use this if want to return NJets found
  void FillJetConstituents(StJet *jet, std::vector<fastjet::PseudoJet> &constituents,
                           std::vector<fastjet::PseudoJet> &constituents_sub, Int_t flag = 0, TString particlesSubName = "");
  Bool_t AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert);    // track accept cuts function
  Bool_t AcceptJetTrack(StPicoTrack *trk, Float_t B, TVector3 Vert); // jet track accept cuts function
  Bool_t AcceptJetTower(StPicoBTowHit *tower, Int_t towerID);        // jet tower accept cuts function
  Int_t GetCentBin(Int_t cent, Int_t nBin) const;                    // centrality bin
  Bool_t GetMomentum(TVector3 &mom, const StPicoBTowHit *tower, Double_t mass, Int_t towerID, Double_t CorrectedEnergy) const;
  void FillEmcTriggersArr();
  Double_t GetMaxTrackPt(); // find max track pt in event
  Double_t GetMaxTowerEt(); // find max tower Et in event
  void RunEventQA();        // function to fill some event QA plots
  void SetSumw2();          // set errors weights

  // may not need any of these except fill jet branch if I want 2 different functions
  void FillJetBranch(int iter);
  void FillJetBGBranch(int iter);
  void InitUtilities();
  void PrepareUtilities();
  void ExecuteUtilities(StJet *jet, Int_t ij);
  void TerminateUtilities();

  Bool_t GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

  double dPhi(double phi1, double phi2);

  // switches
  Bool_t fD0Analysis;                        // Switch on d0analysis
  Bool_t doWriteHistos;                      // write QA histos
  Bool_t doUsePrimTracks;                    // primary track switch
  Int_t fDebugLevel;                         // debug printout level
  Int_t fRunFlag;                            // Run Flag numerator value
  Bool_t doppAnalysis;                       // use pp analysis data
  Bool_t fRequireCentSelection;              // require particular centrality bin
  Bool_t doConstituentSubtr;                 // run constituent subtractor
  Bool_t fDoEffCorr;                         // efficiency correction to tracks
  Bool_t doCorrectTracksforEffBeforeJetReco; // efficiency correction to tracks - not used yet (using 'fDoEffCorr')
  Int_t fTrackEfficiencyType;                // track efficiency type: pt-eta, pt, eta
  Bool_t fMCEventsWithoutCent;               // MC Simulation Without Centrality Or Trigger

  Bool_t fBackground;
  Bool_t fPrintLevel;
  // event cuts
  Double_t fEventZVtxMinCut;     // min event z-vertex cut
  Double_t fEventZVtxMaxCut;     // max event z-vertex cut
  Int_t fCentralitySelectionCut; // centrality selection cut
  Double_t fMaxEventTrackPt;     // max track pt in the event (to cut on)
  Double_t fMaxEventTowerEt;     // max tower Et in the event (to cut on)
  Bool_t doRejectBadRuns;        // switch to reject bad runs and thus skip from analysis

  // event variables
  Double_t Bfield;  // event Bfield
  TVector3 mVertex; // event vertex 3-vector
  Double_t zVtx;    // z-vertex component
  Int_t fRunNumber; // run number

  // event selection types
  UInt_t fEmcTriggerEventType; // Physics selection of event used for signal - HT or JP
  UInt_t fMBEventType;         // MB selection
  UInt_t fTriggerToUse;        // trigger to use for analysis
  Int_t fEmcTriggerArr[8];     // EMCal triggers array: used to select signal and do QA

  // tower to firing trigger type matched array
  Bool_t fTowerToTriggerTypeHT1[4800]; // Tower with corresponding HT1 trigger type array
  Bool_t fTowerToTriggerTypeHT2[4800]; // Tower with corresponding HT2 trigger type array
  Bool_t fTowerToTriggerTypeHT3[4800]; // Tower with corresponding HT3 trigger type array

  // centrality
  Double_t fCentralityScaled; // scaled by 5% centrality
  Int_t ref16;                // multiplicity bin (16)
  Int_t ref9;                 // multiplicity bin (9)

  // output file name string
  TString mOutName;
  TString fTracksName; // name of track collection
  TString fCaloName;   // name of calo cluster collection
  TString fJetsName;   // name of jet collection

  Int_t fJetAlgo;      // jet algorithm (kt, akt, etc)
  Int_t fJetType;      // jet type (full, charged, neutral)
  Int_t fRecombScheme; // recombination scheme used by fastjet

  StFJWrapper fjw; //! fastjet wrapper

  // jet attributes
  Double_t fRadius;     // jet radius
  TString fJetsTag;     // tag of jet collection (usually = "Jets")
  Double_t fMinJetArea; // min area to keep jet in output
  Double_t fMinJetPt;   // min jet pt to keep jet in output
  Double_t fJetPhiMin;  // minimum phi to keep jet in output
  Double_t fJetPhiMax;  // maximum phi to keep jet in output
  Double_t fJetEtaMin;  // minimum eta to keep jet in output
  Double_t fJetEtaMax;  // maximum eta to keep jet in output
  Double_t fGhostArea;  // ghost area

  // set Jets for D0s, D0s Bg US, D0s Bg LS
  Double_t fD0Kind; // D0s (negative sign), Bg D0s LS (positive) 1->Standard cut, 2->Tight, 3->Loose

  // track attributes
  Double_t fMinJetTrackPt;      // min jet track transverse momentum cut
  Double_t fMaxJetTrackPt;      // max jet track transverse momentum cut
  Double_t fMinJetClusPt;       // min jet cluster transverse momentum cut
  Double_t fMinJetClusE;        // min jet cluster energy cut
  Double_t fMinJetTowerE;       // min jet tower energy cut - not used (use mTowerEnergyMin)
  Double_t fTrackEtaMin;        // min track eta cut
  Double_t fTrackEtaMax;        // max track eta cut
  Double_t fTrackPhiMin;        // min track phi cut
  Double_t fTrackPhiMax;        // max track phi cut
  Double_t fJetTrackEtaMin;     // min jet track eta cut
  Double_t fJetTrackEtaMax;     // max jet track eta cut
  Double_t fJetTrackPhiMin;     // min jet track phi cut
  Double_t fJetTrackPhiMax;     // max jet track phi cut
  Double_t fJetTrackDCAcut;     // max jet track dca cut
  Int_t fJetTracknHitsFit;      // requirement for track hits
  Double_t fJetTracknHitsRatio; // requirement for nHitsFit / nHitsMax
  Double_t fTrackEfficiency;    // artificial tracking inefficiency (0...1) - this is NOT used

  // tower attributes
  Double_t fJetTowerEMin;    // min jet tower energy cut
  Double_t fJetTowerEMax;    // max jet tower energy cut
  Double_t fJetTowerEtaMin;  // min jet tower eta cut
  Double_t fJetTowerEtaMax;  // max jet tower eta cut
  Double_t fJetTowerPhiMin;  // min jet tower phi cut
  Double_t fJetTowerPhiMax;  // max jet tower phi cut
  Double_t mTowerEnergyMin;  // min jet tower energy cut
  Float_t mHadronicCorrFrac; // hadronic correction fraction from 0.0 to 1.0
  Int_t fJetHadCorrType;     // hadronic correction type to be used

  // may not need some of next bools
  TObjArray *fUtilities; // jet utilities (gen subtractor, constituent subtractor etc.)
  Bool_t fLocked;        // true if lock is set
  Bool_t fIsInit;        //!=true if already initialized
  Bool_t fLegacyMode;    //!=true to enable FJ 2.x behavior
  Bool_t fFillGhost;     //!=true ghost particles will be filled in StJet obj

  // jet and jet constituent objects
  TClonesArray *fJets;    //! jet collection
  TClonesArray *fJetsCS1; //! jet background subtracted collection
  TClonesArray *fJetsCS2; //! jet background subtracted collection

  std::vector<TClonesArray *> fJetsAreaArr;
  std::vector<TClonesArray *> fJetsCS1Arr;
  std::vector<TClonesArray *> fJetsCS2Arr;

  std::vector<ObjectsToPass> fObjectsToPass;
  std::vector<double> fD0CorrectedPt;
  std::vector<double> fRhoValue;

  vector<fastjet::PseudoJet> fFull_Event;      //! jet input vectors
  vector<fastjet::PseudoJet> fConstituents;    //! jet constituents
  vector<fastjet::PseudoJet> fConstituentsCS1; //! jet constituents
  vector<fastjet::PseudoJet> fConstituentsCS2; //! jet constituents

  vector<TLorentzVector> *fRecoTracks;
  vector<TLorentzVector> *fRecoTowers;

  pair<int, int> *fMCEventInfo;

  Double_t fMaxTrackPt;
  Double_t fMaxTowerEtBeforeHC;
  Double_t fMaxTowerEtAfterHC;

  // Emc geometry
  StEmcGeom *mGeom;

  static const Int_t fgkConstIndexShift; //! contituent index shift

  // track efficiency file
  TFile *fEfficiencyInputFile;

private:
  StMuDst *mu;                        // muDst object
  StPicoDstMaker *mPicoDstMaker;      // PicoDstMaker object
  StPicoDst *mPicoDst;                // PicoDst object
  StPicoEvent *mPicoEvent;            // PicoEvent object
  StCentMaker *mCentMaker;            // Centrality maker object
  StJetFrameworkPicoBase *mBaseMaker; // Base class object

  // position objection
  StEmcGeom *mBemcGeom;
  StEmcPosition2 *mEmcPosition; // Emc position object

  // centrality objects
  StRefMultCorr *grefmultCorr;

  // D0 tagger objects
  StHIOverlay_Test *mHIOverlay; // D0 Tagger object

  std::vector<std::vector<double>> d0TrackIndices;

  Int_t id0candidate;
  Int_t iTrack1;
  Int_t iTrack2;

  //  Double_t                mTowerMatchTrkIndex[4800];
  //  Bool_t                 mTowerStatusArr[4800];
  Double_t mTowerMatchTrkIndex[4800][7];
  Int_t mTowerStatusArr[4800];

  // bad and dead tower list
  std::set<Int_t> badTowers;
  std::set<Int_t> deadTowers;

  // bad run list
  std::set<Int_t> badRuns;

  StHIRecoJetsCombined(const StHIRecoJetsCombined &);            // not implemented
  StHIRecoJetsCombined &operator=(const StHIRecoJetsCombined &); // not implemented

  ClassDef(StHIRecoJetsCombined, 4) // Jet producing task
};
#endif
