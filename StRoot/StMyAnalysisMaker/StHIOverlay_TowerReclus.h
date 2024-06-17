#ifndef StHIOverlay_TowerReclus_h
#define StHIOverlay_TowerReclus_h

#include "StJetFrameworkPicoBase.h"

#include "StEmcUtil/geometry/StEmcGeom.h"

#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TF1.h"
#include "TGraph.h"

#include "StFJWrapper.h"
#include "FJ_includes.h"
#include "StJet.h"

// old file kept
#include "StPicoConstants.h"

class StFJWrapper;

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

class StHIOverlay_TowerReclus : public StJetFrameworkPicoBase {
  public:

    StHIOverlay_TowerReclus(const char *name, StPicoDstMaker *picoMaker, const char *outName, const char* filename);
    virtual ~StHIOverlay_TowerReclus();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    void                    SetNumberOfEventsToOverLay(int a){fSetNumberOfEvents = a;}
    void SetMCJetPtCut(double a) {fMCJetPtCut = a;}


    virtual void            GetFirstNEvents(){randomevents = kFALSE;}
    virtual void            SetACentralityBinForTest(int i) {if (i >= 0 && i <= 2) fCentBin = i; }
    virtual void            SetPrintLevel(int i)                {fPrintLevel = i;}


    TVector3      *GetMCOrigin(){return fOrigin;}
    vector<TLorentzVector> *GetMCTracks(){return fMCEventTracks;}
    vector<TLorentzVector> *GetMCTowers(){return fMCEventTowers;}
    pair<TVector3, TVector3> *GetMCD0()  {return fMCD0Information;}

    vector<TLorentzVector> *GetRecoTracks(){return fRecoEventTracks;}
    vector<TLorentzVector> *GetRecoTowers(){return fRecoEventTowers;}
    pair<TVector3, TVector3> *GetRecoD0()  {return fRecoD0Information;}
    pair<int, int>        *GetMCEventInfo(){return fMCEventInfo;}

    int GetTheNumberOfEventsToOverLay() {return fNumberofeventsoverlayed;}

    std::vector<TClonesArray *> GetJets()                        { return fJetsArr; }
    std::vector<TClonesArray *> GetRecoJets()                        { return fRecoJetsArr; }


    //// Tree Variables

    // Fixed size dimensions of array or collections stored in the TTree if any.
    static const Int_t kMaxEvent = 1;
    static const Int_t kMaxTrack = 1000;
    static const Int_t kMaxEmcTrigger = 89;
    static const Int_t kMaxMtdTrigger = 1;
    static const Int_t kMaxBTowHit = 4800;
    static const Int_t kMaxBTofHit = 83;
    static const Int_t kMaxMtdHit = 4;
    static const Int_t kMaxBbcHit = 1;
    static const Int_t kMaxEpdHit = 1;
    static const Int_t kMaxFmsHit = 1;
    static const Int_t kMaxEmcPidTraits = 7;
    static const Int_t kMaxBTofPidTraits = 31;
    static const Int_t kMaxMtdPidTraits = 1;
    static const Int_t kMaxTrackCovMatrix = 107;
    static const Int_t kMaxBEmcSmdEHit = 1;
    static const Int_t kMaxBEmcSmdPHit = 1;
    static const Int_t kMaxETofHit = 1;
    static const Int_t kMaxETofPidTraits = 1;
    static const Int_t kMaxMcVertex = 1000;
    static const Int_t kMaxMcTrack = 1000;

    // Declaration of leaf types
    Int_t           Event_;
    Int_t           Event_mRunId[kMaxEvent];   //[Event_]
    Int_t           Event_mEventId[kMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexX[kMaxEvent];
    Float_t         Event_mPrimaryVertexY[kMaxEvent];
    Float_t         Event_mPrimaryVertexZ[kMaxEvent];

    Int_t           Track_;
    Float_t         Track_mGMomentumX[kMaxTrack];   //[Track_]
    Float_t         Track_mGMomentumY[kMaxTrack];   //[Track_]
    Float_t         Track_mGMomentumZ[kMaxTrack];   //[Track_]
    Float_t         Track_mOriginX[kMaxTrack];   //[Track_]
    Float_t         Track_mOriginY[kMaxTrack];   //[Track_]
    Float_t         Track_mOriginZ[kMaxTrack];   //[Track_]
   
    Char_t          Track_mNHitsFit[kMaxTrack];   //[Track_]
    UChar_t         Track_mNHitsMax[kMaxTrack];   //[Track_]
    
    Short_t         Track_mNSigmaPion[kMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaKaon[kMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaProton[kMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaElectron[kMaxTrack];   //[Track_]
    Short_t         Track_mBEmcMatchedTowerIndex[kMaxTrack];   //[Track_]
    UShort_t        Track_mIdTruth[kMaxTrack];   //[Track_]
    UShort_t        Track_mQATruth[kMaxTrack];   //[Track_]

    Int_t           BTowHit_;
    Short_t         BTowHit_mE[kMaxBTowHit];   //[BTowHit_]
    
    Int_t           McVertex_;
    Int_t           McVertex_mId[kMaxMcVertex];   //[McVertex_]
    UShort_t        McVertex_mNoDaughters[kMaxMcVertex];   //[McVertex_]
    Int_t           McVertex_mIdParTrk[kMaxMcVertex];   //[McVertex_]
    Int_t           McVertex_mIsInterm[kMaxMcVertex];   //[McVertex_]
    Float_t         McVertex_mTime[kMaxMcVertex];   //[McVertex_]
    Float_t         McVertex_mVx[kMaxMcVertex];   //[McVertex_]
    Float_t         McVertex_mVy[kMaxMcVertex];   //[McVertex_]
    Float_t         McVertex_mVz[kMaxMcVertex];   //[McVertex_]
    Int_t           McTrack_;
    UShort_t        McTrack_mId[kMaxMcTrack];   //[McTrack_]
    Int_t           McTrack_mGePid[kMaxMcTrack];   //[McTrack_]
    Char_t          McTrack_mCharge[kMaxMcTrack];   //[McTrack_]
    UChar_t         McTrack_mHits[kMaxMcTrack][22];   //[McTrack_]
    Float_t         McTrack_mPx[kMaxMcTrack];   //[McTrack_]
    Float_t         McTrack_mPy[kMaxMcTrack];   //[McTrack_]
    Float_t         McTrack_mPz[kMaxMcTrack];   //[McTrack_]
    Float_t         McTrack_mE[kMaxMcTrack];   //[McTrack_]
    Bool_t          McTrack_mIsFromShower[kMaxMcTrack];   //[McTrack_]
    Short_t         McTrack_mIdVtxStart[kMaxMcTrack];   //[McTrack_]
    Short_t         McTrack_mIdVtxStop[kMaxMcTrack];   //[McTrack_]
    Short_t         McTrack_mIdVtxItrmd[kMaxMcTrack];   //[McTrack_]

    // List of branches
    TBranch        *b_Event_;   //!
    TBranch        *b_Event_mRunId;   //!
    TBranch        *b_Event_mEventId;   //!
    TBranch        *b_Event_mPrimaryVertexX;   //!
    TBranch        *b_Event_mPrimaryVertexY;   //!
    TBranch        *b_Event_mPrimaryVertexZ;   //!
    
    TBranch        *b_Track_;   //!
    TBranch        *b_Track_mGMomentumX;   //!
    TBranch        *b_Track_mGMomentumY;   //!
    TBranch        *b_Track_mGMomentumZ;   //!
    TBranch        *b_Track_mOriginX;   //!
    TBranch        *b_Track_mOriginY;   //!
    TBranch        *b_Track_mOriginZ;   //!
    
    TBranch        *b_Track_mNHitsFit;   //!
    TBranch        *b_Track_mNHitsMax;   //!
    TBranch        *b_Track_mNSigmaPion;   //!
    TBranch        *b_Track_mNSigmaKaon;   //!
    TBranch        *b_Track_mNSigmaProton;   //!
    TBranch        *b_Track_mNSigmaElectron;   //!
    TBranch        *b_Track_mTopologyMap;   //!
    TBranch        *b_Track_mBEmcMatchedTowerIndex;   //!
    TBranch        *b_Track_mIdTruth;   //!
    TBranch        *b_Track_mQATruth;   //!

    TBranch        *b_BTowHit_;   //!
    TBranch        *b_BTowHit_mE;   //!

    TBranch        *b_McVertex_;   //!
    TBranch        *b_McVertex_mId;   //!
    TBranch        *b_McVertex_mNoDaughters;   //!
    TBranch        *b_McVertex_mIdParTrk;   //!
    TBranch        *b_McVertex_mIsInterm;   //!
    TBranch        *b_McVertex_mTime;   //!
    TBranch        *b_McVertex_mVx;   //!
    TBranch        *b_McVertex_mVy;   //!
    TBranch        *b_McVertex_mVz;   //!
    TBranch        *b_McTrack_;   //!
    TBranch        *b_McTrack_mId;   //!
    TBranch        *b_McTrack_mGePid;   //!
    TBranch        *b_McTrack_mCharge;   //!
    TBranch        *b_McTrack_mHits;   //!
    TBranch        *b_McTrack_mPx;   //!
    TBranch        *b_McTrack_mPy;   //!
    TBranch        *b_McTrack_mPz;   //!
    TBranch        *b_McTrack_mE;   //!
    TBranch        *b_McTrack_mIsFromShower;   //!
    TBranch        *b_McTrack_mIdVtxStart;   //!
    TBranch        *b_McTrack_mIdVtxStop;   //!
    TBranch        *b_McTrack_mIdVtxItrmd;   //!


    Bool_t IsItASoftEvent();
    void AssignTreeVariables();
    void GetD0AndDaughters();
    vector <int> *GetMCTrackListForVertex();
    void  GetAllTracksFromVertex(int vertexid, vector <int> &trackvec);
    void  MCTracksToDiscard(int D0num);
    void  PrepareSetOfMCTracks(int eventnum, int D0num);
    Bool_t DoesItHaveAGoodMCJet(int eventnum);
    Bool_t DoesItHaveAGoodRecoJet(int eventnum);
    Bool_t IsTrackADescendantOfD0Daughters(int trackid);
    void   PrepareSetOfRecoInput(int eventnum, int D0num);
    int GetMatchedRecoTrackFromMCTrack(int mctrkid);

  protected:
    void                    RunTracks();
    void                    RunTowers();
    void                    RunJets();
    Double_t                RelativePhi(Double_t mphi, Double_t vphi) const;      // relative jet track angle
    Double_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const; // relative jet event plane angle
    void                    FillEmcTriggers();                          // EmcTrigger counter histo
    Bool_t                  DoComparison(int myarr[], int elems);
    void                    SetSumw2(); // set errors weights

    void                    PickARandomFile(); 
    void                    SampleMCEvents();
    void                    SampleMCEventsForTest();
    int                     GetMatchedBtowID(int trkbemcid, TVector3 gMom, TVector3 org, int charge);
    TVector3                FastSimMom(TVector3 p, int pid);
    Bool_t                  KeepTrack(int particleid, int centralitybin, double pt);

    Bool_t                  GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

    int                    fPrintLevel;

    // switches
    Bool_t                  doPrintEventCounter;     // print event # switch
    Bool_t                  fDoEffCorr;              // efficiency correction to tracks

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

    Int_t                   fCentBin;

    double                  fMCJetPtCut;
    
    TClonesArray             *fJets;                   //!jet collection
    std::vector<TClonesArray *> fJetsArr;
    std::vector<TClonesArray *> fRecoJetsArr;




  private:

    double pi0mass = Pico::mMass[0]; // GeV

    // variables
    Int_t                   fRunNumber;

    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);

    int centralitybinforefficiency;

    // histos
    TH1F *hCentrality;//!
    TH1F *hMultiplicity;//!
 
    // jet histos
    TH1F *hJetPt;//!
    TH1F *hJetCorrPt;//!
    TH1F *hNEF;

    TFile *f;
    TTree *fMCPico;


    //This is where I will store MC event by event information. All the processing will be done in this class, and once finished, we will have no access to the picodsts we got the files from.
    //All the track selection cuts for jets are done here.
    vector<TLorentzVector> fMCEventTracks[100]; // For charged particles
    vector<TLorentzVector> fMCEventTowers[100]; // For neutral particles
    vector<TLorentzVector> fRecoEventTracks[100]; // For charged tracks
    vector<TLorentzVector> fRecoEventTowers[100]; // For neutral towers
    pair<TVector3, TVector3>         fMCD0Information[100]; // Pion momenta, Kaon momenta (3 components each)
    pair<TVector3, TVector3>         fRecoD0Information[100]; // Pion momenta, Kaon momenta (3 components each)
    TVector3               fOrigin[100]; // MC Event Origin Information
    pair<int, int>         fMCEventInfo[100]; //RunID, EventID

    vector<int> vertexids;
    vector<int> pionids;
    vector<int> kaonids;

    vector<int> matchedpionids;
    vector<int> matchedkaonids;

    vector <int> *fvertextotrack; //Array of vectors which saves the MC track IDs matched to each vertex (START) // This can be private as I am not calling this function outside this class
    // The dropped MC tracks will be different for each D0 because we only want to drop final state tracks which are from the KPi from D0 decaying.
    vector<int>     fDroppedMCTracks; // Dropped MC Tracks (This will be tracks which came from the KPi from D0 decaying. We don't want them.)

    int centraleventsfound;
    // vector<int> randomlisttoevents;

    int fSetNumberOfEvents;
    int fNumberofeventsoverlayed;

    bool randomevents;

    // bad and dead tower list
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list
    std::set<Int_t>        badRuns;

    // base class pointer object
    StJetFrameworkPicoBase *mBaseMaker;

    // position objection
    StEmcGeom             *mBemcGeom;

    // maker names
    TString                fAnalysisMakerName;
    TString                fEventMixerMakerName;

    TString                fMCFileListName;

    vector<TString>         filenamesforHIOverlay;

    TF1* fKaonMomResolution;
    TF1* fPionMomResolution;
    TF1* fProtonMomResolution;

    TGraph* fPionWeight[3];
    TGraph* fKaonWeight[3];
    TGraph* fProtonWeight[3];
    TGraph* fAProtonWeight[3];

    Double_t               fMinJetTrackPt;          // min jet track transverse momentum cut
    Double_t               fMaxJetTrackPt;          // max jet track transverse momentum cut

    Double_t               fJetTrackEtaMin;         // min jet track eta cut
    Double_t               fJetTrackEtaMax;         // max jet track eta cut
    Double_t               fJetTrackPhiMin;         // min jet track phi cut
    Double_t               fJetTrackPhiMax;         // max jet track phi cut
    Double_t               fJetTrackDCAcut;         // max jet track dca cut
    Int_t                  fJetTracknHitsFit;       // requirement for track hits
    Double_t               fJetTracknHitsRatio;     // requirement for nHitsFit / nHitsMax

    Double_t               fJetTowerEMin;           // min jet tower energy cut
    Double_t               fJetTowerEMax;           // max jet tower energy cut
    Double_t               fJetTowerEtaMin;         // min jet tower eta cut
    Double_t               fJetTowerEtaMax;         // max jet tower eta cut
    Double_t               fJetTowerPhiMin;         // min jet tower phi cut
    Double_t               fJetTowerPhiMax;         // max jet tower phi cut
    Double_t               mTowerEnergyMin;         // min jet tower energy cut


    ClassDef(StHIOverlay_TowerReclus, 1)
};
#endif
