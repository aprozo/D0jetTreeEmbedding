#ifndef StHIOverlay_h
#define StHIOverlay_h

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

class StHIOverlay : public StJetFrameworkPicoBase {
  public:

    StHIOverlay(const char *name, StPicoDstMaker *picoMaker, const char *outName, const char* filename);
    virtual ~StHIOverlay();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    void                    SetNumberOfEventsToOverLay(int a){fSetNumberOfEvents = a;}


    virtual void            GetFirstNEvents(){randomevents = kFALSE;}
    virtual void            SetPrintLevel(int i)                {fPrintLevel = i;}


    TVector3      *GetMCOrigin(){return fOrigin;}
    vector<TLorentzVector> *GetMCTracks(){return fMCEventTracks;}
    vector<TLorentzVector> *GetMCTowers(){return fMCEventTowers;}
    pair<TVector3, TVector3> *GetMCD0()  {return fMCD0Information;}

    vector<TLorentzVector> *GetRecoTracks(){return fRecoEventTracks;}
    vector<TLorentzVector> *GetRecoTowers(){return fRecoEventTowers;}
    pair<TVector3, TVector3> *GetRecoD0()  {return fRecoD0Information;}

    int GetTheNumberOfEventsToOverLay() {return fNumberofeventsoverlayed;}

    std::vector<TClonesArray *> GetJets()                        { return fJetsArr; }

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


    
    TClonesArray             *fJets;                   //!jet collection
    std::vector<TClonesArray *> fJetsArr;




  private:
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


    ClassDef(StHIOverlay, 1)
};
#endif
