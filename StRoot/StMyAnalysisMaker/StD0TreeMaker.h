#ifndef StD0TreeMaker_h
#define StD0TreeMaker_h

#include "StJetFrameworkPicoBase.h"
#include "StD0TreeStruct.h"
#include <vector>

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
class TTree;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

// jet-framework classes
class StCentMaker;
class StEmcPosition2;
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;
class StTagD0Events;

class StD0TreeMaker : public StJetFrameworkPicoBase {
  public:

    StD0TreeMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StD0TreeMaker();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareTree();
    void    BookTree();
    void    WriteTree();

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp;}
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }

    virtual void            SetMinMaxCentrality(Int_t c1, Int_t c2) {fMinCentrality = c1; fMaxCentrality = c2;}

    void                    IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e);
    Int_t                   IsWhatParticle(StPicoTrack *trk);

    // jet setters
    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetMaxTowerEt(Double_t t)       { fTowerBias        = t; }    // tower bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    
    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

    // track setters
    virtual void            SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
    virtual void            SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;} // max track cut
    virtual void            SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
    virtual void            SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
    virtual void            SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
    virtual void            SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
    virtual void            SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }

    // tower setters
    virtual void            SetTowerERange(Double_t enmi, Double_t enmx) { fTowerEMinCut = enmi; fTowerEMaxCut = enmx; }
    virtual void            SetTowerEtaRange(Double_t temi, Double_t temx) { fTowerEtaMinCut = temi; fTowerEtaMaxCut = temx; }
    virtual void            SetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fTowerPhiMinCut = tpmi; fTowerPhiMaxCut = tpmx; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)    { fEmcTriggerEventType = te;  }
    virtual void            SetMBEventType(UInt_t mbe)           { fMBEventType = mbe; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Bool_t effcorr)          { fDoEffCorr = effcorr; }

    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)          { fCorrJetPt = cpt; }

    void                    SetD0Kind(Int_t kind)                { fD0Kind = kind  ; }


  protected:
    void                    RunTracks();
    void                    RunTowers();
    void                    RunJets();
    Double_t                RelativePhi(Double_t mphi, Double_t vphi) const;      // relative jet track angle
    Double_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const; // relative jet event plane angle
    void                    FillEmcTriggers();                          // EmcTrigger counter histo
    Bool_t                  DoComparison(int myarr[], int elems);
    void                    SetSumw2(); // set errors weights
    int                     GetJetPtBin(Double_t pt) ;
    int                     GetD0PtBin(Double_t pt) ;

    Double_t                standardPhi(Double_t phi);


    // switches
    Bool_t                  doPrintEventCounter;     // print event # switch
    Bool_t                  fDoEffCorr;              // efficiency correction to tracks

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

    // set Jets for D0s, D0s Bg US, D0s Bg LS
    Double_t                fD0Kind;                 // D0s, Bg D0s US, LS

  private:

    const double Mpion = 0.139570;
    const double Mkaon = 0.493677;
    const double Mproton = 0.938272;
    // variables
    Int_t                   fRunNumber;
    Int_t                   fMinCentrality;
    Int_t                   fMaxCentrality;
    
    double                  frefCorr2;
    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);

    // position object
    StEmcPosition2         *mEmcPosition;

    // D0 tagger objects
    StTagD0Events         *mD0Tagger;               // D0 Tagger object


    std::vector<std::vector<double> >  d0TrackIndices; 

    Int_t                 id0candidate;
    Int_t                 iTrack1;
    Int_t                 iTrack2;

    TClonesArray     *fJets;//!jet collection
    std::vector<TClonesArray *> fJetsArr;

    std::vector<double> rhovalue;

    StD0TreeStruct fJetTree;

    TTree *jettree;

    // bad and dead tower list
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list
    std::set<Int_t>        badRuns;

    // base class pointer object
    StJetFrameworkPicoBase *mBaseMaker;


    // maker names
    TString                fAnalysisMakerName;
    TString                fEventMixerMakerName;

    ClassDef(StD0TreeMaker, 1)
};
#endif
