#ifndef StSetupMaker_h
#define StSetupMaker_h

#include "StJetFrameworkPicoBase.h"
class StJetFrameworkPicoBase;

#include <set>

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TString;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;

class StSetupMaker : public StJetFrameworkPicoBase {
  public:

    // debug flags for specifics
    enum fDebugFlagEnum {
      kDebugNothing, // don't want lowest elements to be used
      kDebugEmcTrigger,
      kDebugGeneralEvt,
      kDebugCentrality,
    };

    StSetupMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName, bool mDoComments);
    virtual ~StSetupMaker();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();

    // switches
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }

    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
    virtual void            SetBadTowerListVers(UInt_t ibt)    { fBadTowerListVers = ibt; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }
    virtual void            SetBadRunListVers(Int_t i)         { fBadRunListVers = i; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }

    // return setup object for use in other classes
    std::set<Int_t>         GetBadTowers()                   { return badTowers          ; }
    std::set<Int_t>         GetDeadTowers()                  { return deadTowers         ; }
    std::set<Int_t>         GetBadRuns()                     { return badRuns            ; }

    void                    ResetBadTowerList( );
    void                    ResetDeadTowerList( );
    Bool_t                  AddBadTowers(TString csvfile);
    Bool_t                  AddDeadTowers(TString csvfile);
    Bool_t                  IsTowerOK( Int_t mTowId );
    Bool_t                  IsTowerDead( Int_t mTowId );

    void                    ResetBadRunList( );
    Bool_t                  AddBadRuns(TString csvfile);
    Bool_t                  IsRunOK( Int_t mRunId );

  protected:
    // functions

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB

  private:
    // switches
    bool                    doComments;

    // bad and dead tower list functions and arrays
//    void                   ResetBadTowerList( );
//    void                   ResetDeadTowerList( );
//    Bool_t                 AddBadTowers(TString csvfile);
//    Bool_t                 AddDeadTowers(TString csvfile);
//    Bool_t                 IsTowerOK( Int_t mTowId );
//    Bool_t                 IsTowerDead( Int_t mTowId );
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list 
//    void                   ResetBadRunList( );
//    Bool_t                 AddBadRuns(TString csvfile);
//    Bool_t                 IsRunOK( Int_t mRunId );
    std::set<Int_t>        badRuns;

    // maker names
    TString                fAnalysisMakerName;

    ClassDef(StSetupMaker, 1)
};
#endif
