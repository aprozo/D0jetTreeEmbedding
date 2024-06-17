#ifndef StThrowARandomTrack_H
#define StThrowARandomTrack_H

#include <iostream>
#include <fstream>
#include "StJetFrameworkPicoBase.h"
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
class TGraph;
class THnSparse;
class TProfile;
class TString;
class TVector3;
class TLorentzVector;

// STAR classes
//class StarRoot;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;

// class StPicoTrackCovMatrix;
class StRefMultCorr;

// jet-framework classes
class StCentMaker;
class StEmcPosition2;

class StThrowARandomTrack : public StJetFrameworkPicoBase {
  public:

    bool testrun = kTRUE;

    StThrowARandomTrack(const char* name, StPicoDstMaker *picoMaker, const char* outName);
    virtual ~StThrowARandomTrack();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();

    // particle identification
    Bool_t    IsAnAcceptableTrack(StPicoTrack *trk, bool dohistograms);
    Bool_t    IsTrackWithinJet(StJet *jet, StPicoTrack *trk);
    void      IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e);
    Int_t     IsWhatParticle(StPicoTrack *trk);
    void      InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum);
    Double_t  InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2);
    void      FillPidHistograms(StPicoTrack *trk);
    Int_t     TopologicalCuts(StPicoTrack *trk1, StPicoTrack *trk2, bool dohistograms);
    Bool_t    KeepTrack(int trackid);
    void      FindTracksToDrop();

    void      GetD0TracksPair();
    void      ReplacePoorlySmearedTracks();
    int       GetMatchedRecoTrackFromMCTrack(int mctrkid);


    TVector3  FastSimMom(TVector3 p, int pid);

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp;}
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }
    virtual void            SetTestRun(Bool_t trun)            { testrun = trun; }
    
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

    // topological cut setter

    virtual void            SetTopoCutsLevel(int l)             {fTopoLevel = l         ; }

    // tower setters
    virtual void            SetTowerERange(Double_t enmi, Double_t enmx) { fTowerEMinCut = enmi; fTowerEMaxCut = enmx; }
    virtual void            SetTowerEtaRange(Double_t temi, Double_t temx) { fTowerEtaMinCut = temi; fTowerEtaMaxCut = temx; }
    virtual void            SetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fTowerPhiMinCut = tpmi; fTowerPhiMaxCut = tpmx; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)    { fEmcTriggerEventType = te;  }
    virtual void            SetMBEventType(UInt_t mbe)           { fMBEventType = mbe; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Bool_t effcorr)          { fDoEffCorr = effcorr; }


    virtual void            SetAnalysisCentralityBin(Int_t centbin) {fAnaCentBin = centbin;} //Bin 0 = 0-10; Bin 1 = 10-40; Bin 2 = 40-80

    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)          { fCorrJetPt = cpt; }


    //Centrality bin setter
    Int_t                   GetFourCentBin(Double_t centbin) const;

    // event getters
    Bool_t                  DoesEventHaveD0()                  {return fd0;}
    Bool_t                  DoesEventHaveD0BgUS()              {return fd0BgUS;}
    Bool_t                  DoesEventHaveD0BgLS()              {return fd0BgLS;}

    Bool_t                  FixPoorlySmearedTracks() {fReplaceBadTracks = kTRUE;}


    TVector3                GetRandomHardTrack()       {return fHardTrack;}

  protected:
    void                    RunTracks();
    void                    RunTowers();
    void                    TestTracks();
    void                    RunJets();
    Double_t                RelativePhi(Double_t mphi, Double_t vphi) const;      // relative jet track angle
    Double_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const; // relative jet event plane angle
    void                    FillEmcTriggers();                          // EmcTrigger counter histo
    Bool_t                  DoComparison(int myarr[], int elems);
    void                    SetSumw2(); // set errors weights 

    Bool_t                  IsD0MCTrack(int trackid);


    // Topo Cut level
    int                     fTopoLevel;             // Topo cut level (0 is null, 1 is tight, 2 is less tight, 3 is lax)

    int                     fAnaCentBin;

    // switches
    Bool_t                  doPrintEventCounter;     // print event # switch
    Bool_t                  fDoEffCorr;              // efficiency correction to tracks

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

    Bool_t                  fd0;
    Bool_t                  fd0BgUS;
    Bool_t                  fd0BgLS;

    Bool_t                  fReplaceBadTracks;

    TVector3                fHardTrack;

  private:

    // funcs
    // function to calculate relative phi between 2 objects and shift between 0 and 2pi 
    //___________________________________________________________________________________________
    Double_t dPhi(Double_t phi1, Double_t phi2) {
      Double_t deltaPhi;
      deltaPhi = abs(phi1 - phi2); //TODO absolute values
      if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
      if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

      if (deltaPhi > TMath::Pi()) deltaPhi=2*(TMath::Pi()) - deltaPhi;
      return deltaPhi;   // dphi in [0, 2Pi]
    }

    // Standardise phi to 0 to 2 Pi
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

    // variables
    Int_t                   fRunNumber;

    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);

    // position object
    StEmcPosition2         *mEmcPosition;

    TH1F *hCentrality;
    TH1F *hMultiplicity;

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

    ClassDef(StThrowARandomTrack, 1)
};
#endif
