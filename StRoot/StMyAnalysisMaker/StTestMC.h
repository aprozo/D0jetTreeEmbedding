#ifndef STTESTMC_h
#define STTESTMC_h

#include "StJetFrameworkPicoBase.h"
class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3F;
class THnSparse;
class TProfile;
class TString;
class TVector3;

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

class StTestMC : public StJetFrameworkPicoBase {
  public:

    StTestMC(const char *name, StPicoDstMaker *picoMaker, const char *outName, const char *jetMakerName, const char *rhoMakerName);
    virtual ~StTestMC();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();

    void      IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e);
    Int_t     IsWhatParticle(StPicoTrack *trk);
    Bool_t    IsD0MCTrack(int trackid);
    // void    WriteVertexAndEventInfo();

    bool      IsTpcPion(StPicoTrack *trk);
    bool      IsTpcKaon(StPicoTrack *trk);
    bool      IsPion(StPicoTrack *trk);
    bool      IsKaon(StPicoTrack *trk);

    float     GetTofBeta(StPicoTrack *trk);

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp;}
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }

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

  protected:
    void                    RunTracks();
    void                    RunTowers();
    void                    RunJets();
    Double_t                RelativePhi(Double_t mphi, Double_t vphi) const;      // relative jet track angle
    Double_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const; // relative jet event plane angle
    void                    FillEmcTriggers();                          // EmcTrigger counter histo
    Bool_t                  DoComparison(int myarr[], int elems);
    void                    SetSumw2(); // set errors weights 

    Bool_t                  IsD0MCVertex(vector<int> arr[], int idvx);

    // switches
    Bool_t                  doPrintEventCounter;     // print event # switch
    Bool_t                  fDoEffCorr;              // efficiency correction to tracks

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

  private:

    const double Mpion = 0.139570;
    const double Mkaon = 0.493677;
    const double Mproton = 0.938272;
    // variables
    Int_t                   fRunNumber;

    int D0Counter = 0;

    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);

    // position object
    StEmcPosition2         *mEmcPosition;
    StPicoTrack            *SmearMom(TVector3 , int);
   
    // histos
    TH1F *hCentrality;//!
    TH1F *hMultiplicity;//!
 
    // jet histos
    TH1F *hJetPt;//!
    TH1F *hJetCorrPt;//!

    TH1F *hMCD0Eta;

    TH1F *hD0MCMass;

    //____________ TEST Code #FIXME _________________
    TH2F *hTrackPt;//!
    TH2F *hTrackEt;//!

    TH1F *hEventZVertex;
    TH1F *hDCA;

    TH1F *hTPCPt;
    TH1F *hTPCEta;
    TH1F *hTPCPhi;

    TH1F *hHFTPt;
    TH1F *hHFTEta;
    TH1F *hHFTPhi;

    TH1F *hMCPt;
    TH1F *hMCEta;
    TH1F *hMCPhi;

    TH1F *hGenPt;
    TH1F *hGenEta;
    TH1F *hGenPhi;

    TH1F *hMatchedRecoTrack;
    TH1F *hRecoTrack;

    ///////////////////////// D0 Kaon Pion Histograms ///////////////////////////////////////////////

    TH1F *hMCD0KaonPt;
    TH1F *hMCD0PionPt;

    TH1F *hRecoD0KaonPt;
    TH1F *hRecoD0PionPt;

    TH2F *hRecoVMCD0KaonPt;
    TH2F *hRecoVMCD0PionPt;

    TH1F *hMCD0KaonPhi;
    TH1F *hMCD0PionPhi;

    TH1F *hRecoD0KaonPhi;
    TH1F *hRecoD0PionPhi;

    TH2F *hRecoVMCD0KaonPhi;
    TH2F *hRecoVMCD0PionPhi;

    TH1F *hMCD0KaonEta;
    TH1F *hMCD0PionEta;

    TH1F *hRecoD0KaonEta;
    TH1F *hRecoD0PionEta;

    TH2F *hRecoVMCD0KaonEta;
    TH2F *hRecoVMCD0PionEta;

    ////////////////////////// All Kaon Pion Histograms //////////////////////////////////////////

    TH1F *hMCKaonPt;
    TH1F *hMCPionPt;

    TH1F *hRecoKaonPt;
    TH1F *hRecoPionPt;

    TH2F *hRecoVMCKaonPt;
    TH2F *hRecoVMCPionPt;

    TH1F *hMCKaonPhi;
    TH1F *hMCPionPhi;

    TH1F *hRecoKaonPhi;
    TH1F *hRecoPionPhi;

    TH2F *hRecoVMCKaonPhi;
    TH2F *hRecoVMCPionPhi;

    TH1F *hMCKaonEta;
    TH1F *hMCPionEta;

    TH1F *hRecoKaonEta;
    TH1F *hRecoPionEta;

    TH2F *hRecoVMCKaonEta;
    TH2F *hRecoVMCPionEta;

    ////////////////////////// All Kaon Pion (TPC) PID Histograms //////////////////////////////////////////

    TH1F *hRecoTPCKaonPt;
    TH1F *hRecoTPCPionPt;

    TH2F *hRecoVMCTPCKaonPt;
    TH2F *hRecoVMCTPCPionPt;

    TH1F *hRecoTPCKaonPhi;
    TH1F *hRecoTPCPionPhi;

    TH2F *hRecoVMCTPCKaonPhi;
    TH2F *hRecoVMCTPCPionPhi;

    TH1F *hRecoTPCKaonEta;
    TH1F *hRecoTPCPionEta;

    TH2F *hRecoVMCTPCKaonEta;
    TH2F *hRecoVMCTPCPionEta;

    ////////////////////////// All Kaon Pion (TPC + TOF) PID Histograms //////////////////////////////////////////

    TH1F *hRecoTOFKaonPt;
    TH1F *hRecoTOFPionPt;

    TH2F *hRecoVMCTOFKaonPt;
    TH2F *hRecoVMCTOFPionPt;

    TH1F *hRecoTOFKaonPhi;
    TH1F *hRecoTOFPionPhi;

    TH2F *hRecoVMCTOFKaonPhi;
    TH2F *hRecoVMCTOFPionPhi;

    TH1F *hRecoTOFKaonEta;
    TH1F *hRecoTOFPionEta;

    TH2F *hRecoVMCTOFKaonEta;
    TH2F *hRecoVMCTOFPionEta;

    ////////////////////////// All Kaon Pion (TPC + HFT) PID Histograms //////////////////////////////////////////

    TH1F *hRecoHFTKaonPt;
    TH1F *hRecoHFTPionPt;

    TH2F *hRecoVMCHFTKaonPt;
    TH2F *hRecoVMCHFTPionPt;

    TH1F *hRecoHFTKaonPhi;
    TH1F *hRecoHFTPionPhi;

    TH2F *hRecoVMCHFTKaonPhi;
    TH2F *hRecoVMCHFTPionPhi;

    TH1F *hRecoHFTKaonEta;
    TH1F *hRecoHFTPionEta;

    TH2F *hRecoVMCHFTKaonEta;
    TH2F *hRecoVMCHFTPionEta;

    TH1F *hD0Mass;

    TH1F *hD0Daughters;
    TH1F *hD0DecayLength;

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


    TF1* fKaonMomResolution;
    TF1* fPionMomResolution;
    TF1* fProtonMomResolution;

    // std::ofstream vertexfile;

    ClassDef(StTestMC, 1)
};
#endif
