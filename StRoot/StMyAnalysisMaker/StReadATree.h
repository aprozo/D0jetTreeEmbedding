#ifndef StReadATree_h
#define StReadATree_h


#include "StMaker.h"
#include "StJetFrameworkPicoBase.h"
#include <vector>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <utility>
#include "TTree.h"
#include "TBranch.h"

class StJetFrameworkPicoBase;

// ROOT classes
class TTree;
class TBranch;
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH2D;
class TH3F;
class THnSparse;
class TProfile;
class TString;
class TVector3;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoDstReader;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;
class StEmcGeom;
// class StEmcPosition;

// jet-framework classes
class StCentMaker;
class StEmcPosition2;
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxEvent = 1;
const Int_t kMaxTrack = 5000;
const Int_t kMaxEmcTrigger = 5000;
const Int_t kMaxMtdTrigger = 5000;
const Int_t kMaxBTowHit = 5000;
const Int_t kMaxBTofHit = 5000;
const Int_t kMaxMtdHit = 5000;
const Int_t kMaxBbcHit = 5000;
const Int_t kMaxEpdHit = 5000;
const Int_t kMaxFmsHit = 5000;
const Int_t kMaxEmcPidTraits = 5000;
const Int_t kMaxBTofPidTraits = 5000;
const Int_t kMaxMtdPidTraits = 5000;
const Int_t kMaxTrackCovMatrix = 5000;


class StReadATree : public StJetFrameworkPicoBase  {
    
  public:

    StReadATree(const char *name, const char* filename, const char* outName);
    virtual ~StReadATree();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    void    FillTheVectorThatReplacesThePicoDst();
    void    GetTOFMatching();
    void    DeclareHistograms();
    void    WriteHistograms();

    const char *mFilename;

    TFile *f;

    TTree *fVariablesAssignedToTree;

    // Declaration of leaf types
    Int_t           Event_;
    Int_t           Event_mRunId[kMaxEvent];   //[Event_]
    Int_t           Event_mEventId[kMaxEvent];   //[Event_]
    UShort_t        Event_mFillId[kMaxEvent];   //[Event_]
    Float_t         Event_mBField[kMaxEvent];   //[Event_]
    Int_t           Event_mTime[kMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexX[kMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexY[kMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexZ[kMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexErrorX[kMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexErrorY[kMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexErrorZ[kMaxEvent];   //[Event_]
    Float_t         Event_mRanking[kMaxEvent];   //[Event_]
    UShort_t        Event_mNBEMCMatch[kMaxEvent];   //[Event_]
    UShort_t        Event_mNBTOFMatch[kMaxEvent];   //[Event_]
    vector<unsigned int> Event_mTriggerIds[kMaxEvent];
    UShort_t        Event_mRefMultFtpcEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultFtpcWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultNeg[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultPos[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult2NegEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult2PosEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult2NegWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult2PosWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult3NegEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult3PosEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult3NegWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult3PosWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult4NegEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult4PosEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult4NegWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult4PosWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultHalfNegEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultHalfPosEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultHalfNegWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultHalfPosWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mGRefMult[kMaxEvent];   //[Event_]
    UShort_t        Event_mNumberOfGlobalTracks[kMaxEvent];   //[Event_]
    UShort_t        Event_mbTofTrayMultiplicity[kMaxEvent];   //[Event_]
    UShort_t        Event_mNHitsHFT[kMaxEvent][4];   //[Event_]
    UChar_t         Event_mNVpdHitsEast[kMaxEvent];   //[Event_]
    UChar_t         Event_mNVpdHitsWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mNTofT0[kMaxEvent];   //[Event_]
    Float_t         Event_mVzVpd[kMaxEvent];   //[Event_]
    UInt_t          Event_mZDCx[kMaxEvent];   //[Event_]
    UInt_t          Event_mBBCx[kMaxEvent];   //[Event_]
    Float_t         Event_mBackgroundRate[kMaxEvent];   //[Event_]
    Float_t         Event_mBbcBlueBackgroundRate[kMaxEvent];   //[Event_]
    Float_t         Event_mBbcYellowBackgroundRate[kMaxEvent];   //[Event_]
    Float_t         Event_mBbcEastRate[kMaxEvent];   //[Event_]
    Float_t         Event_mBbcWestRate[kMaxEvent];   //[Event_]
    Float_t         Event_mZdcEastRate[kMaxEvent];   //[Event_]
    Float_t         Event_mZdcWestRate[kMaxEvent];   //[Event_]
    UShort_t        Event_mZdcSumAdcEast[kMaxEvent];   //[Event_]
    UShort_t        Event_mZdcSumAdcWest[kMaxEvent];   //[Event_]
    UShort_t        Event_mZdcSmdEastHorizontal[kMaxEvent][8];   //[Event_]
    UShort_t        Event_mZdcSmdEastVertical[kMaxEvent][8];   //[Event_]
    UShort_t        Event_mZdcSmdWestHorizontal[kMaxEvent][8];   //[Event_]
    UShort_t        Event_mZdcSmdWestVertical[kMaxEvent][8];   //[Event_]
    UShort_t        Event_mBbcAdcEast[kMaxEvent][24];   //[Event_]
    UShort_t        Event_mBbcAdcWest[kMaxEvent][24];   //[Event_]
    UChar_t         Event_mHighTowerThreshold[kMaxEvent][4];   //[Event_]
    UChar_t         Event_mJetPatchThreshold[kMaxEvent][4];   //[Event_]
    Int_t           Track_;
    UShort_t        Track_mId[kMaxTrack];   //[Track_]
    UShort_t        Track_mChi2[kMaxTrack];   //[Track_]
    Float_t         Track_mPMomentumX[kMaxTrack];   //[Track_]
    Float_t         Track_mPMomentumY[kMaxTrack];   //[Track_]
    Float_t         Track_mPMomentumZ[kMaxTrack];   //[Track_]
    Float_t         Track_mGMomentumX[kMaxTrack];   //[Track_]
    Float_t         Track_mGMomentumY[kMaxTrack];   //[Track_]
    Float_t         Track_mGMomentumZ[kMaxTrack];   //[Track_]
    Float_t         Track_mOriginX[kMaxTrack];   //[Track_]
    Float_t         Track_mOriginY[kMaxTrack];   //[Track_]
    Float_t         Track_mOriginZ[kMaxTrack];   //[Track_]
    Float16_t       Track_mDedx[kMaxTrack];   //[Track_]
    Float16_t       Track_mDedxError[kMaxTrack];   //[Track_]
    Char_t          Track_mNHitsFit[kMaxTrack];   //[Track_]
    UChar_t         Track_mNHitsMax[kMaxTrack];   //[Track_]
    UChar_t         Track_mNHitsDedx[kMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaPion[kMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaKaon[kMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaProton[kMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaElectron[kMaxTrack];   //[Track_]
    UInt_t          Track_mTopologyMap[kMaxTrack][2];   //[Track_]
    Short_t         Track_mBEmcPidTraitsIndex[kMaxTrack];   //[Track_]
    Short_t         Track_mBTofPidTraitsIndex[kMaxTrack];   //[Track_]
    Short_t         Track_mMtdPidTraitsIndex[kMaxTrack];   //[Track_]
    Int_t           EmcTrigger_;
    UChar_t         EmcTrigger_mFlag[kMaxEmcTrigger];   //[EmcTrigger_]
    UShort_t        EmcTrigger_mId[kMaxEmcTrigger];   //[EmcTrigger_]
    UShort_t        EmcTrigger_mAdc[kMaxEmcTrigger];   //[EmcTrigger_]
    Int_t           MtdTrigger_;
    UShort_t        MtdTrigger_mVpdTacSum[kMaxMtdTrigger];   //[MtdTrigger_]
    UInt_t          MtdTrigger_mTHUBtime[kMaxMtdTrigger][2];   //[MtdTrigger_]
    UShort_t        MtdTrigger_mQTtacSum[kMaxMtdTrigger][8][8];   //[MtdTrigger_]
    UShort_t        MtdTrigger_mMT101Tac[kMaxMtdTrigger][8][2];   //[MtdTrigger_]
    UChar_t         MtdTrigger_mMT101Id[kMaxMtdTrigger][8][2];   //[MtdTrigger_]
    UInt_t          MtdTrigger_mTF201TriggerBit[kMaxMtdTrigger];   //[MtdTrigger_]
    Char_t          MtdTrigger_mShouldHaveRejectEvent[kMaxMtdTrigger];   //[MtdTrigger_]
    Int_t           BTowHit_;
    UShort_t        BTowHit_mAdc[kMaxBTowHit];   //[BTowHit_]
    Short_t         BTowHit_mE[kMaxBTowHit];   //[BTowHit_]
    Int_t           BTofHit_;
    Short_t         BTofHit_mId[kMaxBTofHit];   //[BTofHit_]
    Int_t           MtdHit_;
    Short_t         MtdHit_mgChannel[kMaxMtdHit];   //[MtdHit_]
    UChar_t         MtdHit_mTriggerFlag[kMaxMtdHit];   //[MtdHit_]
    Float_t         MtdHit_mLeadingEdgeTime_first[kMaxMtdHit];   //[MtdHit_]
    Float_t         MtdHit_mLeadingEdgeTime_second[kMaxMtdHit];   //[MtdHit_]
    Float_t         MtdHit_mTrailingEdgeTime_first[kMaxMtdHit];   //[MtdHit_]
    Float_t         MtdHit_mTrailingEdgeTime_second[kMaxMtdHit];   //[MtdHit_]
    Int_t           BbcHit_;
    Short_t         BbcHit_mId[kMaxBbcHit];   //[BbcHit_]
    Int_t           BbcHit_mQTdata[kMaxBbcHit];   //[BbcHit_]
    Int_t           EpdHit_;
    Short_t         EpdHit_mId[kMaxEpdHit];   //[EpdHit_]
    Int_t           EpdHit_mQTdata[kMaxEpdHit];   //[EpdHit_]
    Float_t         EpdHit_mnMIP[kMaxEpdHit];   //[EpdHit_]
    Int_t           FmsHit_;
    UShort_t        FmsHit_mChannelDetectorId[kMaxFmsHit];   //[FmsHit_]
    UShort_t        FmsHit_mAdc[kMaxFmsHit];   //[FmsHit_]
    Int_t           EmcPidTraits_;
    Short_t         EmcPidTraits_mTrackIndex[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcId[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcAdc0[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcE0[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcE[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcZDist[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcPhiDist[kMaxEmcPidTraits];   //[EmcPidTraits_]
    UChar_t         EmcPidTraits_mBemcSmdNEta[kMaxEmcPidTraits];   //[EmcPidTraits_]
    UChar_t         EmcPidTraits_mBemcSmdNPhi[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowId[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Char_t          EmcPidTraits_mBtowId23[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowE[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowE2[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowE3[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowEtaDist[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowPhiDist[kMaxEmcPidTraits];   //[EmcPidTraits_]
    Int_t           BTofPidTraits_;
    Short_t         BTofPidTraits_mTrackIndex[kMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofCellId[kMaxBTofPidTraits];   //[BTofPidTraits_]
    UChar_t         BTofPidTraits_mBTofMatchFlag[kMaxBTofPidTraits];   //[BTofPidTraits_]
    Float_t         BTofPidTraits_mBTof[kMaxBTofPidTraits];   //[BTofPidTraits_]
    UShort_t        BTofPidTraits_mBTofBeta[kMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofYLocal[kMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofZLocal[kMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofHitPosX[kMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofHitPosY[kMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofHitPosZ[kMaxBTofPidTraits];   //[BTofPidTraits_]
    Int_t           MtdPidTraits_;
    Short_t         MtdPidTraits_mTrackIndex[kMaxMtdPidTraits];   //[MtdPidTraits_]
    Short_t         MtdPidTraits_mMtdHitIndex[kMaxMtdPidTraits];   //[MtdPidTraits_]
    Char_t          MtdPidTraits_mMatchFlag[kMaxMtdPidTraits];   //[MtdPidTraits_]
    Short_t         MtdPidTraits_mDeltaY[kMaxMtdPidTraits];   //[MtdPidTraits_]
    Short_t         MtdPidTraits_mDeltaZ[kMaxMtdPidTraits];   //[MtdPidTraits_]
    Float_t         MtdPidTraits_mDeltaTimeOfFlight[kMaxMtdPidTraits];   //[MtdPidTraits_]
    UShort_t        MtdPidTraits_mBeta[kMaxMtdPidTraits];   //[MtdPidTraits_]
    Short_t         MtdPidTraits_mMtdHitChan[kMaxMtdPidTraits];   //[MtdPidTraits_]
    Int_t           TrackCovMatrix_;
    Float16_t       TrackCovMatrix_mImp[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mZ[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mPsi[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mPti[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mTan[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mCurv[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mSigma[kMaxTrackCovMatrix][5];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mCorr[kMaxTrackCovMatrix][10];   //[TrackCovMatrix_]

    //setters

    void         SetMinJetTrackPt(Double_t min)             { fMinJetTrackPt = min; }
    void         SetMaxJetTrackPt(Double_t max)             { fMaxJetTrackPt = max; }
    void         SetJetTrackEtaRange(Double_t etmi, Double_t etma) { fJetTrackEtaMin = etmi; fJetTrackEtaMax = etma; }
    void         SetJetTrackPhiRange(Double_t ptmi, Double_t ptma) { fJetTrackPhiMax = ptmi; fJetTrackPhiMax = ptma; }
    void         SetJetTrackDCAcut(Double_t d)              { fJetTrackDCAcut   = d     ; }
    void         SetJetTracknHitsFit(Double_t h)            { fJetTracknHitsFit = h     ; }
    void         SetJetTracknHitsRatio(Double_t r)          { fJetTracknHitsRatio = r   ; }
    void         SetMinJetTowerE(Double_t min)              { mTowerEnergyMin = min; }
    void         SetJetTowerERange(Double_t enmi, Double_t enmx) { fJetTowerEMin = enmi; fJetTowerEMax = enmx; }
    void         SetJetTowerEtaRange(Double_t temi, Double_t temx) { fJetTowerEtaMin = temi; fJetTowerEtaMax = temx; }
    void         SetJetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fJetTowerPhiMin = tpmi; fJetTowerPhiMax = tpmx; }
    void         SetMinJetClusPt(Double_t min)              { fMinJetClusPt  = min; }
    void         SetMinJetClusE(Double_t min)               { fMinJetClusE   = min; }

    std::vector<std::vector<double>> EventTreeToVector() {return fEventVectorToFill;}
    std::vector<std::vector<double>> TrackTreeToVector() {return fTrackVectorToFill;}
    std::vector<std::vector<double>> TowerTreeToVector() {return fTowerVectorToFill;}

  protected:
    std::vector<std::vector<double>> fEventVectorToFill;
    std::vector<std::vector<double>> fTrackVectorToFill;
    std::vector<std::vector<double>> fTowerVectorToFill;

    // track attributes
    Double_t               fMinJetTrackPt;          // min jet track transverse momentum cut
    Double_t               fMaxJetTrackPt;          // max jet track transverse momentum cut
    Double_t               fMinJetClusPt;           // min jet cluster transverse momentum cut
    Double_t               fMinJetClusE;            // min jet cluster energy cut
    Double_t               fMinJetTowerE;           // min jet tower energy cut - not used (use mTowerEnergyMin)
    Double_t               fTrackEtaMin;            // min track eta cut
    Double_t               fTrackEtaMax;            // max track eta cut
    Double_t               fTrackPhiMin;            // min track phi cut
    Double_t               fTrackPhiMax;            // max track phi cut
    Double_t               fJetTrackEtaMin;         // min jet track eta cut
    Double_t               fJetTrackEtaMax;         // max jet track eta cut
    Double_t               fJetTrackPhiMin;         // min jet track phi cut
    Double_t               fJetTrackPhiMax;         // max jet track phi cut
    Double_t               fJetTrackDCAcut;         // max jet track dca cut
    Int_t                  fJetTracknHitsFit;       // requirement for track hits
    Double_t               fJetTracknHitsRatio;     // requirement for nHitsFit / nHitsMax
    Double_t               fTrackEfficiency;        // artificial tracking inefficiency (0...1) - this is NOT used

    // tower attributes
    Double_t               fJetTowerEMin;           // min jet tower energy cut
    Double_t               fJetTowerEMax;           // max jet tower energy cut
    Double_t               fJetTowerEtaMin;         // min jet tower eta cut
    Double_t               fJetTowerEtaMax;         // max jet tower eta cut
    Double_t               fJetTowerPhiMin;         // min jet tower phi cut
    Double_t               fJetTowerPhiMax;         // max jet tower phi cut
    Double_t               mTowerEnergyMin;         // min jet tower energy cut
    Float_t                mHadronicCorrFrac;       // hadronic correction fraction from 0.0 to 1.0

    int     GetMatchedBtowID(StPicoTrack *trk);

    vector<double> ReadNthLine(TString filename, int N)
    {
      std::ifstream in(filename.Data());

      std::string s;
      //for performance
      s.reserve(100);    

      //skip N lines
      for(int i = 0; i < N; ++i)
         std::getline(in, s);

      std::getline(in,s);

      std::istringstream vertexforthisevent(s); 

      int RunID, EventID;
      double vx, vy, vz;

      vertexforthisevent >> RunID >> EventID >> vx >> vy >> vz;

      std::vector<double> v;
      v.clear();

      v.push_back(RunID);
      v.push_back(EventID);
      v.push_back(vx);
      v.push_back(vy);
      v.push_back(vz);

      return v; 
    }

    
  private:

    StPicoDstMaker *mPicoDstMaker;
    StPicoDst *mPicoDst;
    StPicoEvent *mPicoEvent;

    // position object
    StEmcGeom              *mBemcGeom;
    StEmcPosition2         *mEmcPosition;
    // StCentMaker            *mCentMaker;
    // base class pointer object
    StJetFrameworkPicoBase *mBaseMaker;

    Double_t               mTowerMatchTrkIndex[4800][7];
    Int_t                  mTowerStatusArr[4800];

    TString                fAnalysisMakerName;
    TString                vertexfilename;
    TString                mOutName;

    TH1D *hdEta;
    TH1D *hdPhi;
    TH2D *hdEtadPhi;
    TH1D *hdTowerIndex;
    TH1D *hPercentTracks;

    ClassDef(StReadATree, 1)
};


#endif
