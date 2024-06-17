#ifndef StHIRho_H
#define StHIRho_H

// $Id$
// adapted from the AliROOT class AliAnalysisTaskRho.h

#include "StRhoBase.h"

// additional includes
#include "StMaker.h"

// ROOT classes
class TH2;
class TH2F;
class TClonesArray;

// STAR classes
class StMaker;
class StHIRecoJets;

class StHIRho : public StRhoBase {

 public:
  StHIRho();
  StHIRho(const char *name, Bool_t histo=kFALSE, const char* outName="", const char* jetMakerName="");
  //virtual ~StHIRho() {}
  virtual ~StHIRho();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *opt="");
  virtual Int_t Finish();

  std::vector<double> GetArrayofRhoValue() {return rhovalue;}

  void    SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }

  // set Jets for D0s, D0s Bg US, D0s Bg LS

 protected:
  UInt_t            fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation

  TClonesArray     *fJets;//!jet collection
  std::vector<TClonesArray *> fJetsArr;

  // set Jets for D0s, D0s Bg US, D0s Bg LS
 private:

  StHIRecoJets *mHIRecoJets;//!

  std::vector<double> rhovalue; //!

  StHIRho(const StHIRho&);             // not implemented
  StHIRho& operator=(const StHIRho&);  // not implemented
  
  ClassDef(StHIRho, 2); // Rho task
};
#endif
