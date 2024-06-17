#ifndef STD0RHO_H
#define STD0RHO_H

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
class StD0EventsJetMaker;

class StD0Rho : public StRhoBase {

 public:
  StD0Rho();
  StD0Rho(const char *name, Bool_t histo=kFALSE, const char* outName="", const char* jetMakerName="");
  //virtual ~StD0Rho() {}
  virtual ~StD0Rho();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();

  std::vector<double> GetArrayofRhoValue() {return rhovalue;}

  void    SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }

  // set Jets for D0s, D0s Bg US, D0s Bg LS
  void         SetD0Kind(Int_t kind)                      { fD0Kind           = kind  ; }

 protected:
  UInt_t            fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation

  TClonesArray     *fJets;//!jet collection
  std::vector<TClonesArray *> fJetsArr;

  // set Jets for D0s, D0s Bg US, D0s Bg LS
  Double_t               fD0Kind;                 // D0s, Bg D0s US, LS

 private:
  TH2F             *fHistMultvsRho;//!

  StD0EventsJetMaker *D0JetMaker;//!

  std::vector<double> rhovalue; //!

  StD0Rho(const StD0Rho&);             // not implemented
  StD0Rho& operator=(const StD0Rho&);  // not implemented
  
  ClassDef(StD0Rho, 2); // Rho task
};
#endif
