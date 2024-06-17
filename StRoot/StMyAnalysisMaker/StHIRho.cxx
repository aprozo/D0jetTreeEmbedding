// $Id$
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fOutRhoName".Apppend("_Scaled").
//

#include "StHIRho.h"

// ROOT includes
#include <TClonesArray.h>
#include <TMath.h>
#include "TH2.h"
#include "TH2F.h"
#include "TVector3.h"

// jet-framework includes
#include "StJet.h"
#include "StRhoParameter.h"
// #include "StJetMakerTask.h"
#include "StHIRecoJets.h"
#include "StCentMaker.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

// STAR centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"


// ROOT classes
class TH2;
class TH2F;

ClassImp(StHIRho)

//________________________________________________________________________
StHIRho::StHIRho() : StRhoBase("")
{
  fNExclLeadJets = 0;
  fJets = 0x0;
  // fJetsArr = 0x0;
  // rhovalue = 0x0;
  mOutName = ""; 
  fJetMakerName = "";
  fRhoMakerName = "";
}

//________________________________________________________________________
StHIRho::StHIRho(const char *name, Bool_t histo, const char *outName, const char *jetMakerName) :
  StRhoBase(name, histo, jetMakerName)
{
  fNExclLeadJets = 0;
  fJets = 0x0;
  // fJetsArr = 0x0;
  // rhovalue = 0x0;
  mBaseMaker = 0x0;
  mOutName = outName;
  fJetMakerName = jetMakerName;
  fRhoMakerName = name;
  // Constructor.
  if (!name) return;
  SetName(name);
}

//________________________________________________________________________
StHIRho::~StHIRho()
{ /*  */

}

//________________________________________________________________________
Int_t StHIRho::Init()
{
  // nothing done - base class should take care of that
  // this in effect inherits from StJetFrameworkPicoBase - check it out!
  StRhoBase::Init();


  // Create user objects.
  fJets = new TClonesArray("StJet");
  //fJets->SetName(fJetsName);   

  fJetsArr.clear();
  rhovalue.clear();   

  return kStOk;
}

//________________________________________________________________________
Int_t StHIRho::Finish() {
  //  Write histos to file and close it.
  return kStOK;
}


//________________________________________________________________________
void StHIRho::Clear(Option_t *opt) {
//  StHIRhoBase::Clear();
//  fJets->Clear();
}

//
// Function that runs the analysis for each event
//________________________________________________________________________
Int_t StHIRho::Make() 
{
  fJetsArr.clear();
  rhovalue.clear();   

  // get PicoDstMaker
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();

  // get run number, check bad runs list if desired (kFALSE if bad)
  int fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  // cut event on max track pt > 30.0 GeV

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get vertex 3 vector and declare variables
  TVector3 mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // z-vertex cut - per the Aj analysis (-40, 40) for reference

  // get JetMaker
  mHIRecoJets = static_cast<StHIRecoJets*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!mHIRecoJets) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // if we have JetMaker, get jet collection associated with it
  if(mHIRecoJets) {
    fJetsArr =  mHIRecoJets->GetJets();
    //fJets->SetName("BGJetsRho");  // name is set by Maker who created it
  }
  if(fJetsArr.size()==0) return kStOk;

  // initialize Rho and scaled Rho
  fOutRho->SetVal(0);
  if(fOutRhoScaled) fOutRhoScaled->SetVal(0);

  // cout << "Number of Jet Collections = " << fJetsArr.size() << endl;
  for (int jetcollection = 0; jetcollection < fJetsArr.size(); jetcollection++){

    // cout << "Jet Collection # " << jetcollection << endl;
    fJets->Clear();
    // get number of jets, initialize arrays
    fJets = fJetsArr[jetcollection];
    
    const Int_t Njets = fJetsArr[jetcollection]->GetEntries();

    Int_t maxJetIds[]   = {-1, -1};
    Float_t maxJetPts[] = { 0,  0};

    // exclude leading jets
    if(fNExclLeadJets > 0) {
      // loop over jets
      for(Int_t ij = 0; ij < Njets; ++ij) {
        // get jet pointer
        
        StJet *jet = static_cast<StJet*>(fJets->At(ij));
        if(!jet) { continue; } 

        //if (!AcceptJet(jet)) continue; // FIXME if wanting additional cuts
        // get some jet parameters
        double jetPt = jet->Pt();
        // some threshold cuts for tests
        if(jetPt < 0) continue;

        // get ID and pt of the leading and sub-leading jet   
        if(jetPt > maxJetPts[0]) {
  	maxJetPts[1] = maxJetPts[0];
  	maxJetIds[1] = maxJetIds[0];
  	maxJetPts[0] = jetPt;
  	maxJetIds[0] = ij;
        } else if (jetPt > maxJetPts[1]) {
  	maxJetPts[1] = jetPt;
  	maxJetIds[1] = ij;
        }
      }

      // only set to remove leading jet 
      if(fNExclLeadJets < 2) {
        maxJetIds[1] = -1;
        maxJetPts[1] = 0;
      }
    }

    static Double_t rhovec[999];
    Int_t NjetAcc = 0;

    // push all jets within selected acceptance into stack
    for(Int_t iJets = 0; iJets < Njets; ++iJets) {
      // excluding lead jets
      if(iJets == maxJetIds[0] || iJets == maxJetIds[1])
        continue;

      // get jet pointer
      StJet *jet = static_cast<StJet*>(fJets->At(iJets));
      if(!jet) { continue; } 

      // NEED TO CHECK FOR DEFAULTS - cuts are done at the jet finder level
      //if(!AcceptJet(jet)) continue; //FIXME if wanting additional cuts
      // get some get parameters
      double jetPt = jet->Pt();
      double jetArea = jet->Area();
      // some threshold cuts for tests
      if(jetPt < 0) continue;
   
      rhovec[NjetAcc] = jetPt / jetArea;
      ++NjetAcc;
    }

    // cout << "Number of accepted jets = " << NjetAcc << endl;

    // when we have accepted Jets - calculate and set rho
    if(NjetAcc > 0) {
      // find median value
      Double_t rho = TMath::Median(NjetAcc, rhovec);
      fOutRho->SetVal(rho);

      rhovalue.push_back(fOutRho->GetVal());

      // if we want scaled Rho from charged -> ch+ne
      if(fOutRhoScaled) {
        //Double_t rhoScaled = rho * GetScaleFactor(fCent); //Don't need yet, fCent is in 5% bins
        Double_t rhoScaled = rho * 1.0;
        fOutRhoScaled->SetVal(rhoScaled);
      }
    }

    // StRhoBase::FillHistograms();

  }
  return kStOk;
} 
