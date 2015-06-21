#ifndef SkimLeptonStudies_h
#define SkimLeptonStudies_h

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include <iostream>
#include <memory>

#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"

//
// class declaration
//

class SkimLeptonStudies : public edm::EDFilter
{
 public:
  
  explicit SkimLeptonStudies(const edm::ParameterSet&);
  ~SkimLeptonStudies();
  
  virtual bool filter(edm::Event&, const edm::EventSetup& );
      
 private:

  edm::InputTag electronCollection_;
  edm::InputTag muonCollection_;
  
	bool isEleID_;
	bool isMuonID_;
	double lep_ptLow_;
	double lep_ptHigh_;
	int nLep_ptLow_;
	int nLep_ptHigh_;
	/*
	double ele_ptLow_;
	double ele_ptHigh_;
	int nEle_ptLowMIN_;
	int nEle_ptHighMIN_;
	double muon_ptLow_;
	double muon_ptHigh_;
	int nMuon_ptLowMIN_;
	int nMuon_ptHighMIN_;
	*/
};

#endif
