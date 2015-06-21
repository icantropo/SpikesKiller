#ifndef SkimElectronStudies_h
#define SkimElectronStudies_h

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

#include <iostream>
#include <memory>

#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"

//
// class declaration
//

class SkimElectronStudies : public edm::EDFilter
{
 public:
  
  explicit SkimElectronStudies(const edm::ParameterSet&);
  ~SkimElectronStudies();
  
  virtual bool filter(edm::Event&, const edm::EventSetup& );
      
 private:

  edm::InputTag electronCollection_;
  //edm::InputTag clusterCollection_;
  //edm::InputTag gsftrackCollection_;
  
	double ele_ptLow_;
	double ele_ptHigh_;
	double nEle_ptLowMIN_;
	double nEle_ptHighMIN_;
	double sc_EtLow_;
	double sc_EtHigh_;
	double nSC_EtLowMIN_;
	double nSC_EtHighMIN_;
	
};

#endif
