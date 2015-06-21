// -*- C++ -*-
//
// Package:    EmulSpikesFilter
// Class:      EmulSpikesFilter
// 
/**\class EmulSpikesFilter EmulSpikesFilter.cc EmulSpikesFilter/EmulSpikesFilter/plugins/EmulSpikesFilter.cc

 Description: Search RecHit collection for presence of the spikes (sev=3 or sev=4).
 Keep only events with the spikes.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Iurii Antropov
//         Created:  Fri, 13 Feb 2015 17:39:12 GMT
// $Id$
//


// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>

//includes root
#include "TTree.h"

//// CMSSW
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

//file service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//ECAL includes
#include "CondFormats/EcalObjects/interface/EcalTPGTowerStatus.h" // for the EcalTPGTowerStatusMapIterator
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

// TPG (Nadir study)
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalTrigPrimCompactColl.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveSample.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"

//Clusters
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

// For H/E - Iso on SC
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"

//RECHITS
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

// severity for spikes
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

using namespace std;
using namespace edm;

//
// class declaration
//

class EmulSpikesFilter : public edm::EDFilter {
   public:
      explicit EmulSpikesFilter(const edm::ParameterSet&);
      ~EmulSpikesFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


      // ----------member data ---------------------------
      edm::InputTag EcalRecHitCollectionEB_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EmulSpikesFilter::EmulSpikesFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
	EcalRecHitCollectionEB_ = (iConfig.getParameter<edm::InputTag>("EcalRecHitCollectionEB"));
}


EmulSpikesFilter::~EmulSpikesFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EmulSpikesFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   	// Get severity
   	unsigned long long cacheSevLevel = 0;
   	edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
   	if (cacheSevLevel != iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()) {
   		cacheSevLevel =	iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
   		iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
   	}
   	const EcalSeverityLevelAlgo* sl = sevLevel.product();
   	// Get EB RecHits
   	edm::Handle<EcalRecHitCollection> rechitsEB;
   	EcalRecHitCollection::const_iterator rechitItr;

   	uint32_t sev = 0;
	EBDetId id;

	// Loop over RecHit collection only in the Barrel
	if (iEvent.getByLabel(EcalRecHitCollectionEB_, rechitsEB)) {
   		for (rechitItr = rechitsEB->begin(); rechitItr != rechitsEB->end();	++rechitItr) {
   			id = rechitItr->id();
   			//sev = EcalSeverityLevelAlgo::severityLevel( id, *rechitsEB, *chaStatus );
   			//sev = EcalSeverityLevelAlgo::severityLevel(id,*rechitsEB,*chaStatus, 5.,
   			//EcalSeverityLevelAlgo::kSwissCross,0.95) ;
   			sev = 0;
   			sev = sl->severityLevel(id, *rechitsEB);
   			//if(sev>0) cout << "got sev=" << sev << endl;
   			//if (sev >= 1) {
   			if ((sev == 3) || (sev == 4)) {  // if it's a spike
//   				cout << "Spike Found:" << endl;
//   				cout << id << endl;
   				return true;
   			} // end if sev=3 || sev=4
   		}	//loop RecHits
   	} // if Barrel RecHit

   	return false;

}

// ------------ method called once each job just before starting event loop  ------------
void 
EmulSpikesFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EmulSpikesFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
EmulSpikesFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
EmulSpikesFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
EmulSpikesFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
EmulSpikesFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EmulSpikesFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(EmulSpikesFilter);
