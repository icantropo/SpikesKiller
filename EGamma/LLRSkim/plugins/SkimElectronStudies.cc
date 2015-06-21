#include "EGamma/LLRSkim/plugins/SkimElectronStudies.h"

#include "TMath.h"

//! ctor
SkimElectronStudies::SkimElectronStudies(const edm::ParameterSet& iConfig):
electronCollection_               (iConfig.getParameter<edm::InputTag>("electronCollection")),
ele_ptLow_                        (iConfig.getParameter<double>("ele_ptLow")),
ele_ptHigh_                       (iConfig.getParameter<double>("ele_ptHigh")),
nEle_ptLowMIN_                    (iConfig.getParameter<int>("nEle_ptLowMIN")),
nEle_ptHighMIN_                   (iConfig.getParameter<int>("nEle_ptHighMIN")),
sc_EtLow_                        (iConfig.getParameter<double>("sc_EtLow")),
sc_EtHigh_                       (iConfig.getParameter<double>("sc_EtHigh")),
nSC_EtLowMIN_                    (iConfig.getParameter<int>("nSC_EtLowMIN")),
nSC_EtHighMIN_                   (iConfig.getParameter<int>("nSC_EtHighMIN"))
{}

//! dtor
SkimElectronStudies::~SkimElectronStudies()
{}


bool SkimElectronStudies::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	///int runId = iEvent.id().run();
	///int lumiId = iEvent.luminosityBlock();
	///int eventId = iEvent.id().event();
	
	//GSF Electron	
	edm::Handle<reco::GsfElectronCollection> EleHandle;
	iEvent.getByLabel(electronCollection_.label(), EleHandle);
	
	//EB SuperCluster
	edm::Handle<reco::SuperClusterCollection> sc_coll_EB;
	iEvent.getByLabel(edm::InputTag("correctedHybridSuperClusters"), sc_coll_EB);
	
	//EE SuperCluster
	edm::Handle<reco::SuperClusterCollection> sc_coll_EE;
	iEvent.getByLabel(edm::InputTag("correctedMulti5x5SuperClustersWithPreshower"), sc_coll_EE);
	
	// loop on electrons
	int nEle_ptLow = 0;
	int nEle_ptHigh = 0;
	for(unsigned int eleIt = 0; eleIt < EleHandle->size(); ++eleIt)
	{
		
		reco::GsfElectronRef eleRef(EleHandle, eleIt);
		
		// count electrons
		if( eleRef->pt() > ele_ptLow_) ++nEle_ptLow;
		if( eleRef->pt() > ele_ptHigh_) ++nEle_ptHigh;
		
	} // end loop on electrons
	
	// loop on EB superclusters
	int nSC_EtLow = 0;
	int nSC_EtHigh = 0;
	for( reco::SuperClusterCollection::const_iterator isc=sc_coll_EB->begin(); 
		isc!=sc_coll_EB->end(); 
		++isc) {
		
		double R  = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y() +isc->z()*isc->z());
		double Rt = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y());
		
		
		if( isc->energy()*(Rt/R) > sc_EtLow_) ++nSC_EtLow;
		if( isc->energy()*(Rt/R) > sc_EtHigh_) ++nSC_EtHigh;
	} // end loop on EB superclusters
	
	// loop on EE superclusters
	for( reco::SuperClusterCollection::const_iterator isc=sc_coll_EE->begin(); 
		isc!=sc_coll_EE->end(); 
		++isc) {
		
		double R  = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y() +isc->z()*isc->z());
		double Rt = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y());
		
		
		if( isc->energy()*(Rt/R) > sc_EtLow_) ++nSC_EtLow;
		if( isc->energy()*(Rt/R) > sc_EtHigh_) ++nSC_EtHigh;
	} // end loop on EE superclusters
	
	
	// print out results
	if( nEle_ptLow >= nEle_ptLowMIN_ ) {
		//std::cout << ">>>SkimElectronStudies::run=" << runId << "::lumi=" 
		//<< lumiId << "::eventId=" 
		//<< eventId << "::Found at least "
		//<< nEle_ptLowMIN_ << " electrons with pt > " << ele_ptLow_ << std::endl; 
		return true;
	}
	if( nEle_ptHigh >= nEle_ptHighMIN_ ) {
		//std::cout << ">>>SkimElectronStudies::run=" << runId << "::lumi=" 
		//<< lumiId << "::eventId=" 
		//<< eventId << "::Found at least "
		//<< nEle_ptHighMIN_ << " electrons with pt > " << ele_ptHigh_ << std::endl; 
		return true;
	}
	if( nSC_EtLow >= nSC_EtLowMIN_ ) {
		//std::cout << ">>>SkimElectronStudies::run=" << runId << "::lumi=" 
		//<< lumiId << "::eventId=" 
		//<< eventId << "::Found at least "
		//<< nSC_EtLowMIN_ << " supercluster with Et > " << sc_EtLow_ << std::endl; 
		return true;
	}
	if( nSC_EtHigh >= nSC_EtHighMIN_ ) {
		//std::cout << ">>>SkimElectronStudies::run=" << runId << "::lumi=" 
		//<< lumiId << "::eventId=" 
		//<< eventId << "::Found at least "
		//<< nSC_EtHighMIN_ << " supercluster with Et > " << sc_EtHigh_ << std::endl; 
		return true;
	}
	
	return false;
}
