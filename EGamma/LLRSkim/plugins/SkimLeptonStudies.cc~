#include "EGamma/LLRSkim/plugins/SkimLeptonStudies.h"

#include "TMath.h"

//! ctor
SkimLeptonStudies::SkimLeptonStudies(const edm::ParameterSet& iConfig):
electronCollection_               (iConfig.getParameter<edm::InputTag>("electronCollection")),
muonCollection_                   (iConfig.getParameter<edm::InputTag>("muonCollection")),
isEleID_                          (iConfig.getParameter<bool>("isEleID")),
isMuonID_                         (iConfig.getParameter<bool>("isMuonID")),
lep_ptLow_                        (iConfig.getParameter<double>("lep_ptLow")),
lep_ptHigh_                       (iConfig.getParameter<double>("lep_ptHigh")),
nLep_ptLow_                    (iConfig.getParameter<int>("nLep_ptLow")),
nLep_ptHigh_                   (iConfig.getParameter<int>("nLep_ptHigh"))
/*
ele_ptLow_                        (iConfig.getParameter<double>("ele_ptLow")),
ele_ptHigh_                       (iConfig.getParameter<double>("ele_ptHigh")),
nEle_ptLowMIN_                    (iConfig.getParameter<int>("nEle_ptLowMIN")),
nEle_ptHighMIN_                   (iConfig.getParameter<int>("nEle_ptHighMIN")),
muonCollection_                   (iConfig.getParameter<edm::InputTag>("muonCollection")),
muon_ptLow_                       (iConfig.getParameter<double>("muon_ptLow")),
muon_ptHigh_                      (iConfig.getParameter<double>("muon_ptHigh")),
nMuon_ptLowMIN_                   (iConfig.getParameter<int>("nMuon_ptLowMIN")),
nMuon_ptHighMIN_                  (iConfig.getParameter<int>("nMuon_ptHighMIN"))
*/
{}

//! dtor
SkimLeptonStudies::~SkimLeptonStudies()
{}


bool SkimLeptonStudies::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	/**int runId = iEvent.id().run();
	int lumiId = iEvent.luminosityBlock();
        int eventId = iEvent.id().event();*/
	
	//GSF Electron	
	edm::Handle<reco::GsfElectronCollection> EleHandle;
	iEvent.getByLabel(electronCollection_.label(), EleHandle);
	
	// loop on electrons
	int nEle_ptLow = 0;
	int nEle_ptHigh = 0;
	for(unsigned int eleIt = 0; eleIt < EleHandle->size(); ++eleIt)
	{
		
		reco::GsfElectronRef eleRef(EleHandle, eleIt);
		
		//electron ID
		bool eleIDBool = true ;
		if (isEleID_) {
		  //barrel
		  if (eleRef->isEB()) {
		    if ( fabs(eleRef->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) eleIDBool = false ;
		    if ( fabs(eleRef->deltaPhiSuperClusterTrackAtVtx()) > 0.08 )  eleIDBool = false ;
		    if ( eleRef->hadronicOverEm()  > 0.15 ) eleIDBool = false ;
		    if ( eleRef->sigmaIetaIeta() > 0.01 ) eleIDBool = false ;
		    
		  }
		  //endcap
		  if (eleRef->isEE()) {
		    if ( fabs(eleRef->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) eleIDBool = false ;
		    if ( fabs(eleRef->deltaPhiSuperClusterTrackAtVtx()) > 0.07)  eleIDBool = false ;
		    if ( eleRef->hadronicOverEm()  > 0.07 ) eleIDBool = false ;
		    if ( eleRef->sigmaIetaIeta() > 0.03 ) eleIDBool = false ;
		    
		  } // if endcaps
		} // if Ele ID
		
		// count electrons
		if( eleIDBool && eleRef->pt() > lep_ptLow_) ++nEle_ptLow;
		if( eleIDBool && eleRef->pt() > lep_ptHigh_) ++nEle_ptHigh;
		
	} // end loop on electrons
	

	edm::Handle<reco::MuonCollection> MuonHandle;
	iEvent.getByLabel(muonCollection_, MuonHandle);
	
	// loop on muons
	int nMuon_ptLow = 0;
	int nMuon_ptHigh = 0;
	for(unsigned int muonIt = 0; muonIt < MuonHandle->size(); ++muonIt)
	{ 	
		edm::Ref<reco::MuonCollection> muonRef(MuonHandle,muonIt);
	
		//muon ID
		bool muonIDBool = true ;
		if (isMuonID_) {
			if ( muonRef->isGlobalMuon() == 0 ) muonIDBool = false ;
		}
		// count muons
		if( muonIDBool && muonRef->pt() > lep_ptLow_) ++nMuon_ptLow;
		if( muonIDBool && muonRef->pt() > lep_ptHigh_) ++nMuon_ptHigh;
		
	} // end loop on muons
	

	// print out results
	if ( ( nEle_ptLow + nMuon_ptLow >= nLep_ptLow_ ) &&  
	     ( nEle_ptHigh + nMuon_ptHigh >= nLep_ptHigh_ )  ){
		//std::cout << ">>>SkimLeptonStudies::run=" << runId << "::lumi=" 
		//<< lumiId << "::eventId=" 
		//<< eventId << "::Found at least "
		//<< nLep_ptLow_ << " electrons with pt > " << lep_ptLow_ << " AND " 
		//<< nLep_ptHigh_ << " electrons with pt > " << lep_ptHigh_ << std::endl; 
		return true;
	}
	
	return false;
}
