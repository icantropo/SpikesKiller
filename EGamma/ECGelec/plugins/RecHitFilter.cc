#include "EGamma/ECGelec/plugins/RecHitFilter.h"




 

//! ctor
RecHitFilter::RecHitFilter(const edm::ParameterSet& iConfig):
  aod_ (iConfig.getUntrackedParameter<bool>("AOD"))
{}

// ----------------------------------------------------------------






//! dtor
RecHitFilter::~RecHitFilter()
{}

// ----------------------------------------------------------------






void RecHitFilter::beginJob() 
{}

// ----------------------------------------------------------------






void RecHitFilter::endJob() 
{}

// ----------------------------------------------------------------


//! check the conversionsCollection size
bool RecHitFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  //   int run = iEvent.id().run();
  //   if (run != 124009 &&  run != 124020 && run != 124022 && run != 124023 && run != 124024) 
  //     return false;

  edm::Handle<EcalRecHitCollection> pBarrelEcalRecHits ;
  edm::Handle<EcalRecHitCollection> pEndcapEcalRecHits ;
  if(!aod_){
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEB"), pBarrelEcalRecHits );
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEE"), pEndcapEcalRecHits ) ;
  }
  else{
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEB"), pBarrelEcalRecHits );
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEE"), pEndcapEcalRecHits ) ;
  }

  const EcalRecHitCollection* theBarrelEcalRecHits = pBarrelEcalRecHits.product () ;
  const EcalRecHitCollection* theEndcapEcalRecHits = pEndcapEcalRecHits.product () ;


  if((theBarrelEcalRecHits->size() == 0 && theEndcapEcalRecHits->size() != 0) ||
     (theBarrelEcalRecHits->size() != 0 && theEndcapEcalRecHits->size() == 0))
    {
      //std::cout << "false" << std::endl;
      return false;
    }
  
  else
    {
      //std::cout << "true" << std::endl;
      return true;  
    }
}
