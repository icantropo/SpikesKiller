#ifndef RecHitFilter_h
#define RecHitFilter_h

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"



class RecHitFilter : public edm::EDFilter
{
  
 public:
  
  //! ctor
  explicit RecHitFilter (const edm::ParameterSet&);
  
  //! dtor 
  ~RecHitFilter();
  
  
  
 private:
  
  void beginJob();
  void endJob();
  
  //! the actual filter method 
  bool filter(edm::Event&, const edm::EventSetup&);
  
  
  
 private:
  bool aod_;
  
};

#endif
