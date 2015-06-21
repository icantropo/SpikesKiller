#ifndef CandidateHadIsolation_h
#define CandidateHadIsolation_h

//C++ includes
#include <vector>
#include <functional>

//CMSSW includes
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"

class CandidateHadIsolation {
 public:
  
  //constructors
  CandidateHadIsolation ( ) ;
  CandidateHadIsolation (edm::ESHandle<CaloGeometry> ,
                HBHERecHitMetaCollection*  ,
		const reco::GsfElectron* ) ;
  CandidateHadIsolation (edm::ESHandle<CaloGeometry> , 
                HBHERecHitMetaCollection*  ,
		const reco::GsfElectron* ,
		const edm::View<reco::GsfElectron>* ) ;
  
  CandidateHadIsolation (edm::ESHandle<CaloGeometry> ,
                HBHERecHitMetaCollection*  ,
		const math::XYZPoint* ) ;
  //methods
  void setExtRadius (double extRadius) ;
  void setIntRadius (double intRadius) ;
  void setEtLow (double etLow) ;

  double getEtHadClusters () const ;
  double getHoE () const ;

  //destructor 
  ~CandidateHadIsolation() ;
  
 private:
  
  edm::ESHandle<CaloGeometry>  theCaloGeom_ ;
  HBHERecHitMetaCollection* mhbhe_ ;
  const reco::GsfElectron*  electron_ ;
  const edm::View<reco::GsfElectron> *electronCollection_ ;
  
  double extRadius_ ;
  double intRadius_ ;
  double etLow_ ;
  const math::XYZPoint* CaloPosition_;
  math::XYZPoint theCaloPosition;
};

#endif
