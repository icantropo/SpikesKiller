//C++ includes
#include <vector>
#include <functional>

//ROOT includes
#include <Math/VectorUtil.h>

//CMSSW includes
#include "EGamma/ECGelec/interface/CandidateHadIsolation.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoCaloTools/Selectors/interface/CaloDualConeSelector.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollectionV.h"

using namespace std;

CandidateHadIsolation::CandidateHadIsolation ()
{
}

CandidateHadIsolation::CandidateHadIsolation ( edm::ESHandle<CaloGeometry> theCaloGeom ,
                             HBHERecHitMetaCollection*  mhbhe,
                             const reco::GsfElectron* electron ) :
   theCaloGeom_(theCaloGeom) ,  
   mhbhe_(mhbhe) ,
   electron_(electron) 
{
  electronCollection_ = 0 ;
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  etLow_ = 1.5 ; 
  theCaloPosition = electron_->caloPosition() ;
}

CandidateHadIsolation::CandidateHadIsolation (edm::ESHandle<CaloGeometry> theCaloGeom ,
                            HBHERecHitMetaCollection*  mhbhe,
			    const reco::GsfElectron* electron , 
			    const edm::View<reco::GsfElectron>* electronCollection ) : 
  theCaloGeom_(theCaloGeom) , 
  mhbhe_(mhbhe) ,
  electron_(electron) ,
  electronCollection_(electronCollection)  
{
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  etLow_ = 1.5 ; 
}  

CandidateHadIsolation::CandidateHadIsolation (edm::ESHandle<CaloGeometry> theCaloGeom ,
                             HBHERecHitMetaCollection*  mhbhe,
                             const math::XYZPoint* CaloPosition ) :
   theCaloGeom_(theCaloGeom) ,  
   mhbhe_(mhbhe) ,
   CaloPosition_(CaloPosition) 
{
  electronCollection_ = 0 ;
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  etLow_ = 1.5 ; 
  theCaloPosition = *CaloPosition_;
}

CandidateHadIsolation::~CandidateHadIsolation ()
{
}

void CandidateHadIsolation::setExtRadius (double extRadius)
{
  extRadius_ = extRadius ;
}

void CandidateHadIsolation::setIntRadius (double intRadius)
{  
  intRadius_ = intRadius ;
}

void CandidateHadIsolation::setEtLow (double etLow)
{  
  etLow_ = etLow ;
}

double CandidateHadIsolation::getEtHadClusters () const
{
  double hcalEt = 0.;
  if (mhbhe_) 
   {
      //Take the SC position
      const CaloGeometry* caloGeom = theCaloGeom_.product();
      CaloDualConeSelector * sel = new CaloDualConeSelector(intRadius_ ,extRadius_, caloGeom, DetId::Hcal);
//       math::XYZPoint theCaloPosition = electron_->caloPosition () ;
      GlobalPoint pclu(	theCaloPosition.x () ,
                		theCaloPosition.y () ,
						theCaloPosition.z () );
      
      //Compute the HCAL energy behind ECAL
      std::auto_ptr<CaloRecHitMetaCollectionV> chosen = sel->select(pclu,*mhbhe_);
      for (CaloRecHitMetaCollectionV::const_iterator i = chosen->begin () ; i!= chosen->end () ; ++i)      {
			double hcalHit_eta = caloGeom->getPosition(i->detid()).eta();
			double hcalHit_Et = i->energy()*sin(2*atan(exp(-hcalHit_eta)));
			if ( hcalHit_Et > etLow_)
			hcalEt += hcalHit_Et;
// 			std::cout << "[CandidateHadIsolation] Hit value = " <<  hcalHit_Et << std::endl;
      }
    }
// 	std::cout << "[CandidateHadIsolation] HcalIsoET value = " <<  hcalEt<< std::endl;
  return hcalEt ;
}

double CandidateHadIsolation::getHoE () const
{
  double HoE = 0 ;
  if (mhbhe_) 
   {
     //Take the SC position
     const CaloGeometry* caloGeom = theCaloGeom_.product();
     CaloDualConeSelector * sel = new CaloDualConeSelector(intRadius_ ,extRadius_, caloGeom, DetId::Hcal);
     math::XYZPoint theCaloPosition = electron_->caloPosition () ;
     GlobalPoint pclu (theCaloPosition.x () ,
                       theCaloPosition.y () ,
		       theCaloPosition.z () );
     //Compute the HCAL energy behind ECAL
     double hcalEnergy = 0. ;
     std::auto_ptr<CaloRecHitMetaCollectionV> chosen = sel->select(pclu,*mhbhe_);
     for (CaloRecHitMetaCollectionV::const_iterator i = chosen->begin () ; 
                                                    i!= chosen->end () ; 
						    ++i) 
     {
       hcalEnergy += i->energy();
     }
     //Take the SC energy
     double ecalEnergy = electron_->caloEnergy () ;
     //Compute HoE
     HoE = hcalEnergy/ecalEnergy ;
   } 
  return HoE ;
}

