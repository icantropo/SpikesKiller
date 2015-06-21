
#include "EGamma/ECGelec/plugins/ConversionRemovalSelector.h"

//#include "RecoEgamma/Examples/plugins/ConversionRemovalSelector.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

bool ConversionRemovalSelector::selection(edm::Ref<reco::GsfElectronCollection> eleRef)
{
  bool isAmbiguous = true, isNotFromPixel = true;
  if (eleRef->ambiguousGsfTracksSize() == 0 ) {isAmbiguous = false;}
  
  /*
  TrackingRecHitRef rhit =  eleRef->gsfTrack()->extra()->recHit(0);
  int subdetId = rhit->geographicalId().subdetId();
  int layerId  = 0;
  DetId id = rhit->geographicalId();
  if (id.subdetId()==3) layerId = ((TIBDetId)(id)).layer();
  if (id.subdetId()==5) layerId = ((TOBDetId)(id)).layer();
  if (id.subdetId()==1) layerId = ((PXBDetId)(id)).layer();
  if (id.subdetId()==4) layerId = ((TIDDetId)(id)).wheel();
  if (id.subdetId()==6) layerId = ((TECDetId)(id)).wheel();
  if (id.subdetId()==2) layerId = ((PXFDetId)(id)).disk();
  //std::cout << " subdetIdele layerIdele = " << id.subdetId() << "   " << layerId << std::endl;

  if ((id.subdetId()==1 && layerId == 1) || (id.subdetId()==2 && layerId == 1)) {isNotFromPixel = false;}
  */


  int  mishits = eleRef->gsfTrack()->trackerExpectedHitsInner().numberOfHits(); 
  //std::cout << "mishits = " << mishits << std::endl;
  if (mishits == 0){isNotFromPixel = false;}

  if(isAmbiguous || isNotFromPixel) return false;
  
  return true;

  

}
