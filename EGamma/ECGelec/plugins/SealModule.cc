#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "EGamma/ECGelec/plugins/SimpleNtple.h"
#include "EGamma/ECGelec/plugins/SimpleNtpleCustom.h"
#include "EGamma/ECGelec/plugins/SimpleNtpleSpike.h"
#include "EGamma/ECGelec/plugins/StdPreselectionSelector.h"
#include "EGamma/ECGelec/plugins/RunSelect.h"
#include "EGamma/ECGelec/plugins/RecHitFilter.h"
#include "EGamma/ECGelec/plugins/SpikeRemoval.h"
#include "EGamma/ECGelec/plugins/ConversionRemovalSelector.h"
#include "EGamma/ECGelec/plugins/TdrHzzIsolationProducer.h"
#include "EGamma/ECGelec/plugins/MuonHzzIsolationProducer.h"


typedef ObjectSelector<
  StdPreselectionSelector, 
  edm::RefVector<reco::GsfElectronCollection> 
  > StdardPreselectionSelectorRef ;

typedef ObjectSelector<
  StdPreselectionSelector
  > StdardPreselectionSelector;

typedef ObjectSelector<
  SpikeRemoval
  > SpikeRemovalSelector;

typedef ObjectSelector<
	   ConversionRemovalSelector, 
           edm::RefVector<reco::GsfElectronCollection> 
          > ConvRemovSelectorRef ;

typedef ObjectSelector<
           ConversionRemovalSelector 
          > ConvRemovSelector ;

DEFINE_FWK_MODULE(SimpleNtple);
DEFINE_FWK_MODULE(SimpleNtpleCustom);
DEFINE_FWK_MODULE(SimpleNtpleSpike);

DEFINE_FWK_MODULE(StdardPreselectionSelector);
DEFINE_FWK_MODULE(StdardPreselectionSelectorRef);
DEFINE_FWK_MODULE(RunSelect);
DEFINE_FWK_MODULE(RecHitFilter);
DEFINE_FWK_MODULE(SpikeRemovalSelector);
DEFINE_FWK_MODULE(ConvRemovSelector);
DEFINE_FWK_MODULE(ConvRemovSelectorRef);
DEFINE_FWK_MODULE(TdrHzzIsolationProducer);
DEFINE_FWK_MODULE(MuonHzzIsolationProducer);
