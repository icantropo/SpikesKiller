#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GtTechnicalTriggerRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtTechnicalTrigger.h"



#include "JetMETCorrections/MCJet/plugins/PFDataTreeProducer.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"


DEFINE_FWK_MODULE(PFDataTreeProducer);
