#include <FWCore/PluginManager/interface/ModuleDef.h>


#include "EGamma/LLRSkim/plugins/SkimElectronStudies.h"
DEFINE_FWK_MODULE(SkimElectronStudies);

#include "EGamma/LLRSkim/plugins/SkimLeptonStudies.h"
DEFINE_FWK_MODULE(SkimLeptonStudies);

#include "EGamma/LLRSkim/plugins/SkimAllPathsFilter.h"
DEFINE_FWK_MODULE(SkimAllPathsFilter);

