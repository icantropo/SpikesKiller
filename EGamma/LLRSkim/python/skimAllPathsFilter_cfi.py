import FWCore.ParameterSet.Config as cms

skimAllPathsFilter = cms.EDFilter(
    "SkimAllPathsFilter",
    electronCollection = cms.InputTag("gedGsfElectrons"),
    muonCollection = cms.InputTag("muons"),
    mode = cms.string("ML"),
    eleID = cms.string("VBTF95")
)

