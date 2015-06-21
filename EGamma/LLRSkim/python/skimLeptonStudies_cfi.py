import FWCore.ParameterSet.Config as cms

skimLeptonStudies = cms.EDFilter(
    "SkimLeptonStudies",
    electronCollection = cms.InputTag("gedGsfElectrons"),
    muonCollection = cms.InputTag("muons"),
    isEleID = cms.bool(False),
    isMuonID = cms.bool(False),
    lep_ptLow = cms.double(15.),
    lep_ptHigh = cms.double(15.),
    nLep_ptLow = cms.int32(1),
    nLep_ptHigh = cms.int32(1)
)

