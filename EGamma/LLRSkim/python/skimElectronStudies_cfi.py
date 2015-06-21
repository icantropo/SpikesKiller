import FWCore.ParameterSet.Config as cms

skimElectronStudies = cms.EDFilter(
    "SkimElectronStudies",
    electronCollection = cms.InputTag("gedGsfElectrons"),
    ele_ptLow = cms.double(15.),
    ele_ptHigh = cms.double(15.),
    nEle_ptLowMIN = cms.int32(1),
    nEle_ptHighMIN = cms.int32(1),
    sc_EtLow = cms.double(15.),
    sc_EtHigh = cms.double(15.),
    nSC_EtLowMIN = cms.int32(1),
    nSC_EtHighMIN = cms.int32(1)
)

