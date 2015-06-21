import FWCore.ParameterSet.Config as cms

spikeRemovalSelector = cms.EDFilter("SpikeRemovalSelector",
                                    src = cms.InputTag("gedGsfElectrons")
                                    )
