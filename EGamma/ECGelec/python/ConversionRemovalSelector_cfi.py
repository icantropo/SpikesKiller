
import FWCore.ParameterSet.Config as cms

ConvRemovSelector = cms.EDFilter("ConvRemovSelector",
                                    src = cms.InputTag("gedGsfElectrons")
                                    )
