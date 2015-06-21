import FWCore.ParameterSet.Config as cms

stdardPreselectionSelector = cms.EDFilter("StdardPreselectionSelector",
                                          src = cms.InputTag("gedGsfElectrons")
                                       )
