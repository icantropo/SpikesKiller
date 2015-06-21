import FWCore.ParameterSet.Config as cms

produceNtupleCustom = cms.EDAnalyzer("TPSpikesMiss",
                                     ## Custom Nadir

                                     # bolean parameters
                                     NadL1M = cms.untracked.bool(False),
                                     NadTP = cms.untracked.bool(False),
                                     NadTPmodif = cms.untracked.bool(False),
                                     NadTPemul = cms.untracked.bool(False),
                                     PrintDebug = cms.untracked.bool(False),

                                     # to get the trigger primitives
                                     TPCollectionNormal = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
                                     TPCollectionModif  = cms.InputTag("zeroedEcalTrigPrimDigis"),
                                     TPEmulatorCollection  = cms.InputTag("simEcalTriggerPrimitiveDigis"),
                                     
                                     hcalTowers = cms.InputTag("towerMaker"),
                                     hOverEPtMin = cms.double(0.), 
                                     hOverEConeSize = cms.double(0.15)
       )
