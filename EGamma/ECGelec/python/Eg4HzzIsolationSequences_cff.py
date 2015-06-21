import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *

eleIsoDepositTk.src = "gedGsfElectrons"

eleIsoDepositEcalFromHits.src = "gedGsfElectrons"
eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = "reducedEcalRecHitsEB"
eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = "reducedEcalRecHitsEE"

eleIsoDepositHcalFromTowers.src = "gedGsfElectrons"

eleIsoFromDepsTkOptimized = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositTk"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('muons:0.01', 
                           'gedGsfElectrons:0.015', 
                           'gedGsfElectrons:RectangularEtaPhiVeto(-0.005,0.005,-0.3,0.3)', 
                           'EcalEndcaps:RectangularEtaPhiVeto(-0.005,0.005,-0.5,0.5)',
                           'Threshold(0.7)'),
       skipDefaultVeto = cms.bool(True)
   ))
)

eleIsoFromDepsEcalFromHitsByCrystalOptimized = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositEcalFromHits"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('NumCrystalVeto(3.0)',
                           'EcalBarrel:NumCrystalEtaPhiVeto(1.0,9999.0)',
                           'EcalEndcaps:NumCrystalEtaPhiVeto(1.5,9999.0)',
                           'EcalBarrel:AbsThresholdFromTransverse(0.08)',
                           'EcalEndcaps:AbsThreshold(0.20)',
                           'muons:0.05',
                           'gedGsfElectrons:NumCrystalVeto(3.0)',
                           'gedGsfElectrons:NumCrystalEtaPhiVeto(1.5,15.0)'),
       skipDefaultVeto = cms.bool(True)
   ))
)

eleIsoFromDepsHcalFromTowersOptimized = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       src = cms.InputTag("eleIsoDepositHcalFromTowers"),
       deltaR = cms.double(0.4),
       weight = cms.string('1'),
       vetos = cms.vstring('0.15', 'muons:0.05'),
       skipDefaultVeto = cms.bool(True),
       mode = cms.string('sum')
   ))
) 

from Configuration.StandardSequences.GeometryIdeal_cff import *

ElectronIsolationMakeIsoDeposits = cms.Sequence(
    eleIsoDepositTk+
    eleIsoDepositEcalFromHits#+
    #eleIsoDepositHcalFromTowers
)

ElectronIsolationMakeIsoValueMaps = cms.Sequence(
    eleIsoFromDepsTkOptimized+
    eleIsoFromDepsEcalFromHitsByCrystalOptimized#+
    #eleIsoFromDepsHcalFromTowersOptimized
)

Eg4HzzElectronIsolationDepositSequence = cms.Sequence(
    ElectronIsolationMakeIsoDeposits*
    ElectronIsolationMakeIsoValueMaps
)
