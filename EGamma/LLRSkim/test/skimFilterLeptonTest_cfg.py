import FWCore.ParameterSet.Config as cms

process = cms.Process("SkimAnalyzer")


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.MessageLogger.cerr.threshold = 'INFO'


# source
process.source = cms.Source(
    "PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(  
        #RUN 139459 EG SD
        #'/store/data/Run2010A/EG/RECO/v4/000/139/459/ECA0E034-D088-DF11-B34E-0030487A18A4.root',
        #'/store/data/Run2010A/EG/RECO/v4/000/139/459/CA9AF9F7-CB88-DF11-BACD-001D09F29524.root',
        #'/store/data/Run2010A/EG/RECO/v4/000/139/459/8C531733-D088-DF11-B6D2-0030487C6090.root',
        #RUN 139459 MU SD
        '/store/data/Run2010A/Mu/RECO/v4/000/139/459/CC6B9944-CB88-DF11-853C-0030487C5CE2.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/139/459/848B09F4-C488-DF11-9A73-001617C3B70E.root',

        #'file:/data_CMS/cms/ochando/DATA/run138921_EGrecov4/54B89BD2-CA83-DF11-B2B7-001617C3B6CC.root'
            )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *
process.hltFilter = triggerResultsFilter.clone()
process.hltFilter.hltResults = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.l1tResults = cms.InputTag('')
process.hltFilter.throw = cms.bool(False)
process.hltFilter.daqPartitions = cms.uint32( 1 )
process.hltFilter.triggerConditions = cms.vstring('NOT HLT_Photon10_Cleaned_L1R OR NOT HLT_Ele10_SW_L1R OR NOT HLT_Mu9')


# Skim analyzer
from EGamma.LLRSkim.skimLeptonStudies_cfi import *
process.SingleLeptonSkim = skimLeptonStudies.clone()
process.DoubleLeptonSkim = skimLeptonStudies.clone()
process.TripleLeptonSkim = skimLeptonStudies.clone()
##Single
process.SingleLeptonSkim.isEleID = cms.bool(True)
process.SingleLeptonSkim.isMuonID = cms.bool(True)
process.SingleLeptonSkim.lep_ptLow = cms.double(20.)
process.SingleLeptonSkim.lep_ptHigh = cms.double(20.)
process.SingleLeptonSkim.nLep_ptLow = cms.int32(1)
process.SingleLeptonSkim.nLep_ptHigh = cms.int32(1)
##Double
process.DoubleLeptonSkim.isEleID = cms.bool(False)
process.DoubleLeptonSkim.isMuonID = cms.bool(False)
process.DoubleLeptonSkim.lep_ptLow = cms.double(10.)
process.DoubleLeptonSkim.lep_ptHigh = cms.double(15.)
process.DoubleLeptonSkim.nLep_ptLow = cms.int32(2)
process.DoubleLeptonSkim.nLep_ptHigh = cms.int32(1)
##Triple
process.TripleLeptonSkim.isEleID = cms.bool(False)
process.TripleLeptonSkim.isMuonID = cms.bool(False)
process.TripleLeptonSkim.lep_ptLow = cms.double(5.)
process.TripleLeptonSkim.lep_ptHigh = cms.double(10.)
process.TripleLeptonSkim.nLep_ptLow = cms.int32(3)
process.TripleLeptonSkim.nLep_ptHigh = cms.int32(2)

process.out = cms.OutputModule(
    "PoolOutputModule",
    verbose = cms.untracked.bool(True),
    outputCommands = cms.untracked.vstring('keep *_*_*_*'),
    fileName = cms.untracked.string('LepSkimMuMiniTest.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1','p2','p3')
        )
    )
process.out1 = cms.OutputModule(
    "PoolOutputModule",
    verbose = cms.untracked.bool(True),
    outputCommands = cms.untracked.vstring('keep *_*_*_*'),
    fileName = cms.untracked.string('SingleLepSkim.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
        )
    )
process.out2 = cms.OutputModule(
    "PoolOutputModule",
    verbose = cms.untracked.bool(True),
    outputCommands = cms.untracked.vstring('keep *_*_*_*'),
    fileName = cms.untracked.string('DoubleLepSkim.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p2')
        )
    )
process.out3 = cms.OutputModule(
    "PoolOutputModule",
    verbose = cms.untracked.bool(True),
    outputCommands = cms.untracked.vstring('keep *_*_*_*'),
    fileName = cms.untracked.string('TripleLepSkim.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p3')
        )
    )

# Paths
process.p1 = cms.Path( process.SingleLeptonSkim)
process.p2 = cms.Path( process.DoubleLeptonSkim)
process.p3 = cms.Path( process.TripleLeptonSkim)
#process.p1 = cms.Path( process.hltFilter + process.SingleLeptonSkim)
#process.p2 = cms.Path( process.hltFilter + process.DoubleLeptonSkim)
#process.p3 = cms.Path( process.hltFilter + process.TripleLeptonSkim)

process.o = cms.EndPath(process.out)
#process.o = cms.EndPath(process.out+process.out1+process.out2+process.out3)