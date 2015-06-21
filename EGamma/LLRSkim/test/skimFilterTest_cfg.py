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
        'file:/data_CMS/cms/charlot/FirstBeam7TeV/ZeeEvents/Zee_run136100.root'
            )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# Skim analyzer
from EGamma.LLRSkim.skimElectronStudies_cfi import *
process.mySkimElectronStudies = skimElectronStudies.clone()

process.out = cms.OutputModule(
    "PoolOutputModule",
    verbose = cms.untracked.bool(True),
    outputCommands = cms.untracked.vstring('keep *_*_*_*'),
    fileName = cms.untracked.string('test.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
        )
    )

# Paths
process.p = cms.Path(
    process.mySkimElectronStudies
    )

process.o = cms.EndPath(process.out)