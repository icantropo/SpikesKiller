import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012C/ZeroBias1/RAW/v1/000/197/923/2A0B6407-09C3-E111-928B-001D09F29533.root'
    )
)

process.demo = cms.EDAnalyzer('CodeTesting'
)


process.p = cms.Path(process.demo)
