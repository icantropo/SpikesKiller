# Auto generated configuration file
# using: 
# Revision: 1.381.2.18 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: data -s RAW2DIGI,L1Reco,RECO --conditions auto:com10 -n 500 --filein=file:/build/argiro/data/MinBias-Run2011B-RAW.root --data --scenario pp --process RERECO
import FWCore.ParameterSet.Config as cms

process = cms.Process('RERECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(500)
)

# Input source
process.source = cms.Source("PoolSource",
   secondaryFileNames = cms.untracked.vstring(),
   fileNames = cms.untracked.vstring('/store/data/Run2012C/ZeroBias1/RAW/v1/000/197/923/2A0B6407-09C3-E111-928B-001D09F29533.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
   version = cms.untracked.string('$Revision: 1.381.2.18 $'),
   annotation = cms.untracked.string('data nevts:500'),
   name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
   splitLevel = cms.untracked.int32(0),
   eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
   outputCommands = process.RECOSIMEventContent.outputCommands,
   fileName = cms.untracked.string('data_RAW2DIGI_L1Reco_RECO.root'),
   dataset = cms.untracked.PSet(
       filterName = cms.untracked.string(''),
       dataTier = cms.untracked.string('')
   )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:com10', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.ecallocalreco_step= cms.Path(process.ecalLocalRecoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.ecallocalreco_step,process.endjob_step,process.RECOSIMoutput_step)