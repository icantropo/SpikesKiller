import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.EcalRecAlgos.ecalCleaningAlgo import cleaningAlgoConfig 

process = cms.Process("RecHitSpikesFilter")

process.EmulSpikesFilter = cms.EDFilter('EmulSpikesFilter',
                                       EcalRecHitCollectionEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"), 
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
from RecoEcal.Configuration.RecoEcal_cff import *

process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi");
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi");


process.ecalEBunpacker.InputLabel = cms.InputTag('rawDataCollector');
#process.ecalEBunpacker.InputLabel = cms.InputTag('source');

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")

process.load('Configuration.StandardSequences.RawToDigi_Data_cff') #modif-alex
process.load('Configuration.StandardSequences.SimL1Emulator_cff') 
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('L1Trigger.Configuration.L1Trigger_EventContent_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

# global tag for data
process.GlobalTag.globaltag = 'GR_R_53_V18::All'

# Output definition
process.SpecialEventContent = cms.PSet(
#         outputCommands = cms.untracked.vstring('drop *'),
        outputCommands = cms.untracked.vstring('keep *'),
        splitLevel = cms.untracked.int32(0)
       )

#Some default 'keep' options
# process.SpecialEventContent.outputCommands.extend(process.RECOEventContent.outputCommands)
# process.SpecialEventContent.outputCommands.extend(process.L1TriggerFEVTDEBUG.outputCommands)
#process.SpecialEventContent.outputCommands.append('keep *_l1extraParticlesOnline_*_*')
#process.SpecialEventContent.outputCommands.append('keep *_zeroedEcalTrigPrimDigis_*_*')
# process.SpecialEventContent.outputCommands.append('keep *_ecalDigis_*_*')
# process.SpecialEventContent.outputCommands.append('keep *_hcalDigis_*_*') #modif-new
#process.SpecialEventContent.outputCommands.append('keep *_simEcalTriggerPrimitiveDigis_*_*')
#process.SpecialEventContent.outputCommands.append('keep *_SimGtDigis_*_*') #modif

# Keep RAW data
# process.SpecialEventContent.outputCommands.append('keep *FEDRawDataCollection_*_*_*')

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
                                              splitLevel = cms.untracked.int32(0),
                                              #outputCommands = cms.untracked.vstring('keep *'),
                                              #outputCommands = process.RECOEventContent.outputCommands,
                                              outputCommands = process.SpecialEventContent.outputCommands,
                                              fileName = cms.untracked.string('RecHitSpikesZeroBias1_Run2012C_RAW.root'),
                                              SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring ('filter') ),
                                              dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('EmulSpikesFilter'),
        dataTier = cms.untracked.string('DIGI-RECO')
        )
)
process.output_step = cms.EndPath(process.FEVTDEBUGHLToutput)

#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

#############   Define the source file ###############
process.source = cms.Source("PoolSource",
 fileNames=cms.untracked.vstring(
        # Data
#         'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/azabi/F0B823DE-763B-E211-A312-0025905964B2.root')
        '/store/data/Run2012C/ZeroBias1/RAW/v1/000/197/923/2A0B6407-09C3-E111-928B-001D09F29533.root'
        )
)

# Other statements  (fixes some LaserCorrectionsInfo errors. Antropov)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:com10', '')

#############   Path       ###########################
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.ecallocalreco_step= cms.Path(process.ecalLocalRecoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)

# data - RAW
process.filter = cms.Path(process.EmulSpikesFilter)

# process.schedule = cms.Schedule(process.raw2digi_step, process.L1Reco_step, process.ecallocalreco_step, process.endjob_step, process.filter, process.output_step)
process.schedule = cms.Schedule(process.raw2digi_step, process.ecallocalreco_step, process.endjob_step, process.filter, process.output_step)

# mc - RAW

#process.p = cms.Path(process.noscraping * process.primaryVertexFilter * process.HBHENoiseFilter * process.eeBadScFilter * process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.pfDataTree)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100


