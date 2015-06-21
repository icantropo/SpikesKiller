import FWCore.ParameterSet.Config as cms

process = cms.Process("OfflineSpikeCrystalToOnlineMatch")

# Produce Ntuple Module
#process.produceNtuple.functionName = cms.string("Some Function Name")

# process.load('Configuration.StandardSequences.RawToDigi_Data_cff')

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

process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi");
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi");
process.ecalEBunpacker.InputLabel = cms.InputTag('rawDataCollector');
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
   version = cms.untracked.string('$Revision: 1.381.2.18 $'),
   annotation = cms.untracked.string('data nevts:500'),
   name = cms.untracked.string('PyReleaseValidation')
)

# Other statements  (fixes some LaserCorrectionsInfo errors. Antropov)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:com10', '')

# # global tag for data
# process.GlobalTag.globaltag = 'GR_R_53_V18::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/azabi/F0B823DE-763B-E211-A312-0025905964B2.root'
        #'file:/home/llr/cms/antropov/Run2012C/ZeroBias1/RAW/v1/000/197/923/2A0B6407-09C3-E111-928B-001D09F29533.root'
        #'root://xrootd.unl.edu//store/data/Run2012D/ZeroBias1/RECO/PromptReco-v1/000/204/840/78A494B6-8015-E211-B959-001D09F26C5C.root'
#         '/store/data/Run2012C/ZeroBias1/RAW/v1/000/197/923/80AC1792-00C3-E111-B1DD-001D09F24FEC.root',
#         '/store/data/Run2012C/ZeroBias1/RAW/v1/000/197/923/848BBF77-58C3-E111-AC7C-001D09F2437B.root',
#         '/store/data/Run2012C/ZeroBias1/RAW/v1/000/197/923/2A0B6407-09C3-E111-928B-001D09F29533.root'
        'file:/home/llr/cms/antropov/CMSSW_6_2_12/src/EmulSpikesFilter/EmulSpikesFilter/l1EmulatorFromRaw_RAW2DIGI_L1_pRECO.root'
#         'file:/data_CMS/cms/antropov/Run2012C/ZeroBias1/RAW/v1/000/197/923/2A0B6407-09C3-E111-928B-001D09F29533.root'
        #'/store/data/Run2012C/ZeroBias1/RAW-RECO/25Feb2013-v1/10000/04C1CA46-BB7F-E211-82DA-003048678E24.root'
        #'file:/home/llr/cms/antropov/Run2012C/04C1CA46-BB7F-E211-82DA-003048678E24.root'
    )
)

# to emulate with default parameters
# process.load("SimCalorimetry.EcalTrigPrimProducers.ecalTriggerPrimitiveDigis_cff")

# to simulate with config file
process.load('SimCalorimetry.EcalTrigPrimProducers.ecalTrigPrimESProducer_cff')
process.EcalTrigPrimESProducer.DatabaseFile = 'TPG_beamv6_trans_spikekill.txt.gz' 

# process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
process.simEcalTriggerPrimitiveDigis.Label = 'ecalEBunpacker'
process.simEcalTriggerPrimitiveDigis.InstanceEB =  'ebDigis'
process.simEcalTriggerPrimitiveDigis.InstanceEE =  'eeDigis'
process.simEcalTriggerPrimitiveDigis.BarrelOnly = True

process.OfflineSpikeCrystalToOnlineMatch = cms.EDAnalyzer('OfflineSpikeCrystalToOnlineMatch',
                              histogramFile         = cms.string('iurii_TPFile_OffSearch_hists.root'),
                              TPEmulatorCollection  = cms.InputTag("simEcalTriggerPrimitiveDigis"),
                              TPOnlineCollection    = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
                              # to get the RecHits --> spikes
                              EcalRecHitCollectionEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                              do_reconstruct_amplitude_allTTs=cms.bool(True),
                              do_rechit_spikes_search=cms.bool(True),
                              do_reconstruct_amplitudes_spikes=cms.bool(True), #overrides do_rechit_spikes_search
                              do_3DSpike_plots=cms.bool(True)
                              )

process.TFileService = cms.Service ("TFileService",
                                    fileName = cms.string ("iurii_TPFile_OffSearch.root")
                                    )

# # Output definition
# 
# process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    outputCommands = process.RECOSIMEventContent.outputCommands,
#    fileName = cms.untracked.string('data_RAW2DIGI_L1Reco_RECO.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('')
#    )
# )

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.ecallocalreco_step= cms.Path(process.ecalLocalRecoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)  # write output

# Schedule definition
process.p = cms.Path(process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.OfflineSpikeCrystalToOnlineMatch)
# Run the unpacker, UncalibRecHit producer and the RecHit producer + p process.
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.ecallocalreco_step,process.endjob_step, process.p)#,process.RECOSIMoutput_step)

# process.p = cms.Path(process.RawToDigi * process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.OfflineSpikeCrystalToOnlineMatch)
