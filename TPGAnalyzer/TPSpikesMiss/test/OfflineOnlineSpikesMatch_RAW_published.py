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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

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
#         '/store/data/Run2012C/ZeroBias1/RAW/v1/000/197/923/2A0B6407-09C3-E111-928B-001D09F29533.root'
#         '/store/user/iantropo/ZeroBias1/Run2012C-v1_RAW_Spikes_publish_test/371aa169538416419fdcbb2bba27730c/RecHitSpikesZeroBias1_Run2012C_RAW_1_1_0ZL.root'
#         '/store/user/iantropo/ZeroBias1/Run2012C-v1_RAW_Spikes_publish_test_2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_5_1_Mq0.root'
#         '/store/user/iantropo/ZeroBias1/Run2012C-v1_RAW_Spikes_publish_test_2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_1_1_qXs.rootcd vb'
        '/store/user/iantropo/ZeroBias1/Run2012C_v1_RAW_Spikes_v2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_97_1_kiA.root',
        '/store/user/iantropo/ZeroBias1/Run2012C_v1_RAW_Spikes_v2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_96_1_DEb.root',
        '/store/user/iantropo/ZeroBias1/Run2012C_v1_RAW_Spikes_v2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_95_1_L1w.root',
        '/store/user/iantropo/ZeroBias1/Run2012C_v1_RAW_Spikes_v2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_2746_1_RyX.root',
        '/store/user/iantropo/ZeroBias1/Run2012C_v1_RAW_Spikes_v2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_2747_1_xaA.root',
        '/store/user/iantropo/ZeroBias1/Run2012C_v1_RAW_Spikes_v2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_2765_1_y6E.root',
        '/store/user/iantropo/ZeroBias1/Run2012C_v1_RAW_Spikes_v2/71fb363099bf4cf205945e9168ab2c25/RecHitSpikesZeroBias1_Run2012C_RAW_2744_1_zyj.root'
    )
)

# to emulate with default parameters
# process.load("SimCalorimetry.EcalTrigPrimProducers.ecalTriggerPrimitiveDigis_cff")

# to simulate with config file
process.load('SimCalorimetry.EcalTrigPrimProducers.ecalTrigPrimESProducer_cff')
process.EcalTrigPrimESProducer.DatabaseFile = 'TPG_beamv6_trans_spikekill.txt.gz' 

process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
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
                              do_3DSpike_plots=cms.bool(False)
                              )

process.TFileService = cms.Service ("TFileService",
                                    fileName = cms.string ("SpikesEmulation.root")
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

process.p = cms.Path(process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.OfflineSpikeCrystalToOnlineMatch)

# Schedule definition
# Run the unpacker, UncalibRecHit producer and the RecHit producer + p process.
# p process - unpack EB, emulate Ecal Trigger Primitives
#             + TPG Analizer - collects some TPG Online and TPG Emulated information,
#             makes some cuts on it with Offline spikes. L1 electrons has been added.

# Assume presence of RAW, RawToDigi output and ecalLocalRecoSequence output in the .root file.
# faster
process.schedule = cms.Schedule(process.p)#,process.RECOSIMoutput_step)

# Run same directly on the RAW data. (slowerß)
# process.schedule = cms.Schedule(process.L1Reco_step,process.endjob_step, process.p)#,process.RECOSIMoutput_step)

process.MessageLogger.cerr.FwkReport.reportEvery = 100

# process.p = cms.Path(process.RawToDigi * process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.OfflineSpikeCrystalToOnlineMatch)
