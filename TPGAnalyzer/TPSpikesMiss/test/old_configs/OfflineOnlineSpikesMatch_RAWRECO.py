import FWCore.ParameterSet.Config as cms

process = cms.Process("OfflineSpikeCrystalToOnlineMatch")

# ---------------------------------------------------------------------
# Produce Ntuple Module
# ---------------------------------------------------------------------
#
## standard parameters
#process.produceNtuple.functionName = cms.string("Some Function Name")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')


process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi");
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi");
process.ecalEBunpacker.InputLabel = cms.InputTag('rawDataCollector');
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

# global tag for data
process.GlobalTag.globaltag = 'GR_R_53_V18::All'
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'root:///afs/cern.ch/work/i/iantropo/161E7F57-74CA-E111-AC13-0025901D631E.root'
        #'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/azabi/F0B823DE-763B-E211-A312-0025905964B2.root'
        #'file:/home/llr/cms/antropov/Run2012C/ZeroBias1/RAW/v1/000/197/923/2A0B6407-09C3-E111-928B-001D09F29533.root'
        #'root://xrootd.unl.edu//store/data/Run2012D/ZeroBias1/RECO/PromptReco-v1/000/204/840/78A494B6-8015-E211-B959-001D09F26C5C.root'
        #root://xrootd.unl.edu//store/data/Run2012C/ZeroBias1/RAW-RECO/25Feb2013-v1/10000/04C1CA46-BB7F-E211-82DA-003048678E24.root'
        'file:/home/llr/cms/antropov/Run2012C/04C1CA46-BB7F-E211-82DA-003048678E24.root'
        #'/store/data/Run2012C/ZeroBias1/RAW-RECO/25Feb2013-v1/10000/04C1CA46-BB7F-E211-82DA-003048678E24.root'
       
    )
)

process.load("SimCalorimetry.EcalTrigPrimProducers.ecalTriggerPrimitiveDigis_cff")
process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
#process.simEcalTriggerPrimitiveDigis.Label = 'ecalEBunpacker'
process.simEcalTriggerPrimitiveDigis.InstanceEB =  'ebDigis'
process.simEcalTriggerPrimitiveDigis.InstanceEE =  'eeDigis'
#process.simEcalTriggerPrimitiveDigis.BarrelOnly = True


process.OfflineSpikeCrystalToOnlineMatch = cms.EDAnalyzer('OfflineSpikeCrystalToOnlineMatch',
                              histogramFile         = cms.string('iurii_TPFile_OffSearch_hists.root'),
                              TPEmulatorCollection  = cms.InputTag("simEcalTriggerPrimitiveDigis"),
                              TPOnlineCollection    = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
                              # to get the RecHits --> spikes
                              EcalRecHitCollectionEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                              do_reconstruct_amplitude_allTTs=cms.bool(True),
                              do_rechit_spikes_search=cms.bool(True),
                              do_reconstruct_amplitudes_spikes=cms.bool(True), #overrides do_rechit_spikes_search
                              do_3DSpike_plots=cms.bool(True),
                              do_l1extraparticles=cms.bool(True)
                              )

process.TFileService = cms.Service ("TFileService",
                                    fileName = cms.string ("iurii_TPFile_OffSearch.root")
                                    )

process.p = cms.Path(process.RawToDigi * process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.OfflineSpikeCrystalToOnlineMatch)
#process.p = cms.Path(process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.OfflineSpikeCrystalToOnlineMatch)
