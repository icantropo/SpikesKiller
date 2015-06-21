import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.EcalRecAlgos.ecalCleaningAlgo import cleaningAlgoConfig 
#from RecoEcal.EgammaClusterProducers.ecalRecHitFlags_cfi import *
#from RecoEcal.EgammaClusterProducers.ecalSeverityLevelAlgos_cfi import *
#from RecoEcal.EgammaClusterProducers.ecalSeverityLevelFlags_cfi import *


process = cms.Process("Ana")


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
process.load('L1Trigger.Configuration.L1Trigger_EventContent_cff')

# global tag for run-dependent MC
#process.GlobalTag.globaltag = 'START53_V7N::All'

# global tag for data
process.GlobalTag.globaltag = 'GR_R_53_V18::All'


# Output definition
process.SpecialEventContent = cms.PSet(
       outputCommands = cms.untracked.vstring('drop *'),
          splitLevel = cms.untracked.int32(0)
       )

#process.SpecialEventContent.outputCommands.extend(process.RECOEventContent.outputCommands)
process.SpecialEventContent.outputCommands.extend(process.L1TriggerFEVTDEBUG.outputCommands)
process.SpecialEventContent.outputCommands.append('keep *_l1extraParticlesOnline_*_*')
process.SpecialEventContent.outputCommands.append('keep *_zeroedEcalTrigPrimDigis_*_*')
process.SpecialEventContent.outputCommands.append('keep *_ecalDigis_*_*')
process.SpecialEventContent.outputCommands.append('keep *_hcalDigis_*_*') #modif-new
process.SpecialEventContent.outputCommands.append('keep *_simEcalTriggerPrimitiveDigis_*_*')
process.SpecialEventContent.outputCommands.append('keep *_SimGtDigis_*_*') #modif

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
                                              splitLevel = cms.untracked.int32(0),
                                              #outputCommands = cms.untracked.vstring('keep *'),
                                              #outputCommands = process.RECOEventContent.outputCommands,
                                              outputCommands = process.SpecialEventContent.outputCommands,
                                              fileName = cms.untracked.string('l1EmulatorFromRaw_RAW2DIGI_L1_pRECO.root'),
                                              dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DIGI-RECO')
        )
)
process.output_step = cms.EndPath(process.FEVTDEBUGHLToutput)

#
# customization fragment to run L1 emulator starting from a RAW file
#
# V.M. Ghete 2010-06-09
Def customise(process):
    
    #
    # (re-)run the  L1 emulator starting from a RAW file
    #
    from L1Trigger.Configuration.L1Trigger_custom import customiseL1EmulatorFromRaw
    process=customiseL1EmulatorFromRaw(process)
    
    #
    # special configuration cases (change to desired configuration in customize_l1TriggerConfiguration)
    #
    from L1Trigger.Configuration.customise_l1TriggerConfiguration import customiseL1TriggerConfiguration
    process=customiseL1TriggerConfiguration(process)
    
    #
    # customization of output commands
    #
    from L1Trigger.Configuration.L1Trigger_custom import customiseOutputCommands
    process=customiseOutputCommands(process)
    
    return (process)

# Customise the process as-is
process = customise(process)


#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
 fileNames=cms.untracked.vstring(
        # Data
        'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/azabi/F0B823DE-763B-E211-A312-0025905964B2.root')
)



#############   Geometry  ###############

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");


#############    Analyser options   ############################
process.pfDataTree = cms.EDAnalyzer("PFDataTreeProducer",
jets                 = cms.string('ak5PFJets'),
histogramFile        = cms.string('ak5PFDataTree_data_50ns_tpgen_rechitmatch_tpoot.root'),
tracks               = cms.string('generalTracks'),
vertex               = cms.string('offlinePrimaryVertices'),
JetCorrectionService = cms.string('L2L3JetCorrectorAK5PF'),
EBRecHitCollection   = cms.string('ReducedEcalRecHitsEB'),
EERecHitCollection   = cms.string('ReducedEcalRecHitsEE'),
IsMC                 = cms.bool(False),  # set to True if using MC
OnlineTPs            = cms.InputTag("ecalDigis", "EcalTriggerPrimitives"), #modif-alex
cleaningConfig       = cleaningAlgoConfig,
badsc_coordinatesEE  = cms.vint32(-1023023,1048098,-1078063),
phoProducer          = cms.string('photons'),                             
photonCollection     = cms.string(''),

# bunch structure, run 200473 - 50ns
bunchstartbx         = cms.vint32(66,146,226,306,413,493,573,653,773,853,960,1040,1120,1200,1307,1387,1467,1547,1667,1747,1854,1934,2014,2094,2201,2281,2361,2441,2549,2629,2736,2816,2896,2976,3083,3163,3243,3323)
)

process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.HLTPaths = cms.vstring("HLT_ZeroBias*")

#modif-alex
process.load('SimCalorimetry.EcalTrigPrimProducers.ecalTrigPrimESProducer_cff')
#process.EcalTrigPrimESProducer.DatabaseFile = 'TPG_beamv5_mc_ideal.txt.gz'
process.EcalTrigPrimESProducer.DatabaseFile = 'TPG_beamv6_trans_spikekill.txt.gz'
#process.EcalTrigPrimESProducer.DatabaseFile = 'TPG_beamv6_notrans_spikekill_newscale.txt.gz'
#process.EcalTrigPrimESProducer.DatabaseFile = 'TPG_beamv6_notrans_spikekill_nonlinear.txt.gz'
#process.EcalTrigPrimESProducer.DatabaseFile = 'TPG_beamv6_notrans_spikekill_newscale_nonlinear.txt.gz'

process.load("SimCalorimetry.EcalTrigPrimProducers.ecalTriggerPrimitiveDigis_cff")
process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
#process.simEcalTriggerPrimitiveDigis.Label = 'ecalEBunpacker'
process.simEcalTriggerPrimitiveDigis.InstanceEB =  'ebDigis'
process.simEcalTriggerPrimitiveDigis.InstanceEE =  'eeDigis'
#process.simEcalTriggerPrimitiveDigis.BarrelOnly = True

# emulator trigger
process.simRctDigis.ecalDigis = cms.VInputTag(cms.InputTag("simEcalTriggerPrimitiveDigis"))
process.simRctDigis.hcalDigis = cms.VInputTag(cms.InputTag("hcalDigis"))
process.simGctDigis.inputLabel = cms.InputTag("simRctDigis")

# L1 extra for the re-simulated candidates
process.l1extraParticles = cms.EDProducer("L1ExtraParticlesProd",
                                          muonSource = cms.InputTag("gtDigis"),
                                          etTotalSource = cms.InputTag("simGctDigis"),
                                          nonIsolatedEmSource = cms.InputTag("simGctDigis","nonIsoEm"),
                                          etMissSource = cms.InputTag("simGctDigis"),
                                          htMissSource = cms.InputTag("simGctDigis"),
                                          produceMuonParticles = cms.bool(False),
                                          forwardJetSource = cms.InputTag("simGctDigis","forJets"),
                                          centralJetSource = cms.InputTag("simGctDigis","cenJets"),
                                          produceCaloParticles = cms.bool(True),
                                          tauJetSource = cms.InputTag("simGctDigis","tauJets"),
                                          isolatedEmSource = cms.InputTag("simGctDigis","isoEm"),
                                          etHadSource = cms.InputTag("simGctDigis"),
                                          hfRingEtSumsSource = cms.InputTag("simGctDigis"),
                                          hfRingBitCountsSource = cms.InputTag("simGctDigis"),
                                          centralBxOnly = cms.bool(True),
                                          ignoreHtMiss = cms.bool(False)
                                          )


# L1 extra for the online candidates
process.l1extraParticlesOnline = cms.EDProducer("L1ExtraParticlesProd",
                                                muonSource = cms.InputTag("gtDigis"),
                                                etTotalSource = cms.InputTag("gctDigis"),
                                                nonIsolatedEmSource = cms.InputTag("gctDigis","nonIsoEm"),
                                                etMissSource = cms.InputTag("gctDigis"),
                                                htMissSource = cms.InputTag("gctDigis"),
                                                produceMuonParticles = cms.bool(False),
                                                forwardJetSource = cms.InputTag("gctDigis","forJets"),
                                                centralJetSource = cms.InputTag("gctDigis","cenJets"),
                                                produceCaloParticles = cms.bool(True),
                                                tauJetSource = cms.InputTag("gctDigis","tauJets"),
                                                isolatedEmSource = cms.InputTag("gctDigis","isoEm"),
                                                etHadSource = cms.InputTag("gctDigis"),
                                                hfRingEtSumsSource = cms.InputTag("gctDigis"),
                                                hfRingBitCountsSource = cms.InputTag("gctDigis"),
                                                centralBxOnly = cms.bool(True),
                                                ignoreHtMiss = cms.bool(False)
                                                )


process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger.HLTfilters.hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND NOT (36 OR 37 OR 38 OR 39)')



process.load('RecoMET.METFilters.eeBadScFilter_cfi')



process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )


process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


#############   Path       ###########################



# data - RAW
#process.p = cms.Path(process.hltHighLevel * process.noscraping * process.primaryVertexFilter * process.HBHENoiseFilter * process.eeBadScFilter * process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.pfDataTree)
#process.p = cms.Path(process.RawToDigi * process.noscraping * process.primaryVertexFilter * process.HBHENoiseFilter * process.eeBadScFilter * process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.pfDataTree)
process.p = cms.Path(process.RawToDigi * process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.simRctDigis * process.simGctDigis * process.simGtDigis * process.l1extraParticles * process.l1extraParticlesOnline * process.pfDataTree)

process.schedule = cms.Schedule(process.p,process.output_step)

# mc - RAW

#process.p = cms.Path(process.noscraping * process.primaryVertexFilter * process.HBHENoiseFilter * process.eeBadScFilter * process.ecalEBunpacker * process.simEcalTriggerPrimitiveDigis * process.pfDataTree)


#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100


