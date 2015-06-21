import FWCore.ParameterSet.Config as cms

process = cms.Process("electronTreeProducer")

# import of standard configurations

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

# 4_4_0
#process.GlobalTag.globaltag = 'GR_R_44_V1::All'

# Florian FastSim 441
process.GlobalTag.globaltag = 'START44_V6::All'

HLT_name = 'HLT'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1)
                                        #SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )

# ---------------------------------------------------------------------
# Input Files
# ---------------------------------------------------------------------
process.source = cms.Source("PoolSource",
                            #debugFlag = cms.untracked.bool(True),
                            #debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring(
                            'root:///afs/cern.ch/work/t/thennequ/private/0C2732B7-C4C5-E311-A8E5-00266CFFA120.root'
                            #'castor:/castor/cern.ch/cms/store/mc/Fall13dr/DYToEE_M-50_Tune4C_13TeV-pythia8/AODSIM/castor_tsg_PU1bx50_POSTLS162_V1-v1/00000/04662F02-D8EC-E311-AA38-0025905A6132.root'
                            #'rfio:///store/mc/Fall13dr/DYToEE_M-50_Tune4C_13TeV-pythia8/AODSIM/castor_tsg_PU1bx50_POSTLS162_V1-v1/00000/04662F02-D8EC-E311-AA38-0025905A6132.root'
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/020931C6-CC0B-E211-8819-001D09F24D8A.root',
                            #'/store/mc/Summer13/DYToEE_M-20_TuneZ2star_14TeV-pythia6/GEN-SIM-RECO/UpgradePhase1Age3H_DR61SLHCx_PU140Bx25_STAR17_61_V6A-v1/10000/EE36C273-21CC-E211-8133-002618FDA204.root',
                            #'/store/relval/CMSSW_5_2_7-START52_V10_FastSim/RelValZEE/GEN-SIM-DIGI-RECO/v1/00000/0841244A-4A06-E211-A5AA-003048FFD730.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/10942A9F-CA0B-E211-9361-003048D2C0F0.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/24931AB6-C90B-E211-842A-003048D37456.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/24F74F8C-C80B-E211-B21A-0030486780A8.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/4881FB57-CB0B-E211-8711-003048D37538.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/5476662C-C70B-E211-80D1-003048D373F6.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/741A79EE-C90B-E211-BFA8-0030486780B8.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/7C78A145-D00B-E211-B02C-001D09F24303.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/8C357D1E-CC0B-E211-919A-003048D373AE.root',
                            #'store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/92CB7B80-D40B-E211-A1AA-001D09F24D8A.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/9A7CD4A0-CA0B-E211-BB53-003048673374.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/A454A6FF-D00B-E211-AAE1-001D09F253D4.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/CA30DD6F-D20B-E211-901C-001D09F25109.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/CC80BA0A-D10B-E211-9136-001D09F253C0.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/D2D031EF-C90B-E211-B42A-0030486780E6.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/FAF33CB6-C50B-E211-91E0-003048D2C01E.root',
                            #'/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/853/FE540270-D20B-E211-B0CE-001D09F2441B.root'
    ),                         
                            )

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ---------------------------------------------------------------------
# Ouptut File
# ---------------------------------------------------------------------
process.TFileService = cms.Service ("TFileService", 
                                    fileName = cms.string ("tree_testFastSim_TPG.root")
                                    )

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

                                       
process.runSelection = cms.EDFilter("RunSelect",
    requireNoTimeScan = cms.untracked.bool(True) ,
    requireCollidingBX = cms.untracked.bool(False),
    requireNoLumiScan = cms.untracked.bool(False),
    debug = cms.untracked.bool(False)
    )

# ---------------------------------------------------------------------
# Skim ALL Path Filter
# ---------------------------------------------------------------------
#load the EDfilter to select just skim data
process.load("EGamma.LLRSkim.skimAllPathsFilter_cfi")
from EGamma.LLRSkim.skimAllPathsFilter_cfi import *
process.skimAllPathsFilter = skimAllPathsFilter.clone()

# Nadir TagAndProbe : at least 2 ele with eT>5 and at least 1 ele passing eleID==VBTF95
process.skimAllPathsFilter.mode = "TP_nadir"
process.skimAllPathsFilter.eleID= "VBTF95"


# ---------------------------------------------------------------------
# JETS
# ---------------------------------------------------------------------
# JPT
#process.load('RecoJets.Configuration.RecoJPTJets_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
#JEC Corrections... to come !
# for 360: create colection of L2L3 corrected JPT jets: ak5JPTJetsL2L3  
# one need set of tags will be provided be JES
# process.p1 = cms.Path(process.ak5JPTJetsL2L3*process.dump)

# ---------------------------------------------------------------------
# Fast Jet Rho Correction
# ---------------------------------------------------------------------
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

# ---------------------------------------------------------------------
# Vertexing DA
4# ---------------------------------------------------------------------
#process.load("RecoVertex.Configuration.RecoVertex_cff")
from RecoVertex.Configuration.RecoVertex_cff import *
process.vertexreco = cms.Sequence(offlinePrimaryVertices*offlinePrimaryVerticesWithBS)

# ---------------------------------------------------------------------
# Produce eID infos
# ---------------------------------------------------------------------
###process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentification_cfi")
###New optimization
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
# ---------------------------------------------------------------------
# Produce eIso infos from HZZ official package
# ---------------------------------------------------------------------
process.load("EGamma.ECGelec.HzzIsolationSequences_cff")
# ---------------------------------------------------------------------
# Produce muIso infos from HZZ official package
# ---------------------------------------------------------------------
#process.load("EGamma.ECGelec.muonHzzIsolationSequences_cff")
# ---------------------------------------------------------------------
# Produce Ntuple Module
# ---------------------------------------------------------------------

#process.load("EGamma.ECGelec.NtupleProducer_cfi")
#from EGamma.ECGelec.NtupleProducer_cfi import *
#process.produceNtuple = produceNtuple.clone()

process.load("EGamma.ECGelec.NtupleProducer_custom_cfi")
from EGamma.ECGelec.NtupleProducer_custom_cfi import *
process.produceNtuple = produceNtupleCustom.clone()

## Nadir's parameters
process.produceNtuple.NadL1M = cms.untracked.bool(False)
process.produceNtuple.NadTP = cms.untracked.bool(False)
process.produceNtuple.NadTPemul = cms.untracked.bool(True)
process.produceNtuple.NadTPmodif = cms.untracked.bool(False)
process.produceNtuple.PrintDebug = cms.untracked.bool(False)
#process.produceNtuple.PrintDebug_HLT = cms.untracked.bool(False)

## standard parameters
process.produceNtuple.type = 'MC'
process.produceNtuple.AOD = cms.untracked.bool(True)
process.produceNtuple.FillSC = cms.untracked.bool(True)
process.produceNtuple.functionName = cms.string("EcalClusterEnergyUncertainty")
# Trigger Stuff
process.produceNtuple.HLTTag          = 'TriggerResults::' + HLT_name
process.produceNtuple.TriggerEventTag = 'hltTriggerSummaryAOD::' + HLT_name
process.produceNtuple.HLTElePaths     = cms.vstring(
    'HLT_Ele17_SW_TighterEleIdIsol_L1R_v3', 'HLT_Ele17_SW_TighterEleIdIsol_L1R_v2', 'HLT_Ele17_SW_TighterEleIdIsol_L1R_v1',
    'HLT_Ele17_SW_TightEleIdIsol_L1R', 'HLT_DoubleEle17_SW_L1R_v1', 'HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2',
    'HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1')
process.produceNtuple.HLTMuonPaths    = cms.vstring('HLT_Mu9')
process.produceNtuple.HLTFilters      = cms.VInputTag('hltL1NonIsoHLTNonIsoSingleElectronEt17TighterEleIdIsolTrackIsolFilter::'+HLT_name,
                                                      'hltL1NonIsoHLTNonIsoDoubleElectronEt17PixelMatchFilter::'+HLT_name,
                                                      'hltL1NonIsoHLTNonIsoSingleElectronEt17TightCaloEleIdEle8HEPixelMatchFilter::'+HLT_name,
                                                      'hltL1NonIsoHLTNonIsoSingleElectronEt17TightCaloEleIdEle8HEDoublePixelMatchFilter::'+HLT_name,
                                                      # Muon Trigger
                                                      'hltSingleMu9L3Filtered9')

#hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter::'+HLT_name)
#should add one for the Cleaned trigger?!

## HLT Filter from S. Beauceron
import HLTrigger.HLTfilters.hltHighLevel_cfi

process.MyHLTSelection = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    # SingleElectron paths
    #HLTPaths = [ 'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
    #             'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',
    #             'HLT_Ele45_CaloIdVT_TrkIdT_v3'
    #             ]
    # DoubleElectron paths
##     HLTPaths = [ 'HLT_Ele17_CaloIdL_CaloIsoVL_v3',
##                  'HLT_Ele8_CaloIdL_CaloIsoVL_v3',
##                  'HLT_Ele8_CaloIdL_TrkIdVL_v3',
##                  'HLT_Ele8_v3'
##                  'HLT_Ele17_CaloIdL_CaloIsoVL_v2',
##                  'HLT_Ele8_CaloIdL_CaloIsoVL_v2',
##                  'HLT_Ele8_CaloIdL_TrkIdVL_v2',
##                  'HLT_Ele8_v2'
##                  'HLT_Ele17_CaloIdL_CaloIsoVL_v1',
##                  'HLT_Ele8_CaloIdL_CaloIsoVL_v1',
##                  'HLT_Ele8_CaloIdL_TrkIdVL_v1',
##                  'HLT_Ele8_v1'
##                  ],

    # to get the spikes
    #HLTPaths = [ 'HLT_Activity_Ecal_SC*_*' ],
    #HLTPaths = [ 'HLT_Activity_Ecal_SC*_*',
    #             'HLT_L1SingleEG5_*' ],
    HLTPaths = [ 'HLT_L1SingleEG*' ],
    
    throw = False
    #dont throw except on unknown path name
    )
#process.HLTfilter = cms.Path( process.MyHLTSelection )

# ---------------------------------------------------------------------
# Sequence PATH
# ---------------------------------------------------------------------
process.p = cms.Path (
    #process.MyHLTSelection +
    process.vertexreco + 
    process.skimAllPathsFilter +   
    process.kt6PFJets + 
    
    #process.runSelection +

    #produce the eID CiC value maps
    process.eidVeryLoose+
    process.eidLoose+
    process.eidMedium+
    process.eidTight+
    process.eidSuperTight+
    process.eidHyperTight1+
    #process.eidHyperTight2+
    #process.eidHyperTight3+
    #process.eidHyperTight4+
    process.HzzIsolationSequence +
    #process.MuonHZZIsolationSequence +

    process.produceNtuple 
    )

#process.schedule = cms.Schedule( process.p )


