// -*- C++ -*-
//
// Package:    OfflineSpikeCrystalToOnlineMatch
// Class:      OfflineSpikeCrystalToOnlineMatch
//
/**\class OfflineSpikeCrystalToOnlineMatch OfflineSpikeCrystalToOnlineMatch.cc TPGAnalyzer/TPSpikesMiss/plugins/OfflineSpikeCrystalToOnlineMatch.cc

 Description: [one line class summary]

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Iurii Antropov
//         Created:  Fri, 31 Oct 2014 16:36:59 GMT
// $Id$
//
//
// system include files
#include <memory>
#include <iostream>
#include <cmath>

//includes root
#include "TTree.h"
//#include "TNtuple.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

//// CMSSW
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

//file service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//ECAL includes
#include "CondFormats/EcalObjects/interface/EcalTPGTowerStatus.h" // for the EcalTPGTowerStatusMapIterator
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

// TPG (Nadir study)
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalTrigPrimCompactColl.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveSample.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"

//inter-calibration//calibration
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"

//Clusters
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

// For H/E - Iso on SC
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"

//RECHITS
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

//ECAL TT Masked
//#include "DQM/EcalCommon/interface/Numbers.h"
#include "CondFormats/DataRecord/interface/EcalTPGTowerStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTPGTowerStatus.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
//#include "DataFormats/DetId/interface/DetId.h"

// severity for spikes
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

//for L1ExtraParticles collection
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

// For RCT Em candidates
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"

// from trigger twiki
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#Trigger
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//l1trigger
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMask.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskTechTrigRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h"

//
// class declaration
//

using namespace std;
//using namespace reco;
using namespace edm;
//using namespace IPTools;

class OfflineSpikeCrystalToOnlineMatch: public edm::EDAnalyzer {
public:
	explicit OfflineSpikeCrystalToOnlineMatch(const edm::ParameterSet&);
	~OfflineSpikeCrystalToOnlineMatch();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;

	void fireL1(int L1noniso, int L1iso, vector<int> & firedEG,
			vector<int> EG);
	void globalFireL1_Normal(double eleRCT_eT, int noniso, int iso,
				 vector<int> & firedEG_N, vector<int> menu);

	void cout_L1Extra_candidate(l1extra::L1EmParticleCollection::const_iterator & it_l1extra);
	void cout_L1Extra_candidate(l1extra::L1EmParticle & it_l1extra);
	//	void cout_all_L1Extra_candidates(map<double, l1extra::L1EmParticleCollection::const_iterator> & map_candidates);
	void cout_all_L1Extra_candidates(edm::Handle< l1extra::L1EmParticleCollection > OnlineIso,
																	   edm::Handle< l1extra::L1EmParticleCollection > OnlineNonIso,
																	   edm::Handle< l1extra::L1EmParticleCollection > EmulatedIso,
																	   edm::Handle< l1extra::L1EmParticleCollection > EmulatedNonIso);

	// handles to get TPGs collections
	edm::Handle<EcalTrigPrimDigiCollection> * ecal_tp_;
	edm::Handle<EcalTrigPrimDigiCollection> * emulatorTP;

	std::string histogramFile_;
	edm::InputTag tpEmulatorCollection_;
	edm::InputTag tpOnlineCollection_;
	edm::InputTag EcalRecHitCollectionEB_;

	//spikes
	std::vector<EcalTrigTowerDetId> TrigTowersWithRecHitSpikes; // not used currently

	//Things defined in python configuration file
	const bool do_reconstruct_amplitude_allTTs_;
	const bool do_rechit_spikes_search_;
	const bool do_reconstruct_amplitudes_spikes_; //overrides do_rechit_spikes_search
	const bool do_3DSpike_plots_;
	const bool do_l1extraparticles_;
	const bool do_l1EG5Cut_;

	//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	// ----------member data ---------------------------

//	TTree * mytree;
//	//for the tree branches.
//	double ttonline_et;
//	double ttemulator_et;
//	double ttemulator_spikes_to_ttonline_match_et;
//	double rechits_et;
//	double rechit_spikes_et;
//	double rechit_spikes_match_to_ttonline_et;
//	double rechit_spikes_to_ttonline_et_sFGVB0On;
//	double rechit_spikes_to_ttonline_et_sFGVB1On;

//	double num_ttonline_et;
//	double num_ttemulator_et;
//	double num_ttemulator_spikes_to_ttonline_match_et;
//	double num_rechits_et;
//	double num_rechit_spikes_et;
//	double num_rechit_spikes_match_to_ttonline_et;
//	double num_rechit_spikes_to_ttonline_et_sFGVB0On;
//	double num_rechit_spikes_to_ttonline_et_sFGVB1On;
//
//	double ttonline_et;
//	double ttemulator_et;
//	double ttemulator_spikes_to_ttonline_match_et;
//	double rechits_et;
//	double rechit_spikes_et;
//	double rechit_spikes_match_to_ttonline_et;
//	double rechit_spikes_to_ttonline_et_sFGVB0On;
//	double rechit_spikes_to_ttonline_et_sFGVB1On;

//	double tt_amplitude_online;
//	double tt_spikes_amplitude_online;
//  double tt_amplitude_calibrated_online;
//  double tt_spikes_amplitude_calibrated_online;

//	int sev3_sev4_count=0;
//	int sev3_sev4_with_sFGVB1=0;
//	int sev3_sev4_with_sFGVB1_with_dataframe=0;

	//output file   (not used)
	TFile* file_;
	//histograms

//	TH2F * h_TPGsaturated;

	// hitograms to make LEGO2 distributions of energy inside Trigger Tower with offline spike.
	std::vector<TH2D*> hvect_Spikes; // ADC counts on z axis
	std::vector<TH2D*> hvect_Spikes_icalib; // (ADC counts)*icalib constant; on z axis
	std::vector<TH2D*> hvect_Spikes_icalib_ADCtoGeV; // GeV on Z axis
	int spiked_hists_number = 0;

	TH1F *h_nofevents;

	TH1F *h_ttonline_et;
	TH1F *h_ttonline_barrel_et;
	TH1F *h_ttonline_barrel_iphi;
	TH1F *h_ttonline_barrel_ieta;
	TH1F *h_ttonline_barrel_saturated_iphi;
	TH1F *h_ttonline_barrel_saturated_ieta;
	TH1F *h_ttonline_et_nozerocut;
	TH1F *h_ttonline_barrel_et_nozerocut;
	TH1F *h_ttonline_raw;
	TH1F *h_ttonline_et_sFGVB0On;
	TH1F *h_ttonline_et_sFGVB1On;
	TH1F *h_ttonline_et_sFGVB0Em;
	TH1F *h_ttonline_et_sFGVB1Em;
	TH2F *h_ttonline_et_against_ttemulator_et;
	TH2F *h_ttonline_et_against_ttemulator_et_nozero;
	TH2F *h_ttonline_et_against_ttemulator_et_noOnlineZero;
	TH2F *h_ttonline_et_against_ttemulator_et_noEmulatorZero;
	TH1F *h_ttonline_et_zeroInemulator;
	TH1F *h_ttonline_et_zeroInemulator_sFGVB0Em;
	TH1F *h_ttonline_et_zeroInemulator_sFGVB1Em;
	TH1F *h_ttemulator_et;
	TH1F *h_ttemulator_barrel_et;
	TH1F *h_ttemulator_barrel_iphi;
	TH1F *h_ttemulator_barrel_ieta;
	TH1F *h_ttemulator_barrel_saturated_iphi;
	TH1F *h_ttemulator_barrel_saturated_ieta;
	TH1F *h_ttemulator_et_nozerocut;
	TH1F *h_ttemulator_et_sFGVB0On;
	TH1F *h_ttemulator_et_sFGVB1On;
	TH1F *h_ttemulator_et_sFGVB0Em;
	TH1F *h_ttemulator_et_sFGVB1Em;
	TH1F *h_ttemulator_masked_et;
	TH1F *h_ttemulator_masked_et_nozerocut;
	TH1F *h_ttemulator_masked_et_sFGVB0On;
	TH1F *h_ttemulator_masked_et_sFGVB1On;
	TH1F *h_ttemulator_masked_et_sFGVB0Em;
	TH1F *h_ttemulator_masked_et_sFGVB1Em;
	TH1F *h_ttonline_saturated_to_ttoffline_match_et;

	TH1F *h_rechits_et;
	TH1F *h_rechit_spikes_et;
	TH1F *h_rechit_spikes_et_FoundInOnlineTTs;
	TH1F *h_rechit_spikes_et_FoundInEmulatorTTs;
	TH1F *h_rechit_spikes_et_sFGVB0On;
	TH1F *h_rechit_spikes_et_sFGVB1On;
	TH1F *h_rechit_spikes_et_sFGVB0Em;
	TH1F *h_rechit_spikes_et_sFGVB1Em;

	TH1F *h_rechits_energy;
	TH1F *h_rechit_spikes_energy;
	TH1F *h_rechit_spikes_energy_FoundInOnlineTTs;
	TH1F *h_rechit_spikes_energy_FoundInEmulatorTTs;
	TH1F *h_rechit_spikes_energy_sFGVB0On;
	TH1F *h_rechit_spikes_energy_sFGVB1On;
	TH1F *h_rechit_spikes_energy_sFGVB0Em;
	TH1F *h_rechit_spikes_energy_sFGVB1Em;

	TH1F *h_rechit_spikes_match_to_ttonline_et;
	TH1F *h_rechit_spikes_to_ttonline_et_sFGVB0On;
	TH1F *h_rechit_spikes_to_ttonline_et_sFGVB1On;
	TH1F *h_rechit_spikes_match_to_ttemulator_et;
	TH1F *h_rechit_spikes_to_ttemulator_et_sFGVB0Em;
	TH1F *h_rechit_spikes_to_ttemulator_et_sFGVB1Em;

	TH1F *h_rechit_spikes_match_to_ttonline_et_nozerocut;
	TH1F *h_rechit_spikes_to_ttonline_et_sFGVB0On_nozerocut;
	TH1F *h_rechit_spikes_to_ttonline_et_sFGVB1On_nozerocut;
	TH1F *h_rechit_spikes_match_to_ttemulator_et_nozerocut;
	TH1F *h_rechit_spikes_to_ttemulator_et_sFGVB0Em_nozerocut;
	TH1F *h_rechit_spikes_to_ttemulator_et_sFGVB1Em_nozerocut;

	TH1F *h_rechit_spikes_ttonline_et_zeroInEmulator;
	TH1F *h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0On;
	TH1F *h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1On;
	TH1F *h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0Em;
	TH1F *h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1Em;

	TH1F *h_tt_amplitude_online;
	TH1F *h_tt_amplitude_calibrated_online;

	TH1F *h_tt_amplitude_online_above_zero;
	TH1F *h_tt_amplitude_icalib_online_above_zero;
	TH1F *h_tt_amplitude_icalib_ADCtoGeV_online_above_zero;

	TH1F *h_tt_spikes_amplitude_online;
	TH1F *h_tt_spikes_amplitude_icalib_online;
	TH1F *h_tt_spikes_amplitude_icalib_ADCtoGeV_online;

	TH1F *h_tt_amplitude_icalib_ADCtoGeV_online_above_zero_sintheta;
	TH1F *h_tt_spikes_amplitude_icalib_ADCtoGeV_online_sintheta;

	TH2D *h_crystalIds_gain2_3;
	TH2D *h_masked_area;
	TH2D *h_TPGs_unmatch_online_offline;
	TH2D *h_TPGs_unmatch_online_offline_nonzero;
	TH2D *h_ttonline_barrel_map;
	TH2D *h_ttemulator_barrel_map;
	TH2D *h_ttonline_barrel_saturated_map;
	TH2D *h_ttemulator_barrel_saturated_map;
	TH2D *h_l1em_Online_Particles_map;
	TH2D *h_l1em_Online_Particles_index_map;
	TH2D *h_l1em_Online_Particles_Isolated_map;
	TH2D *h_l1em_Online_Particles_Isolated_index_map;
	TH2D *h_l1em_Online_Particles_NonIsolated_map;
	TH2D *h_l1em_Online_Particles_NonIsolated_index_map;
	TH2D *h_l1em_Emul_Particles_map;
	TH2D *h_l1em_Emul_Particles_index_map;
	TH2D *h_l1em_Emul_Particles_Isolated_map;
	TH2D *h_l1em_Emul_Particles_Isolated_index_map;
	TH2D *h_l1em_Emul_Particles_NonIsolated_map;
	TH2D *h_l1em_Emul_Particles_NonIsolated_index_map;

	TH2D *h_l1em_saturated_Online_Particles_index_map;
	TH2D *h_l1em_saturated_Emul_Particles_index_map;

	//L1 EM Particles
	TH1F *h_l1em_Emul_Particles_Isolated_et;
	TH1F *h_l1em_Emul_Particles_NonIsolated_et;
	TH1F *h_l1em_Online_Particles_Isolated_et;
	TH1F *h_l1em_Online_Particles_NonIsolated_et;
	TH1F *h_l1em_Emul_Particles_et;
	TH1F *h_l1em_Online_Particles_et;
	TH1F *h_l1em_Emul_Particles_eta1567_et;
	TH1F *h_l1em_Online_Particles_eta1567_et;
	TH1F *h_l1em_Emul_Particles_barrel_et;
	TH1F *h_l1em_Online_Particles_barrel_et;

	TH1F *h_l1em_Emul_Particles_noBadTowerCut_et;
	TH1F *h_l1em_Online_Particles_noBadTowerCut_et;
	TH1F *h_l1em_Emul_Particles_noBadTowerCut_eta1567_et;
	TH1F *h_l1em_Online_Particles_noBadTowerCut_eta1567_et;
	TH1F *h_l1em_Emul_Particles_noBadTowerCut_barrel_et;
	TH1F *h_l1em_Online_Particles_noBadTowerCut_barrel_et;

	TH1F *h_l1em_Emul_Particles_Highest_et;
	TH1F *h_l1em_Online_Particles_Highest_et;

	TH1F *h_l1em_Emul_Particles_Isolated_energy;
	TH1F *h_l1em_Emul_Particles_NonIsolated_energy;
	TH1F *h_l1em_Online_Particles_Isolated_energy;
	TH1F *h_l1em_Online_Particles_NonIsolated_energy;
	TH1F *h_l1em_Emul_Particles_energy;
	TH1F *h_l1em_Online_Particles_energy;

	TH1F *h_l1em_Emul_Particles_Highest_energy;
	TH1F *h_l1em_Online_Particles_Highest_energy;

	TH1F *h_l1em_highest_energies_ratio_not_equal;
	TH1F *h_l1em_highest_energies_ratio_equal;
	TH1F *h_l1em_highest_ets_ratio_not_equal;
	TH1F *h_l1em_highest_ets_ratio_equal;
	TH1F *h_l1em_highest_energies_ratio;
	TH1F *h_l1em_highest_ets_ratio;

	TH1F *h_l1Candides_N_Online;
	TH1F *h_l1Candides_N_Emul;
	TH1F *h_l1Candides_N_difference;

//	TH1F * h_temp_hist;

	//Creating ntuple
	edm::Service<TFileService> fs;

	//some geometry variables
	const CaloSubdetectorGeometry * theEndcapGeometry_;
	const CaloSubdetectorGeometry * theBarrelGeometry_;

	int spike_N;

	//weights for reconstruction of the signal amplitude in crystal from the ADC counts.
	double weight[10]; 	//online
	double weight_offline[10];


	//L1TriggerBit
	  edm::InputTag m_l1GTReadoutRecTag_;
//	  std::vector<int> l1Accepts_;
//	  std::vector<std::string> l1Names_;

	  //Trigger bits
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
//    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
//    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;


    // map to find RCT region ieta iphi (DetId) knowing TPG ieta iphi (DetId).
    map<EBDetId, L1CaloRegionDetId> map_DetId_RegionID;

	//Map to store info about the spikes
	//Map covers only barrel regions
	map<L1CaloRegionDetId, bool> map_RCT_regs_with_offline_spikes;


};

OfflineSpikeCrystalToOnlineMatch::OfflineSpikeCrystalToOnlineMatch(	const edm::ParameterSet& iConfig) :
				histogramFile_(iConfig.getParameter<std::string>("histogramFile")),
				tpEmulatorCollection_(iConfig.getParameter<edm::InputTag>("TPEmulatorCollection")),
				tpOnlineCollection_(iConfig.getParameter<edm::InputTag>("TPOnlineCollection")),
				do_reconstruct_amplitude_allTTs_(iConfig.getParameter<bool>("do_reconstruct_amplitude_allTTs")),
				do_rechit_spikes_search_(iConfig.getParameter<bool>("do_rechit_spikes_search")),
				do_reconstruct_amplitudes_spikes_(iConfig.getParameter<bool>("do_reconstruct_amplitudes_spikes")), //overrides do_rechit_spikes_search
				do_3DSpike_plots_(iConfig.getParameter<bool>("do_3DSpike_plots")),
				do_l1extraparticles_(iConfig.getParameter<bool>("do_l1extraparticles")),
				do_l1EG5Cut_(iConfig.getParameter<bool>("do_l1EG5Cut")),
				triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits")))
//				triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
//				triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
{
	//now do what ever initialization is needed
	EcalRecHitCollectionEB_ = (iConfig.getParameter<edm::InputTag>(	"EcalRecHitCollectionEB" ));

	//Histograms
	h_nofevents = fs->make<TH1F>("h_nofevents", "h_nofevents", 2, 0, 2);

	h_ttonline_et = fs->make<TH1F>("h_ttonline_et", "h_ttonline_et", 256, -0.5, 255.5);
	h_ttonline_et_nozerocut = fs->make<TH1F>("h_ttonline_et_nozerocut", "h_ttonline_et_nozerocut", 306, -50, 256);
	h_ttonline_barrel_et = fs->make<TH1F>("h_ttonline_barrel_et", "h_ttonline_barrel_et", 256, -0.5, 255.5);
	h_ttonline_barrel_ieta = fs->make<TH1F>("h_ttonline_barrel_ieta", "h_ttonline_barrel_ieta", 37,-18.5,18.5);
	h_ttonline_barrel_iphi = fs->make<TH1F>("h_ttonline_barrel_iphi", "h_ttonline_barrel_iphi", 72, 0.5, 72.5);
	h_ttonline_barrel_saturated_ieta = fs->make<TH1F>("h_ttonline_barrel_saturated_ieta", "h_ttonline_barrel_saturated_ieta", 37,-18.5,18.5);
	h_ttonline_barrel_saturated_iphi = fs->make<TH1F>("h_ttonline_barrel_saturated_iphi", "h_ttonline_barrel_saturated_iphi", 72, 0.5, 72.5);
	h_ttonline_barrel_et_nozerocut = fs->make<TH1F>("h_ttonline_barrel_et_nozerocut", "h_ttonline_barrel_et_nozerocut", 256, -0.5, 255.5);
	h_ttonline_et_sFGVB0On = fs->make<TH1F>("h_ttonline_et_sFGVB0On", "h_ttonline_et_sFGVB0On", 256, -0.5, 255.5);
	h_ttonline_et_sFGVB1On = fs->make<TH1F>("h_ttonline_et_sFGVB1On", "h_ttonline_et_sFGVB1On", 256, -0.5, 255.5);
	h_ttonline_et_sFGVB0Em = fs->make<TH1F>("h_ttonline_et_sFGVB0Em", "h_ttonline_et_sFGVB0Em", 256, -0.5, 255.5);
	h_ttonline_et_sFGVB1Em = fs->make<TH1F>("h_ttonline_et_sFGVB1Em", "h_ttonline_et_sFGVB1Em", 256, -0.5, 255.5);
	h_ttonline_raw = fs->make<TH1F>("h_ttonline_raw", "h_ttonline_raw", 256, -0.5, 255.5);
	h_ttonline_et_against_ttemulator_et = fs->make<TH2F>("h_ttonline_et_against_ttemulator_et", "h_ttonline_et_against_ttemulator_et", 256, -0.5, 255.5, 256, -0.5, 255.5);
	h_ttonline_et_against_ttemulator_et->GetXaxis()->SetTitle("onlineTT_et");
	h_ttonline_et_against_ttemulator_et->GetYaxis()->SetTitle("emulatedTT_et");

	h_ttonline_et_against_ttemulator_et_nozero = fs->make<TH2F>("h_ttonline_et_against_ttemulator_et_nozero", "h_ttonline_et_against_ttemulator_et_nozero", 256, -0.5, 255.5, 256, -0.5, 255.5);
	h_ttonline_et_against_ttemulator_et_nozero->GetXaxis()->SetTitle("onlineTT_et");
	h_ttonline_et_against_ttemulator_et_nozero->GetYaxis()->SetTitle("emulatedTT_et");

	h_ttonline_et_against_ttemulator_et_noOnlineZero = fs->make<TH2F>("h_ttonline_et_against_ttemulator_et_noOnlineZero", "h_ttonline_et_against_ttemulator_et_noOnlineZero", 256, -0.5, 255.5, 256, -0.5, 255.5);
	h_ttonline_et_against_ttemulator_et_noOnlineZero->GetXaxis()->SetTitle("onlineTT_et");
	h_ttonline_et_against_ttemulator_et_noOnlineZero->GetYaxis()->SetTitle("emulatedTT_et");

	h_ttonline_et_against_ttemulator_et_noEmulatorZero = fs->make<TH2F>("h_ttonline_et_against_ttemulator_et_noEmulatorZero", "h_ttonline_et_against_ttemulator_et_noEmulatorZero", 256, -0.5, 255.5, 256, -0.5, 255.5);
	h_ttonline_et_against_ttemulator_et_noEmulatorZero->GetXaxis()->SetTitle("onlineTT_et");
	h_ttonline_et_against_ttemulator_et_noEmulatorZero->GetYaxis()->SetTitle("emulatedTT_et");

	h_ttonline_et_zeroInemulator = fs->make<TH1F>("h_ttonline_et_zeroInemulator", "h_ttonline_et_zeroInemulator", 256, -0.5, 255.5);
	h_ttonline_et_zeroInemulator_sFGVB0Em = fs->make<TH1F>("h_ttonline_et_zeroInemulator_sFGVB0Em", "h_ttonline_et_zeroInemulator_sFGVB0Em", 256, -0.5, 255.5);
	h_ttonline_et_zeroInemulator_sFGVB1Em = fs->make<TH1F>("h_ttonline_et_zeroInemulator_sFGVB1Em", "h_ttonline_et_zeroInemulator_sFGVB1Em", 256, -0.5, 255.5);

	h_ttemulator_et = fs->make<TH1F>("h_ttemulator_et", "h_ttemulator_et", 256, -0.5, 255.5);
	h_ttemulator_barrel_et = fs->make<TH1F>("h_ttemulator_barrel_et", "h_ttemulator_barrel_et", 256, -0.5, 255.5);
	h_ttemulator_barrel_ieta = fs->make<TH1F>("h_ttemulator_barrel_ieta", "h_ttemulator_barrel_ieta", 37,-18.5,18.5);
	h_ttemulator_barrel_iphi = fs->make<TH1F>("h_ttemulator_barrel_iphi", "h_ttemulator_barrel_iphi", 72, 0.5, 72.5);
	h_ttemulator_barrel_saturated_ieta = fs->make<TH1F>("h_ttemulator_barrel_saturated_ieta", "h_ttemulator_barrel_saturated_ieta", 37,-18.5,18.5);
	h_ttemulator_barrel_saturated_iphi = fs->make<TH1F>("h_ttemulator_barrel_saturated_iphi", "h_ttemulator_barrel_saturated_iphi", 72, 0.5, 72.5);
	h_ttemulator_et_nozerocut = fs->make<TH1F>("h_ttemulator_et_nozerocut", "h_ttemulator_et_nozerocut", 256, -0.5, 255.5);
	h_ttemulator_et_sFGVB0On = fs->make<TH1F>("h_ttemulator_et_sFGVB0On", "h_ttemulator_et_sFGVB0On", 256, -0.5, 255.5);
	h_ttemulator_et_sFGVB1On = fs->make<TH1F>("h_ttemulator_et_sFGVB1On", "h_ttemulator_et_sFGVB1On", 256, -0.5, 255.5);
	h_ttemulator_et_sFGVB0Em = fs->make<TH1F>("h_ttemulator_et_sFGVB0Em", "h_ttemulator_et_sFGVB0Em", 256, -0.5, 255.5);
	h_ttemulator_et_sFGVB1Em = fs->make<TH1F>("h_ttemulator_et_sFGVB1Em", "h_ttemulator_et_sFGVB1Em", 256, -0.5, 255.5);
	h_ttemulator_masked_et = fs->make<TH1F>("h_ttemulator_masked_et", "h_ttemulator_masked_et", 256, -0.5, 255.5);
	h_ttemulator_masked_et_nozerocut = fs->make<TH1F>("h_ttemulator_masked_et_nozerocut", "h_ttemulator_masked_et_nozerocut", 256, -0.5, 255.5);
	h_ttemulator_masked_et_sFGVB0On = fs->make<TH1F>("h_ttemulator_masked_et_sFGVB0On", "h_ttemulator_masked_et_sFGVB0On", 256, -0.5, 255.5);
	h_ttemulator_masked_et_sFGVB1On = fs->make<TH1F>("h_ttemulator_masked_et_sFGVB1On", "h_ttemulator_masked_et_sFGVB1On", 256, -0.5, 255.5);
	h_ttemulator_masked_et_sFGVB0Em = fs->make<TH1F>("h_ttemulator_masked_et_sFGVB0Em", "h_ttemulator_masked_et_sFGVB0Em", 256, -0.5, 255.5);
	h_ttemulator_masked_et_sFGVB1Em = fs->make<TH1F>("h_ttemulator_masked_et_sFGVB1Em", "h_ttemulator_masked_et_sFGVB1Em", 256, -0.5, 255.5);
	h_ttonline_saturated_to_ttoffline_match_et = fs->make<TH1F>("h_ttonline_saturated_to_ttoffline_match_et",	"h_ttonline_saturated_to_ttoffline_match_et", 256, -0.5, 255.5);

	h_rechits_et = fs->make<TH1F>("h_rechits_et", "h_rechits_et", 256, -0.5, 255.5);
	h_rechit_spikes_et = fs->make<TH1F>("h_rechit_spikes_et", "h_rechit_spikes_et", 256, -0.5, 255.5);
	h_rechit_spikes_et_FoundInOnlineTTs = fs->make<TH1F>("h_rechit_spikes_et_FoundInOnlineTTs", "h_rechit_spikes_et_FoundInOnlineTTs", 256, -0.5, 255.5);
	h_rechit_spikes_et_sFGVB0On = fs->make<TH1F>("h_rechit_spikes_et_sFGVB0On", "h_rechit_spikes_et_sFGVB0On", 256, -0.5, 255.5);
	h_rechit_spikes_et_sFGVB1On = fs->make<TH1F>("h_rechit_spikes_et_sFGVB1On", "h_rechit_spikes_et_sFGVB1On", 256, -0.5, 255.5);
	h_rechit_spikes_et_FoundInEmulatorTTs = fs->make<TH1F>("h_rechit_spikes_et_FoundInEmulatorTTs", "h_rechit_spikes_et_FoundInEmulatorTTs", 256, -0.5, 255.5);
	h_rechit_spikes_et_sFGVB0Em = fs->make<TH1F>("h_rechit_spikes_et_sFGVB0Em", "h_rechit_spikes_et_sFGVB0Em", 256, -0.5, 255.5);
	h_rechit_spikes_et_sFGVB1Em = fs->make<TH1F>("h_rechit_spikes_et_sFGVB1Em", "h_rechit_spikes_et_sFGVB1Em", 256, -0.5, 255.5);

	h_rechits_energy = fs->make<TH1F>("h_rechits_energy", "h_rechits_energy", 256, -0.5, 255.5);
	h_rechit_spikes_energy = fs->make<TH1F>("h_rechit_spikes_energy", "h_rechit_spikes_energy", 256, -0.5, 255.5);
	h_rechit_spikes_energy_FoundInOnlineTTs = fs->make<TH1F>("h_rechit_spikes_energy_FoundInOnlineTTs", "h_rechit_spikes_energy_FoundInOnlineTTs", 256, -0.5, 255.5);
	h_rechit_spikes_energy_sFGVB0On = fs->make<TH1F>("h_rechit_spikes_energy_sFGVB0On", "h_rechit_spikes_energy_sFGVB0On", 256, -0.5, 255.5);
	h_rechit_spikes_energy_sFGVB1On = fs->make<TH1F>("h_rechit_spikes_energy_sFGVB1On", "h_rechit_spikes_energy_sFGVB1On", 256, -0.5, 255.5);
	h_rechit_spikes_energy_FoundInEmulatorTTs = fs->make<TH1F>("h_rechit_spikes_energy_FoundInEmulatorTTs", "h_rechit_spikes_energy_FoundInEmulatorTTs", 256, -0.5, 255.5);
	h_rechit_spikes_energy_sFGVB0Em = fs->make<TH1F>("h_rechit_spikes_energy_sFGVB0Em", "h_rechit_spikes_energy_sFGVB0Em", 256, -0.5, 255.5);
	h_rechit_spikes_energy_sFGVB1Em = fs->make<TH1F>("h_rechit_spikes_energy_sFGVB1Em", "h_rechit_spikes_energy_sFGVB1Em", 256, -0.5, 255.5);

	h_rechit_spikes_match_to_ttonline_et = fs->make<TH1F>("h_rechit_spikes_match_to_ttonline_et", "h_rechit_spikes_match_to_ttonline_et", 256, -0.5, 255.5);
	h_rechit_spikes_to_ttonline_et_sFGVB0On = fs->make<TH1F>("h_rechit_spikes_to_ttonline_et_sFGVB0On", "h_rechit_spikes_to_ttonline_et_sFGVB0On", 256, -0.5, 255.5);
	h_rechit_spikes_to_ttonline_et_sFGVB1On = fs->make<TH1F>("h_rechit_spikes_to_ttonline_et_sFGVB1On", "h_rechit_spikes_to_ttonline_et_sFGVB1On", 256, -0.5, 255.5);
	h_rechit_spikes_match_to_ttemulator_et = fs->make<TH1F>("h_rechit_spikes_match_to_ttemulator_et", "h_rechit_spikes_match_to_ttemulator_et", 256, -0.5, 255.5);
	h_rechit_spikes_to_ttemulator_et_sFGVB0Em = fs->make<TH1F>("h_rechit_spikes_to_ttemulator_et_sFGVB0Em", "h_rechit_spikes_to_ttemulator_et_sFGVB0Em", 256, -0.5, 255.5);
	h_rechit_spikes_to_ttemulator_et_sFGVB1Em = fs->make<TH1F>("h_rechit_spikes_to_ttemulator_et_sFGVB1Em", "h_rechit_spikes_to_ttemulator_et_sFGVB1Em", 256, -0.5, 255.5);
	h_rechit_spikes_match_to_ttonline_et_nozerocut = fs->make<TH1F>("h_rechit_spikes_match_to_ttonline_et_nozerocut", "h_rechit_spikes_match_to_ttonline_et_nozerocut", 256, -0.5, 255.5);
	h_rechit_spikes_to_ttonline_et_sFGVB0On_nozerocut = fs->make<TH1F>("h_rechit_spikes_to_ttonline_et_sFGVB0On_nozerocut", "h_rechit_spikes_to_ttonline_et_sFGVB0On_nozerocut", 256, -0.5, 255.5);
	h_rechit_spikes_to_ttonline_et_sFGVB1On_nozerocut = fs->make<TH1F>("h_rechit_spikes_to_ttonline_et_sFGVB1On_nozerocut", "h_rechit_spikes_to_ttonline_et_sFGVB1On_nozerocut", 256, -0.5, 255.5);
	h_rechit_spikes_match_to_ttemulator_et_nozerocut = fs->make<TH1F>("h_rechit_spikes_match_to_ttemulator_et_nozerocut", "h_rechit_spikes_match_to_ttemulator_et_nozerocut", 256, -0.5, 255.5);
	h_rechit_spikes_to_ttemulator_et_sFGVB0Em_nozerocut = fs->make<TH1F>("h_rechit_spikes_to_ttemulator_et_sFGVB0Em_nozerocut", "h_rechit_spikes_to_ttemulator_et_sFGVB0Em_nozerocut", 256, -0.5, 255.5);
	h_rechit_spikes_to_ttemulator_et_sFGVB1Em_nozerocut = fs->make<TH1F>("h_rechit_spikes_to_ttemulator_et_sFGVB1Em_nozerocut", "h_rechit_spikes_to_ttemulator_et_sFGVB1Em_nozerocut", 256, -0.5, 255.5);
//	h_TPGsaturated = new TH2F("h_TPGsaturated", "h_TPGsaturated", 60, -30, 30, 72, 0, 72);

	h_rechit_spikes_ttonline_et_zeroInEmulator = fs->make<TH1F>("h_rechit_spikes_ttonline_et_zeroInEmulator", "h_rechit_spikes_ttonline_et_zeroInEmulator", 256, -0.5, 255.5);
	h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0On = fs->make<TH1F>("h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0On", "h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0On", 256, -0.5, 255.5);
	h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1On = fs->make<TH1F>("h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1On", "h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1On", 256, -0.5, 255.5);
	h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0Em = fs->make<TH1F>("h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0Em", "h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0Em", 256, -0.5, 255.5);
	h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1Em = fs->make<TH1F>("h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1Em", "h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1Em", 256, -0.5, 255.5);

	h_tt_amplitude_online = fs->make<TH1F>("h_tt_amplitude_online", "h_tt_amplitude_online", 256, 0, 256);
	h_tt_amplitude_calibrated_online = fs->make<TH1F>("h_tt_amplitude_calibrated_online", "h_tt_amplitude_calibrated_online", 256, 0, 256);

	h_tt_amplitude_online_above_zero = fs->make<TH1F>("h_tt_amplitude_online_above_zero", "h_tt_amplitude_online_above_zero", 256, 0, 256);
	h_tt_amplitude_icalib_online_above_zero = fs->make<TH1F>("h_tt_amplitude_icalib_online_above_zero", "h_tt_amplitude_icalib_online_above_zero", 256, 0, 256);
	h_tt_amplitude_icalib_ADCtoGeV_online_above_zero = fs->make<TH1F>("h_tt_amplitude_icalib_ADCtoGeV_online_above_zero", "h_tt_amplitude_icalib_ADCtoGeV_online_above_zero", 256, 0, 256);

	h_tt_spikes_amplitude_online = fs->make<TH1F>("h_tt_spikes_amplitude_online", "h_tt_spikes_amplitude_online", 256, 0, 256);
	h_tt_spikes_amplitude_icalib_online = fs->make<TH1F>("h_tt_spikes_amplitude_icalib_online", "h_tt_spikes_amplitude_icalib_online", 256, 0, 256);
	h_tt_spikes_amplitude_icalib_ADCtoGeV_online = fs->make<TH1F>("h_tt_spikes_amplitude_icalib_ADCtoGeV_online", "h_tt_spikes_amplitude_icalib_ADCtoGeV_online", 256, 0, 256);

	h_tt_amplitude_icalib_ADCtoGeV_online_above_zero_sintheta = fs->make<TH1F>("h_tt_amplitude_icalib_ADCtoGeV_online_above_zero_sintheta", "h_tt_amplitude_icalib_ADCtoGeV_online_above_zero_sintheta", 256, 0, 256);
	h_tt_spikes_amplitude_icalib_ADCtoGeV_online_sintheta = fs->make<TH1F>("h_tt_spikes_amplitude_icalib_ADCtoGeV_online_sintheta", "h_tt_spikes_amplitude_icalib_ADCtoGeV_online_sintheta", 256, 0, 256);

	h_crystalIds_gain2_3 = fs->make<TH2D>("h_crystalIds_gain2_3","h_crystalIds_gain2_3", 171, -85.5, 85.5, 360, -0.5, 359.5);
	h_crystalIds_gain2_3->GetXaxis()->SetTitle("ieta");
	h_crystalIds_gain2_3->GetYaxis()->SetTitle("iphi");

	h_masked_area=fs->make<TH2D>("h_masked_area","h_masked_area", 37,-18.5,18.5,72, 0.5, 72.5); //TT masked area
	h_masked_area->GetXaxis()->SetTitle("ieta");
	h_masked_area->GetYaxis()->SetTitle("iphi");

	h_TPGs_unmatch_online_offline=fs->make<TH2D>("h_TPGs_unmatch_online_offline","h_TPGs_unmatch_online_offline", 37,-18.5,18.5,72, 0.5, 72.5); //TT masked area
	h_TPGs_unmatch_online_offline->GetXaxis()->SetTitle("ieta");
	h_TPGs_unmatch_online_offline->GetYaxis()->SetTitle("iphi");

	h_TPGs_unmatch_online_offline_nonzero=fs->make<TH2D>("h_TPGs_unmatch_online_offline_nonzero","h_TPGs_unmatch_online_offline_nonzero", 37,-18.5,18.5,72, 0.5, 72.5);
	h_TPGs_unmatch_online_offline_nonzero->GetXaxis()->SetTitle("ieta");
	h_TPGs_unmatch_online_offline_nonzero->GetYaxis()->SetTitle("iphi");

	h_ttonline_barrel_map=fs->make<TH2D>("h_ttonline_barrel_map","h_ttonline_barrel_map", 37,-18.5,18.5,72, 0.5, 72.5);
	h_ttonline_barrel_map->GetXaxis()->SetTitle("ieta");
	h_ttonline_barrel_map->GetYaxis()->SetTitle("iphi");

	h_ttemulator_barrel_map=fs->make<TH2D>("h_ttemulator_barrel_map","h_ttemulator_barrel_map", 37,-18.5,18.5,72, 0.5, 72.5);
	h_ttemulator_barrel_map->GetXaxis()->SetTitle("ieta");
	h_ttemulator_barrel_map->GetYaxis()->SetTitle("iphi");

	h_ttonline_barrel_saturated_map=fs->make<TH2D>("h_ttonline_barrel_saturated_map","h_ttonline_barrel_saturated_map", 37,-18.5,18.5,72, 0.5, 72.5);
	h_ttonline_barrel_saturated_map->GetXaxis()->SetTitle("ieta");
	h_ttonline_barrel_saturated_map->GetYaxis()->SetTitle("iphi");

	h_ttemulator_barrel_saturated_map=fs->make<TH2D>("h_ttemulator_barrel_saturated_map","h_ttemulator_barrel_saturated_map", 37,-18.5,18.5,72, 0.5, 72.5);
	h_ttemulator_barrel_saturated_map->GetXaxis()->SetTitle("ieta");
	h_ttemulator_barrel_saturated_map->GetYaxis()->SetTitle("iphi");

	h_l1em_Online_Particles_map=fs->make<TH2D>("h_l1em_Online_Particles_map","h_l1em_Online_Particles_map", 28,-2.7707,2.7707,36,-(TMath::Pi()), TMath::Pi() );
	h_l1em_Online_Particles_map->GetXaxis()->SetTitle("eta");
	h_l1em_Online_Particles_map->GetYaxis()->SetTitle("phi");

	h_l1em_Online_Particles_index_map=fs->make<TH2D>("h_l1em_Online_Particles_index_map","h_l1em_Online_Particles_index_map", 14,3.5,17.5, 18,-0.5,17.5);
	h_l1em_Online_Particles_index_map->GetXaxis()->SetTitle("ieta");
	h_l1em_Online_Particles_index_map->GetYaxis()->SetTitle("iphi");

	h_l1em_Online_Particles_Isolated_map=fs->make<TH2D>("h_l1em_Online_Particles_Isolated_map","h_l1em_Online_Particles_Isolated_map", 28,-2.7707,2.7707,36,-(TMath::Pi()), TMath::Pi() );
	h_l1em_Online_Particles_Isolated_map->GetXaxis()->SetTitle("eta");
	h_l1em_Online_Particles_Isolated_map->GetYaxis()->SetTitle("phi");

	h_l1em_Online_Particles_Isolated_index_map=fs->make<TH2D>("h_l1em_Online_Particles_Isolated_index_map","h_l1em_Online_Particles_Isolated_index_map", 14,3.5,17.5, 18,-0.5,17.5);
	h_l1em_Online_Particles_Isolated_index_map->GetXaxis()->SetTitle("ieta");
	h_l1em_Online_Particles_Isolated_index_map->GetYaxis()->SetTitle("iphi");

	h_l1em_Online_Particles_NonIsolated_map=fs->make<TH2D>("h_l1em_Online_Particles_NonIsolated_map","h_l1em_Online_Particles_NonIsolated_map", 28,-2.7707,2.7707,36,-(TMath::Pi()), TMath::Pi() );
	h_l1em_Online_Particles_NonIsolated_map->GetXaxis()->SetTitle("eta");
	h_l1em_Online_Particles_NonIsolated_map->GetYaxis()->SetTitle("phi");

	h_l1em_Online_Particles_NonIsolated_index_map=fs->make<TH2D>("h_l1em_Online_Particles_NonIsolated_index_map","h_l1em_Online_Particles_NonIsolated_index_map", 14,3.5,17.5, 18,-0.5,17.5);
	h_l1em_Online_Particles_NonIsolated_index_map->GetXaxis()->SetTitle("ieta");
	h_l1em_Online_Particles_NonIsolated_index_map->GetYaxis()->SetTitle("iphi");

	h_l1em_Emul_Particles_map=fs->make<TH2D>("h_l1em_Emul_Particles_map","h_l1em_Emul_Particles_map", 28,-2.7707,2.7707,36,-(TMath::Pi()), TMath::Pi() );
	h_l1em_Emul_Particles_map->GetXaxis()->SetTitle("eta");
	h_l1em_Emul_Particles_map->GetYaxis()->SetTitle("phi");

	h_l1em_Emul_Particles_index_map=fs->make<TH2D>("h_l1em_Emul_Particles_index_map","h_l1em_Emul_Particles_index_map", 14,3.5,17.5, 18,-0.5,17.5);
	h_l1em_Emul_Particles_index_map->GetXaxis()->SetTitle("ieta");
	h_l1em_Emul_Particles_index_map->GetYaxis()->SetTitle("iphi");

	h_l1em_Emul_Particles_Isolated_map=fs->make<TH2D>("h_l1em_Emul_Particles_Isolated_map","h_l1em_Emul_Particles_Isolated_map", 28,-2.7707,2.7707,36,-(TMath::Pi()), TMath::Pi() );
	h_l1em_Emul_Particles_Isolated_map->GetXaxis()->SetTitle("eta");
	h_l1em_Emul_Particles_Isolated_map->GetYaxis()->SetTitle("phi");

	h_l1em_Emul_Particles_Isolated_index_map=fs->make<TH2D>("h_l1em_Emul_Particles_Isolated_index_map","h_l1em_Emul_Particles_Isolated_index_map", 14,3.5,17.5, 18,-0.5,17.5);
	h_l1em_Emul_Particles_Isolated_index_map->GetXaxis()->SetTitle("ieta");
	h_l1em_Emul_Particles_Isolated_index_map->GetYaxis()->SetTitle("iphi");

	h_l1em_Emul_Particles_NonIsolated_map=fs->make<TH2D>("h_l1em_Emul_Particles_NonIsolated_map","h_l1em_Emul_Particles_NonIsolated_map", 28,-2.7707,2.7707,36,-(TMath::Pi()), TMath::Pi() );
	h_l1em_Emul_Particles_NonIsolated_map->GetXaxis()->SetTitle("eta");
	h_l1em_Emul_Particles_NonIsolated_map->GetYaxis()->SetTitle("phi");

	h_l1em_Emul_Particles_NonIsolated_index_map=fs->make<TH2D>("h_l1em_Emul_Particles_NonIsolated_index_map","h_l1em_Emul_Particles_NonIsolated_index_map", 14,3.5,17.5, 18,-0.5,17.5);
	h_l1em_Emul_Particles_NonIsolated_index_map->GetXaxis()->SetTitle("ieta");
	h_l1em_Emul_Particles_NonIsolated_index_map->GetYaxis()->SetTitle("iphi");

	h_l1em_saturated_Online_Particles_index_map=fs->make<TH2D>("h_l1em_saturated_Online_Particles_index_map","h_l1em_saturated_Online_Particles_index_map", 14,3.5,17.5, 18,-0.5,17.5);
	h_l1em_saturated_Online_Particles_index_map->GetXaxis()->SetTitle("ieta");
	h_l1em_saturated_Online_Particles_index_map->GetYaxis()->SetTitle("iphi");

	h_l1em_saturated_Emul_Particles_index_map=fs->make<TH2D>("h_l1em_saturated_Emul_Particles_index_map","h_l1em_saturated_Emul_Particles_index_map", 14,3.5,17.5, 18,-0.5,17.5);
	h_l1em_saturated_Emul_Particles_index_map->GetXaxis()->SetTitle("ieta");
	h_l1em_saturated_Emul_Particles_index_map->GetYaxis()->SetTitle("iphi");

	//L1 EM Particles
	h_l1em_Emul_Particles_Isolated_et = fs->make<TH1F>("h_l1em_Emul_Particles_Isolated_et", "h_l1em_Emul_Particles_Isolated_et", 256, 0, 256);
	h_l1em_Emul_Particles_NonIsolated_et = fs->make<TH1F>("h_l1em_Emul_Particles_NonIsolated_et", "h_l1em_Emul_Particles_NonIsolated_et", 256, 0, 256);
	h_l1em_Online_Particles_Isolated_et = fs->make<TH1F>("h_l1em_Online_Particles_Isolated_et", "h_l1em_Online_Particles_Isolated_et", 256, 0, 256);
	h_l1em_Online_Particles_NonIsolated_et = fs->make<TH1F>("h_l1em_Online_Particles_NonIsolated_et", "h_l1em_Online_Particles_NonIsolated_et", 256, 0, 256);
	h_l1em_Emul_Particles_et = fs->make<TH1F>("h_l1em_Emul_Particles_et", "h_l1em_Emul_Particles_et", 256, 0, 256);
	h_l1em_Online_Particles_et = fs->make<TH1F>("h_l1em_Online_Particles_et", "h_l1em_Online_Particles_et", 256, 0, 256);
	h_l1em_Emul_Particles_eta1567_et = fs->make<TH1F>("h_l1em_Emul_Particles_eta1567_et", "h_l1em_Emul_Particles_eta1567_et", 256, 0, 256);
	h_l1em_Online_Particles_eta1567_et = fs->make<TH1F>("h_l1em_Online_Particles_eta1567_et", "h_l1em_Online_Particles_eta1567_et", 256, 0, 256);
	h_l1em_Emul_Particles_barrel_et = fs->make<TH1F>("h_l1em_Emul_Particles_barrel_et", "h_l1em_Emul_Particles_barrel_et", 256, 0, 256);
	h_l1em_Online_Particles_barrel_et = fs->make<TH1F>("h_l1em_Online_Particles_barrel_et", "h_l1em_Online_Particles_barrel_et", 256, 0, 256);
	h_l1em_Emul_Particles_Highest_et = fs->make<TH1F>("h_l1em_Emul_Particles_Highest_et", "h_l1em_Emul_Particles_Highest_et", 256, 0, 256);
	h_l1em_Online_Particles_Highest_et = fs->make<TH1F>("h_l1em_Online_Particles_Highest_et", "h_l1em_Online_Particles_Highest_et", 256, 0, 256);

	h_l1em_Emul_Particles_noBadTowerCut_et = fs->make<TH1F>("h_l1em_Emul_Particles_noBadTowerCut_et", "h_l1em_Emul_Particles_noBadTowerCut_et", 256, 0, 256);
	h_l1em_Online_Particles_noBadTowerCut_et = fs->make<TH1F>("h_l1em_Online_Particles_noBadTowerCut_et", "h_l1em_Online_Particles_noBadTowerCut_et", 256, 0, 256);
	h_l1em_Emul_Particles_noBadTowerCut_eta1567_et = fs->make<TH1F>("h_l1em_Emul_Particles_noBadTowerCut_eta1567_et", "h_l1em_Emul_Particles_noBadTowerCut_eta1567_et", 256, 0, 256);
	h_l1em_Online_Particles_noBadTowerCut_eta1567_et = fs->make<TH1F>("h_l1em_Online_Particles_noBadTowerCut_eta1567_et", "h_l1em_Online_Particles_noBadTowerCut_eta1567_et", 256, 0, 256);
	h_l1em_Emul_Particles_noBadTowerCut_barrel_et = fs->make<TH1F>("h_l1em_Emul_Particles_noBadTowerCut_barrel_et", "h_l1em_Emul_Particles_noBadTowerCut_barrel_et", 256, 0, 256);
	h_l1em_Online_Particles_noBadTowerCut_barrel_et = fs->make<TH1F>("h_l1em_Online_Particles_noBadTowerCut_barrel_et", "h_l1em_Online_Particles_noBadTowerCut_barrel_et", 256, 0, 256);


	//L1 EM Particles
	h_l1em_Emul_Particles_Isolated_energy = fs->make<TH1F>("h_l1em_Emul_Particles_Isolated_energy", "h_l1em_Emul_Particles_Isolated_energy", 256, 0, 256);
	h_l1em_Emul_Particles_NonIsolated_energy = fs->make<TH1F>("h_l1em_Emul_Particles_NonIsolated_energy", "h_l1em_Emul_Particles_NonIsolated_energy", 256, 0, 256);
	h_l1em_Online_Particles_Isolated_energy = fs->make<TH1F>("h_l1em_Online_Particles_Isolated_energy", "h_l1em_Online_Particles_Isolated_energy", 256, 0, 256);
	h_l1em_Online_Particles_NonIsolated_energy = fs->make<TH1F>("h_l1em_Online_Particles_NonIsolated_energy", "h_l1em_Online_Particles_NonIsolated_energy", 256, 0, 256);
	h_l1em_Emul_Particles_energy = fs->make<TH1F>("h_l1em_Emul_Particles_energy", "h_l1em_Emul_Particles_energy", 256, 0, 256);
	h_l1em_Online_Particles_energy = fs->make<TH1F>("h_l1em_Online_Particles_energy", "h_l1em_Online_Particles_energy", 256, 0, 256);

	h_l1em_Emul_Particles_Highest_energy = fs->make<TH1F>("h_l1em_Emul_Particles_Highest_energy", "h_l1em_Emul_Particles_Highest_energy", 256, 0, 256);
	h_l1em_Online_Particles_Highest_energy = fs->make<TH1F>("h_l1em_Online_Particles_Highest_energy", "h_l1em_Online_Particles_Highest_energy", 256, 0, 256);

	h_l1em_highest_ets_ratio_not_equal = fs->make<TH1F>("h_l1em_highest_ets_ratio_not_equal", "h_l1em_highest_ets_ratio_not_equal", 256, 0-1/512, 1+1/512);
	h_l1em_highest_ets_ratio_equal = fs->make<TH1F>("h_l1em_highest_ets_ratio_equal", "h_l1em_highest_ets_ratio_equal", 256, 0-1/512, 1+1/512);
	h_l1em_highest_energies_ratio_not_equal = fs->make<TH1F>("h_l1em_highest_energies_ratio_not_equal", "h_l1em_highest_energies_ratio_not_equal", 256, 0-1/512, 1+1/512);
	h_l1em_highest_energies_ratio_equal = fs->make<TH1F>("h_l1em_highest_energies_ratio_equal", "h_l1em_highest_energies_ratio_equal", 256, 0-1/512, 1+1/512);
	h_l1em_highest_energies_ratio = fs->make<TH1F>("h_l1em_highest_energies_ratio", "h_l1em_highest_energies_ratio", 256, 0-1/512, 1+1/512);
	h_l1em_highest_ets_ratio = fs->make<TH1F>("h_l1em_highest_ets_ratio", "h_l1em_highest_ets_ratio", 256, 0-1/512, 1+1/512);


	h_l1Candides_N_Online = fs->make<TH1F>("h_l1Candides_N_Online", "h_l1Candides_N_Online", 9, -0.5, 8.5);
	h_l1Candides_N_Emul = fs->make<TH1F>("h_l1Candides_N_Emul", "h_l1Candides_N_Emul", 9, -0.5, 8.5);
	h_l1Candides_N_difference = fs->make<TH1F>("h_l1Candides_N_difference", "h_l1Candides_N_difference", 17, -8.5, 8.5);

	//L1TrigerBit
	m_l1GTReadoutRecTag_     = iConfig.getParameter<edm::InputTag>("L1GlobalReadoutRecord");

	//weights for reconstruction of the signal amplitude in crystal from the ADC counts.
	//online
	weight[0] = 0.0;
	weight[1] = 0.0;
	weight[2] = -0.560599;
	weight[3] = -0.54671;
	weight[4] = 0.246449;
	weight[5] = 0.491237;
	weight[6] = 0.369624;
	weight[7] = 0.0;
	weight[8] = 0.0;
	weight[9] = 0.0;
	//offline
	weight_offline[0] = -0.3812788;
	weight_offline[1] = -0.3812788;
	weight_offline[2] = -0.3812788;
	weight_offline[3] = 0.0;
	weight_offline[4] = 0.235699;
	weight_offline[5] = 0.4228363;
	weight_offline[6] = 0.3298652;
	weight_offline[7] = 0.1575187;
	weight_offline[8] = -0.002082776;
	weight_offline[9] = 0.0;

	// Create map EBDetID -> RCT_regionId
	// Code tested only for Barrel TPGs   18 < ieta < 18
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/RCTMap

	for (int tow_iEta=-32; tow_iEta < 33; tow_iEta++){
			if (tow_iEta==0) continue;
			for (int tow_iPhi=1; tow_iPhi < 73; tow_iPhi++){
				EBDetId tow_id(tow_iEta, tow_iPhi);

				unsigned int reg_ieta, reg_iphi;

				reg_ieta = floor ( (double) tow_iEta / 4) + 11;

				if ( (tow_iPhi == 71) || (tow_iPhi == 72) ){
					reg_iphi = 0;
				} else {
					reg_iphi = floor ( (double) (tow_iPhi + 1) / 4);
				}

				L1CaloRegionDetId reg_id = L1CaloRegionDetId(reg_ieta, reg_iphi);
				map_DetId_RegionID[tow_id] = reg_id;
			}
		}



//	h_temp_hist = fs->make<TH1F>("h_temp_hist", "debug_hist", 200, 0, -1);

//	mytree = fs->make<TTree>("SpikeHistsTree", "SpikeHistsTree");
//
//	mytree->Branch("ttonline_et", &ttonline_et, "ttonline_et/D");
//	mytree->Branch("ttemulator_et", &ttemulator_et, "ttemulator_et/D");
//	mytree->Branch("ttemulator_spikes_to_ttonline_match_et",
//			&ttemulator_spikes_to_ttonline_match_et,
//			"ttemulator_spikes_to_ttonline_match_et/D");
//	mytree->Branch("rechits_et", &rechits_et, "rechits_et/D");
//	mytree->Branch("rechit_spikes_et", &rechit_spikes_et, "rechit_spikes_et/D");
//	mytree->Branch("rechit_spikes_match_to_ttonline_et",
//			&rechit_spikes_match_to_ttonline_et,
//			"rechit_spikes_match_to_ttonline_et/D");
//	mytree->Branch("rechit_spikes_to_ttonline_et_sFGVB0On",
//			&rechit_spikes_to_ttonline_et_sFGVB0On,
//			"rechit_spikes_to_ttonline_et_sFGVB0On/D");
//	mytree->Branch("rechit_spikes_to_ttonline_et_sFGVB1On",
//			&rechit_spikes_to_ttonline_et_sFGVB1On,
//			"rechit_spikes_to_ttonline_et_sFGVB1On/D");
//	mytree->Branch("rechit_spikes_to_ttonline_et_sFGVB1On",
//				&rechit_spikes_to_ttonline_et_sFGVB1On,
//				"rechit_spikes_to_ttonline_et_sFGVB1On/D");

//	mytree->Branch("tt_amplitude_online",
//					&tt_amplitude_online,
//					"tt_amplitude_online/D");
//	mytree->Branch("tt_spikes_amplitude_online",
//						&tt_spikes_amplitude_online,
//						"tt_spikes_amplitude_online/D");
//	mytree->Branch("tt_amplitude_calibrated_online",
//						&tt_amplitude_calibrated_online,
//						"tt_amplitude_calibrated_online/D");
//	mytree->Branch("tt_spikes_amplitude_calibrated_online",
//						&tt_spikes_amplitude_calibrated_online,
//						"tt_spikes_amplitude_calibrated_online/D");

//	mytree->Branch("sev3_sev4_count", &sev3_sev4_count,	"sev3_sev4_count/I");
//	mytree->Branch("sev3_sev4_with_sFGVB1",	&sev3_sev4_with_sFGVB1, "sev3_sev4_with_sFGVB1/I");
//	mytree->Branch("sev3_sev4_with_sFGVB1_with_dataframe",	&sev3_sev4_with_sFGVB1_with_dataframe, "sev3_sev4_with_sFGVB1_with_dataframe/I");
	//mytree->Branch("rechit_spikes_to_ttonline_et_sFGVB1On",&rechit_spikes_to_ttonline_et_sFGVB1On,"rechit_spikes_to_ttonline_et_sFGVB1On/D");

}

OfflineSpikeCrystalToOnlineMatch::~OfflineSpikeCrystalToOnlineMatch() {
	// do anything here that needs to be done at destruction time
	// (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void OfflineSpikeCrystalToOnlineMatch::analyze(const edm::Event& iEvent,
		const edm::EventSetup& iSetup) {

	h_nofevents->Fill(1);

	//ONLINE TPs (Alex's code for Iurii the best)
	//Get handles for online TPGs
	edm::Handle<EcalTrigPrimDigiCollection> onlineTP;
	iEvent.getByLabel(tpOnlineCollection_, onlineTP);
	EcalTrigPrimDigiCollection const &onlineTPs = *(onlineTP.product());
	EcalTrigPrimDigiCollection::const_iterator onTP_it;

	//Get handles for emulated TPGs
	edm::Handle<EcalTrigPrimDigiCollection> emulatorTP;
	iEvent.getByLabel(tpEmulatorCollection_, emulatorTP);
	EcalTrigPrimDigiCollection const &emulatorTPs = *(emulatorTP.product());
	EcalTrigPrimDigiCollection::const_iterator ofTP_it;

	//ECAL Masked TriggerTowers
	vector<int> MaskEta;
	vector<int> MaskPhi;

	edm::ESHandle<EcalTPGTowerStatus> theEcalTPGTowerStatus_handle;
	iSetup.get<EcalTPGTowerStatusRcd>().get(theEcalTPGTowerStatus_handle);
	const EcalTPGTowerStatus * ecaltpgTowerStatus=theEcalTPGTowerStatus_handle.product();

	const EcalTPGTowerStatusMap &towerMap=ecaltpgTowerStatus->getMap();
	EcalTPGTowerStatusMapIterator  it;

	uint nMaskedChannels = 0;

	//Create map of masked trigger tower
	int MaskedFlagTTs[towerMap.size()];
	for (it=towerMap.begin(); it!=towerMap.end(); ++it) {
	//const EcalTrigTowerDetId  towerId = (*it).first;
	//cout << towerId << " " << (*it).second << endl;
	if ((*it).second > 0)
	  {
		const EcalTrigTowerDetId  ttId((*it).first);
//		std::cout << " iphi " <<  ttId.iphi()  << std::endl;
//		std::cout <<  "masked towers iEta==" <<  ttId.ieta()<< " iPhi=" << ttId.iphi()   << std::endl;
		MaskEta.push_back(ttId.ieta());
		MaskPhi.push_back(ttId.iphi());
		h_masked_area->Fill(ttId.ieta(),ttId.iphi());
		nMaskedChannels++;
//		 EcalAux_.nMasked++;

	  }//masked
	}//loop towers

	  //// Loop through emulated TPGs
	  //if (print_) std::cout<<"TPEmulator collection size="<<tpEmul.product()->size()<<std::endl ;

	int i_emTPs=0;
	for ( auto & emTP_it : emulatorTPs) {
		double EmulatorCompressedEt;
//		EmulatorCompressedEt = (emTP_it)[2].raw() & 0xff;
		EmulatorCompressedEt = emTP_it.compressedEt();
		int emulatorsFGVB = emTP_it.sFGVB();
		EcalTrigTowerDetId towerid = emTP_it.id();
		int towerieta( towerid.ieta() );
		int toweriphi( towerid.iphi() );

		//Get rid of noisy towers. TEMPORARY!!!
		if ( (towerieta == 14)  &&  ( toweriphi == 28 ) ) continue;
		if ( (towerieta == -8)  &&  ( toweriphi == 66 ) ) continue;

		//Get masked flag

		int MaskedFlag;
		auto pair = towerMap.find(emTP_it.id());
		if (pair!=towerMap.end()){
			MaskedFlag=pair->second;
		} else{
			MaskedFlag=0;
		}
		MaskedFlagTTs[i_emTPs]=MaskedFlag;
		i_emTPs++;

		h_ttemulator_et_nozerocut->Fill(EmulatorCompressedEt);
		if (MaskedFlag<=0) h_ttemulator_masked_et_nozerocut->Fill(EmulatorCompressedEt);

		if (EmulatorCompressedEt > 0) {
//			ttemulator_et = EmulatorCompressedEt;
			h_ttemulator_et->Fill(EmulatorCompressedEt);
			if (MaskedFlag<=0) h_ttemulator_masked_et->Fill(EmulatorCompressedEt);
			//find matching tower in the onlineTPGs collection
			EcalTrigPrimDigiCollection::const_iterator onlineMatch_it;
			onlineMatch_it = onlineTPs.find(emTP_it.id());

			if ( ( -18 < towerieta) && (towerieta < 18) ) {
							h_ttemulator_barrel_et->Fill(EmulatorCompressedEt);
							h_ttemulator_barrel_ieta->Fill(towerieta);
							h_ttemulator_barrel_iphi->Fill(toweriphi);
							h_ttemulator_barrel_map->Fill(towerieta, toweriphi);

							if (EmulatorCompressedEt == 255) {
								h_ttemulator_barrel_saturated_ieta->Fill(towerieta);
								h_ttemulator_barrel_saturated_iphi->Fill(toweriphi);
								h_ttemulator_barrel_saturated_map->Fill(towerieta, toweriphi);
							}
			}

			if (emulatorsFGVB == 0) {
				h_ttemulator_et_sFGVB0Em->Fill(EmulatorCompressedEt);
				if (MaskedFlag<=0) h_ttemulator_masked_et_sFGVB0Em->Fill(EmulatorCompressedEt);
			}
			if (emulatorsFGVB == 1)	{
				h_ttemulator_et_sFGVB1Em->Fill(EmulatorCompressedEt);
				if (MaskedFlag<=0) h_ttemulator_masked_et_sFGVB1Em->Fill(EmulatorCompressedEt);
			}
			if (onlineMatch_it!=onlineTPs.end()){
				double OnlineMatchCompressedEt;
				OnlineMatchCompressedEt = onlineMatch_it->compressedEt();

				if (OnlineMatchCompressedEt > 0) {
					if (emulatorsFGVB == 0) h_ttonline_et_sFGVB0Em->Fill(OnlineMatchCompressedEt);
					if (emulatorsFGVB == 1)	h_ttonline_et_sFGVB1Em->Fill(OnlineMatchCompressedEt);
				}

			} // if online match for the TT was found

		}// if EmulatorCompressedEt > 0

		// Hunt for hot tower
		if ( (EmulatorCompressedEt == 132) ||
			 (EmulatorCompressedEt == 48)  ||
			 (EmulatorCompressedEt == 192)  ||
			 (EmulatorCompressedEt == 240)  ||
			 (EmulatorCompressedEt == 248)  ||
			 (EmulatorCompressedEt == 255)
				) {
			cout << emTP_it;

			EcalTrigPrimDigiCollection::const_iterator onlineMatch_it;
			onlineMatch_it = onlineTPs.find(emTP_it.id());

			cout << "Online:" << endl;
			cout << *onlineMatch_it;
			cout <<"Masked flag: " << MaskedFlag << endl << endl;

		}
		}//loop over emulator TPGs


	//Looping over Online TPGs collection
	//std::cout << "looping on online TPG" << std::endl;

	for (onTP_it = onlineTPs.begin(); onTP_it != onlineTPs.end(); ++onTP_it) {
		double OnlineCompressedEt=-1;
		OnlineCompressedEt = onTP_it->compressedEt();

		int online_sFGVB;
		online_sFGVB = onTP_it->sFGVB();

		EcalTrigTowerDetId towerid = onTP_it->id();
		int online_tower_ieta( towerid.ieta() );
		int online_tower_iphi( towerid.iphi() );

		//Get rid of noisy towers. TEMPORARY!!!
		//(To match cuts which were applied to ttemulator)
		if ( (online_tower_ieta == 14)  &&  ( online_tower_iphi == 28 ) ) continue;
		if ( (online_tower_ieta == -8)  &&  ( online_tower_iphi == 66 ) ) continue;

		if (OnlineCompressedEt > 0) {
			h_ttonline_et->Fill(OnlineCompressedEt);
			h_ttonline_raw->Fill((*onTP_it)[0].raw() & 0xff);
			switch (online_sFGVB)
			{
				case 0: h_ttonline_et_sFGVB0On->Fill(OnlineCompressedEt);
				break;
				case 1: h_ttonline_et_sFGVB1On->Fill(OnlineCompressedEt);
				break;
			}

			// Store only barrel online TPGs
			if ( ( -18 < online_tower_ieta) && (online_tower_ieta < 18) )
				h_ttonline_barrel_et->Fill(OnlineCompressedEt);
				h_ttonline_barrel_ieta->Fill(online_tower_ieta);
				h_ttonline_barrel_iphi->Fill(online_tower_iphi);
				h_ttonline_barrel_map->Fill(online_tower_ieta, online_tower_iphi);

				if (OnlineCompressedEt == 255) {
					h_ttonline_barrel_saturated_ieta->Fill(online_tower_ieta);
					h_ttonline_barrel_saturated_iphi->Fill(online_tower_iphi);
					h_ttonline_barrel_saturated_map->Fill(online_tower_ieta, online_tower_iphi);
					}

		}  //non zero TP
		h_ttonline_et_nozerocut->Fill(OnlineCompressedEt);
		if ( ( -18 < online_tower_ieta) && (online_tower_ieta < 18) )
			h_ttonline_barrel_et_nozerocut->Fill(OnlineCompressedEt);

		//find matching tower in the emulator collection
		EcalTrigPrimDigiCollection::const_iterator emulatorMatch_it;
		emulatorMatch_it = emulatorTPs.find(onTP_it->id());

		if (emulatorMatch_it != emulatorTPs.end()) {
			double EmulatorMatchCompressedEt;
			EmulatorMatchCompressedEt = emulatorMatch_it->compressedEt();
			int emulatorMatch_sFGVB;
			emulatorMatch_sFGVB=emulatorMatch_it->sFGVB();

			EcalTrigTowerDetId emul_match_towerid = emulatorMatch_it->id();
			int emul_match_tower_ieta( emul_match_towerid.ieta() );
			int emul_match_tower_iphi( emul_match_towerid.iphi() );

//			if ( (online_tower_ieta != emul_match_tower_ieta) ||
//				 (online_tower_iphi != emul_match_tower_iphi)
//			){
//				cout << "Unmatched towers: " << endl;
//				cout << "Online:" << endl;
//				cout << "ieta: " << online_tower_ieta << " iphi: " << online_tower_iphi << endl;
//				cout << "ADC: " << ( (*onTP_it)[0].raw()&0xff )<< " TTF: " << ( (*onTP_it)[0].ttFlag() ) << endl;
//				cout << endl;
//				cout << "Emulated:" << endl;
//				cout << "ieta: " << emul_match_tower_ieta << " iphi: " << emul_match_tower_iphi << endl;
//				cout << "ADC: " << ( (*emulatorMatch_it)[0].raw()&0xff )<< " TTF: " <<( (*emulatorMatch_it)[0].ttFlag() ) << endl;
//				cout << endl;
//
//			}


			int MaskedFlag;
			auto pair = towerMap.find(emulatorMatch_it->id());
			if (pair!=towerMap.end()){
				MaskedFlag=pair->second;
			} else{
				MaskedFlag=0;
			}

			//different online sFGVB
			if ( EmulatorMatchCompressedEt>0 ){
				switch (online_sFGVB)
					{
						case 0: h_ttemulator_et_sFGVB0On->Fill(EmulatorMatchCompressedEt);
								if (MaskedFlag<=0) h_ttemulator_masked_et_sFGVB0On->Fill(EmulatorMatchCompressedEt);
						break;
						case 1: h_ttemulator_et_sFGVB1On->Fill(EmulatorMatchCompressedEt);
								if (MaskedFlag<=0) h_ttemulator_masked_et_sFGVB1On->Fill(EmulatorMatchCompressedEt);
						break;
					}
					h_ttonline_et_zeroInemulator->Fill(OnlineCompressedEt);
					if (emulatorMatch_sFGVB==0) h_ttonline_et_zeroInemulator_sFGVB0Em->Fill(OnlineCompressedEt);
					if (emulatorMatch_sFGVB==1) h_ttonline_et_zeroInemulator_sFGVB1Em->Fill(OnlineCompressedEt);
			}
			//correlation between online and emulated TPGs
			h_ttonline_et_against_ttemulator_et->Fill(OnlineCompressedEt, EmulatorMatchCompressedEt);

			if ( OnlineCompressedEt>0 ) {
				h_ttonline_et_against_ttemulator_et_noOnlineZero->Fill(OnlineCompressedEt, EmulatorMatchCompressedEt);
				if (EmulatorMatchCompressedEt>0){
					h_ttonline_et_against_ttemulator_et_nozero->Fill(OnlineCompressedEt, EmulatorMatchCompressedEt);
				}
			}

			if (EmulatorMatchCompressedEt>0)
					h_ttonline_et_against_ttemulator_et_noEmulatorZero->Fill(OnlineCompressedEt, EmulatorMatchCompressedEt);


			if ( OnlineCompressedEt != EmulatorMatchCompressedEt) {
				if ( ( -18 < online_tower_ieta) && (online_tower_ieta < 18) ){
//					cout << "Unmatch towers : " << endl;
//					cout << "Online et: " << OnlineCompressedEt << endl;
//					cout << "Emulated et: " << EmulatorMatchCompressedEt << endl;
//					cout << "ieta: " << online_tower_ieta << " iphi: " << online_tower_iphi << endl;
					h_TPGs_unmatch_online_offline->Fill(online_tower_ieta, online_tower_iphi);

					if ( (OnlineCompressedEt > 0) && (EmulatorMatchCompressedEt > 0) ){
//						cout << "Unmatch nonzero towers : " << endl;
//						cout << "Online et: " << OnlineCompressedEt << endl;
//						cout << "Emulated et: " << EmulatorMatchCompressedEt << endl;
//						cout << "ieta: " << online_tower_ieta << " iphi: " << online_tower_iphi << endl;
						h_TPGs_unmatch_online_offline_nonzero->Fill(online_tower_ieta, online_tower_iphi);

//						if ( (online_tower_ieta == 14) && (online_tower_iphi == 28) ) {
//							cout << "Huge peak: " << endl;
//							cout << "Online:" << endl;
//							cout << "ieta: " << online_tower_ieta << " iphi: " << online_tower_iphi << endl;
//							cout << "ADC: " << ( (*onTP_it)[0].raw()&0xff )<< " TTF: " << ( (*onTP_it)[0].ttFlag() ) << endl;
//
//							cout << "Emulated:" << endl;
//							cout << "ieta: " << emul_match_towerid.ieta() << " iphi: " << emul_match_towerid.iphi() << endl;
//							cout << "ADC: " << ( (*emulatorMatch_it)[0].raw()&0xff )<< " TTF: " <<( (*emulatorMatch_it)[0].ttFlag() ) << endl;
//
//						}

					}

//					//Printing out info for specific huge peak
//					if ( (online_tower_ieta == 14) && (online_tower_iphi == 28) ) {
//						cout << "Huge peak: " << endl;
//						cout << "Online:" << endl;
//						cout << "ieta: " << online_tower_ieta << " iphi: " << online_tower_iphi << endl;
//						cout << "ADC: " << ( (*onTP_it)[0].raw()&0xff )<< " TTF: " << ( (*onTP_it)[0].ttFlag() ) << endl;
//
//						cout << "Emulated:" << endl;
//						cout << "ieta: " << emul_match_towerid.ieta() << " iphi: " << emul_match_towerid.iphi() << endl;
////						cout << "ADC: " << ( (*emulatorMatch_it)[0].raw()&0xff )<< " TTF: " <<( (*emulatorMatch_it)[0].ttFlag() ) << endl;
////						cout << "ADC1: " << ( (*emulatorMatch_it)[1].raw()&0xff )<< " TTF1: " <<( (*emulatorMatch_it)[1].ttFlag() ) << endl;
//						cout << "ADC: " << ( (*emulatorMatch_it)[2].raw()&0xff )<< " TTF: " <<( (*emulatorMatch_it)[2].ttFlag() ) << endl;
//
//					}
				}
			}


			//std::cout << "Offline match"<< std::endl;
			//std::cout << "ET="<<emulatorMatch_it->compressedEt()<<"  iphi="<< emulatorMatch_it->id().iphi() << "  ieta=" << offlineMatchForSpike_it->id().ieta()<< std::endl;

			//if online TT is saturated
			if (onTP_it->compressedEt() == 255) {
				//std::cout << "SaturatedTP ------------" << std::endl;
				//std::cout << "ET=255  ieta=" << onTP_it->id().ieta()
				//		<< "  iphi=" << onTP_it->id().iphi()
				//		<< " sFVGB=" << onTP_it->sFGVB()
				//	<< std::endl;

				//h_TPGsaturated->Fill(onTP_it->id().ieta(), onTP_it->id().iphi());

				h_ttonline_saturated_to_ttoffline_match_et->Fill( (*emulatorMatch_it)[2].raw() & 0xff );

			} //  end if saturated

		} else {
			//std::cout << "No match tower for saturated tower was found in the emulator TT collection:"<< std::endl;
			//std::cout << "ET=255  iphi="<< onTP_it->id().iphi() << "  ieta=" << onTP_it->id().ieta()<< std::endl;
		} // end else


	} //loop online TPG



	//------------------------//
	// GET RecHit Spikes	  //
	//------------------------//

	// geometry (used for L1 trigger)
	//   //cout << "get the geometry" << endl;
	edm::ESHandle<CaloSubdetectorGeometry> theBarrelGeometry_handle;
	iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",
			theBarrelGeometry_handle);
	//iSetup.get<IdealGeometryRecord>().get(eTTmap_);
	theBarrelGeometry_ = &(*theBarrelGeometry_handle);

	//cout << "starting getting the spike-like rechits " << endl;

	// channel status : old way to do it
	/*
	 edm::ESHandle<EcalChannelStatus> pChannelStatus;
	 iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);
	 const EcalChannelStatus *chaStatus = pChannelStatus.product();
	 */
	// for 42x
	unsigned long long cacheSevLevel = 0;
	edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
	if (cacheSevLevel
			!= iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()) {
		cacheSevLevel =
				iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
		iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
	}
	const EcalSeverityLevelAlgo* sl = sevLevel.product();
	// Get EB rechits
	edm::Handle<EcalRecHitCollection> rechitsEB;
	EcalRecHitCollection::const_iterator rechitItr;

	EBDetId id;

	int flag = 0;
	uint32_t sev = 0;
	double thetaHit, etaHit, phiHit;
	int rhit_ieta, rhit_iphi;

	int i = 0;
	//cout << "looking" << endl;

	std::map<EBDetId, EBDataFrame> mapdigis_;

	// Create digis map for the amplitude reconstruction.
	//modif-digis-beg
	edm::Handle<EBDigiCollection> digiEB;
	iEvent.getByLabel("ecalEBunpacker", "ebDigis", digiEB);
	const EBDigiCollection *theBarrelEcalDigis = digiEB.product();
	//cout << "There are " << theBarrelEcalDigis->size() << " digis" << endl;
	mapdigis_.clear();
	for (unsigned int i = 0; i < theBarrelEcalDigis->size(); i++) {
		EBDataFrame df = (*(theBarrelEcalDigis))[i];
		EBDetId id = df.id();
		mapdigis_.insert(make_pair(id, df));
//		theBarrelEcalDigis->find()
	}

  //intercalibration constants
  edm::ESHandle<EcalIntercalibConstants> ical;
  iSetup.get<EcalIntercalibConstantsRcd>().get(ical);
  const EcalIntercalibConstantMap& icalMap = ical->getMap();

  //calibration
  edm::ESHandle<EcalADCToGeVConstant> agc;
  iSetup.get<EcalADCToGeVConstantRcd>().get(agc);
  float calibAdcToGeV = float(agc->getEBValue());
//  cout << "calib=" << float(agc->getEBValue()) << endl;
//  cout << "calibAdcToGeV = " << calibAdcToGeV;



////// RECONSTRUCT AMPLITUDES OF TRIGGER TOWERS FROM CRYSTALs ADCs

if (do_reconstruct_amplitude_allTTs_==true){
//  Reconstruct amplitudes for all trigger towers
  for (onTP_it = onlineTPs.begin(); onTP_it != onlineTPs.end(); ++onTP_it) {
	  	const EcalTrigTowerDetId towid = onTP_it->id();

	  	// calculate crystal ieta,iphy indexes from tower indexes.
	  	int cr_start_ieta, cr_start_iphi;
		cr_start_ieta = ( signbit(towid.ieta()) ) ? ( towid.ieta()*5 ) : ( towid.ieta()*5-4 );
		cr_start_iphi = towid.iphi()*5-4;
		//cr_start_iphi = (int)(rhit_iphi / 5) * 5 + 1;

		double TriggerTowerAmplitude=0;
		double TriggerTowerAmplitude_icalib=0;
		double TriggerTowerAmplitude_icalib_ADCtoGev=0;
		double TriggerTowerAmplitude_icalib_ADCtoGev_sintheta=0;

		bool there_are_digis_for_tower=false;
		for (int ieta_it = cr_start_ieta; ieta_it < cr_start_ieta + 5; ieta_it++) {
			for (int iphi_it = cr_start_iphi; iphi_it < cr_start_iphi + 5; iphi_it++) {
				EBDetId neighbour_id(ieta_it, iphi_it);

				map<EBDetId, EBDataFrame>::const_iterator idigis = 	mapdigis_.find(neighbour_id);
				if (idigis == mapdigis_.end()) { continue ;}
				there_are_digis_for_tower=true;
				EBDataFrame df = idigis->second;

				int gain = 0;
				int adc = 0;
				double mean = 0;
				int maxADC = 0;

				double amplitudetpg = 0.0;
				double amplitudeoffline = 0.0;

				bool shouldstoregain = false;
				for (int samp = 0; samp < 10; samp++) {
					adc = df[samp].adc();
					gain = df[samp].gainId();
					if (samp < 3)
						mean += adc;
					if (adc > maxADC) {
						maxADC = adc;
					}	// maxsamp= samp;} //steph
					

					//temporary
					if (gain>1)  {shouldstoregain=true;}


					amplitudetpg += weight[samp] * adc;
					amplitudeoffline += weight_offline[samp] * adc;

				}	//loop samples
				mean /= 3.0;

				if (shouldstoregain) h_crystalIds_gain2_3->Fill(ieta_it, iphi_it,1); //cout <<df<<endl;

//				if (shouldstoregain) {
//					cout << "High gain : " << *onTP_it << endl;
//					cout << "Emulator match : " << endl;
//					EcalTrigPrimDigiCollection::const_iterator emulatorMatch_it;
//					emulatorMatch_it = emulatorTPs.find(onTP_it->id());
//					//find matching tower in the emulator collection
//					if (emulatorMatch_it!=emulatorTPs.end()){
//						cout << *emulatorMatch_it << endl;
//					}
//				}

				// first intercalibration constants
				EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(neighbour_id);
				EcalIntercalibConstant icalconst = 1;
				if( icalit!=icalMap.end() ) {
					icalconst = (*icalit);
				} else {
					cout << "No intercalib const found for xtal "
					<< neighbour_id.rawId()
					<< "! something wrong with EcalIntercalibConstants in your DB? ";
				}

				//h_temp_hist->Fill(icalconst);

				double amplitudetpg_icalib = 0.0;
				double amplitudetpg_icalib_ADCtoGeV = 0.0;
				if (amplitudetpg < 0 ) {
					amplitudetpg=0;
				} else {
					amplitudetpg_icalib = amplitudetpg*icalconst;
					amplitudetpg_icalib_ADCtoGeV = amplitudetpg_icalib*calibAdcToGeV;
				}
				//									  cout << "           rechit e=" << myhit.energy() << " time=" << myhit.time() << " " << towid
//										   << " amplitudetpg=" << amplitudetpg << " offline=" << amplitudeoffline
//										   << " amplitudetpg_icalib_ADCtoGeV=" << amplitudetpg*calibAdcToGeV*icalconst << " " << amplitudeoffline*calibAdcToGeV*icalconst
//										   << " iconst=" << icalconst << endl; //<< " uncalibamplitude=" << myuhit.amplitude() << " " << maxADC - mean << endl;

				TriggerTowerAmplitude += amplitudetpg;
				TriggerTowerAmplitude_icalib+=amplitudetpg_icalib;
				TriggerTowerAmplitude_icalib_ADCtoGev += amplitudetpg_icalib_ADCtoGeV;
				double thetacrystal;
				thetacrystal = (theBarrelGeometry_->getGeometry(neighbour_id)->getPosition()).theta();
				TriggerTowerAmplitude_icalib_ADCtoGev_sintheta += amplitudetpg_icalib_ADCtoGeV * sin(thetacrystal);

				//cout<<"amplitudetpg_icalib_ADCtoGeV = " << amplitudetpg_icalib_ADCtoGeV <<  " amplitudetpg_icalib_ADCtoGeV_sintheta = " <<  amplitudetpg_icalib_ADCtoGeV * sin(thetacrystal) << " sin = " << sin(thetacrystal) << endl;
				//if (sin(thetacrystal)<0) cout << "!! sin = " << sin(thetacrystal) << endl;
			} // end iphi loop
		}  // end ieta loop

		if (there_are_digis_for_tower) { // only saving data if at least for one crystals digis are stored
			h_tt_amplitude_online->Fill(TriggerTowerAmplitude);
			h_tt_amplitude_calibrated_online->Fill(TriggerTowerAmplitude_icalib_ADCtoGev);
			if (TriggerTowerAmplitude>0){ h_tt_amplitude_online_above_zero->Fill(TriggerTowerAmplitude); }
			if (TriggerTowerAmplitude_icalib>0){h_tt_amplitude_icalib_online_above_zero->Fill(TriggerTowerAmplitude_icalib);}
			if (TriggerTowerAmplitude_icalib_ADCtoGev>0){h_tt_amplitude_icalib_ADCtoGeV_online_above_zero->Fill(TriggerTowerAmplitude_icalib_ADCtoGev);}
			if (TriggerTowerAmplitude_icalib_ADCtoGev_sintheta>0) h_tt_amplitude_icalib_ADCtoGeV_online_above_zero_sintheta->Fill(TriggerTowerAmplitude_icalib_ADCtoGev_sintheta);
		}
	}
}// if do_reconstruct_amplitude_allTTs_


//////  BROWSE RecHit COLLECTION //////////////

//Map to store info about the spikes
//Map covers only barrel regions
//Clearing map

map_RCT_regs_with_offline_spikes.clear();
for (unsigned int reg_ieta=6; reg_ieta < 15; reg_ieta++){
	for (unsigned int reg_iphi=0; reg_iphi < 18; reg_iphi++){
		L1CaloRegionDetId RCT_region_id(reg_ieta, reg_iphi);
		map_RCT_regs_with_offline_spikes[RCT_region_id]=false;
	}
}


// starting rechit search
if (do_rechit_spikes_search_ || do_reconstruct_amplitudes_spikes_) {
//Find spikes in the RecHit collection and reconstruct amplitudes of crystals inside the trigger tower
//	cout << "Inside RecHit Loop." << endl;
	if (iEvent.getByLabel(EcalRecHitCollectionEB_, rechitsEB)) {
		//cout << "looping" << endl;
		for (rechitItr = rechitsEB->begin(); rechitItr != rechitsEB->end();
				++rechitItr) {
			//cout << "i=" << i << endl;
			if (i >= 5000) {
				cout << "more than 5000 spikes in run " << iEvent.id().run()
						<< " , event " << iEvent.id().event() << endl;
				//break;
			}
			id = rechitItr->id();
			//sev = EcalSeverityLevelAlgo::severityLevel( id, *rechitsEB, *chaStatus );
			//sev = EcalSeverityLevelAlgo::severityLevel(id,*rechitsEB,*chaStatus, 5.,
			//EcalSeverityLevelAlgo::kSwissCross,0.95) ;
			sev = 0;
			sev = sl->severityLevel(id, *rechitsEB);
			//if(sev>0) cout << "got sev=" << sev << endl;
			//if (sev >= 1) {
			thetaHit =
					(theBarrelGeometry_->getGeometry(id)->getPosition()).theta();
//			double rechits_et;
//			rechits_et = (rechitItr->energy()) * sin(thetaHit);
			h_rechits_et->Fill( (rechitItr->energy()) * sin(thetaHit) );
			h_rechits_energy->Fill(rechitItr->energy());

			double TriggerTowerAmplitude=0;
			double TriggerTowerAmplitude_icalib=0;
			double TriggerTowerAmplitude_icalib_ADCtoGev=0;
			double TriggerTowerAmplitude_icalib_ADCtoGev_sintheta=0;

//			cout << id << endl;
//			cout << "sev = " << sev<< endl;
			if ((sev == 3) || (sev == 4)) {  // if it's a spike
//				sev3_sev4_count++;
//				cout << "Inside Spikes Loop." << endl;
				//thetaHit = (theBarrelGeometry_->getGeometry(id)->getPosition()).theta();
//				etaHit = (theBarrelGeometry_->getGeometry(id)->getPosition()).eta();
//				phiHit = (theBarrelGeometry_->getGeometry(id)->getPosition()).phi();
//				rhit_ieta = id.ieta();
				rhit_iphi = id.iphi();

//				cout << endl << endl << "Inside loop" << endl << endl;

				const EcalTrigTowerDetId towid = id.tower();

				TrigTowersWithRecHitSpikes.push_back(towid);
				flag = 0;
				flag = rechitItr->recoFlag();

				// Fill map of ECT regions with spikes
				EBDetId tower_detId(towid.ieta(), towid.iphi());
				L1CaloRegionDetId & RCT_region_match = map_DetId_RegionID[tower_detId];
				map_RCT_regs_with_offline_spikes[RCT_region_match]=true;
//				cout << towid << endl;
//				cout << "Reg ieta : " << RCT_region_match.ieta() << endl;
//				cout << "Reg iphi : " << RCT_region_match.iphi() << endl;

				if (do_reconstruct_amplitudes_spikes_==true) {
					EcalTrigPrimDigiCollection::const_iterator TrigPrimWithReghitSpike_it;
					TrigPrimWithReghitSpike_it = onlineTPs.find(towid);

					if (TrigPrimWithReghitSpike_it != onlineTPs.end()){ // if towid exists in towid collection
						if (TrigPrimWithReghitSpike_it->sFGVB()==1) { // if spike wasn't recognized like spike (or TrigTower energy to low)
	//						sev3_sev4_with_sFGVB1++;
							int cr_start_ieta, cr_start_iphi;
							cr_start_ieta = ( signbit(towid.ieta()) ) ? ( towid.ieta()*5 ) : ( towid.ieta()*5-4 );
							//cr_start_iphi = towid.iphi()*5-4;
							cr_start_iphi = (int)(rhit_iphi / 5) * 5 + 1;

							/*cout << "-- RecHit : " << "flag=" << flag << " | sev="
									<< sev << " | thetaHit=" << thetaHit << " | phiHit="
									<< phiHit << " | etaHit=" << etaHit << " | ieta="
									<< rhit_ieta << " | iphi=" << rhit_iphi << " | E="
									<< rechitItr->energy() << " | Et=" << rechits_et
									<< " | TT=" << towid.ieta() << " " << towid.iphi()
									//<< " | ADC=" << rechitItr->adc();
									<< endl;
							*/
							//cout << id << endl;

							//cout << "-- RecHit and neighbors: " << endl;

							TH2D * newSpike2DHist;
							TH2D * newSpike2DHist_icalib;
							TH2D * newSpike2DHist_icalib_ADCtoGeV;

							if (do_3DSpike_plots_) {
								spiked_hists_number++;
								std::string hist_name;
		//						hist_name = "h_SpikedTT_"
		//								+ std::to_string(spiked_hists_number) + "_Ev" + std::to_string((int)(iEvent.eventAuxiliary().event()));
								hist_name = (std::string)"h_SpikedTT"
										+ "_Run" + std::to_string( ( iEvent.id().run() ))
										+ "_Ev" + std::to_string( ( iEvent.id().event() ))
										+ "_Lum" + std::to_string( ( iEvent.id().luminosityBlock() ));
		//						TH2D * newSpike2DHist = new TH2D(hist_name.c_str(),
		//								hist_name.c_str(), 5, cr_start_ieta, cr_start_ieta + 4, 5, cr_start_iphi, cr_start_iphi + 4);
								newSpike2DHist = fs->make<TH2D>(hist_name.c_str(),
										hist_name.c_str(), 5, cr_start_ieta, cr_start_ieta + 4, 5, cr_start_iphi, cr_start_iphi + 4);
								hvect_Spikes.push_back(newSpike2DHist);
								newSpike2DHist->GetZaxis()->SetTitle("amplitude");
								newSpike2DHist->GetYaxis()->SetTitle("iphi");
								newSpike2DHist->GetXaxis()->SetTitle("ieta");

								std::string icalib_hist_name;
								icalib_hist_name = hist_name + "_icalib";
								newSpike2DHist_icalib = fs->make<TH2D>(icalib_hist_name.c_str(),
										icalib_hist_name.c_str(), 5, cr_start_ieta, cr_start_ieta + 4, 5, cr_start_iphi, cr_start_iphi + 4);
								hvect_Spikes_icalib.push_back(newSpike2DHist_icalib);
								newSpike2DHist_icalib->GetZaxis()->SetTitle("amplitude");
								newSpike2DHist_icalib->GetYaxis()->SetTitle("iphi");
								newSpike2DHist_icalib->GetXaxis()->SetTitle("ieta");

								std::string icalib_ADCtoGeV_hist_name;
								icalib_ADCtoGeV_hist_name = hist_name + "_icalib_ADCtoGeV";
								newSpike2DHist_icalib_ADCtoGeV = fs->make<TH2D>(icalib_ADCtoGeV_hist_name.c_str(),
										icalib_ADCtoGeV_hist_name.c_str(), 5, cr_start_ieta, cr_start_ieta + 4, 5, cr_start_iphi, cr_start_iphi + 4);
								hvect_Spikes_icalib_ADCtoGeV.push_back(newSpike2DHist_icalib_ADCtoGeV);
								newSpike2DHist_icalib_ADCtoGeV->GetZaxis()->SetTitle("GeV");
								newSpike2DHist_icalib_ADCtoGeV->GetYaxis()->SetTitle("iphi");
								newSpike2DHist_icalib_ADCtoGeV->GetXaxis()->SetTitle("ieta");
							}

							for (int ieta_it = cr_start_ieta;
									ieta_it < cr_start_ieta + 5; ieta_it++) {
								for (int iphi_it = cr_start_iphi;
										iphi_it < cr_start_iphi + 5; iphi_it++) {
									EBDetId neighbour_id(ieta_it, iphi_it);
									//cout << neighbour_id << endl;

									map<EBDetId, EBDataFrame>::const_iterator idigis =
											mapdigis_.find(neighbour_id);
										  if (idigis == mapdigis_.end()) { continue ;}
									EBDataFrame df = idigis->second;

									int gain = 0;
									int adc = 0;
									double mean = 0;
									int maxADC = 0;

									double amplitudetpg = 0.0;
									double amplitudeoffline = 0.0;
									for (int samp = 0; samp < 10; samp++) {
										adc = df[samp].adc();
										gain = df[samp].gainId();
										if (samp < 3)
											mean += adc;
										if (adc > maxADC) {
											maxADC = adc;
										}	// maxsamp= samp;} //steph

										amplitudetpg += weight[samp] * adc;
										amplitudeoffline += weight_offline[samp] * adc;

									}	//loop samples
									mean /= 3.0;

									// first intercalibration constants
									EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(neighbour_id);
									EcalIntercalibConstant icalconst = 1;
									if( icalit!=icalMap.end() ) {
										icalconst = (*icalit);
									} else {
										cout << "No intercalib const found for xtal "
										<< neighbour_id.rawId()
										<< "! something wrong with EcalIntercalibConstants in your DB? ";
									}



									double amplitudetpg_icalib = 0.0;
									double amplitudetpg_icalib_ADCtoGeV = 0.0;
									if (amplitudetpg < 0 ) {
										amplitudetpg=0;
									} else {
										amplitudetpg_icalib = amplitudetpg*icalconst;
										amplitudetpg_icalib_ADCtoGeV = amplitudetpg_icalib*calibAdcToGeV;
									}
	//									  cout << "           rechit e=" << myhit.energy() << " time=" << myhit.time() << " " << towid
	//										   << " amplitudetpg=" << amplitudetpg << " offline=" << amplitudeoffline
	//										   << " amplitudetpg_icalib_ADCtoGeV=" << amplitudetpg*calibAdcToGeV*icalconst << " " << amplitudeoffline*calibAdcToGeV*icalconst
	//										   << " iconst=" << icalconst << endl; //<< " uncalibamplitude=" << myuhit.amplitude() << " " << maxADC - mean << endl;

									TriggerTowerAmplitude += amplitudetpg;
									TriggerTowerAmplitude_icalib+=amplitudetpg_icalib;
									TriggerTowerAmplitude_icalib_ADCtoGev += amplitudetpg_icalib_ADCtoGeV;

									double theta_cryst;
									theta_cryst = (theBarrelGeometry_->getGeometry(neighbour_id)->getPosition()).theta();
									TriggerTowerAmplitude_icalib_ADCtoGev_sintheta+=amplitudetpg_icalib_ADCtoGeV * sin(theta_cryst);

									if (do_3DSpike_plots_){
										newSpike2DHist->Fill(ieta_it, iphi_it, amplitudetpg);
										newSpike2DHist_icalib->Fill(ieta_it, iphi_it, amplitudetpg_icalib);
										newSpike2DHist_icalib_ADCtoGeV->Fill(ieta_it, iphi_it, amplitudetpg_icalib_ADCtoGeV);
									}
								} // end iphi loop
							}  // end ieta loop

							if (TriggerTowerAmplitude > 0) h_tt_spikes_amplitude_online->Fill(TriggerTowerAmplitude);
							if (TriggerTowerAmplitude_icalib > 0) h_tt_spikes_amplitude_icalib_online->Fill(TriggerTowerAmplitude_icalib);
							if (TriggerTowerAmplitude_icalib_ADCtoGev > 0) h_tt_spikes_amplitude_icalib_ADCtoGeV_online->Fill(TriggerTowerAmplitude_icalib_ADCtoGev);
							if (TriggerTowerAmplitude_icalib_ADCtoGev_sintheta > 0) h_tt_spikes_amplitude_icalib_ADCtoGeV_online_sintheta->Fill(TriggerTowerAmplitude_icalib_ADCtoGev_sintheta );

						} // end if towid sFGVB
					} //end if towid found
				} //if do_reconstruct_amplitudes_spikes_

			    // fill histograms (Spikes, matches in online, matches in online with sFGVB)
				double rechit_spikes_et;
				double rechit_spikes_energy;
				rechit_spikes_et = (rechitItr->energy()) * sin(thetaHit);
				rechit_spikes_energy = rechitItr->energy();

				if (rechit_spikes_et > 0) h_rechit_spikes_et->Fill( rechit_spikes_et );
				if (rechit_spikes_energy > 0) h_rechit_spikes_energy->Fill( rechit_spikes_energy );

//				if (rechit_spikes_et > 30)
//					cout << "BIG SPIKE-----=" << towid << endl;

				//find matching TT in the onilne collection
				EcalTrigPrimDigiCollection::const_iterator onlineMatchForSpike_it;
				onlineMatchForSpike_it = onlineTPs.find(towid);
				if (onlineMatchForSpike_it != onlineTPs.end()) {
					double rechit_spikes_match_to_ttonline_et;
					rechit_spikes_match_to_ttonline_et = onlineMatchForSpike_it->compressedEt();

					h_rechit_spikes_et_FoundInOnlineTTs->Fill(rechit_spikes_et);
					h_rechit_spikes_energy_FoundInOnlineTTs->Fill(rechit_spikes_energy);

					int sFGVB = onlineMatchForSpike_it->sFGVB();
					if (sFGVB == 0)	h_rechit_spikes_et_sFGVB0On->Fill( rechit_spikes_et );
					if (sFGVB == 1) h_rechit_spikes_et_sFGVB1On->Fill( rechit_spikes_et );
					if (sFGVB == 0)	h_rechit_spikes_energy_sFGVB0On->Fill( rechit_spikes_energy );
					if (sFGVB == 1) h_rechit_spikes_energy_sFGVB1On->Fill( rechit_spikes_energy );

					h_rechit_spikes_match_to_ttonline_et_nozerocut->Fill(	rechit_spikes_match_to_ttonline_et );
					if (sFGVB == 0) {
						h_rechit_spikes_to_ttonline_et_sFGVB0On_nozerocut->Fill(onlineMatchForSpike_it->compressedEt() );
					}
					if (sFGVB == 1) {
						h_rechit_spikes_to_ttonline_et_sFGVB1On_nozerocut->Fill(onlineMatchForSpike_it->compressedEt() );
					}

					if (rechit_spikes_match_to_ttonline_et > 0 ) {
						h_rechit_spikes_match_to_ttonline_et->Fill(	rechit_spikes_match_to_ttonline_et );

						if (sFGVB == 0) {
							//std::cout << "sFGVB=0" << std::endl;
//							rechit_spikes_to_ttonline_et_sFGVB0On =
//									onlineMatchForSpike_it->compressedEt();
							h_rechit_spikes_to_ttonline_et_sFGVB0On->Fill( onlineMatchForSpike_it->compressedEt() );
						}
						if (sFGVB == 1) {
							//std::cout << "sFGVB=1" << std::endl;
//							rechit_spikes_to_ttonline_et_sFGVB1On =
//									onlineMatchForSpike_it->compressedEt();
							h_rechit_spikes_to_ttonline_et_sFGVB1On->Fill( onlineMatchForSpike_it->compressedEt() );
						}
					}
				} else {
					cout << "Spike TT wasn't found in online collection!"
							<< std::endl;
				} // end if spike not found online

				//search for spike TT in the emulator collection
				EcalTrigPrimDigiCollection::const_iterator emulatorMatchForSpike_it;
				emulatorMatchForSpike_it = emulatorTPs.find(towid);
				if (emulatorMatchForSpike_it != emulatorTPs.end()) {
						double rechit_spikes_match_to_ttemulator_et;
						rechit_spikes_match_to_ttemulator_et = emulatorMatchForSpike_it->compressedEt();

						h_rechit_spikes_et_FoundInEmulatorTTs->Fill(rechit_spikes_et);
						h_rechit_spikes_energy_FoundInEmulatorTTs->Fill(rechit_spikes_energy);

						int sFGVB = emulatorMatchForSpike_it->sFGVB();
						if (sFGVB == 0)	h_rechit_spikes_et_sFGVB0Em->Fill( rechit_spikes_et );
						if (sFGVB == 1) h_rechit_spikes_et_sFGVB1Em->Fill( rechit_spikes_et );
						if (sFGVB == 0)	h_rechit_spikes_energy_sFGVB0Em->Fill( rechit_spikes_energy );
						if (sFGVB == 1) h_rechit_spikes_energy_sFGVB1Em->Fill( rechit_spikes_energy );

						h_rechit_spikes_match_to_ttemulator_et_nozerocut->Fill(	rechit_spikes_match_to_ttemulator_et );
						if (sFGVB == 0)	h_rechit_spikes_to_ttemulator_et_sFGVB0Em_nozerocut->
											Fill( emulatorMatchForSpike_it->compressedEt() );
						if (sFGVB == 1) h_rechit_spikes_to_ttemulator_et_sFGVB1Em_nozerocut->
											Fill( emulatorMatchForSpike_it->compressedEt() );

						if (rechit_spikes_match_to_ttemulator_et > 0 ) {
							h_rechit_spikes_match_to_ttemulator_et->Fill(	rechit_spikes_match_to_ttemulator_et );
							if (sFGVB == 0)	h_rechit_spikes_to_ttemulator_et_sFGVB0Em->
												Fill( emulatorMatchForSpike_it->compressedEt() );
							if (sFGVB == 1) h_rechit_spikes_to_ttemulator_et_sFGVB1Em->
												Fill(emulatorMatchForSpike_it->compressedEt() );
						}
					} else {
						cout << "Spike TT wasn't found in emulator collection!"
								<< std::endl;
					} // end if spike not found in emulator collection

				if (emulatorMatchForSpike_it != emulatorTPs.end() and
							onlineMatchForSpike_it != onlineTPs.end()){
					double rechit_spikes_match_to_ttonline_et;
					rechit_spikes_match_to_ttonline_et = onlineMatchForSpike_it->compressedEt();

					double rechit_spikes_match_to_ttemulator_et;
					rechit_spikes_match_to_ttemulator_et = emulatorMatchForSpike_it->compressedEt();

					int onlineMatch_sFGVB;
					onlineMatch_sFGVB=onlineMatchForSpike_it->sFGVB();
					int emulatorMatch_sFGVB;
					emulatorMatch_sFGVB=emulatorMatchForSpike_it->sFGVB();

					if (  rechit_spikes_match_to_ttonline_et>0   &&
						!(rechit_spikes_match_to_ttemulator_et>0)	){
							h_rechit_spikes_ttonline_et_zeroInEmulator->Fill(rechit_spikes_match_to_ttonline_et);
							if (onlineMatch_sFGVB==0) h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0On->Fill(rechit_spikes_match_to_ttonline_et);
							if (onlineMatch_sFGVB==1) h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1On->Fill(rechit_spikes_match_to_ttonline_et);
							if (emulatorMatch_sFGVB==0) h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB0Em->Fill(rechit_spikes_match_to_ttonline_et);
							if (emulatorMatch_sFGVB==1) h_rechit_spikes_ttonline_et_zeroInEmulator_sFGVB1Em->Fill(rechit_spikes_match_to_ttonline_et);

					}
				}

				//cout << "is gonna fill the variables" << endl;
				if (flag == EcalRecHit::kOutOfTime){;}
					 //spike_outOfTime[i] = 1;
				//else spike_outOfTime[i] = 0;
				/*spike_severityLevel[i] = sev ;
				 spike_time[i] = rechitItr->time();
				 spike_Et[i] = (rechitItr->energy())*sin(thetaHit);
				 spike_phi[i] = phiHit;
				 spike_eta[i] = etaHit;
				 spike_theta[i] = thetaHit;
				 spike_TTiphi[i] = towid.iphi();
				 spike_TTieta[i] = towid.ieta();
				 spike_Riphi[i] = getGCTRegionPhi(towid.iphi());
				 spike_Rieta[i] = getGCTRegionEta(towid.ieta());
				 */
				i++;
			} // end if sev=3 || sev=4
		}	//loop rechit

	}	//rechits EB
} // if do_rechit_spikes_search_

////Loop through L1EmParticles
if (do_l1extraparticles_){
	 	   //L1 trigger menu
		   edm::ESHandle<L1GtTriggerMenu> menuRcd;
		   iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
		   const L1GtTriggerMenu* menu = menuRcd.product();

		   //l1trigger masked
		   edm::ESHandle< L1GtTriggerMask > l1GtTmAlgo;
		   iSetup.get< L1GtTriggerMaskAlgoTrigRcd >().get( l1GtTmAlgo );
		   std::vector<unsigned int> triggerMaskAlgoTrig = l1GtTmAlgo.product()->gtTriggerMask();

		   //BX
		   int bxNb = iEvent.bunchCrossing();
//		   cout << "BXnumber is " << bxNb << endl;
		   edm::Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
		   iEvent.getByLabel(m_l1GTReadoutRecTag_,L1GTRR);

		   bool isEcalL1 = false;
		   const unsigned int sizeOfDecisionWord(L1GTRR->decisionWord().size());
//		   cout << "size of L1 decision word is " << sizeOfDecisionWord << endl;

		   const  DecisionWord& gtDecisionWordBeforeMask = L1GTRR->decisionWord();
		   DecisionWord dWord2 = L1GTRR->decisionWord();

		   bool l1SingleEG5 = menu->gtAlgorithmResult( "L1_SingleEG5", gtDecisionWordBeforeMask);

		   bool processL1extraParticles=true;
		   if (do_l1EG5Cut_) {processL1extraParticles=l1SingleEG5;}

		   if(processL1extraParticles){

				//Handles for l1extraparticles
			   	//Online collection
				edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl ;
				edm::Handle< l1extra::L1EmParticleCollection > emIsolColl ;
				iEvent.getByLabel("l1extraParticlesOnline","NonIsolated", emNonisolColl ) ;
				iEvent.getByLabel("l1extraParticlesOnline","Isolated", emIsolColl ) ;
				//Emulated collection
				edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl_M ;
				edm::Handle< l1extra::L1EmParticleCollection > emIsolColl_M ;
				iEvent.getByLabel("l1extraParticles","NonIsolated", emNonisolColl_M ) ;
				iEvent.getByLabel("l1extraParticles","Isolated", emIsolColl_M ) ;

				//Energies and Ets of all candidates
				std::vector<double> l1CandidatesEnergies_Online;
				std::vector<double> l1CandidatesEnergies_Emulated;
				std::vector<double> l1CandidatesEts_Online;
				std::vector<double> l1CandidatesEts_Emulated;

				// Maps: Ets -> reference to L1 candidate
				map<double, l1extra::L1EmParticleCollection::const_iterator> map_Ets_candidates_Online;
				map<double, l1extra::L1EmParticleCollection::const_iterator> map_Ets_candidates_Emul;

				//// Fill histograms and Ets vectors
				//Browse online Isolated collection
				for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl->begin(); emItr != emIsolColl->end() ;++emItr) {

						h_l1em_Online_Particles_noBadTowerCut_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.479) && (emItr->eta() > -1.479) )
							h_l1em_Online_Particles_noBadTowerCut_barrel_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.567) && (emItr->eta() > -1.567) )
							h_l1em_Online_Particles_noBadTowerCut_eta1567_et->Fill(emItr->et() );

						// Cut out regions with hot towers
						if ( (emItr->gctEmCand()->regionId().ieta() == 14) && (emItr->gctEmCand()->regionId().iphi() == 7) )
							continue;
						if ( (emItr->gctEmCand()->regionId().ieta() == 9) && (emItr->gctEmCand()->regionId().iphi() == 16) )
													continue;

						h_l1em_Online_Particles_Isolated_et->Fill( emItr->et() );
						h_l1em_Online_Particles_Isolated_energy->Fill( emItr->energy() );
						h_l1em_Online_Particles_et->Fill( emItr->et() );
						// In the barrel
						if ( (emItr->eta() < 1.479) && (emItr->eta() > -1.479) ) {
								h_l1em_Online_Particles_barrel_et->Fill( emItr->et() );
								l1CandidatesEnergies_Online.push_back(emItr->energy());
								l1CandidatesEts_Online.push_back(emItr->et());
								map_Ets_candidates_Online[emItr->et()]=emItr;
						}
						if ( (emItr->eta() < 1.567) && (emItr->eta() > -1.567) )
								h_l1em_Online_Particles_eta1567_et->Fill( emItr->et() );
						h_l1em_Online_Particles_energy->Fill( emItr->energy() );
						h_l1em_Online_Particles_Isolated_map->Fill(emItr->eta(), emItr->phi());
						h_l1em_Online_Particles_Isolated_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																		emItr->gctEmCand()->regionId().iphi());
						h_l1em_Online_Particles_map->Fill(emItr->eta(), emItr->phi());
						h_l1em_Online_Particles_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																		emItr->gctEmCand()->regionId().iphi());

						if (emItr->et() == 63) h_l1em_saturated_Online_Particles_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																				emItr->gctEmCand()->regionId().iphi());
					}
				// Online non-isolated
				for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl->begin(); emItr != emNonisolColl->end() ;++emItr) {

						h_l1em_Online_Particles_noBadTowerCut_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.479) && (emItr->eta() > -1.479) )
							h_l1em_Online_Particles_noBadTowerCut_barrel_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.567) && (emItr->eta() > -1.567) )
							h_l1em_Online_Particles_noBadTowerCut_eta1567_et->Fill(emItr->et() );

						if ( (emItr->gctEmCand()->regionId().ieta() == 14) && (emItr->gctEmCand()->regionId().iphi() == 7) )
							continue;
						if ( (emItr->gctEmCand()->regionId().ieta() == 9) && (emItr->gctEmCand()->regionId().iphi() == 16) )
													continue;

						h_l1em_Online_Particles_NonIsolated_et->Fill( emItr->et() );
						h_l1em_Online_Particles_NonIsolated_energy->Fill( emItr->energy() );
						h_l1em_Online_Particles_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.479) && (emItr->eta() > -1.479) ) {
										h_l1em_Online_Particles_barrel_et->Fill( emItr->et() );
										l1CandidatesEnergies_Online.push_back(emItr->energy());
										l1CandidatesEts_Online.push_back(emItr->et());
										map_Ets_candidates_Online[emItr->et()]=emItr;
						}
						if ( (emItr->eta() < 1.567) && (emItr->eta() > -1.567) )
								h_l1em_Online_Particles_eta1567_et->Fill( emItr->et() );
						h_l1em_Online_Particles_energy->Fill( emItr->energy() );
						h_l1em_Online_Particles_NonIsolated_map->Fill(emItr->eta(), emItr->phi());
						h_l1em_Online_Particles_NonIsolated_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																		emItr->gctEmCand()->regionId().iphi());
						h_l1em_Online_Particles_map->Fill(emItr->eta(), emItr->phi());
						h_l1em_Online_Particles_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																		emItr->gctEmCand()->regionId().iphi());

						if (emItr->et() == 63) h_l1em_saturated_Online_Particles_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																				emItr->gctEmCand()->regionId().iphi());
					}
				// Emulated isolated
				for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl_M->begin(); emItr != emIsolColl_M->end() ;++emItr) {
						h_l1em_Emul_Particles_noBadTowerCut_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.479) && (emItr->eta() > -1.479) )
							h_l1em_Emul_Particles_noBadTowerCut_barrel_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.567) && (emItr->eta() > -1.567) )
							h_l1em_Emul_Particles_noBadTowerCut_eta1567_et->Fill(emItr->et() );

						if ( (emItr->gctEmCand()->regionId().ieta() == 14) && (emItr->gctEmCand()->regionId().iphi() == 7) )
							continue;
						if ( (emItr->gctEmCand()->regionId().ieta() == 9) && (emItr->gctEmCand()->regionId().iphi() == 16) )
													continue;


						h_l1em_Emul_Particles_Isolated_et->Fill( emItr->et() );
						h_l1em_Emul_Particles_Isolated_energy->Fill( emItr->energy() );
						h_l1em_Emul_Particles_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.479) && (emItr->eta() > -1.479) ){
										h_l1em_Emul_Particles_barrel_et->Fill( emItr->et() );
										l1CandidatesEnergies_Emulated.push_back(emItr->energy());
										l1CandidatesEts_Emulated.push_back(emItr->et());
										map_Ets_candidates_Emul[emItr->et()]=emItr;
						}
						if ( (emItr->eta() < 1.567) && (emItr->eta() > -1.567) )
								h_l1em_Emul_Particles_eta1567_et->Fill( emItr->et() );
						h_l1em_Emul_Particles_energy->Fill( emItr->energy() );
			//			l1CandidatesEnergies_Emulated.push_back(emItr->energy());
			//			l1CandidatesEts_Emulated.push_back(emItr->et());

						h_l1em_Emul_Particles_Isolated_map->Fill(emItr->eta(), emItr->phi());
						h_l1em_Emul_Particles_Isolated_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																		emItr->gctEmCand()->regionId().iphi());
						h_l1em_Emul_Particles_map->Fill(emItr->eta(), emItr->phi());
						h_l1em_Emul_Particles_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																		emItr->gctEmCand()->regionId().iphi());

						if (emItr->et() == 63) h_l1em_saturated_Emul_Particles_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																				emItr->gctEmCand()->regionId().iphi());

					}
				// Emulated non-isolated
				for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl_M->begin(); emItr != emNonisolColl_M->end() ;++emItr) {
						h_l1em_Emul_Particles_noBadTowerCut_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.479) && (emItr->eta() > -1.479) )
							h_l1em_Emul_Particles_noBadTowerCut_barrel_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.567) && (emItr->eta() > -1.567) )
							h_l1em_Emul_Particles_noBadTowerCut_eta1567_et->Fill(emItr->et() );

						if ( (emItr->gctEmCand()->regionId().ieta() == 14) && (emItr->gctEmCand()->regionId().iphi() == 7) )
							continue;
						if ( (emItr->gctEmCand()->regionId().ieta() == 9) && (emItr->gctEmCand()->regionId().iphi() == 16) )
													continue;

						h_l1em_Emul_Particles_NonIsolated_et->Fill( emItr->et() );
						h_l1em_Emul_Particles_NonIsolated_energy->Fill( emItr->energy() );
						h_l1em_Emul_Particles_et->Fill( emItr->et() );
						if ( (emItr->eta() < 1.479) && (emItr->eta() > -1.479) ){
										h_l1em_Emul_Particles_barrel_et->Fill( emItr->et() );
										l1CandidatesEnergies_Emulated.push_back(emItr->energy());
										l1CandidatesEts_Emulated.push_back(emItr->et());
										map_Ets_candidates_Emul[emItr->et()]=emItr;
						}
						if ( (emItr->eta() < 1.567) && (emItr->eta() > -1.567) )
								h_l1em_Emul_Particles_eta1567_et->Fill( emItr->et() );
						h_l1em_Emul_Particles_energy->Fill( emItr->energy() );
								h_l1em_Emul_Particles_NonIsolated_map->Fill(emItr->eta(), emItr->phi());
								h_l1em_Emul_Particles_NonIsolated_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																				emItr->gctEmCand()->regionId().iphi());
								h_l1em_Emul_Particles_map->Fill(emItr->eta(), emItr->phi());
								h_l1em_Emul_Particles_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																	  emItr->gctEmCand()->regionId().iphi());
								if (emItr->et() == 63) h_l1em_saturated_Emul_Particles_index_map->Fill(emItr->gctEmCand()->regionId().ieta(),
																						emItr->gctEmCand()->regionId().iphi());
					}

				// Sort collection of L1 candidates Ets and fill highest Et histograms
				if (l1CandidatesEnergies_Online.size() > 0){
					std::sort(l1CandidatesEnergies_Online.begin(), l1CandidatesEnergies_Online.end());
				//	cout<<"L1 Online energy highest: "<< l1CandidatesEnergies_Online.back()<<endl;
					h_l1em_Online_Particles_Highest_energy->Fill(l1CandidatesEnergies_Online.back());
				}
				if (l1CandidatesEnergies_Emulated.size() > 0){
					std::sort(l1CandidatesEnergies_Emulated.begin(), l1CandidatesEnergies_Emulated.end());
				//	cout<<"L1 Emulated energy highest: "<< l1CandidatesEnergies_Emulated.back()<<endl;
					h_l1em_Emul_Particles_Highest_energy->Fill(l1CandidatesEnergies_Emulated.back());
				}

				if (l1CandidatesEts_Online.size() > 0){
					std::sort(l1CandidatesEts_Online.begin(), l1CandidatesEts_Online.end());
				//	cout<<"L1 Online et highest: "<< l1CandidatesEts_Online.back()<<endl;
					h_l1em_Online_Particles_Highest_et->Fill(l1CandidatesEts_Online.back());
				}
				if (l1CandidatesEts_Emulated.size() > 0){
					std::sort(l1CandidatesEts_Emulated.begin(), l1CandidatesEts_Emulated.end());
				//	cout<<"L1 Emulated et highest: "<< l1CandidatesEts_Emulated.back()<<endl;
					h_l1em_Emul_Particles_Highest_et->Fill(l1CandidatesEts_Emulated.back());
				}

				h_l1Candides_N_Online->Fill(l1CandidatesEts_Online.size());
				h_l1Candides_N_Emul->Fill(l1CandidatesEts_Emulated.size());
				h_l1Candides_N_difference->Fill(l1CandidatesEts_Online.size()-l1CandidatesEts_Emulated.size());

//				if (l1CandidatesEts_Online.size() > 0){
//					if (l1CandidatesEts_Emulated.size()==0){
				if (l1CandidatesEts_Emulated.size()!=0){
					if ( (l1CandidatesEts_Emulated.back() >= 12 ) && ( l1CandidatesEts_Emulated.back() <= 18 ) ){
						bool cout_event = false;
						if ( l1CandidatesEts_Online.size()!=0 ) {
							if ( (l1CandidatesEts_Online.back() < 12 ) || (l1CandidatesEts_Online.back() > 18 ) ) cout_event = true;
						}
						if ( l1CandidatesEts_Online.size()==0 )	cout_event = true;
						if (cout_event){
							cout << "Emul et [12,18]; online et <12 || >18 :" << endl;
							cout << "EventN : "<< iEvent.id().event() << endl << endl;
							cout << "online highest Et:" << l1CandidatesEts_Online.back() << endl;
							cout << "emulated highest Et:" << l1CandidatesEts_Emulated.back() << endl;
							for (auto it: l1CandidatesEts_Online){
								cout << "OnlineEt : " << it << endl;
							}
//							cout << "Online candidates:" << endl;
//							cout_all_L1Extra_candidates(map_Ets_candidates_Online);
//							cout << endl;
//							cout << "Emul candidates:" << endl;
//							cout_all_L1Extra_candidates(map_Ets_candidates_Emul);
//							cout << endl;
							cout_all_L1Extra_candidates(emIsolColl, emNonisolColl, emIsolColl_M, emNonisolColl_M);

//							cout << " Online : "<< endl;
//							for (auto it : l1CandidatesEts_Online){
//								cout << " Et : " << it << endl;
//							}
//							cout << " Emulated : "<< endl;
//							for (auto it : l1CandidatesEts_Emulated){
//								cout << " Et : " << it << endl;
						}
					}
				}

											//HLT TriggerPath
												edm::Handle<edm::TriggerResults> triggerBits;
											//				edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
											//				edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

												iEvent.getByToken(triggerBits_, triggerBits);
											//				iEvent.getByToken(triggerObjects_, triggerObjects);
											//				iEvent.getByToken(triggerPrescales_, triggerPrescales);

												const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
												std::cout << "\n === TRIGGER PATHS === " << std::endl;
												for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
													std::cout << "Trigger " << names.triggerName(i) <<
											//			                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
															": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
															<< std::endl;
												}


//				if ( (l1CandidatesEts_Online.size() == 1 ) &&
//						l1CandidatesEts_Online.size()!=l1CandidatesEts_Emulated.size()){
//							cout << "l1CandidatesEnergies_Online.size() :"<< l1CandidatesEnergies_Online.size() << endl;
//							cout << "l1CandidatesEnergies_Emulated.size() :"<< l1CandidatesEnergies_Emulated.size() << endl;
//							cout << "l1CandidatesEts_Online.size() :"<< l1CandidatesEts_Online.size() << endl;
//							cout << "l1CandidatesEts_Emulated.size() :"<< l1CandidatesEts_Emulated.size() << endl<< endl;
//

//
//				}


			//	//HLT TriggerPath
			//	edm::Handle<edm::TriggerResults> triggerBits;
			////				edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
			////				edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
			//
			//	iEvent.getByToken(triggerBits_, triggerBits);
			////				iEvent.getByToken(triggerObjects_, triggerObjects);
			////				iEvent.getByToken(triggerPrescales_, triggerPrescales);
			//
			//	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
			//    std::cout << "\n === TRIGGER PATHS === " << std::endl;
			//    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
			//        std::cout << "Trigger " << names.triggerName(i) <<
			////			                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
			//                ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
			//                << std::endl;
			//    }

				// Make ratios of Ets from Highest Et candidates.
				if (l1CandidatesEnergies_Online.size() > 0
					&& l1CandidatesEnergies_Emulated.size() > 0
					&& l1CandidatesEts_Online.size() > 0
					&& l1CandidatesEts_Emulated.size() > 0
					){
						if (l1CandidatesEnergies_Online.back() > 0
								&& l1CandidatesEts_Online.back() > 0){
							double highest_energies_ratio = l1CandidatesEnergies_Emulated.back() / l1CandidatesEnergies_Online.back();
							double highest_ets_ratio = l1CandidatesEts_Emulated.back() / l1CandidatesEts_Online.back();

				//			cout << "highest_energies_ratio = " << highest_energies_ratio << endl;
				//			cout << "highest_ets_ratio = " << highest_ets_ratio << endl;

							h_l1em_highest_energies_ratio->Fill(highest_energies_ratio);
							h_l1em_highest_ets_ratio->Fill(highest_ets_ratio);

							if (l1CandidatesEnergies_Online.back() == l1CandidatesEnergies_Emulated.back()){
				//				cout << "equal" << endl;
								h_l1em_highest_energies_ratio_equal->Fill(highest_energies_ratio);
								h_l1em_highest_ets_ratio_equal->Fill(highest_ets_ratio);
							}
							else{
				//				cout << "not equal"<< endl;
								h_l1em_highest_energies_ratio_not_equal->Fill(highest_energies_ratio);
								h_l1em_highest_ets_ratio_not_equal->Fill(highest_ets_ratio);
							}
						}

				}

//				if (iEvent.id().event() == 13077624){
//					cout_all_L1Extra_candidates(emIsolColl, emNonisolColl, emIsolColl_M, emNonisolColl_M);
//				}

		   };//if processL1extraParticles
}//if do_l1extraparticles

//// cout RCT regions with offline spikes (RecHit sev=3 or sev=4)
//for (auto RCT_it : map_RCT_regs_with_offline_spikes){
//	if (RCT_it.second == true){
//		cout <<"Spike in region : " << endl;
//		cout << "reg ieta : " << RCT_it.first.ieta() << " ; "
//				<< "reg iphi : " << RCT_it.first.iphi() << endl;
//	}
//}


//	mytree->Fill();


}

void OfflineSpikeCrystalToOnlineMatch::cout_L1Extra_candidate(l1extra::L1EmParticleCollection::const_iterator & it_l1extra){
	cout << "Et : " << it_l1extra->et()
		 << " ; energy :" << it_l1extra->energy() << endl
		 << "iEta : " << it_l1extra->gctEmCand()->regionId().ieta()
		 << " ; iPhi : "<< it_l1extra->gctEmCand()->regionId().iphi() << endl
		 << "eta : " << it_l1extra->eta()
		 << " ; phi : "<< it_l1extra->phi() << endl
		 << "Offline spike in region : " <<  map_RCT_regs_with_offline_spikes[it_l1extra->gctEmCand()->regionId()]
		 << endl << endl;

};

void OfflineSpikeCrystalToOnlineMatch::cout_L1Extra_candidate(l1extra::L1EmParticle & it_l1extra){
	cout << "Et : " << it_l1extra.et()
		 << " ; energy :" << it_l1extra.energy() << endl
		 << "iEta : " << it_l1extra.gctEmCand()->regionId().ieta()
		 << " ; iPhi : "<< it_l1extra.gctEmCand()->regionId().iphi() << endl
		 << "eta : " << it_l1extra.eta()
		 << " ; phi : "<< it_l1extra.phi() << endl
		 << "Offline spike in region : " <<  map_RCT_regs_with_offline_spikes[it_l1extra.gctEmCand()->regionId()]
		 << endl << endl;

};

//void OfflineSpikeCrystalToOnlineMatch::cout_all_L1Extra_candidates(map<double, l1extra::L1EmParticleCollection::const_iterator> & map_candidates){
//	for (auto it_map: map_candidates){
//		cout_L1Extra_candidate(it_map.second);
//	}
//
//};

void OfflineSpikeCrystalToOnlineMatch::cout_all_L1Extra_candidates(edm::Handle< l1extra::L1EmParticleCollection > OnlineIso,
																   edm::Handle< l1extra::L1EmParticleCollection > OnlineNonIso,
																   edm::Handle< l1extra::L1EmParticleCollection > EmulatedIso,
																   edm::Handle< l1extra::L1EmParticleCollection > EmulatedNonIso){
	cout << "Online Isolated:" << endl << endl;
	for (l1extra::L1EmParticleCollection::const_iterator emItr = OnlineIso->begin(); emItr != OnlineIso->end() ;++emItr) cout_L1Extra_candidate(emItr);
	cout << "Online NonIsolated:" << endl << endl;
	for (l1extra::L1EmParticleCollection::const_iterator emItr = OnlineNonIso->begin(); emItr != OnlineNonIso->end() ;++emItr) cout_L1Extra_candidate(emItr);
	cout << "Emulated Isolated:" << endl << endl;
	for (l1extra::L1EmParticleCollection::const_iterator emItr = EmulatedIso->begin(); emItr != EmulatedIso->end() ;++emItr) cout_L1Extra_candidate(emItr);
	cout << "Emulated NonIsolated:" << endl << endl;
	for (l1extra::L1EmParticleCollection::const_iterator emItr = EmulatedNonIso->begin(); emItr != EmulatedNonIso->end() ;++emItr) cout_L1Extra_candidate(emItr);

};

//Emulate L1 trigger decision. Taken from L1Studies/EGamma/Macro/Efficiency/selectPairs.C
void OfflineSpikeCrystalToOnlineMatch::fireL1(int L1noniso, int L1iso, vector<int> & firedEG, vector<int> EG) {

  const int n=EG.size();

  int fired = -1;

  for(int i=0 ; i<n ; i++) {
    if( L1iso >= EG[n-1-i] || L1noniso >= EG[n-1-i] ) {
    //if(i>3) {
      fired = n-i-1;
      //cout << "L1iso=" << L1iso << " | L1noniso=" << L1noniso << " | fired EG : " << EG[n-1-i] << endl;
      break;
    }
  }

  if(fired>=0)
    for(int i=0 ; i<fired+1 ; i++)
      firedEG[i]=1;

}

//Emulate L1 trigger decision. Taken from L1Studies/EGamma/Macro/Efficiency/selectPairs.C
void OfflineSpikeCrystalToOnlineMatch::globalFireL1_Normal(double eleRCT_eT, int noniso, int iso,
			 vector<int> & firedEG_N, vector<int> menu) {

  if( eleRCT_eT > 2 ) {// interesting regions <=> et > 2GeV
    fireL1( noniso , iso , firedEG_N, menu );
  }
}

// ------------ method called once each job just before starting event loop  ------------
void OfflineSpikeCrystalToOnlineMatch::beginJob() {
//	file_ = new TFile(histogramFile_.c_str(), "RECREATE");
//	file_->cd();
//	h_ttonline_et = new TH1F("h_ttonline_et", "h_ttonline_et", 256, 0, 256);
//	h_ttonline_et_nozerocut = new TH1F("h_ttonline_et_nozerocut", "h_ttonline_et_nozerocut", 306, -50, 256);
//	h_ttemulator_et = new TH1F("h_ttemulator_et", "h_ttemulator_et", 256, 0, 256);
//	h_ttonline_saturated_to_ttoffline_match_et = new TH1F(
//			"h_ttonline_saturated_to_ttoffline_match_et",
//			"h_ttonline_saturated_to_ttoffline_match_et", 256, 0, 256);
//	h_rechits_et = new TH1F("h_rechits_et", "h_rechits_et", 256, 0, 256);
//	h_rechit_spikes_et = new TH1F("h_rechit_spikes_et", "h_rechit_spikes_et",
//			256, 0, 256);
//	h_rechit_spikes_match_to_ttonline_et = new TH1F(
//			"h_rechit_spikes_match_to_ttonline_et",
//			"h_rechit_spikes_match_to_ttonline_et", 256, 0, 256);
//	h_rechit_spikes_to_ttonline_et_sFGVB0On = new TH1F(
//			"h_rechit_spikes_to_ttonline_et_sFGVB0On",
//			"h_rechit_spikes_to_ttonline_et_sFGVB0On", 256, 0, 256);
//	h_rechit_spikes_to_ttonline_et_sFGVB1On = new TH1F(
//			"h_rechit_spikes_to_ttonline_et_sFGVB1On",
//			"h_rechit_spikes_to_ttonline_et_sFGVB1On", 256, 0, 256);
////	h_TPGsaturated = new TH2F("h_TPGsaturated", "h_TPGsaturated", 60, -30, 30, 72, 0, 72);
//
////	h_tt_amplitude_online = new TH1F("h_tt_amplitude_online", "h_tt_amplitude_online", 306, -50, 256);
////	h_tt_amplitude_calibrated_online = new TH1F("h_tt_amplitude_calibrated_online", "h_tt_amplitude_calibrated_online", 306, -50, 128);
////
////	h_tt_amplitude_online_above_zero = new TH1F("h_tt_amplitude_online_above_zero", "h_tt_amplitude_online_above_zero", 256, 0, 256);
////	h_tt_amplitude_icalib_ADCtoGeV_online_above_zero = new TH1F("h_tt_amplitude_icalib_ADCtoGeV_online_above_zero", "h_tt_amplitude_icalib_ADCtoGeV_online_above_zero", 256, 0, 256);
////
////	h_tt_spikes_amplitude_online = new TH1F("h_tt_spikes_amplitude_online", "h_tt_spikes_amplitude_online", 306, -25, 256);
////	h_tt_spikes_amplitude_icalib_ADCtoGeV_online = new TH1F("h_tt_spikes_amplitude_icalib_ADCtoGeV_online", "h_tt_spikes_amplitude_icalib_ADCtoGeV_online", 306, -25, 128);
////	h_rechit_tt_amplitude_online = new TH1F("h_rechit_tt_amplitude_online", "h_rechit_tt_amplitude_online", 306, -50, 256);
////	h_rechit_tt_amplitude_calibrated_online = new TH1F("h_rechit_tt_amplitude_calibrated_online", "h_rechit_tt_amplitude_calibrated_online", 306, -25, 128);
//
//	h_tt_amplitude_online = new TH1F("h_tt_amplitude_online", "h_tt_amplitude_online", 256, 0, 256);
//	h_tt_amplitude_calibrated_online = new TH1F("h_tt_amplitude_calibrated_online", "h_tt_amplitude_calibrated_online", 256, 0, 256);
//
//	h_tt_amplitude_online_above_zero = new TH1F("h_tt_amplitude_online_above_zero", "h_tt_amplitude_online_above_zero", 256, 0, 256);
//	h_tt_amplitude_icalib_ADCtoGeV_online_above_zero = new TH1F("h_tt_amplitude_icalib_ADCtoGeV_online_above_zero", "h_tt_amplitude_icalib_ADCtoGeV_online_above_zero", 256, 0, 256);
//
//	h_tt_spikes_amplitude_online = new TH1F("h_tt_spikes_amplitude_online", "h_tt_spikes_amplitude_online", 256, 0, 256);
//	h_tt_spikes_amplitude_icalib_online = new TH1F("h_tt_spikes_amplitude_icalib_online", "h_tt_spikes_amplitude_icalib_online", 256, 0, 256);
//	h_tt_spikes_amplitude_icalib_ADCtoGeV_online = new TH1F("h_tt_spikes_amplitude_icalib_ADCtoGeV_online", "h_tt_spikes_amplitude_icalib_ADCtoGeV_online", 256, 0, 256);
//
//	h_tt_amplitude_icalib_ADCtoGeV_online_above_zero_sintheta = new TH1F("h_tt_amplitude_icalib_ADCtoGeV_online_above_zero_sintheta", "h_tt_amplitude_icalib_ADCtoGeV_online_above_zero_sintheta", 256, 0, 256);
//	h_tt_spikes_amplitude_icalib_ADCtoGeV_online_sintheta = new TH1F("h_tt_spikes_amplitude_icalib_ADCtoGeV_online_sintheta", "h_tt_spikes_amplitude_icalib_ADCtoGeV_online_sintheta", 256, 0, 256);

//	h_temp_hist = new TH1F("h_temp_hist", "debug_hist", 200, 0, -1);

}
// ------------ method called once each job just after ending the event loop  ------------
void OfflineSpikeCrystalToOnlineMatch::endJob() {

//	if (file_ != 0) {
//		file_->cd();
////		h_temp_hist->Write();
//
//		for (std::vector<TH2D*>::const_iterator hist_iter =
//				hvect_Spikes.begin(); hist_iter != hvect_Spikes.end();
//				hist_iter++) {
//			(*hist_iter)->Write();
//		}
//
//		std::cout << "Writing plots" << endl;
//		file_->Close();
//	}
//	file_ = 0;
}

// ------------ method called when starting to processes a run  ------------
/*
 void
 OfflineSpikeCrystalToOnlineMatch::beginRun(edm::Run const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when ending the processing of a run  ------------
/*
 void
 OfflineSpikeCrystalToOnlineMatch::endRun(edm::Run const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
 void
 OfflineSpikeCrystalToOnlineMatch::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
 void
 OfflineSpikeCrystalToOnlineMatch::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void OfflineSpikeCrystalToOnlineMatch::fillDescriptions(
		edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OfflineSpikeCrystalToOnlineMatch);
