#ifndef PF_DATA_TREE_PRODUCER_H
#define PF_DATA_TREE_PRODUCER_H

#include "TTree.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TFile.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

class EcalCleaningAlgo;

class PFDataTreeProducer : public edm::EDAnalyzer 
{
  public:
    explicit PFDataTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& e, edm::EventSetup const& iSetup);
    virtual void endJob();
    ~PFDataTreeProducer();
    void scan5x5(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo,const edm::ESHandle<CaloGeometry>  &geometry, int &nHits, float & totEt);


  private:
    edm::InputTag onlineTPTag_; //modif-alex
    edm::InputTag L1GtReadoutRecordTag;
    std::string histogramFile_;
    std::string jets_;
    std::string tracks_;
    std::string vertex_coll_;
    std::string jetcorr_;
    std::string ebhitcoll_;
    std::string eehitcoll_;
    std::string pfhitcoll1_;
    std::string pfhitcoll2_;
    bool isMC_;
    std::string photonProducer_;
    std::string photonCollection_;


    TFile* file_tt;

    TH2F* ietamapTT;
    TH2F* iphimapTT;


    TFile* file_;
    TTree* dataTree_;
    TH1D *histo_event_;

    //modif-alex
    TH1F *h_ttonline_et;
    TH1F *h_ttoffline_et;

    TH1D *eb_rechitenergy_;
    TH1D *ee_rechitenergy_;

    TH1D *ee_rechitenergy_notypeb_;
 
    TH1D *eb_rechitenergy_02;
    TH1D *eb_rechitenergy_04;
    TH1D *eb_rechitenergy_06;
    TH1D *eb_rechitenergy_08;
    TH1D *eb_rechitenergy_10;
    TH1D *eb_rechitenergy_12;
    TH1D *eb_rechitenergy_14;
    TH1D *eb_rechitenergy_148;


    TH1D *ee_rechitenergy_16;
    TH1D *ee_rechitenergy_18;
    TH1D *ee_rechitenergy_20;
    TH1D *ee_rechitenergy_22;
    TH1D *ee_rechitenergy_24;
    TH1D *ee_rechitenergy_26;
    TH1D *ee_rechitenergy_28;
    TH1D *ee_rechitenergy_30;



    TH1D *ee_rechitet_16;
    TH1D *ee_rechitet_18;
    TH1D *ee_rechitet_20;
    TH1D *ee_rechitet_22;
    TH1D *ee_rechitet_24;
    TH1D *ee_rechitet_26;
    TH1D *ee_rechitet_28;
    TH1D *ee_rechitet_30;



    TH1D *eb_rechitet_02;
    TH1D *eb_rechitet_04;
    TH1D *eb_rechitet_06;
    TH1D *eb_rechitet_08;
    TH1D *eb_rechitet_10;
    TH1D *eb_rechitet_12;
    TH1D *eb_rechitet_14;
    TH1D *eb_rechitet_148;



    TH2D *eb_rechitetvspu_05;
    TH2D *eb_rechitetvspu_10;
    TH2D *eb_rechitetvspu_15;

    TH2D *ee_rechitetvspu_20;
    TH2D *ee_rechitetvspu_25;
    TH2D *ee_rechitetvspu_30;


 

 
    TH1D *eb_rechitet_;
    TH1D *ee_rechitet_;


    TH2D *eb_rechiten_vs_eta;

    TH2D *eb_rechitet_vs_eta;



    TH2D *eep_rechitet_vs_eta;
    TH2D *eep_rechitet_vs_phi;

    TH2D *eem_rechitet_vs_eta;
    TH2D *eem_rechitet_vs_phi;



    TH2D *eep_rechiten_vs_eta;
    TH2D *eep_rechiten_vs_phi;

    TH2D *eem_rechiten_vs_eta;
    TH2D *eem_rechiten_vs_phi;


    TH2F *ebocc; 
    TH2F *eboccgt1;
    TH2F *eboccgt1et;
    TH2F *eboccet;
    TH2F *eboccetgt1et;
    TH2F *eboccen;
    TH2F *eboccengt1;

    TH2F *eeocc;
    TH2F *eeoccgt1;
    TH2F *eeoccgt1et;
    TH2F *eeoccet;
    TH2F *eeoccetgt1et;
    TH2F *eeoccen;
    TH2F *eeoccengt1;


    TH1D *numtp_vs_ieta;
    TH1D *numtp_vs_ieta_samp2;
    TH1D *numtp_geq1_vs_ieta;
    TH1D *numtp_etweighted_vs_ieta;



    TH2F *calotower_map;
    TH2F *tp_map;
    TH2F *tp_map_samp2;
    TH2F *tp_map2;
    TH2F *rechit_map;
    TH2F *rechit_map_time;

    TH2F *rechit_map_ixiy;
    TH2F *rechit_map_time_ixiy;


    TH2F *calotower_map_cumul;
    TH2F *tpmap_cumul;
    TH2F *tpmap_cumul2;
    TH2F *rechitmap_cumul;


    TH2D *matched_timevset;
    TH2D *unmatched_timevset;

    TH2D *matched_timevsieta;
    TH2D *unmatched_timevsieta;

    TH2D *matched_timevsieta_etweighted;
    TH2D *unmatched_timevsieta_etweighted;
 
    TH3D *matched_timevsetvsieta;
    TH3D *unmatched_timevsetvsieta;

    TH1D *ietamatched;
    TH1D *ietatotal;




    TH2F *tower_map_cumul_matched;
    TH2F *tower_map_cumul_matched2;

    TH2F *tower_map_cumul_matched_rh;
    TH2F *tower_map_cumul_matched_rh2;


    TH2F *tower_map_cumul_unmatched1;
    TH2F *tower_map_cumul_unmatched2;


    TH1F *rechiteta_all;
    TH1F *rechiteta_gt1et;
    TH1F *rechiteta_etweight;
    TH1F *rechiteta_etweight_gt1et;

    TH1F *rechiteta_gt1et_pu10;
    TH1F *rechiteta_gt1et_pu20;
    TH1F *rechiteta_gt1et_pu30;




    TH1F *calotowereta_all;
    TH1F *calotowereta_gt1et;
    TH1F *calotowereta_etweight;
    TH1F *calotowereta_etweight_gt1et;

    TH1F *calotowereta_gt1et_pu10;
    TH1F *calotowereta_gt1et_pu20;
    TH1F *calotowereta_gt1et_pu30;



    TH1F *sceta_all_pueq01;
    TH1F *sceta_severity0_pueq01;

    TH1F *sceta_all_pueq02;
    TH1F *sceta_severity0_pueq02;

    TH1F *sceta_all_pueq03;
    TH1F *sceta_severity0_pueq03;

    TH1F *sceta_all_pueq04;
    TH1F *sceta_severity0_pueq04;

    TH1F *sceta_all_pueq05;
    TH1F *sceta_severity0_pueq05;

    TH1F *sceta_all_pueq06;
    TH1F *sceta_severity0_pueq06;

    TH1F *sceta_all_pueq07;
    TH1F *sceta_severity0_pueq07;

    TH1F *sceta_all_pueq08;
    TH1F *sceta_severity0_pueq08;

    TH1F *sceta_all_pueq09;
    TH1F *sceta_severity0_pueq09;




    TH1F *sceta_all;
    TH1F *sceta_severity0;
    TH1F *sceta_koot0;


    TH1F *sceta_all_gt2;
    TH1F *sceta_severity0_gt2;
    TH1F *sceta_koot0_gt2;


    TH1F *sceta_all_gt5;
    TH1F *sceta_severity0_gt5;
    TH1F *sceta_koot0_gt5;


    TH1F *sceta_all_gt10;
    TH1F *sceta_severity0_gt10;
    TH1F *sceta_koot0_gt10;


    TH1F *scet_eb_all;
    TH1F *scet_eb_severity0;
    TH1F *scet_eb_koot0;

    TH1F *scet_ee_all;
    TH1F *scet_ee_severity0;
    TH1F *scet_ee_koot0;

    TH1F *scsumet_eb_all;
    TH1F *scsumet_eb_severity0;
    TH1F *scsumet_eb_koot0;

    TH1F *scsumet_ee_all;
    TH1F *scsumet_ee_severity0;
    TH1F *scsumet_ee_koot0;


    TH1F *scsumet_eb_all_gt2;
    TH1F *scsumet_eb_severity0_gt2;
    TH1F *scsumet_eb_koot0_gt2;

    TH1F *scsumet_ee_all_gt2;
    TH1F *scsumet_ee_severity0_gt2;
    TH1F *scsumet_ee_koot0_gt2;

    TH1F *scsumet_eb_all_gt5;
    TH1F *scsumet_eb_severity0_gt5;
    TH1F *scsumet_eb_koot0_gt5;

    TH1F *scsumet_ee_all_gt5;
    TH1F *scsumet_ee_severity0_gt5;
    TH1F *scsumet_ee_koot0_gt5;

    TH1F *scsumet_eb_all_gt10;
    TH1F *scsumet_eb_severity0_gt10;
    TH1F *scsumet_eb_koot0_gt10;

    TH1F *scsumet_ee_all_gt10;
    TH1F *scsumet_ee_severity0_gt10;
    TH1F *scsumet_ee_koot0_gt10;


    TH2F *scocc_eb_gt50;
    TH2F *scocc_ee_gt50;


    TH1F *scet_eb_all_eta15;
    TH1F *scet_eb_all_eta20;
    TH1F *scet_eb_all_eta25;
 
    TH1F *scet_ee_all_eta15;
    TH1F *scet_ee_all_eta20;
    TH1F *scet_ee_all_eta25;



    TH1F *scet_eb_all_eta15_pu10;
    TH1F *scet_eb_all_eta20_pu10;
    TH1F *scet_eb_all_eta25_pu10;
 
    TH1F *scet_ee_all_eta15_pu10;
    TH1F *scet_ee_all_eta20_pu10;
    TH1F *scet_ee_all_eta25_pu10;

    TH1F *scet_eb_all_eta15_pu20;
    TH1F *scet_eb_all_eta20_pu20;
    TH1F *scet_eb_all_eta25_pu20;
 
    TH1F *scet_ee_all_eta15_pu20;
    TH1F *scet_ee_all_eta20_pu20;
    TH1F *scet_ee_all_eta25_pu20;

    TH1F *scet_eb_all_eta15_pu30;
    TH1F *scet_eb_all_eta20_pu30;
    TH1F *scet_eb_all_eta25_pu30;
 
    TH1F *scet_ee_all_eta15_pu30;
    TH1F *scet_ee_all_eta20_pu30;
    TH1F *scet_ee_all_eta25_pu30;


    TH1F *scet_eb_all_eta15_pueq10;
    TH1F *scet_eb_all_eta20_pueq10;
    TH1F *scet_eb_all_eta25_pueq10;
 
    TH1F *scet_ee_all_eta15_pueq10;
    TH1F *scet_ee_all_eta20_pueq10;
    TH1F *scet_ee_all_eta25_pueq10;

    TH1F *scet_eb_all_eta15_pueq20;
    TH1F *scet_eb_all_eta20_pueq20;
    TH1F *scet_eb_all_eta25_pueq20;
 
    TH1F *scet_ee_all_eta15_pueq20;
    TH1F *scet_ee_all_eta20_pueq20;
    TH1F *scet_ee_all_eta25_pueq20;
 


    TH1F *scet_eb_all_eta15_pueq01;
    TH1F *scet_eb_all_eta20_pueq01;
    TH1F *scet_eb_all_eta25_pueq01;
 
    TH1F *scet_ee_all_eta15_pueq01;
    TH1F *scet_ee_all_eta20_pueq01;
    TH1F *scet_ee_all_eta25_pueq01;


    TH1F *scsumet_eb_all_eta15_pueq01;
    TH1F *scsumet_eb_all_eta20_pueq01;
    TH1F *scsumet_eb_all_eta25_pueq01;
 
    TH1F *scsumet_ee_all_eta15_pueq01;
    TH1F *scsumet_ee_all_eta20_pueq01;
    TH1F *scsumet_ee_all_eta25_pueq01;



    TH1F *scet_eb_all_eta15_pueq02;
    TH1F *scet_eb_all_eta20_pueq02;
    TH1F *scet_eb_all_eta25_pueq02;
 
    TH1F *scet_ee_all_eta15_pueq02;
    TH1F *scet_ee_all_eta20_pueq02;
    TH1F *scet_ee_all_eta25_pueq02;


    TH1F *scsumet_eb_all_eta15_pueq02;
    TH1F *scsumet_eb_all_eta20_pueq02;
    TH1F *scsumet_eb_all_eta25_pueq02;
 
    TH1F *scsumet_ee_all_eta15_pueq02;
    TH1F *scsumet_ee_all_eta20_pueq02;
    TH1F *scsumet_ee_all_eta25_pueq02;



    TH1F *scet_eb_all_eta15_pueq03;
    TH1F *scet_eb_all_eta20_pueq03;
    TH1F *scet_eb_all_eta25_pueq03;
 
    TH1F *scet_ee_all_eta15_pueq03;
    TH1F *scet_ee_all_eta20_pueq03;
    TH1F *scet_ee_all_eta25_pueq03;


    TH1F *scsumet_eb_all_eta15_pueq03;
    TH1F *scsumet_eb_all_eta20_pueq03;
    TH1F *scsumet_eb_all_eta25_pueq03;
 
    TH1F *scsumet_ee_all_eta15_pueq03;
    TH1F *scsumet_ee_all_eta20_pueq03;
    TH1F *scsumet_ee_all_eta25_pueq03;



    TH1F *scet_eb_all_eta15_pueq04;
    TH1F *scet_eb_all_eta20_pueq04;
    TH1F *scet_eb_all_eta25_pueq04;
 
    TH1F *scet_ee_all_eta15_pueq04;
    TH1F *scet_ee_all_eta20_pueq04;
    TH1F *scet_ee_all_eta25_pueq04;


    TH1F *scsumet_eb_all_eta15_pueq04;
    TH1F *scsumet_eb_all_eta20_pueq04;
    TH1F *scsumet_eb_all_eta25_pueq04;
 
    TH1F *scsumet_ee_all_eta15_pueq04;
    TH1F *scsumet_ee_all_eta20_pueq04;
    TH1F *scsumet_ee_all_eta25_pueq04;



    TH1F *scet_eb_all_eta15_pueq05;
    TH1F *scet_eb_all_eta20_pueq05;
    TH1F *scet_eb_all_eta25_pueq05;
 
    TH1F *scet_ee_all_eta15_pueq05;
    TH1F *scet_ee_all_eta20_pueq05;
    TH1F *scet_ee_all_eta25_pueq05;


    TH1F *scsumet_eb_all_eta15_pueq05;
    TH1F *scsumet_eb_all_eta20_pueq05;
    TH1F *scsumet_eb_all_eta25_pueq05;
 
    TH1F *scsumet_ee_all_eta15_pueq05;
    TH1F *scsumet_ee_all_eta20_pueq05;
    TH1F *scsumet_ee_all_eta25_pueq05;



    TH1F *scet_eb_all_eta15_pueq06;
    TH1F *scet_eb_all_eta20_pueq06;
    TH1F *scet_eb_all_eta25_pueq06;
 
    TH1F *scet_ee_all_eta15_pueq06;
    TH1F *scet_ee_all_eta20_pueq06;
    TH1F *scet_ee_all_eta25_pueq06;


    TH1F *scsumet_eb_all_eta15_pueq06;
    TH1F *scsumet_eb_all_eta20_pueq06;
    TH1F *scsumet_eb_all_eta25_pueq06;
 
    TH1F *scsumet_ee_all_eta15_pueq06;
    TH1F *scsumet_ee_all_eta20_pueq06;
    TH1F *scsumet_ee_all_eta25_pueq06;



    TH1F *scet_eb_all_eta15_pueq07;
    TH1F *scet_eb_all_eta20_pueq07;
    TH1F *scet_eb_all_eta25_pueq07;
 
    TH1F *scet_ee_all_eta15_pueq07;
    TH1F *scet_ee_all_eta20_pueq07;
    TH1F *scet_ee_all_eta25_pueq07;


    TH1F *scsumet_eb_all_eta15_pueq07;
    TH1F *scsumet_eb_all_eta20_pueq07;
    TH1F *scsumet_eb_all_eta25_pueq07;
 
    TH1F *scsumet_ee_all_eta15_pueq07;
    TH1F *scsumet_ee_all_eta20_pueq07;
    TH1F *scsumet_ee_all_eta25_pueq07;



    TH1F *scet_eb_all_eta15_pueq08;
    TH1F *scet_eb_all_eta20_pueq08;
    TH1F *scet_eb_all_eta25_pueq08;
 
    TH1F *scet_ee_all_eta15_pueq08;
    TH1F *scet_ee_all_eta20_pueq08;
    TH1F *scet_ee_all_eta25_pueq08;


    TH1F *scsumet_eb_all_eta15_pueq08;
    TH1F *scsumet_eb_all_eta20_pueq08;
    TH1F *scsumet_eb_all_eta25_pueq08;
 
    TH1F *scsumet_ee_all_eta15_pueq08;
    TH1F *scsumet_ee_all_eta20_pueq08;
    TH1F *scsumet_ee_all_eta25_pueq08;



    TH1F *scet_eb_all_eta15_pueq09;
    TH1F *scet_eb_all_eta20_pueq09;
    TH1F *scet_eb_all_eta25_pueq09;
 
    TH1F *scet_ee_all_eta15_pueq09;
    TH1F *scet_ee_all_eta20_pueq09;
    TH1F *scet_ee_all_eta25_pueq09;


    TH1F *scsumet_eb_all_eta15_pueq09;
    TH1F *scsumet_eb_all_eta20_pueq09;
    TH1F *scsumet_eb_all_eta25_pueq09;
 
    TH1F *scsumet_ee_all_eta15_pueq09;
    TH1F *scsumet_ee_all_eta20_pueq09;
    TH1F *scsumet_ee_all_eta25_pueq09;







    TH1F *scsumet_eb_all_eta15;
    TH1F *scsumet_eb_all_eta20;
    TH1F *scsumet_eb_all_eta25;
 
    TH1F *scsumet_ee_all_eta15;
    TH1F *scsumet_ee_all_eta20;
    TH1F *scsumet_ee_all_eta25;



    TH1F *scsumet_eb_all_eta15_pu10;
    TH1F *scsumet_eb_all_eta20_pu10;
    TH1F *scsumet_eb_all_eta25_pu10;
 
    TH1F *scsumet_ee_all_eta15_pu10;
    TH1F *scsumet_ee_all_eta20_pu10;
    TH1F *scsumet_ee_all_eta25_pu10;

    TH1F *scsumet_eb_all_eta15_pu20;
    TH1F *scsumet_eb_all_eta20_pu20;
    TH1F *scsumet_eb_all_eta25_pu20;
 
    TH1F *scsumet_ee_all_eta15_pu20;
    TH1F *scsumet_ee_all_eta20_pu20;
    TH1F *scsumet_ee_all_eta25_pu20;

    TH1F *scsumet_eb_all_eta15_pu30;
    TH1F *scsumet_eb_all_eta20_pu30;
    TH1F *scsumet_eb_all_eta25_pu30;
 
    TH1F *scsumet_ee_all_eta15_pu30;
    TH1F *scsumet_ee_all_eta20_pu30;
    TH1F *scsumet_ee_all_eta25_pu30;


    TH1F *scsumet_eb_all_eta15_pueq10;
    TH1F *scsumet_eb_all_eta20_pueq10;
    TH1F *scsumet_eb_all_eta25_pueq10;
 
    TH1F *scsumet_ee_all_eta15_pueq10;
    TH1F *scsumet_ee_all_eta20_pueq10;
    TH1F *scsumet_ee_all_eta25_pueq10;

    TH1F *scsumet_eb_all_eta15_pueq20;
    TH1F *scsumet_eb_all_eta20_pueq20;
    TH1F *scsumet_eb_all_eta25_pueq20;
 
    TH1F *scsumet_ee_all_eta15_pueq20;
    TH1F *scsumet_ee_all_eta20_pueq20;
    TH1F *scsumet_ee_all_eta25_pueq20;
 






    TH1D *eb_timing_0;
    TH1D *eb_timing_200;
    TH1D *eb_timing_400;
    TH1D *eb_timing_600;
    TH1D *eb_timing_800;
    TH1D *eb_timing_1000;
    TH1D *eb_timing_2000;
    TH1D *eb_timing_3000;
    TH1D *eb_timing_5000;

    TH1D *eb_r4_0;
    TH1D *eb_r4_200;
    TH1D *eb_r4_400;
    TH1D *eb_r4_600;
    TH1D *eb_r4_800;
    TH1D *eb_r4_1000;
    TH1D *eb_r4_2000;
    TH1D *eb_r4_3000;
    TH1D *eb_r4_5000;
  
    
    TH1D *eb_timing_r4_0;
    TH1D *eb_timing_r4_200;
    TH1D *eb_timing_r4_400;
    TH1D *eb_timing_r4_600;
    TH1D *eb_timing_r4_800;
    TH1D *eb_timing_r4_1000;
    TH1D *eb_timing_r4_2000;
    TH1D *eb_timing_r4_3000;
    TH1D *eb_timing_r4_5000;

    TH2D *eb_timing_vs_r4_0;
    TH2D *eb_timing_vs_r4_200;
    TH2D *eb_timing_vs_r4_400;
    TH2D *eb_timing_vs_r4_600;
    TH2D *eb_timing_vs_r4_800;
    TH2D *eb_timing_vs_r4_1000;
    TH2D *eb_timing_vs_r4_2000;
    TH2D *eb_timing_vs_r4_3000;
    TH2D *eb_timing_vs_r4_5000;

 
    TH2D *ebtime_vs_bxtrain_01;
    TH2D *ebtime_vs_bxtrain_05;
    TH2D *eetime_vs_bxtrain_01;
    TH2D *eetime_vs_bxtrain_05;


    TH2D *sceta_vs_bxtrain;
    TH2D *rechiteta_vs_bxtrain_01;
    TH2D *rechiteta_vs_bxtrain_05;

    TH2D *eb_digi_01;
    TH2D *ee_digi_01;

    TH2D *eb_digi_05;
    TH2D *ee_digi_05;

    TH2D *eb_digi_30;
    TH2D *ee_digi_30;

    TH2D *eb_digi_0105;
    TH2D *ee_digi_0105;
    
    TH2D *eb_digi_0530;
    TH2D *ee_digi_0530;


    TH2D *eb_digi_0105_vs_time;
    TH2D *ee_digi_0105_vs_time;
    
    TH2D *eb_digi_0530_vs_time;
    TH2D *ee_digi_0530_vs_time;


    TH2D *eb_digi_0105_vs_bxtrain;
    TH2D *ee_digi_0105_vs_bxtrain;
    
    TH2D *eb_digi_0530_vs_bxtrain;
    TH2D *ee_digi_0530_vs_bxtrain;


    TH2D *eb_digi_0105_vs_time_norm;
    TH2D *ee_digi_0105_vs_time_norm;
    
    TH2D *eb_digi_0530_vs_time_norm;
    TH2D *ee_digi_0530_vs_time_norm;


    TH2D *eb_digi_0105_vs_bxtrain_norm;
    TH2D *ee_digi_0105_vs_bxtrain_norm;
    
    TH2D *eb_digi_0530_vs_bxtrain_norm;
    TH2D *ee_digi_0530_vs_bxtrain_norm;






    TH2D *eb_digi_0105_vs_time_eta15;
    TH2D *ee_digi_0105_vs_time_eta15;
    
    TH2D *eb_digi_0530_vs_time_eta15;
    TH2D *ee_digi_0530_vs_time_eta15;


    TH2D *eb_digi_0105_vs_bxtrain_eta15;
    TH2D *ee_digi_0105_vs_bxtrain_eta15;
    
    TH2D *eb_digi_0530_vs_bxtrain_eta15;
    TH2D *ee_digi_0530_vs_bxtrain_eta15;


    TH2D *eb_digi_0105_vs_time_norm_eta15;
    TH2D *ee_digi_0105_vs_time_norm_eta15;
    
    TH2D *eb_digi_0530_vs_time_norm_eta15;
    TH2D *ee_digi_0530_vs_time_norm_eta15;


    TH2D *eb_digi_0105_vs_bxtrain_norm_eta15;
    TH2D *ee_digi_0105_vs_bxtrain_norm_eta15;
    
    TH2D *eb_digi_0530_vs_bxtrain_norm_eta15;
    TH2D *ee_digi_0530_vs_bxtrain_norm_eta15;




    TH2D *eb_digi_0105_vs_time_eta20;
    TH2D *ee_digi_0105_vs_time_eta20;
    
    TH2D *eb_digi_0530_vs_time_eta20;
    TH2D *ee_digi_0530_vs_time_eta20;


    TH2D *eb_digi_0105_vs_bxtrain_eta20;
    TH2D *ee_digi_0105_vs_bxtrain_eta20;
    
    TH2D *eb_digi_0530_vs_bxtrain_eta20;
    TH2D *ee_digi_0530_vs_bxtrain_eta20;


    TH2D *eb_digi_0105_vs_time_norm_eta20;
    TH2D *ee_digi_0105_vs_time_norm_eta20;
    
    TH2D *eb_digi_0530_vs_time_norm_eta20;
    TH2D *ee_digi_0530_vs_time_norm_eta20;


    TH2D *eb_digi_0105_vs_bxtrain_norm_eta20;
    TH2D *ee_digi_0105_vs_bxtrain_norm_eta20;
    
    TH2D *eb_digi_0530_vs_bxtrain_norm_eta20;
    TH2D *ee_digi_0530_vs_bxtrain_norm_eta20;



    TH2D *eb_digi_0105_vs_time_eta25;
    TH2D *ee_digi_0105_vs_time_eta25;
    
    TH2D *eb_digi_0530_vs_time_eta25;
    TH2D *ee_digi_0530_vs_time_eta25;


    TH2D *eb_digi_0105_vs_bxtrain_eta25;
    TH2D *ee_digi_0105_vs_bxtrain_eta25;
    
    TH2D *eb_digi_0530_vs_bxtrain_eta25;
    TH2D *ee_digi_0530_vs_bxtrain_eta25;


    TH2D *eb_digi_0105_vs_time_norm_eta25;
    TH2D *ee_digi_0105_vs_time_norm_eta25;
    
    TH2D *eb_digi_0530_vs_time_norm_eta25;
    TH2D *ee_digi_0530_vs_time_norm_eta25;


    TH2D *eb_digi_0105_vs_bxtrain_norm_eta25;
    TH2D *ee_digi_0105_vs_bxtrain_norm_eta25;
    
    TH2D *eb_digi_0530_vs_bxtrain_norm_eta25;
    TH2D *ee_digi_0530_vs_bxtrain_norm_eta25;


    TH1D *ct_et_all;
    TH1D *ct_et_matched;
    TH1D *ct_et_matched2;
    TH1D *ct_et_unmatched;

    TH1D *rh_et_matched_rh;
    TH1D *rh_et_matched_rh2;



    TH1D *rh_et_all_ee;
    TH1D *rh_et_all_ee24;




    TH1D *tp_et_all;
    TH1D *tp_et_matched;
    TH1D *tp_et_matched2;
 
    TH1D *tp_et_matched_rh;
    TH1D *tp_et_matched_rh2
;
    TH1D *tp_et_unmatched;


    TH1D *ctminustp_et_matched;
    TH1D *ctminustp_et_matched2;

    TH1D *rhminustp_et_matched_rh;
    TH1D *rhminustp_et_matched_rh2;



    TH1D *ct_et_all_eb;
    TH1D *ct_et_matched_eb;
    TH1D *ct_et_matched_eb2;


    TH1D *rh_et_matched_eb_rh;
    TH1D *rh_et_matched_eb_rh2;

    TH1D *ct_et_unmatched_eb;

    TH1D *tp_et_all_eb;
    TH1D *tp_et_matched_eb;
    TH1D *tp_et_matched_eb2;
    TH1D *tp_et_matched_eb_rh;
    TH1D *tp_et_matched_eb_rh2;
    TH1D *tp_et_unmatched_eb;


    TH1D *ctminustp_et_matched_eb;
    TH1D *ctminustp_et_matched_eb2;

    TH1D *rhminustp_et_matched_eb_rh;
    TH1D *rhminustp_et_matched_eb_rh2;


    TH1D *ct_et_all_ee;
    TH1D *ct_et_matched_ee;
    TH1D *ct_et_matched_ee2;
    TH1D *ct_et_unmatched_ee;

    TH1D *rh_et_matched_ee_rh;
    TH1D *rh_et_matched_ee24_rh;
    TH1D *rh_et_matched_ee_rh2;



    TH1D *tp_et_all_ee;
    TH1D *tp_et_all_ee24;
    TH1D *tp_et_matched_ee;
    TH1D *tp_et_matched_ee2;
    TH1D *tp_et_matched_ee_rh;
    TH1D *tp_et_matched_ee24_rh;
    TH1D *tp_et_matched_ee_rh2;
    TH1D *tp_et_unmatched_ee;


    TH1D *ctminustp_et_matched_ee;
    TH1D *ctminustp_et_matched_ee2;

    TH1D *rhminustp_et_matched_ee_rh;
    TH1D *rhminustp_et_matched_ee_rh2;



    TH2D *ct_et_vs_ieta_all;
    TH2D *ct_et_vs_ieta_matched;
    TH2D *ct_et_vs_ieta_matched2;
    TH2D *ct_et_vs_ieta_unmatched;

    TH2D *rh_et_vs_ieta_matched_rh;
    TH2D *rh_et_vs_ieta_matched_rh2;



    TH2D *tp_et_vs_ieta_all;
    TH2D *tp_et_vs_ieta_matched;
    TH2D *tp_et_vs_ieta_matched2;
    TH2D *tp_et_vs_ieta_matched_rh;
    TH2D *tp_et_vs_ieta_matched_rh2;

    TH2D *tp_et_vs_ieta_unmatched;


    TH2D *ctminustp_et_vs_ieta_matched;
    TH2D *ctminustp_et_vs_ieta_matched2;


    TH2D *rhminustp_et_vs_ieta_matched_rh;
    TH2D *rhminustp_et_vs_ieta_matched_rh2;

    TH2D *tp_et_vs_ct_et_matched;
    TH2D *tp_et_vs_ct_et_matched2;

    TH2D *tp_et_vs_ct_et_matched_eb;
    TH2D *tp_et_vs_ct_et_matched_ee;
 
    TH2D *tp_et_vs_ct_et_matched_eb2;
    TH2D *tp_et_vs_ct_et_matched_ee2;




    TH2D *tp_et_vs_rh_et_matched_rh;
    TH2D *tp_et_vs_rh_et_matched_rh2;

    TH2D *tp_et_vs_rh_et_matched_eb_rh;
    TH2D *tp_et_vs_rh_et_matched_ee_rh;
 
    TH2D *tp_et_vs_rh_et_matched_eb_rh2;
    TH2D *tp_et_vs_rh_et_matched_ee_rh2;

    TH2D *tp_et_vs_rh_et_matched_ee_ieta1820_rh;
    TH2D *tp_et_vs_rh_et_matched_ee_ieta2122_rh;
    TH2D *tp_et_vs_rh_et_matched_ee_ieta2324_rh;
    TH2D *tp_et_vs_rh_et_matched_ee_ieta2526_rh;
    TH2D *tp_et_vs_rh_et_matched_ee_ieta2728_rh;





    int   rank_;
    int   ntrk;
    int   goodtrk;

    float vtx_x,vtx_y,vtx_z,vtx_x_err,vtx_y_err,vtx_z_err;
    float vtx_chi2,vtx_ndof;
    int   numvtx,vtx_ntracks,vtx_isfake,numgoodvtx;
    int   vtx_good;
    
    float scale;

    int   bit3,bit4,bit9,bit40,bit41,bit36,bit37,bit38,bit39,run,even,lumi,bx,orbit,bit0;
    int   eg1,eg2,eg5, algo124;


    int algo54;
    int algo55;
    int algo56;
    int algo57;
    int algo58;
    int algo59;
    int algo60;
    int algo61;
    int algo62;

    int algo106;
    int algo107;

    int numtp_eb, numtp_ee;
    int numtp_samp2_eb, numtp_samp2_ee;
    int numtp_geq1_eb, numtp_geq1_ee;

    float tpsumet_eb, tpsumet_ee;
    float tpsumet_samp2_eb, tpsumet_samp2_ee;

    float tpsumet_geq1_eb, tpsumet_geq1_ee;

    float tpsumet_mod1_eb, tpsumet_mod2_eb, tpsumet_mod3_eb, tpsumet_mod4_eb; 
    float tpsumet_geq1_mod1_eb, tpsumet_geq1_mod2_eb, tpsumet_geq1_mod3_eb, tpsumet_geq1_mod4_eb; 

    float tpsumet_ring1_ee, tpsumet_ring2_ee, tpsumet_ring3_ee;
    float tpsumet_geq1_ring1_ee, tpsumet_geq1_ring2_ee, tpsumet_geq1_ring3_ee;




    int   physdeclared;
    float time;
    float ptJet_,phiJet_,etaJet_,chfJet_,nhfJet_,cemfJet_,nemfJet_; 
    float energy_pf_,energyc_pf_,energyn_pf_,energyg_pf_;
    int   cmultiJet_,nmultiJet_,nrjets_,ncr_;   
    float energy_ecal, energy_hcal; 

    float ebmax, eemax, ebtime, eetime;
    int   eb_ieta,eb_iphi,ebhits,ebhits1GeV, ebflags, recoflag, sevlev;
    int   eeix,eeiy,eeiz,eehits,eehits1GeV, eeflags;
    float eb_eta,eb_phi,ebmaxet, eb_r9, eb_r4;
    float ee_eta,ee_phi,eemaxet, ee_r9;
    float ebchi2, ebchi2oot;

    float eb2chi2, eb2chi2oot;


      Int_t ee_kGood;
      Int_t ee_kPoorReco;
      Int_t ee_kOutOfTime;
      Int_t ee_kFaultyHardware;
      Int_t ee_kNoisy;
      Int_t ee_kPoorCalib;
      Int_t ee_kSaturated;
      Int_t ee_kLeadingEdgeRecovered;
      Int_t ee_kNeighboursRecovered;
      Int_t ee_kTowerRecovered;
      Int_t ee_kDead;
      Int_t ee_kKilled;
      Int_t ee_kTPSaturated;
      Int_t ee_kL1SpikeFlag;
      Int_t ee_kWeird;
      Int_t ee_kDiWeird;
      Int_t ee_kUnknown;



      int numcalotower;
   int numtrigprim;
   int numcalotower_eb;
   int numtrigprim_eb;
   int numcalotower_ee;
   int numtrigprim_ee;
   int numtrigprim_ee24;

   int numcalotower_notmatched1;
   int numtrigprim_notmatched1;
   int numcalotower_notmatched1_eb;
   int numtrigprim_notmatched1_eb;
   int numcalotower_notmatched1_ee;
   int numtrigprim_notmatched1_ee;
 

   int numcalotower_notmatched2;
   int numtrigprim_notmatched2;
   int numcalotower_notmatched2_eb;
   int numtrigprim_notmatched2_eb;
   int numcalotower_notmatched2_ee;
   int numtrigprim_notmatched2_ee;
 
   int numcalotower_matched;
   int numtrigprim_matched;
   int numcalotower_matched_eb;
   int numtrigprim_matched_eb;
   int numcalotower_matched_ee;
   int numtrigprim_matched_ee;
 

   int numcalotower_matched2;
   int numtrigprim_matched2;
   int numcalotower_matched_eb2;
   int numtrigprim_matched_eb2;
   int numcalotower_matched_ee2;
   int numtrigprim_matched_ee2;
 

   int numrechittower_matched_rh;
   int numtrigprim_matched_rh;
   int numrechittower_matched_eb_rh;
   int numtrigprim_matched_eb_rh;
   int numrechittower_matched_ee_rh;
   int numrechittower_matched_ee24_rh;
   int numtrigprim_matched_ee_rh;
   int numtrigprim_matched_ee24_rh;
 

   int numrechittower_matched_rh2;
   int numtrigprim_matched_rh2;
   int numrechittower_matched_eb_rh2;
   int numtrigprim_matched_eb_rh2;
   int numrechittower_matched_ee_rh2;
   int numtrigprim_matched_ee_rh2;
 




   float sumcalotower;
   float sumtrigprim;
   float sumcalotower_eb;
   float sumtrigprim_eb;
   float sumcalotower_ee;
   float sumtrigprim_ee;
   float sumtrigprim_ee24;


   float sumtrigprim_samp2_eb;
   float sumtrigprim_samp2_ee;



   float sumcalotower_notmatched1;
   float sumtrigprim_notmatched1;
   float sumcalotower_notmatched1_eb;
   float sumtrigprim_notmatched1_eb;
   float sumcalotower_notmatched1_ee;
   float sumtrigprim_notmatched1_ee;
 

   float sumcalotower_notmatched2;
   float sumtrigprim_notmatched2;
   float sumcalotower_notmatched2_eb;
   float sumtrigprim_notmatched2_eb;
   float sumcalotower_notmatched2_ee;
   float sumtrigprim_notmatched2_ee;
 
   float sumcalotower_matched;
   float sumtrigprim_matched;
   float sumcalotower_matched_eb;
   float sumtrigprim_matched_eb;
   float sumcalotower_matched_ee;
   float sumtrigprim_matched_ee;
 


   float sumcalotower_matched2;
   float sumtrigprim_matched2;
   float sumcalotower_matched_eb2;
   float sumtrigprim_matched_eb2;
   float sumcalotower_matched_ee2;
   float sumtrigprim_matched_ee2;
 

   float sumrechittower_matched_rh;
   float sumtrigprim_matched_rh;
   float sumrechittower_matched_eb_rh;
   float sumtrigprim_matched_eb_rh;
   float sumrechittower_matched_ee_rh;
   float sumtrigprim_matched_ee_rh;
   float sumrechittower_matched_ee24_rh;
   float sumtrigprim_matched_ee24_rh;
 
   float sumrechittower_matched_rh2;
   float sumtrigprim_matched_rh2;
   float sumrechittower_matched_eb_rh2;
   float sumtrigprim_matched_eb_rh2;
   float sumrechittower_matched_ee_rh2;
   float sumtrigprim_matched_ee_rh2;
 

   float sumtrigprim_ee_eta1820, sumtrigprim_ee_eta2122, 
     sumtrigprim_ee_eta2324,
     sumtrigprim_ee_eta2526, sumtrigprim_ee_eta2728;

   float sumtrigprim_ee_bxm2, sumtrigprim_ee_bxm1, sumtrigprim_ee_bx0, 
     sumtrigprim_ee_bxp1, sumtrigprim_ee_bxp2;


   float sumtrigprim_ee_matched_eta1820, sumtrigprim_ee_matched_eta2122, 
     sumtrigprim_ee_matched_eta2324, sumtrigprim_ee_matched_eta2526, 
     sumtrigprim_ee_matched_eta2728;


   float sumrechittower_ee_matched_eta1820, sumrechittower_ee_matched_eta2122, 
     sumrechittower_ee_matched_eta2324, sumrechittower_ee_matched_eta2526, 
     sumrechittower_ee_matched_eta2728;



      
      int ebhits1GeVet, eehits1GeVet;
      int ebhits2GeV, ebhits4GeV, eehits2GeV, eehits4GeV;
      int ebhits2GeVet, ebhits4GeVet, eehits2GeVet, eehits4GeVet;


    float ebmax2, eemax2, ebtime2, eetime2;
    int   eb_ieta2,eb_iphi2,ebhits1GeV2, ebflags2;
    int   eeix2,eeiy2,eeiz2,eehits1GeV2, eeflags2;
    float eb_eta2,eb_phi2,ebmaxet2, eb_r92, eb_r42;
    float ee_eta2,ee_phi2,eemaxet2, ee_r92;

    int adc01,adc02,adc03,adc04,adc05,adc06,adc07,adc08,adc09,adc10;
    int gain01,gain02,gain03,gain04,gain05,gain06,gain07,gain08,gain09,gain10;

    float eb_e9, eb_e25;

    int eephits, eemhits;

    int ebflag_kgood, ebflag_kpoorreco, ebflag_koutoftime, ebflag_kfake;


    float tmean_en,  terr_en;
    float tmean_sig, terr_sig;


    int r4count;

    int e2e9count_thresh0, e2e25count_thresh0;
    int e2e9count_thresh1, e2e25count_thresh1;

    int e2e9count_thresh0_nor4, e2e25count_thresh0_nor4;
    int e2e9count_thresh1_nor4, e2e25count_thresh1_nor4;


    int r4_algo_count;
    int e2e9_algo_count;
    int e2e9_algo_count_5_1;
    int e2e9_algo_count_5_0;

    float swisscross_algo, e2e9_algo;
    EcalCleaningAlgo * cleaningAlgo_;  

    float Emin_;
    int side_;

    const std::vector<int> badsc_;


    const std::vector<int> bunchstartbx_;


   // photon variables

    float phoEta, phoPhi, phoEt;
    float pho_maxen_xtal, pho_e3x3, pho_e5x5, pho_r9;
    float etaPhoSc;

    float pho_ntrksolidconedr4;
    float pho_ntrksumptsolidconedr4;
    float pho_ntrkhollowconedr4;
    float pho_ntrksumpthollowconedr4;

    float pho_ecaletconedr4;
    float pho_hcaletconedr4;
    float pho_sigmaietaieta;

    float pho_hovere;
    int phoBarrel;



    int numrechittower, numrechittower_eb, numrechittower_ee, numrechittower_ee24;
    float sumrechittower, sumrechittower_eb, sumrechittower_ee, sumrechittower_ee24;



    float eesum_gt1, eesum_gt2, eesum_gt4;
    float ebsum_gt1, ebsum_gt2, ebsum_gt4;



    float eesum_gt1et, eesum_gt2et, eesum_gt4et;
    float ebsum_gt1et, ebsum_gt2et, ebsum_gt4et;


    int ncalotower, ncalotowereb, ncalotoweree, ncalotowerhf;

    int ncalotowerebgt1, ncalotowerebgt2, ncalotowerebgt5, ncalotowerebgt10;

    int ncalotowereegt1,  ncalotowereegt2,    ncalotowereegt5,    ncalotowereegt10;

    int ncalotowerhfgt1,    ncalotowerhfgt2,    ncalotowerhfgt5,    ncalotowerhfgt10;



    int ncalotowerebgt1had,    ncalotowerebgt2had,    ncalotowerebgt5had,    ncalotowerebgt10had;
    int ncalotowereegt1had,    ncalotowereegt2had,    ncalotowereegt5had,    ncalotowereegt10had;
    int ncalotowerhfgt1had,    ncalotowerhfgt2had,    ncalotowerhfgt5had,    ncalotowerhfgt10had;

    int ncalotowerebgt1hof,    ncalotowerebgt2hof,    ncalotowerebgt5hof,    ncalotowerebgt10hof;
    int ncalotowereegt1hof,    ncalotowereegt2hof,    ncalotowereegt5hof,    ncalotowereegt10hof;
    int ncalotowerhfgt1hof,    ncalotowerhfgt2hof,    ncalotowerhfgt5hof,    ncalotowerhfgt10hof;

    int ncalotowerebgt1all,    ncalotowerebgt2all,    ncalotowerebgt5all,    ncalotowerebgt10all;
    int ncalotowereegt1all,    ncalotowereegt2all,    ncalotowereegt5all,    ncalotowereegt10all;
    int ncalotowerhfgt1all,    ncalotowerhfgt2all,    ncalotowerhfgt5all,    ncalotowerhfgt10all;


    float ctsumebgt1,    ctsumebgt2,    ctsumebgt5,    ctsumebgt10;
    float ctsumeegt1,    ctsumeegt2,    ctsumeegt5,    ctsumeegt10;
    float ctsumhfgt1,    ctsumhfgt2,    ctsumhfgt5,    ctsumhfgt10;

    float ctsumebgt05, ctsumeegt05, ctsumhfgt05;
    float ctsumebgt05had, ctsumeegt05had, ctsumhfgt05had;
    float ctsumebgt05hof, ctsumeegt05hof, ctsumhfgt05hof;
    float ctsumebgt05all, ctsumeegt05all, ctsumhfgt05all;



    float ctsumebgt1had,    ctsumebgt2had,    ctsumebgt5had,    ctsumebgt10had;
    float ctsumeegt1had,    ctsumeegt2had,    ctsumeegt5had,    ctsumeegt10had;
    float ctsumhfgt1had,    ctsumhfgt2had,    ctsumhfgt5had,    ctsumhfgt10had;

    float ctsumebgt1hof,    ctsumebgt2hof,    ctsumebgt5hof,    ctsumebgt10hof;
    float ctsumeegt1hof,    ctsumeegt2hof,    ctsumeegt5hof,    ctsumeegt10hof;
    float ctsumhfgt1hof,    ctsumhfgt2hof,    ctsumhfgt5hof,    ctsumhfgt10hof;

    float ctsumebgt1all,    ctsumebgt2all,    ctsumebgt5all,    ctsumebgt10all;
    float ctsumeegt1all,    ctsumeegt2all,    ctsumeegt5all,    ctsumeegt10all;
    float ctsumhfgt1all,    ctsumhfgt2all,    ctsumhfgt5all,    ctsumhfgt10all;


    float rechitsumet_eb_all, rechitsumet_eb_01, rechitsumet_eb_05;
    float rechitsumet_ee_all, rechitsumet_ee_01, rechitsumet_ee_05;

    float rechitsumet_eb_0105, rechitsumet_eb_0530;
    float rechitsumet_ee_0105, rechitsumet_ee_0530;


    int bunchintrain;



    float ebscsumet_all, eescsumet_all;
    float ebscsumet_all_eta15, ebscsumet_all_eta20, ebscsumet_all_eta25;
    float eescsumet_all_eta15, eescsumet_all_eta20, eescsumet_all_eta25;

    int ebnumsc_all, eenumsc_all;

    int ebnumrechits_01, ebnumrechits_0105, ebnumrechits_05, ebnumrechits_0530;
    int eenumrechits_01, eenumrechits_0105, eenumrechits_05, eenumrechits_0530;


    int ncalotowerebgt05, ncalotowereegt05, ncalotowerhfgt05;
    int ncalotowerebgt05had, ncalotowereegt05had, ncalotowerhfgt05had;
    int ncalotowerebgt05hof, ncalotowereegt05hof, ncalotowerhfgt05hof;
    int ncalotowerebgt05all, ncalotowereegt05all, ncalotowerhfgt05all;



};
#endif
