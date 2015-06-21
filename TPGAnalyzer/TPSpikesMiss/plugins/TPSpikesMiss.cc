// -*- C++ -*-
//
// Package:    TPSpikesMiss
// Class:      TPSpikesMiss
// 
/**\class TPSpikesMiss TPSpikesMiss.cc TPGAnalyzer/TPSpikesMiss/plugins/TPSpikesMiss.cc

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

//include root
#include "TTree.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TFile.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

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

// CMSSW
/*#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
*/

//
// class declaration
//

using namespace std;
//using namespace reco;
using namespace edm;
//using namespace IPTools;

class TPSpikesMiss : public edm::EDAnalyzer {
   public:
      explicit TPSpikesMiss(const edm::ParameterSet&);
      ~TPSpikesMiss();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // handles to get the TPs
      edm::Handle<EcalTrigPrimDigiCollection> * ecal_tp_;
      edm::Handle<EcalTrigPrimDigiCollection> * offlineTP;
      EcalTrigTowerDetId TPtowid_;
      EcalTrigTowerDetId TPtowidM_;

      // emulated TP
      int _trig_tower_N_emul, _trig_tower_ieta_emul[4032],_trig_tower_iphi_emul[4032],
               _trig_tower_adc_emul[4032][5], _trig_tower_sFGVB_emul[4032][5]; 

      std::string   histogramFile_;  
      edm::InputTag tpEmulatorCollection_ ;
      edm::InputTag tpOnlineCollection_ ;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      TFile* file_;
      TH1F *h_ttonline_et;
      TH1F *h_ttoffline_et;
      TH1F *h_ttoffline_spikes_matches_et;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TPSpikesMiss::TPSpikesMiss(const edm::ParameterSet& iConfig) :
  histogramFile_ (iConfig.getParameter<std::string>("histogramFile")),
  tpEmulatorCollection_ (iConfig.getParameter<edm::InputTag> ("TPEmulatorCollection")),
  tpOnlineCollection_ (iConfig.getParameter<edm::InputTag> ("TPOnlineCollection"))
{
   //now do what ever initialization is needed
	
}


TPSpikesMiss::~TPSpikesMiss()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete file_;
  delete h_ttonline_et;
  delete h_ttoffline_et;
  delete h_ttoffline_spikes_matches_et;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TPSpikesMiss::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   //ONLINE TPs (Alex's code for Iurii the best)
   edm::Handle<EcalTrigPrimDigiCollection> onlineTP;
   iEvent.getByLabel(tpOnlineCollection_, onlineTP);
   EcalTrigPrimDigiCollection const &onlineTPs  = *(onlineTP.product());
   EcalTrigPrimDigiCollection::const_iterator onTP_it;

   edm::Handle<EcalTrigPrimDigiCollection>  offlineTP;
   iEvent.getByLabel(tpEmulatorCollection_, offlineTP);
   EcalTrigPrimDigiCollection const &offlineTPs = *(offlineTP.product());
   EcalTrigPrimDigiCollection::const_iterator ofTP_it;    
  
   //std::cout << "looping on online TPG" << std::endl;
   for(onTP_it = onlineTPs.begin(); onTP_it != onlineTPs.end(); ++onTP_it){
     if(onTP_it->compressedEt() > 0) {
       //std::cout << "ET=" << onTP_it->compressedEt() << " raw=" << ((*onTP_it)[0].raw()&0xff) << std::endl;
       h_ttonline_et->Fill(onTP_it->compressedEt());
     }//non zero TP
      if(onTP_it->compressedEt() == 255){
       //std::cout << "ET=255  iphi="<< onTP_it->id().iphi() << "  ieta=" << onTP_it->id().ieta()<< std::endl; 
       EcalTrigPrimDigiCollection::const_iterator  offlineMatchForSpike_it;
       offlineMatchForSpike_it =  offlineTPs.find(onTP_it->id());
       if (offlineMatchForSpike_it!=offlineTPs.end()){
          h_ttoffline_spikes_matches_et->Fill((*offlineMatchForSpike_it)[2].raw()&0xff);
	  //std::cout << "Offline match"<< std::endl;
          //std::cout << "ET="<<offlineMatchForSpike_it->compressedEt()<<"  iphi="<< offlineMatchForSpike_it->id().iphi() << "  ieta=" << offlineMatchForSpike_it->id().ieta()<< std::endl;
          } else {
                 std::cout << "No match tower for next spike was found offline:"<< std::endl;
                 std::cout << "ET=255  iphi="<< onTP_it->id().iphi() << "  ieta=" << onTP_it->id().ieta()<< std::endl;
	         } //end else 
       
      }//  end if spike

     /*edm::ESHandle<EcalTPGTowerStatus> theEcalTPGTowerStatus_handle;
     iSetup.get<EcalTPGTowerStatusRcd>().get(theEcalTPGTowerStatus_handle);
     const EcalTPGTowerStatus * ecaltpgTowerStatus=theEcalTPGTowerStatus_handle.product();
     const EcalTPGTowerStatusMap &towerMap=ecaltpgTowerStatus->getMap();
     */
     //EcalTPGTowerStatusMapIterator towIt; // = towerMap.find(onIt_d->id().rawId());
     //std::cout << "ET=255  iphi="<< towIt->id().iphi() << "  ieta=" << towIt->id().ieta()<< std::endl; 
     
   }//loop online TPG

   
   // EMULATOR TPs (Nadir's code  modified by me :)
   
   //if (print_) std::cout<<"TPEmulator collection size="<<tpEmul.product()->size()<<std::endl ;
 
   for (ofTP_it = offlineTPs.begin(); ofTP_it != offlineTPs.end() ; ++ofTP_it) {
     //TPtowidM_ = ofTP_it->id();
     //_trig_tower_iphi_emul[i] = TPtowidM_.iphi() ;
     //_trig_tower_ieta_emul[i] = TPtowidM_.ieta() ;
     
     h_ttoffline_et->Fill((*ofTP_it)[2].raw()&0xff);
     
     /*for (int j=0 ; j<5 ; j++) {
       _trig_tower_adc_emul[i][j] = (dM_[j].raw()&0xff) ;
       //_trig_tower_sFGVB_emul[i][j] = d[j].l1aSpike(); 
       //_trig_tower_sFGVB_emul[i][j] = dM_[j].sFGVB(); 
       if(showit)
	 cout << (dM_[j].raw()&0xff) << " " ;
     }*/

     /*onTP_it = onlineTPs.find(ofTP_it->id());
     if(onTP_it == onlineTPs.end())
       std::cout << "ONLINE TP NOT FOUND" << std::endl;
     else {
       cout << "compare values=" << onTP_it->compressedEt() << " " << ofTP_it->compressedEt() << endl;
     }// matching offline TP and online TP 
     */
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
TPSpikesMiss::beginJob()
{

  file_          = new TFile(histogramFile_.c_str(),"RECREATE");
  h_ttonline_et  = new TH1F( "h_ttonline_et", "h_ttonline_et", 256,0,256);
  h_ttoffline_et = new TH1F( "h_ttoffline_et", "h_ttoffline_et", 256,0,256);
  h_ttoffline_spikes_matches_et = new TH1F( "h_ttoffline_spikes_matches_et", "h_ttoffline_et", 256,0,256);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TPSpikesMiss::endJob() 
{

  if (file_ !=0) 
    {
      file_->cd();

      h_ttonline_et->Write();
      h_ttoffline_et->Write();
      h_ttoffline_spikes_matches_et->Write();
    }
  file_ = 0;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
TPSpikesMiss::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TPSpikesMiss::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TPSpikesMiss::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TPSpikesMiss::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TPSpikesMiss::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TPSpikesMiss);
