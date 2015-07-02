// -*- C++ -*-
//
// Package:    CodeTesting
// Class:      CodeTesting
// 
/**\class CodeTesting CodeTesting.cc TPGAnalyzer/CodeTesting/plugins/CodeTesting.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Iurii Antropov
//         Created:  Wed, 01 Jul 2015 14:48:09 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


// TowerID -> RCT Id transition
#include <map>
#include <iostream>
#include <math.h>
//#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"

using namespace std;

//
// class declaration
//

class CodeTesting : public edm::EDAnalyzer {
   public:
      explicit CodeTesting(const edm::ParameterSet&);
      ~CodeTesting();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      map<EBDetId, L1CaloRegionDetId> map_DetId_RegionID;

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
CodeTesting::CodeTesting(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed


	// Create map EBDetID -> RCT_regionId
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
}


CodeTesting::~CodeTesting()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CodeTesting::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   EBDetId tow_id = EBDetId((int)-28,(int)71);
   cout << tow_id << endl
   		 << "reg ieta : " << map_DetId_RegionID[tow_id].ieta() << " ; "
		 << "reg iphi : " << map_DetId_RegionID[tow_id].iphi() << endl;

   tow_id = EBDetId((int)25,(int)3);
      cout << tow_id << endl
      		 << "reg ieta : " << map_DetId_RegionID[tow_id].ieta() << " ; "
   		 << "reg iphi : " << map_DetId_RegionID[tow_id].iphi() << endl;

      tow_id = EBDetId((int)4,(int)2);
            cout << tow_id << endl
            		 << "reg ieta : " << map_DetId_RegionID[tow_id].ieta() << " ; "
         		 << "reg iphi : " << map_DetId_RegionID[tow_id].iphi() << endl;

}


// ------------ method called once each job just before starting event loop  ------------
void 
CodeTesting::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CodeTesting::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
CodeTesting::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
CodeTesting::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
CodeTesting::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
CodeTesting::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CodeTesting::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CodeTesting);
