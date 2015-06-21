// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//
// class declaration
//

class GoodVertexFilter : public edm::EDFilter {
   public:
      explicit GoodVertexFilter(const edm::ParameterSet&);
      ~GoodVertexFilter();

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      edm::InputTag vertexSrc;        
      unsigned int minNumTracks;
      double maxAbsZ;
      double maxd0;
      bool newSelection;
      unsigned int minNdof;
      // ----------member data ---------------------------
};

GoodVertexFilter::GoodVertexFilter(const edm::ParameterSet& iConfig)
{
  vertexSrc = iConfig.getParameter<edm::InputTag>("vertexCollection");
  minNumTracks = iConfig.getParameter<unsigned int>("minimumNumberOfTracks");
  maxAbsZ = iConfig.getParameter<double>("maxAbsZ");
  maxd0 = iConfig.getParameter<double>("maxd0");
  newSelection = iConfig.getParameter<bool>("newSel");
  minNdof = iConfig.getParameter<unsigned int>("minndof");
}


GoodVertexFilter::~GoodVertexFilter()
{
}

bool
GoodVertexFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 bool result = false; 
 edm::Handle<reco::VertexCollection> pvHandle; 
 iEvent.getByLabel(vertexSrc,pvHandle);
 const reco::VertexCollection & vertices = *pvHandle.product();
 for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
  {
      if (newSelection) {
        if (it->ndof() > minNdof)
         result = true ;
      }
      if (!newSelection) {
        if (it->tracksSize() > minNumTracks && 
            ( (maxAbsZ <=0 ) || fabs(it->z()) <= maxAbsZ ) &&
            ( (maxd0 <=0 ) || fabs(it->position().rho()) <= maxd0 )
          ) result = true;
      }
  }

   return result;
}


//define this as a plug-in
DEFINE_FWK_MODULE(GoodVertexFilter);
