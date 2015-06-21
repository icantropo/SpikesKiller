/**\class TdrHzzIsolationProducerProducer
 *
 *
 * Original Author:  Nicola De Filippis
 *
 * based on the HadIsolation module of R. Salerno 
 * 
 * Compute isolation for cones around electron candidates
 */

// #include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/TdrHzzIsolationProducer.h"
#include "EGamma/ECGelec/plugins/TdrHzzIsolationProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"


// Electrons
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>

// Calo
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoCaloTools/Selectors/interface/CaloDualConeSelector.h"

//#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/CandidateHadIsolation.h"
//#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/CandidateTkIsolation.h"
// #include "HiggsAnalysis/HiggsToZZ4Leptons/interface/CandidateHadIsolation.h"
// #include "HiggsAnalysis/HiggsToZZ4Leptons/interface/CandidateTkIsolation.h"
#include "EGamma/ECGelec/interface/CandidateHadIsolation.h"
#include "EGamma/ECGelec/interface/CandidateTkIsolation.h"



#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

// constructor
TdrHzzIsolationProducer::TdrHzzIsolationProducer(const edm::ParameterSet& pset) {

  // for tracker isolation
  tracksTag_        = pset.getParameter<edm::InputTag>("TracksLabel");
  muonsTag_         = pset.getParameter<edm::InputTag>("MuonsLabel");
  isoCone           = pset.getParameter<double>("isolationCone");
  isoVeto           = pset.getParameter<double>("isolationConeVeto");
  //isoVarTag       = pset.getParameter<edm::InputTag>("isoVarTag");
  isoVarCut         = pset.getParameter<double>("isoVarCut");
  
  // for hadronic isolation
  electronTag_      = pset.getParameter<edm::InputTag>("ElectronsLabel");
  electronVetoTag_  = pset.getParameter<edm::InputTag>("ElectronsVetoLabel");
  hcalrhitsLabel_   = pset.getParameter<string>("hcalrhits");
  radiusConeIntHad_ = pset.getParameter<double>("radiusConeIntHad");
  radiusConeExtHad_ = pset.getParameter<double>("radiusConeExtHad");
  eTMinHad_         = pset.getParameter<double>("eTMinHad");

  string alias;
  string iName = "TdrHzzIsolation";

  produces<vector<float> >( alias = iName + "SumPTtrackerele" ).setBranchAlias( alias );
  produces<vector<float> >( alias = iName + "SumEThad" ).setBranchAlias( alias );
  produces<vector<float> >( alias = iName + "SumEThadoverpT" ).setBranchAlias( alias );
  produces<vector<float> >( alias = iName + "SumEThad2" ).setBranchAlias( alias );
  produces<vector<float> >( alias = iName + "Xele" ).setBranchAlias( alias );
  produces<reco::GsfElectronCollection>();
  produces<edm::ValueMap<float> >();
  produces<edm::ValueMap<float> >(alias = "Tk" ).setBranchAlias( alias );
  produces<edm::ValueMap<float> >(alias = "Ecal" ).setBranchAlias( alias );
  produces<edm::ValueMap<float> >(alias = "Hcal" ).setBranchAlias( alias );

}


// destructor
TdrHzzIsolationProducer::~TdrHzzIsolationProducer() {

}


void TdrHzzIsolationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector<float> > electronIsoSumPTtrackerele( new vector<float> );	
  auto_ptr<vector<float> > electronIsoSumEThad( new vector<float> );
  auto_ptr<vector<float> > electronIsoSumEThad_pT( new vector<float> );
  auto_ptr<vector<float> > electronIsoSumEThad2( new vector<float> );
  auto_ptr<vector<float> > electronIsoXele( new vector<float> );
  auto_ptr<reco::GsfElectronCollection> Isolelec( new reco::GsfElectronCollection );
  auto_ptr<edm::ValueMap<float> > IsoEleMap(new edm::ValueMap<float> ());

  auto_ptr<edm::ValueMap<float> > IsoEleTkMap(new edm::ValueMap<float> ());
  auto_ptr<edm::ValueMap<float> > IsoEleEcalMap(new edm::ValueMap<float> ());
  auto_ptr<edm::ValueMap<float> > IsoEleHcalMap(new edm::ValueMap<float> ());
  
  edm::ValueMap<float>::Filler filler(*IsoEleMap);
  edm::ValueMap<float>::Filler fillerTk(*IsoEleTkMap);
  edm::ValueMap<float>::Filler fillerEcal(*IsoEleEcalMap);
  edm::ValueMap<float>::Filler fillerHcal(*IsoEleHcalMap);

  iSetup.get<CaloGeometryRecord>().get(theCaloGeom_);


   // get muons
  Handle<edm::View<reco::Muon> > AllMuons0;
  iEvent.getByLabel(muonsTag_.label(), AllMuons0);
  const edm::View<reco::Muon>* muonCollectionmu =  AllMuons0.product () ;


  //Tracks
  Handle<edm::View<Track> > tracks;
  iEvent.getByLabel(tracksTag_.label(), tracks);
  const edm::View<reco::Track>* trackCollection = tracks.product () ;
  
  //  hcal rechit collection
  edm::Handle<HBHERecHitCollection> hbhe;
  iEvent.getByLabel(hcalrhitsLabel_, hbhe); 
  HBHERecHitMetaCollection mhbhe =  HBHERecHitMetaCollection(*hbhe);    
  
  // Get electron candidates  
  Handle<edm::View<GsfElectron> > EleCandidates;
  iEvent.getByLabel(electronTag_.label(), EleCandidates);
  
  Handle<edm::View<GsfElectron> > EleVetoCandidates;
  iEvent.getByLabel(electronVetoTag_.label(), EleVetoCandidates);
  const edm::View<reco::GsfElectron>* Resolved_Collection = EleVetoCandidates.product () ;

  // First isolation with tracker then with hadronic depositions
  size_t n = EleCandidates->size();
  std::vector<float> iso(n);
  std::vector<float> isoTk(n);
  std::vector<float> isoEcal(n);
  std::vector<float> isoHcal(n);

  unsigned int index=0;
  double PtTracksCorr, EtHadClusters;
 
  for (edm::View<reco::GsfElectron>::const_iterator cand = EleCandidates->begin(); 
       cand != EleCandidates->end(); ++cand) {

    CandidateTkIsolation myTkIsolation(&(*cand),trackCollection,Resolved_Collection,muonCollectionmu,1) ;
    myTkIsolation.setIntRadius(isoVeto);
    myTkIsolation.setExtRadius(isoCone);
    PtTracksCorr = myTkIsolation.getPtTracksCorr()[0];
    electronIsoSumPTtrackerele->push_back(PtTracksCorr / cand->pt());
    float myHadET=0.;

    CandidateHadIsolation myHadIsolation(theCaloGeom_,&mhbhe,&(*cand)) ;
    myHadIsolation.setEtLow (eTMinHad_) ;
    myHadIsolation.setExtRadius(radiusConeExtHad_ );
    myHadIsolation.setIntRadius(radiusConeIntHad_ );
    EtHadClusters = myHadIsolation.getEtHadClusters();
    
    electronIsoSumEThad->push_back(EtHadClusters);
    electronIsoSumEThad_pT->push_back(EtHadClusters / cand->pt());
    electronIsoSumEThad2->push_back(EtHadClusters * EtHadClusters);
    myHadET=EtHadClusters/cand->pt();
    electronIsoXele->push_back(PtTracksCorr/cand->pt()+myHadET);

    iso[index] = float(PtTracksCorr / cand->pt()+myHadET);
    isoTk[index] = float(PtTracksCorr);
    isoEcal[index] = float(-999.);
    isoHcal[index] = float(myHadET);

    Isolelec->push_back( *cand );

    index++;    
    
  }

  // filling map
  filler.insert(EleCandidates, iso.begin(), iso.end());
  filler.fill();

  fillerTk.insert(EleCandidates, isoTk.begin(), isoTk.end());
  fillerTk.fill();

  fillerEcal.insert(EleCandidates, isoEcal.begin(), isoEcal.end());
  fillerEcal.fill();

  fillerHcal.insert(EleCandidates, isoHcal.begin(), isoHcal.end());
  fillerHcal.fill();

  //sorting
  //sort(electronIsoSumPTtrackerele->begin(),electronIsoSumPTtrackerele->end());
  //sort(electronIsoSumEThad->begin(),electronIsoSumEThad->end());
  //sort(electronIsoSumEThad_pT->begin(),electronIsoSumEThad_pT->end());
  //sort(electronIsoSumEThad2->begin(),electronIsoSumEThad2->end());
  //sort(electronIsoXele->begin(),electronIsoXele->end());
    
  const string & isoName = "TdrHzzIsolation";
  iEvent.put( electronIsoSumPTtrackerele, isoName + "SumPTtrackerele" );
  iEvent.put( electronIsoSumEThad, isoName + "SumEThad" );
  iEvent.put( electronIsoSumEThad_pT, isoName + "SumEThadoverpT" );
  iEvent.put( electronIsoSumEThad2, isoName + "SumEThad2" );
  iEvent.put( electronIsoXele, isoName + "Xele" );
  const string iName = "";
  iEvent.put( Isolelec, iName );
  iEvent.put( IsoEleMap, iName );
  iEvent.put( IsoEleTkMap, "Tk" );
  iEvent.put( IsoEleEcalMap, "Ecal" );
  iEvent.put( IsoEleHcalMap, "Hcal" );

}


