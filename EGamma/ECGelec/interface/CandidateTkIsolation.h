#ifndef CandidateTkIsolation_h
#define CandidateTkIsolation_h

//C++ includes
#include <vector>
//#include <functional>
#include <stdio.h>
#include <string>

//CMSSW includes 
// #include "FWCore/Framework/interface/ESHandle.h" //?????
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

class CandidateTkIsolation {
 public:
  
  //constructors
  CandidateTkIsolation () ;
  CandidateTkIsolation (const reco::Candidate * c, 
			const reco::CandidateCollection * c1coll ) ;
  CandidateTkIsolation (const reco::Candidate * c, int trackertrack ) ;
  CandidateTkIsolation (const reco::GsfElectron * electron, 
			const edm::View<reco::Track>* trackCollection,
			const edm::View<reco::GsfElectron>* electronCollection,
			int trackertrack); 

  CandidateTkIsolation (const reco::GsfElectron * electron, 
			const edm::View<reco::Track>* trackCollection,
			const edm::View<reco::GsfElectron>* electronCollection,
			const edm::View<reco::Muon>* muonCollection,
			reco::TrackBase::Point beamPoint,
			int trackertrack); 

  CandidateTkIsolation (const reco::GsfElectron * electron, 
			const edm::View<reco::Track>* trackCollection,
			const edm::View<reco::GsfElectron>* electronCollection,
			const edm::View<reco::Muon>* muonCollection,
			int trackertrack); 

  CandidateTkIsolation (const math::PtEtaPhiMLorentzVector* direction,	
			const edm::View<reco::Track>* trackCollection,
			const edm::View<reco::GsfElectron>* electronCollection,
			const edm::View<reco::Muon> *muonCollection,
			int trackertrack);
			

   //methods
  void setExtRadius (double extRadius) ;
  void setIntRadius (double intRadius) ;
  void setPtLow (double ptLow) ;
  void setLip (double lip) ;
  void setZmassinf ( double Zmassinf);
  void setZmasssup ( double Zmasssup);
  void setJurrasic (double strip); 
  void setDebugMode(bool debug);

  int getNumberTracks() const ;
  double getPtTracks () const ;
  double getPtTracks2() const ;
  std::vector<double> getPtTracksCorr() const ;

  //destructor 
  ~CandidateTkIsolation() ;
  
 private:

  math::XYZVector thisElectronMomentumAtVtx; 
  math::XYZTLorentzVector thisElectronLorentzVector;

  const math::PtEtaPhiMLorentzVector* eleP4;
  const reco::GsfElectron*  electron_ ;
  const edm::View<reco::Track> *trackCollection_ ;
  const edm::View<reco::GsfElectron> *electronCollection_ ;
  const edm::View<reco::Muon> *muonCollection_ ;

  double extRadius_ ;
  double intRadius_ ;
  double stripWidth_;
  double ptLow_ ;
  double lip_ ;
  float ZMASSPDG_;
  double Zmass_inf_;
  double Zmass_sup_;
  reco::TrackBase::Point beamPoint_;
  bool debug_;

  bool isRC;
  bool tracker_track_;
  bool isGoodTrack (edm::View<reco::Track>::const_iterator track) const;
  bool isGoodTrackerTrack(edm::View<reco::Track>::const_iterator track) const;
  bool isGoodTrackerTrack(reco::TrackRef track) const;
  void Initialize();

};

#endif
