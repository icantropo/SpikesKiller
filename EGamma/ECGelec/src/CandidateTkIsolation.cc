
//C++ includes
#include <vector>
#include <functional>

//ROOT includes
#include <Math/VectorUtil.h>

//CMSSW includes
// #include "HiggsAnalysis/HiggsToZZ4Leptons/interface/CandidateTkIsolation.h"
#include "EGamma/ECGelec/interface/CandidateTkIsolation.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleanerBySharedHits.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
//#include "TrackingTools/PatternTools/interface/TrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

using namespace std ;
using namespace reco;
using namespace ROOT::Math::VectorUtil ;

CandidateTkIsolation::CandidateTkIsolation () {}

CandidateTkIsolation::CandidateTkIsolation (const reco::Candidate * electron,
					    const reco::CandidateCollection * trackCollection)  
{
  Initialize();
 }

CandidateTkIsolation::CandidateTkIsolation (const reco::Candidate *electron,  int trackertrack)
{

  //const Candidate * candele;
  //candele = &* electron->masterClone();
  //electron_ = dynamic_cast<const GsfElectron *>(candele);
  //electron_ = dynamic_cast<const GsfElectron& >(electron);
  Initialize();
  if (trackertrack==1) tracker_track_=true;
  else if (trackertrack==0) tracker_track_= false;
  else {
    cout<< "TkIsolation::not a good value for tracker track, chosed trackertrack= false"<<endl;
    tracker_track_= false;
  }
}  


CandidateTkIsolation::CandidateTkIsolation (const reco::GsfElectron* electron, 
				   const edm::View<reco::Track>* trackCollection,
				   const edm::View<reco::GsfElectron>* electronCollection,
				   int trackertrack) : 
  electron_(electron) ,
  trackCollection_(trackCollection) ,
  electronCollection_(electronCollection) 
{
  thisElectronMomentumAtVtx = electron_->gsfTrack()->innerMomentum () ; //(*tmpTrack)
  thisElectronLorentzVector = electron_->p4() ;
  Initialize();
  if (trackertrack==1) tracker_track_=true;
  else if (trackertrack==0) tracker_track_= false;
  else {
    cout<< "TkIsolation::not a good value for tracker track, chosed trackertrack= false"<<endl;
    tracker_track_= false;
  }
}  

//we normally use this constructor
CandidateTkIsolation::CandidateTkIsolation (const reco::GsfElectron* electron, 
				   const edm::View<reco::Track>* trackCollection,
				   const edm::View<reco::GsfElectron>* electronCollection,
				   const edm::View<reco::Muon> *muonCollection,
				   reco::TrackBase::Point beamPoint,
				   int trackertrack) : 
  electron_(electron) ,
//   thisElectronMomentumAtVtx(electron_->gsfTrack()->innerMomentum ()),
//   thisElectronLorentzVector(electron->p4()), 
  trackCollection_(trackCollection) ,
  electronCollection_(electronCollection) ,
  muonCollection_(muonCollection),
  beamPoint_(beamPoint)
{
  thisElectronMomentumAtVtx = electron_->gsfTrack()->innerMomentum () ; //(*tmpTrack)
  thisElectronLorentzVector = electron_->p4() ;
  Initialize();
  if (trackertrack==1) tracker_track_=true;
  else if (trackertrack==0) tracker_track_= false;
  else {
    cout<< "TkIsolation::not a good value for tracker track, chosed trackertrack= false"<<endl;
    tracker_track_= false;
  };
} 


CandidateTkIsolation::CandidateTkIsolation (const reco::GsfElectron* electron, 
				   const edm::View<reco::Track>* trackCollection,
				   const edm::View<reco::GsfElectron>* electronCollection,
				   const edm::View<reco::Muon> *muonCollection,
				   int trackertrack) : 
  electron_(electron) ,
  trackCollection_(trackCollection) ,
  electronCollection_(electronCollection) ,
  muonCollection_(muonCollection)  
{
  thisElectronMomentumAtVtx = electron_->gsfTrack()->innerMomentum () ; //(*tmpTrack)
  thisElectronLorentzVector = electron_->p4() ;

  beamPoint_.SetXYZ(0.,0.,0.);
  Initialize();
  if (trackertrack==1) tracker_track_=true;
  else if (trackertrack==0) tracker_track_= false;
  else {
    cout<< "TkIsolation::not a good value for tracker track, chosed trackertrack= false"<<endl;
    tracker_track_= false;
  };
} 
//Constructor used with random cone technique
CandidateTkIsolation::CandidateTkIsolation (const math::PtEtaPhiMLorentzVector* direction,	//const reco::GsfElectron* electron, 
				   const edm::View<reco::Track>* trackCollection,
				   const edm::View<reco::GsfElectron>* electronCollection,
				   const edm::View<reco::Muon> *muonCollection,
				   int trackertrack) : 
  eleP4(direction) ,
  trackCollection_(trackCollection) ,
  electronCollection_(electronCollection) ,
  muonCollection_(muonCollection)
{
  thisElectronLorentzVector.SetXYZT(eleP4->x(),eleP4->y(),eleP4->z(),eleP4->t()); 
  thisElectronMomentumAtVtx.SetXYZ(eleP4->px(),eleP4->py(),eleP4->pz());
  beamPoint_.SetXYZ(0.,0.,0.);
  Initialize();
  if (trackertrack==1) tracker_track_=true;
  else if (trackertrack==0) tracker_track_= false;
  else {
    cout<< "TkIsolation::not a good value for tracker track, choosed trackertrack= false"<<endl;
    tracker_track_= false;
  };
  isRC=true;
}

CandidateTkIsolation::~CandidateTkIsolation (){}
void CandidateTkIsolation::setExtRadius (double extRadius){ extRadius_ = extRadius ;}
void CandidateTkIsolation::setIntRadius (double intRadius){ intRadius_ = intRadius ;}
void CandidateTkIsolation::setJurrasic(double strip){  stripWidth_ = strip;}
void CandidateTkIsolation::setPtLow (double ptLow){ ptLow_ = ptLow ;}
void CandidateTkIsolation::setLip (double lip){  lip_ = lip ;}
void CandidateTkIsolation::setZmassinf ( double Zmassinf){ Zmass_inf_ = Zmassinf;}
void CandidateTkIsolation::setZmasssup ( double Zmasssup){ Zmass_sup_ = Zmasssup;}
void CandidateTkIsolation::setDebugMode (bool debug) {debug_= debug;}


int CandidateTkIsolation::getNumberTracks () const
{  
  //counter for the tracks in the isolation cone
  int dummyCounter = 0;    

//   if(electron_!=0) {
	for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
		itrTr != (*trackCollection_).end(); ++itrTr){
			math::XYZVector tmpTrackMomentumAtVtx = itrTr->innerMomentum () ; 
			if (tracker_track_ && !isGoodTrackerTrack(itrTr) ) continue;
			if (!tracker_track_ && !isGoodTrack(itrTr)) continue;
			double dr = DeltaR(tmpTrackMomentumAtVtx,thisElectronMomentumAtVtx) ;
			double deta = (*itrTr).eta() - thisElectronMomentumAtVtx.eta();
			//double this_pt  = sqrt( itrTr->innerMomentum().Perp2() );
			if ( fabs(dr) < extRadius_ &&  fabs(dr) > intRadius_ && fabs(deta) >= stripWidth_) 
				++dummyCounter ;
  	}
//   }
  return dummyCounter ;
}

double CandidateTkIsolation::getPtTracks () const
{
  //dummy counter for the pT of tracks inside the cone
  double dummypT = 0 ;
	for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
		itrTr != (*trackCollection_).end(); ++itrTr){
			math::XYZVector tmpTrackMomentumAtVtx = itrTr->innerMomentum () ; 
			if (tracker_track_ && !isGoodTrackerTrack(itrTr) ) continue;
			if (!tracker_track_ && !isGoodTrack(itrTr)) continue;
			double this_pt  = sqrt( itrTr->innerMomentum().Perp2() );
			double dr = DeltaR(tmpTrackMomentumAtVtx,thisElectronMomentumAtVtx) ;
			double deta = (*itrTr).eta() - thisElectronMomentumAtVtx.eta();
			if ( fabs(dr) < extRadius_ &&  fabs(dr) > intRadius_ && fabs(deta) >= stripWidth_){
				dummypT += this_pt;
			}
  	}
  if (debug_) cout  << "[CandidateTkIsolation]::getPtTracks: sumPt = " << dummypT << endl;
  return dummypT ;
}


double CandidateTkIsolation::getPtTracks2 () const
{
  //dummy counter for the pT of tracks inside the cone
  double dummypT = 0 ;
    reco::GsfTrackRef tmpTrack = electron_->gsfTrack() ;
//     math::XYZVector thisElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 
	// thisElectronMomentumAtVtx
	for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
		itrTr != (*trackCollection_).end(); ++itrTr){
			math::XYZVector tmpTrackMomentumAtVtx = itrTr->innerMomentum () ; 
			if (tracker_track_ && !isGoodTrackerTrack(itrTr) ) continue;
			if (!tracker_track_ && !isGoodTrack(itrTr)) continue;			
			double this_pt  = sqrt( itrTr->innerMomentum().Perp2() );			
			double dr = DeltaR(tmpTrackMomentumAtVtx,thisElectronMomentumAtVtx) ;
			double deta = (*itrTr).eta() - thisElectronMomentumAtVtx.eta();
			if ( fabs(dr) < extRadius_ &&  fabs(dr) > intRadius_ && fabs(deta) >= stripWidth_){
			      dummypT += this_pt*this_pt;
			}
  	}
  return dummypT;
}


std::vector<double> CandidateTkIsolation::getPtTracksCorr () const
{
  //dummy counter for the pT of tracks inside the cone
  vector<double> dummypTout;
  double dummypT=0,dummypT1=0,dummypT2=0;
  double tmpElPt1 = 0.,tmpElPt2 =0, drTrackMatch = 0.;

  if (debug_) cout  << "[CandidateTkIsolation]::getPtTracksCorr" << endl;

  //track
  std::vector<math::XYZVector> tracks;
  
//   dummypT=getPtTracks();
	for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
		itrTr != (*trackCollection_).end(); ++itrTr){
			math::XYZVector tmpTrackMomentumAtVtx = itrTr->innerMomentum () ; 
			
			if (tracker_track_ && !isGoodTrackerTrack(itrTr) ) continue;
			if (!tracker_track_ && !isGoodTrack(itrTr)) continue;
			double this_pt  = sqrt( itrTr->innerMomentum().Perp2() );
			double dr = DeltaR(tmpTrackMomentumAtVtx,thisElectronMomentumAtVtx) ;
			double deta = (*itrTr).eta() - thisElectronMomentumAtVtx.eta();
			if ( fabs(dr) < extRadius_ &&  fabs(dr) > intRadius_ && fabs(deta) >= stripWidth_){
				dummypT += this_pt;
				tracks.push_back(tmpTrackMomentumAtVtx);
			}
  	}


  if (debug_) cout  << "[CandidateTkIsolation]::getPtTracksCorr dummypT=  " << dummypT <<  endl;

  if (dummypT==0) {
		if (debug_) cout  << "Electron isolated: no tracks around it.." << endl;
		dummypTout.push_back(dummypT);
  		dummypTout.push_back(dummypT);
  		return dummypTout;	//no need to do any additional subtraction of ele || mu pt
  }
  
//   if (debug_) cout  << "[CandidateTkIsolation]::getPtTracksCorr::LoopOnMuons" << endl;

  for ( edm::View<reco::Muon>::const_iterator itrTr  = muonCollection_->begin() ; 
 	itrTr != muonCollection_->end()   ; 
 	++itrTr) {

// 	check if the muon is tracker muon and global muon 

	if (!(itrTr->isTrackerMuon() || itrTr->isGlobalMuon()) ) continue;
	reco::TrackRef track = itrTr->innerTrack();  
	if (tracker_track_ && !isGoodTrackerTrack(track) ) continue;
	math::XYZVector tmpTrackMomentumAtVtxMu(itrTr->innerTrack()->px(),itrTr->innerTrack()->py(),itrTr->innerTrack()->pz());
	double dr= DeltaR(tmpTrackMomentumAtVtxMu,thisElectronMomentumAtVtx) ;
	double deta = tmpTrackMomentumAtVtxMu.eta() - thisElectronMomentumAtVtx.eta();

      
	if ( fabs(dr) < extRadius_ && fabs(dr) > intRadius_ && fabs(deta) >= stripWidth_){
	    int closestTrackPt2MuPtIndex =-1;
	    int MuTrackIndex=9999;
	    double deltaPt=999.;
	    double thismuon_pt  = sqrt( tmpTrackMomentumAtVtxMu.Perp2 () );
    	    if (debug_) cout  << "Found a muon in the isolation cone of electrons..." << endl;
	    
	    //check if muon track is considered when summing tracks
	    double epsilon1 = 0.05, epsilon2 = 0.1, epsilon3=0.4; //???????
	    double ERRpt, ERReta, ERRphi, ERR, minERR=999.;
	    for (unsigned int i=0; i<tracks.size(); i++){
		if ( fabs(  sqrt( tracks[i].Perp2() ) - thismuon_pt ) < deltaPt ) closestTrackPt2MuPtIndex=i;
		ERRpt  = fabs( (sqrt( tracks[i].Perp2() ) - thismuon_pt)/thismuon_pt );
		ERReta = fabs( (tracks[i].eta() - tmpTrackMomentumAtVtxMu.eta()));
		ERRphi = fabs( (tracks[i].phi() - tmpTrackMomentumAtVtxMu.phi()));
		ERR = ERRpt + ERReta + ERRphi;
		if (ERRpt < epsilon1 && ERReta < epsilon2 && ERRphi < epsilon3 && ERR < minERR) MuTrackIndex=i;
		minERR=ERR;
	    }
	    if ( MuTrackIndex == 9999 && debug_) cout << "No muon track in TrackCollection. " 
						      << "Closest track to muon in cone (pt/eta/phi): "	
						      <<  sqrt( tracks[closestTrackPt2MuPtIndex].Perp2() ) << " "
						      <<  tracks[closestTrackPt2MuPtIndex].eta() << " "
						      <<  tracks[closestTrackPt2MuPtIndex].phi() << endl;

	    else {
		if (debug_) cout  << "Subtracting muon track (pt/eta/phi): "	<<  sqrt( tracks[MuTrackIndex].Perp2() ) << " "
										<<  tracks[MuTrackIndex].eta() << " "
										<<  tracks[MuTrackIndex].phi() << endl;
		dummypT -= sqrt( tracks[MuTrackIndex].Perp2() ) ;
	    }

	    if (debug_) {
	      cout  << "Muon Class: ";
	      if (itrTr->isGlobalMuon()) cout << "G";
	      if (itrTr->isTrackerMuon()) cout << "T";
	      if (itrTr->isStandAloneMuon()) cout << "S";
	      cout << endl;
	      cout << "This muon track: (pt/eta/phi): " << itrTr->innerTrack()->pt()  << " " 
							<< itrTr->innerTrack()->eta() << " "
							<< itrTr->innerTrack()->phi() << endl;
	    }
	}
  }
  

  // correct for close pairs from Z(Z/gamma* decay)
  // idea: look for an electron candidate of oposite sign in the cone, substract its pT from the sum if 
  // invariant mass with original electron is greater then a threshold
  double invMassThr = 12. ;
  
  int iloop=0, iloop2=0;
//   if (debug_) cout  << "[CandidateTkIsolation]::getPtTracksCorr::LoopOnOtherElectrons" << endl;

  for ( edm::View<reco::GsfElectron>::const_iterator ele  = electronCollection_->begin(); 
	ele != electronCollection_->end(); 
	++ele ) {
    //if (debug_) cout  <<"electron number: "<<iloop<<endl;
    bool elveto=false;
    reco::GsfTrackRef tmpTrack = ele->gsfTrack() ;
    math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 

//     if (debug_) cout  << "[CandidateTkIsolation]::getPtTracksCorr: Check if it is the same electron" << endl;
    if (thisElectronMomentumAtVtx == tmpElectronMomentumAtVtx) {
      //   if (debug_) cout <<"same electron"<<endl;
      iloop++;
      continue ;
    }
    
    double dr = DeltaR(thisElectronMomentumAtVtx,tmpElectronMomentumAtVtx);
    double deta = tmpElectronMomentumAtVtx.eta() - thisElectronMomentumAtVtx.eta();

    //   if (debug_) cout <<"dr= "<<dr<<endl;
    if ( fabs(dr) < extRadius_  && fabs(dr) > intRadius_ && fabs(deta) >= stripWidth_ ) { 
		// if (debug_) cout <<"inside cone"<<endl;
		math::XYZTLorentzVector tmpElectronLorentzVector  = ele->p4() ;
		math::XYZTLorentzVector twoElectrons ;
		twoElectrons = thisElectronLorentzVector + tmpElectronLorentzVector ;
	
		//Because of charge misId, maybe we should not reqire opposite charge like here???
		if (!isRC){
		  if (ele->charge() != electron_->charge() && (twoElectrons.mass() > invMassThr)){
			  elveto= true;
		  //	if (debug_) cout <<"e vetoed Z*"<<endl;
		  }
		}
		
		float massmin=100;
		
		if (debug_) cout  << "[CandidateTkIsolation]::getPtTracksCorr::LoopOnOtherElectronsToVeto" << endl;
		for (edm::View<reco::GsfElectron>::const_iterator oele  = electronCollection_->begin(); 
		    oele != electronCollection_->end(); ++oele ) {
			if (debug_) cout <<"checking electrons outside cone : "<<iloop2<<endl;
			math::XYZTLorentzVector otherElectronLorentzVector = oele->p4() ;
			if ((oele->gsfTrack()->innerMomentum() != tmpElectronMomentumAtVtx) 
			&& (oele->gsfTrack()->innerMomentum() != thisElectronMomentumAtVtx)){
				if (oele->charge() != ele->charge()){	//electrons from opposite traveling Z/Z*
					twoElectrons = tmpElectronLorentzVector + otherElectronLorentzVector;
					float mass= twoElectrons.mass();
					// if (debug_) cout <<"mass= "<<mass<<endl;
					if (fabs(mass-ZMASSPDG_) <fabs( massmin))
						massmin= mass-ZMASSPDG_;
				}
				if (!isRC){
				  if (oele->charge() != electron_->charge()){ //electrons from same Z/Z* highly boosted
					  twoElectrons = thisElectronLorentzVector + otherElectronLorentzVector;
					  float mass= twoElectrons.mass();
					  //if (debug_) cout <<"mass= "<<mass<<endl;
					  if (fabs(mass-ZMASSPDG_) <fabs( massmin))
						  massmin= mass-ZMASSPDG_;
				  }
				}
			}
			iloop2++;
		}
		if (debug_) cout <<"massmin= "<<massmin<<endl;
		if (debug_) cout <<"zmassinf= "<<Zmass_inf_<<" Zmasssup= "<<Zmass_sup_<<endl;
		if (massmin > Zmass_inf_ && massmin < Zmass_sup_) elveto=true; //ROKO: ????????Which cuts to put
		
		if (elveto && debug_)cout<<"this electron is vetoed"<<endl;
		else if (debug_) cout<<"this electron is not vetoed"<<endl;


		if (debug_)cout << "[CandidateTkIsolation]::getPtTracksCorr::Search for electron track to remove" << endl;
  		// make a loop on found KF tracks and find if any of them matches this electron
		for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
		itrTr != (*trackCollection_).end()   ; 
		++itrTr) {

			tmpElPt1=0.,tmpElPt2=0.,drTrackMatch = 0.;
			
			//make sure the track that could be substracted are one that are already in the sumPt
			if (tracker_track_ && !isGoodTrackerTrack(itrTr) ) continue;
			if (!tracker_track_ && !isGoodTrack(itrTr)) continue;
			math::XYZVector thisTrackMomentumAtVtx = itrTr->innerMomentum() ;
			drTrackMatch = DeltaR(tmpElectronMomentumAtVtx,thisTrackMomentumAtVtx);
			
			if ( drTrackMatch < intRadius_ ) 	//match track to electron
			{
				dr = DeltaR(thisElectronMomentumAtVtx, thisTrackMomentumAtVtx);
				double deta = thisTrackMomentumAtVtx.eta() - thisElectronMomentumAtVtx.eta();
				//check again if the matched track is in alowed isolation area
				if (fabs(dr) < extRadius_  && fabs(dr) > intRadius_ && fabs(deta) >= stripWidth_){
					  
				      tmpElPt1=sqrt( thisTrackMomentumAtVtx.Perp2 () );
				      if (debug_) cout  << "Found track matched to electron:" << " all ele pt= "<<tmpElPt1<<endl;
				      if (elveto) {
					tmpElPt2=sqrt( thisTrackMomentumAtVtx.Perp2 () );
					if (debug_) cout <<"veto ele pt= "<<tmpElPt2<<endl;
				      }
				      break;
				}
			}
			else if (debug_) cout <<"no KF track match"<<endl;
		}//end loop KF tracks
    }//end if in internal cone
    iloop++;
  }//end loop over electrons
  
  //  }
//   if (debug_) cout  << "[CandidateTkIsolation]::getPtTracksCorr::Sending back values" << endl;

  if (debug_) cout <<"tmpelpt1 (all)= "<<tmpElPt1<<" tmpelpt2 (Z)= "<<tmpElPt2<<endl;
  dummypT1=dummypT;
  dummypT2= dummypT;
  if (tmpElPt1 != 0.)  //Should be > ptLow_ ?
    if (tmpElPt1 <= dummypT) 
      dummypT1=  dummypT - tmpElPt1; 
  
  if (tmpElPt2 != 0.) 
    if (tmpElPt2 <= dummypT) 
      dummypT2=  dummypT - tmpElPt2;
  
  dummypTout.push_back(dummypT1);	//all electrons removed from cone
  dummypTout.push_back(dummypT2);	//only mass vetoed electrons removed from cone
  return dummypTout;

}


bool CandidateTkIsolation::isGoodTrack (edm::View<reco::Track>::const_iterator track) const
{
  double pt = sqrt( track->innerMomentum().Perp2() );  
  double dz;
  if (!isRC) dz = fabs( track->dz() - (electron_->gsfTrack()->dz()) ) ;
  else dz = fabs( track->dz()) ;
  if (pt < ptLow_) return false;
  if (dz > lip_) return false;
  return true;

}


bool CandidateTkIsolation::isGoodTrackerTrack(edm::View<reco::Track>::const_iterator track) const
{

  double pt  = sqrt( track->innerMomentum().Perp2() );
  double d0 = fabs(track->dxy(beamPoint_));
  double ed0 = track->dxyError();
  double edz_tr = track->dzError();
  double edz_ele;
  double dz;
  if (!isRC) {
    edz_ele= electron_->gsfTrack()->dzError();
    dz = fabs( track->dz(beamPoint_) - (electron_->gsfTrack()->dz(beamPoint_)));	    
  }
  else {
    edz_ele=0;
    dz = fabs( track->dz(beamPoint_));	    
  }
  double edz= sqrt(edz_tr*edz_tr + edz_ele*edz_ele);
  int nhits = track->recHitsSize();
  if (pt < 1) return false;
  if ( nhits < 8 ) {
    if ( d0 > 0.04 ) return false;
    if ( dz > 0.50 ) return false;
    if ( d0 / ed0 > 7 ) return false;
    if ( dz / edz > 7 ) return false;
  }
  else if ( nhits < 10 ) {
    if ( d0 > 0.20 ) return false;
    if ( dz > 2.00 ) return false;
    if ( d0 / ed0 > 10. ) return false;
    if ( dz / edz > 10. ) return false;
  }
  else {
    if ( d0 > 1.00 ) return false;
    if ( dz > 5.00 ) return false;
  }

  //    cout<<"quality true"<<endl;
  return true ;

}


bool CandidateTkIsolation::isGoodTrackerTrack(reco::TrackRef track) const
{

  double pt  = sqrt( track->innerMomentum().Perp2() );
  double d0 = fabs(track->dxy(beamPoint_));
  double ed0 = track->dxyError();
  double edz_tr = track->dzError();
  double edz_ele;
  double dz;
  if (!isRC) {
    edz_ele= electron_->gsfTrack()->dzError();
    dz = fabs( track->dz(beamPoint_) - (electron_->gsfTrack()->dz(beamPoint_)));	    
  }
  else {
    edz_ele=0;
    dz = fabs( track->dz(beamPoint_));	    
  }
  double edz= sqrt(edz_tr*edz_tr + edz_ele*edz_ele);
  int nhits = track->recHitsSize();
  if (pt < 1) return false;
  if ( nhits < 8 ) {
    if ( d0 > 0.04 ) return false;
    if ( dz > 0.50 ) return false;
    if ( d0 / ed0 > 7 ) return false;
    if ( dz / edz > 7 ) return false;
  }
  else if ( nhits < 10 ) {
    if ( d0 > 0.20 ) return false;
    if ( dz > 2.00 ) return false;
    if ( d0 / ed0 > 10. ) return false;
    if ( dz / edz > 10. ) return false;
  }
  else {
    if ( d0 > 1.00 ) return false;
    if ( dz > 5.00 ) return false;
  }

  //    cout<<"quality true"<<endl;
  return true ;

}


void CandidateTkIsolation::Initialize()
{
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  stripWidth_ = -999;
  ptLow_ = 1.5 ; 
  lip_ = 0.1 ; 
  Zmass_inf_= 40;	//I think these values should be turned upside down (inf and sup)
  Zmass_sup_= 10;
  ZMASSPDG_= 91.1876;
  isRC=false;
  debug_=false;
  return;
}
