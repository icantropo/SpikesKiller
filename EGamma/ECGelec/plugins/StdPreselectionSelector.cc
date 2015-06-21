#include "EGamma/ECGelec/plugins/StdPreselectionSelector.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

bool StdPreselectionSelector::selection(edm::Ref<reco::GsfElectronCollection> eleRef)
{
  
  bool eg = eleRef->core()->ecalDrivenSeed();
  bool pf = eleRef->core()->trackerDrivenSeed() && !eleRef->core()->ecalDrivenSeed() ;

  // Or MVA
  if ( (eg && (eleRef->mva()>=-0.4)) || (pf && (eleRef->mva()>=-0.4)) ) return true ;

  // Et cut
  if (eg && eleRef->isEB() && (eleRef->superCluster()->energy()/cosh(eleRef->superCluster()->eta()) < 4.)) return false;
  if (eg && eleRef->isEE() && (eleRef->superCluster()->energy()/cosh(eleRef->superCluster()->eta()) < 4.)) return false;
  if (pf && eleRef->isEB() && (eleRef->superCluster()->energy()/cosh(eleRef->superCluster()->eta()) < 0.)) return false;
  if (pf && eleRef->isEE() && (eleRef->superCluster()->energy()/cosh(eleRef->superCluster()->eta()) < 0.)) return false;

  // E/p cut
  if (eg && eleRef->isEB() && (eleRef->eSuperClusterOverP() > 999999999.)) return false;
  if (eg && eleRef->isEE() && (eleRef->eSuperClusterOverP() > 999999999.)) return false;
  if (eg && eleRef->isEB() && (eleRef->eSuperClusterOverP() < 0.)) return false;
  if (eg && eleRef->isEE() && (eleRef->eSuperClusterOverP() < 0.)) return false;
  if (pf && eleRef->isEB() && (eleRef->eSuperClusterOverP() > 999999999.)) return false;
  if (pf && eleRef->isEE() && (eleRef->eSuperClusterOverP() > 999999999.)) return false;
  if (pf && eleRef->isEB() && (eleRef->eSuperClusterOverP() < 0.)) return false;
  if (pf && eleRef->isEE() && (eleRef->eSuperClusterOverP() < 0.)) return false;

  // HoE cuts
  double had = eleRef->hcalOverEcal()*eleRef->superCluster()->energy();
  bool HoEveto = false;
  const reco::CaloCluster & seedCluster = *(eleRef->superCluster()->seed()) ;
  int detector = seedCluster.hitsAndFractions()[0].first.subdetId() ;
  if (detector==EcalBarrel && (had<0. || (had/eleRef->superCluster()->energy())<0.15)) HoEveto=true;
  else if (detector==EcalEndcap && (had<0. || (had/eleRef->superCluster()->energy())<0.15)) HoEveto=true;
  if ( eg && !HoEveto ) return false;
  bool HoEvetoPflow = false;
  if (detector==EcalBarrel && (had<0. || (had/eleRef->superCluster()->energy())<999999999.)) HoEvetoPflow=true;
  else if (detector==EcalEndcap && (had<0. || (had/eleRef->superCluster()->energy())<999999999.)) HoEvetoPflow=true;
  if ( pf && !HoEvetoPflow ) return false;

  // delta eta criteria
  double deta = eleRef->deltaEtaSuperClusterTrackAtVtx();
  if (eg && eleRef->isEB() && (fabs(deta) > 0.02)) return false;
  if (eg && eleRef->isEE() && (fabs(deta) > 0.02)) return false;
  if (pf && eleRef->isEB() && (fabs(deta) > 999999999.)) return false;
  if (pf && eleRef->isEE() && (fabs(deta) > 999999999.)) return false;

  // delta phi criteria
  double dphi = eleRef->deltaPhiSuperClusterTrackAtVtx();
  if (eg && eleRef->isEB() && (fabs(dphi) > 0.15)) return false;
  if (eg && eleRef->isEE() && (fabs(dphi) > 0.15)) return false;
  if (pf && eleRef->isEB() && (fabs(dphi) > 999999999.)) return false;
  if (pf && eleRef->isEE() && (fabs(dphi) > 999999999.)) return false;

  //
  if (eg && eleRef->isEB() && (eleRef->sigmaIetaIeta() > 999999999.)) return false;
  if (eg && eleRef->isEE() && (eleRef->sigmaIetaIeta() > 999999999.)) return false;
  if (pf && eleRef->isEB() && (eleRef->sigmaIetaIeta() > 999999999.)) return false;
  if (pf && eleRef->isEE() && (eleRef->sigmaIetaIeta() > 999999999.)) return false;

  return true ;
    
}
