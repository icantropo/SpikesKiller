#ifndef TdrHzzIsolationProducer_h
#define TdrHzzIsolationProducer_h

/**\class HZZ4LeptonsHadIsolationProducer
 *
 *
 * Original Author: Nicola De Filippis
 *
 * Compute isolation for cones around electron candidates
 */
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h" 
#include "FWCore/Framework/interface/DataKeyTags.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/plugins/CaloGeometryBuilder.h"

// using namespace edm;
 
class TdrHzzIsolationProducer : public edm::EDProducer {

 public:
  explicit TdrHzzIsolationProducer(const edm::ParameterSet&);
  ~TdrHzzIsolationProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag tracksTag_ ;
  double isoCone;
  double isoVeto;
  double isoVarCut;

  edm::InputTag electronTag_,electronVetoTag_,muonsTag_;
  std::string hcalrhitsLabel_;
  double radiusConeIntHad_,radiusConeExtHad_,eTMinHad_;

  edm::ESHandle<CaloGeometry> theCaloGeom_;
};

#endif
