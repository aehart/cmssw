#pragma once

#include <list>
#include <vector>
#include <limits>
#include <string>
#include <atomic>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include "RecoTracker/ExoRegionalTracking/interface/ExoRegion.h"

using namespace std;

struct ExoRegionTFCache {
  ExoRegionTFCache() : graphDef(nullptr) {}

  atomic<tensorflow::GraphDef *> graphDef;
};

class ExoRegionProducer : public edm::stream::EDProducer<edm::GlobalCache<ExoRegionTFCache> > {
public:
  ExoRegionProducer(const edm::ParameterSet &, const ExoRegionTFCache *);
  ~ExoRegionProducer() override;
  static unique_ptr<ExoRegionTFCache> initializeGlobalCache(const edm::ParameterSet &);
  static void globalEndJob(const ExoRegionTFCache *);
  void produce(edm::Event &, const edm::EventSetup &) override;

private:
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<edm::View<reco::VertexCompositeCandidate> > trackClustersToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;

  // clustering parameters
  double rParam_;

  // selection parameters
  double minRadius_;
  double discriminatorCut_;
  string genMatch_;
  vector<string> input_names_;
  vector<string> output_names_;

  tensorflow::Session *session_;

  void getDecayVertices(const edm::View<reco::GenParticle> &, vector<math::XYZVector> &) const;
  const bool roiContainsDecayVertex(const math::XYZVector &, const vector<math::XYZVector> &, const double) const;
  const double getDiscriminatorValue(const ExoRegion &, const reco::BeamSpot &) const;
};
