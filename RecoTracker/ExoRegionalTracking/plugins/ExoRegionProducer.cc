#include "RecoTracker/ExoRegionalTracking/plugins/ExoRegionProducer.h"

ExoRegionProducer::ExoRegionProducer(const edm::ParameterSet &cfg, const ExoRegionTFCache *cache)
    : rParam_(cfg.getParameter<double>("rParam")),

      minRadius_(cfg.getParameter<double>("minRadius")),
      discriminatorCut_(cfg.getParameter<double>("discriminatorCut")),
      genMatch_(cfg.getParameter<string>("genMatch")),
      input_names_(cfg.getParameter<vector<string> >("input_names")),
      output_names_(cfg.getParameter<vector<string> >("output_names")) {

  beamSpotToken_ = consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"));
  trackClustersToken_ =
      consumes<edm::View<reco::VertexCompositeCandidate> >(cfg.getParameter<edm::InputTag>("trackClusters"));
  if (genMatch_ != "None")
    genParticlesToken_ = consumes<edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"));

  if (genMatch_ != "None" && genMatch_ != "HToSSTobbbb" && genMatch_ != "stopToBottom")
    throw cms::Exception("BadConfig") << "Parameter genMatch must be set to \"None\", \"HToSSTobbbb\", or \"stopToBottom\".";

  unsigned nThreads = cfg.getParameter<unsigned>("nThreads");
  string singleThreadPool = cfg.getParameter<string>("singleThreadPool");
  tensorflow::SessionOptions sessionOptions;
  tensorflow::setThreading(sessionOptions, nThreads, singleThreadPool);
  session_ = tensorflow::createSession(cache->graphDef, sessionOptions);

  produces<vector<reco::Vertex> >();
}

ExoRegionProducer::~ExoRegionProducer() {
  if (session_ != nullptr)
    tensorflow::closeSession(session_);
}

unique_ptr<ExoRegionTFCache> ExoRegionProducer::initializeGlobalCache(const edm::ParameterSet &cfg) {
  tensorflow::setLogging("3");
  string pbFile = cfg.getParameter<edm::FileInPath>("graph_path").fullPath();
  ExoRegionTFCache *cache = new ExoRegionTFCache();
  cache->graphDef = tensorflow::loadGraphDef(pbFile);

  return unique_ptr<ExoRegionTFCache>(cache);
}

void ExoRegionProducer::globalEndJob(const ExoRegionTFCache *cache) {
  if (cache->graphDef != nullptr)
    delete cache->graphDef;
}

void ExoRegionProducer::produce(edm::Event &event, const edm::EventSetup &setup) {
  edm::Handle<reco::BeamSpot> beamSpot;
  edm::Handle<edm::View<reco::VertexCompositeCandidate> > trackClusters;
  edm::Handle<edm::View<reco::GenParticle> > genParticles;

  event.getByToken(beamSpotToken_, beamSpot);
  const math::XYZVector bs(beamSpot->x0(), beamSpot->y0(), beamSpot->z0());
  event.getByToken(trackClustersToken_, trackClusters);
  vector<math::XYZVector> decayVertices;
  if (genMatch_ != "None") {
    event.getByToken(genParticlesToken_, genParticles);
    getDecayVertices(*genParticles, decayVertices);
  }

  // Initialize distances.
  list<ExoRegion> pseudoROIs;
  list<Distance> distances;
  const double minTrackClusterRadius = minRadius_ - rParam_;
  for (unsigned i = 0; i < trackClusters->size(); i++) {
    const reco::VertexCompositeCandidate &trackCluster = trackClusters->at(i);
    const math::XYZVector x(trackCluster.vx(), trackCluster.vy(), trackCluster.vz());
    if (minRadius_ < 0.0 || minTrackClusterRadius < 0.0 || (x - bs).rho() > minTrackClusterRadius)
      pseudoROIs.emplace_back(trackClusters, i, rParam_);
  }
  if (pseudoROIs.size() > 1) {
    ExoRegionItr secondToLast = pseudoROIs.end();
    secondToLast--;
    for (ExoRegionItr i = pseudoROIs.begin(); i != secondToLast; i++) {
      ExoRegionItr j = i;
      j++;
      for (; j != pseudoROIs.end(); j++)
        distances.emplace_back(i, j);
    }
  }

  // Track clusters farther apart than 4 times rParam_ cannot wind up in the
  // same ROI, so remove these pairs.
  const auto pred = [&](const Distance &a) { return a.distance() > 4.0 * rParam_; };
  distances.remove_if(pred);

  // Do clustering.
  while (!distances.empty()) {
    const auto comp = [](const Distance &a, const Distance &b) { return a.distance() <= b.distance(); };
    distances.sort(comp);
    DistanceItr d = distances.begin();
    if (d->distance() > rParam_)
      break;

    d->entities().first->merge(*d->entities().second);
    d->entities().second->setInvalid();

    const auto distancePred = [](const Distance &a) {
      return (!a.entities().first->valid() || !a.entities().second->valid());
    };
    const auto pseudoROIPred = [](const ExoRegion &a) { return !a.valid(); };
    distances.remove_if(distancePred);
    pseudoROIs.remove_if(pseudoROIPred);
  }

  // Select valid ROIs.
  const auto roiPred = [&](const ExoRegion &roi) {
    if (!roi.valid()) return true;
    const math::XYZVector x(roi.vx(), roi.vy(), roi.vz());
    if ((minRadius_ >= 0.0 && (x - bs).rho() < minRadius_) ||
        (genMatch_ != "None" && !roiContainsDecayVertex(x, decayVertices, roi.rParam())))
      return true;
    const double discriminatorValue = ((discriminatorCut_ > 0.0) ? getDiscriminatorValue(roi, *beamSpot) : 1.0);
    if (discriminatorValue < discriminatorCut_)
      return true;
    return false;
  };
  pseudoROIs.remove_if(roiPred);

  auto regionsOfInterest = make_unique<vector<reco::Vertex> >();
  regionsOfInterest->clear();

  const vector<double> error = {1.0, 0.0, 1.0, 0.0, 0.0, 1.0};

  for (const auto &roi : pseudoROIs)
    regionsOfInterest->emplace_back(reco::Vertex::Point(roi.vx(), roi.vy(), roi.vz()),
                                    reco::Vertex::Error(error.begin(), error.end(), true, true));

  event.put(move(regionsOfInterest));
}

void ExoRegionProducer::getDecayVertices(const edm::View<reco::GenParticle> &genParticles,
                                                vector<math::XYZVector> &decayVertices) const {
  for (const auto &genParticle : genParticles) {
    if (genMatch_ == "HToSSTobbbb") {
      if (genParticle.pdgId() != 5)
        continue;
      if (!genParticle.isHardProcess())
        continue;
      if (abs(genParticle.mother()->pdgId()) != 9000006)
        continue;
    } else if (genMatch_ == "stopToBottom") {
      if (abs(genParticle.pdgId()) != 5)
        continue;
      if (!genParticle.isHardProcess())
        continue;
      if (abs(genParticle.mother()->pdgId()) < 1000000)
        continue;
    }

    decayVertices.emplace_back(genParticle.vx(), genParticle.vy(), genParticle.vz());
  }
}

const bool ExoRegionProducer::roiContainsDecayVertex(const math::XYZVector &roi,
                                                            const vector<math::XYZVector> &decayVertices,
                                                            const double r) const {
  for (const auto &decayVertex : decayVertices)
    if ((roi - decayVertex).r() < r)
      return true;
  return false;
}

const double ExoRegionProducer::getDiscriminatorValue(const ExoRegion &roi, const reco::BeamSpot &bs) const {
  tensorflow::Tensor vertexTensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 40, 23}));
  auto vertex_map = vertexTensor.tensor<float, 3>();
  tensorflow::Tensor annulusTensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 10, 8}));
  auto annulus_map = annulusTensor.tensor<float, 3>();

  for (int i = 0, map_i = 0; map_i < 40; i++, map_i++) {
    if (i >= static_cast<int>(roi.nConstituents()))
      for (unsigned j = 0; j < 23; j++)
        vertex_map(0, map_i, j) = 0.0;
    else {
      const auto &trackClusterRef = roi.constituent(i);
      const auto &track0 = *trackClusterRef->daughter(0)->bestTrack();
      const auto &track1 = *trackClusterRef->daughter(1)->bestTrack();

      vertex_map(0, map_i, 0) = trackClusterRef->vx() - bs.x0();
      vertex_map(0, map_i, 1) = trackClusterRef->vy() - bs.y0();
      vertex_map(0, map_i, 2) = trackClusterRef->vz() - bs.z0();

      vertex_map(0, map_i, 3) = trackClusterRef->vertexCovariance()(0, 0);
      vertex_map(0, map_i, 4) = trackClusterRef->vertexCovariance()(0, 1);
      vertex_map(0, map_i, 5) = trackClusterRef->vertexCovariance()(0, 2);
      vertex_map(0, map_i, 6) = trackClusterRef->vertexCovariance()(1, 1);
      vertex_map(0, map_i, 7) = trackClusterRef->vertexCovariance()(1, 2);
      vertex_map(0, map_i, 8) = trackClusterRef->vertexCovariance()(2, 2);

      vertex_map(0, map_i, 9) = track0.charge() * track0.pt();
      vertex_map(0, map_i, 10) = track0.eta();
      vertex_map(0, map_i, 11) = track0.phi();
      vertex_map(0, map_i, 12) = track0.dxy(bs);
      vertex_map(0, map_i, 13) = track0.dz(reco::TrackBase::Point(bs.x0(), bs.y0(), bs.z0()));
      vertex_map(0, map_i, 14) = track0.normalizedChi2();
      vertex_map(0, map_i, 15) = track0.quality(reco::Track::highPurity) ? 1 : 0;

      vertex_map(0, map_i, 16) = track1.charge() * track1.pt();
      vertex_map(0, map_i, 17) = track1.eta();
      vertex_map(0, map_i, 18) = track1.phi();
      vertex_map(0, map_i, 19) = track1.dxy(bs);
      vertex_map(0, map_i, 20) = track1.dz(reco::TrackBase::Point(bs.x0(), bs.y0(), bs.z0()));
      vertex_map(0, map_i, 21) = track1.normalizedChi2();
      vertex_map(0, map_i, 22) = track1.quality(reco::Track::highPurity) ? 1 : 0;
    }
  }

  for (int i = 0; i < 10; i++)
    for (unsigned j = 0; j < 8; j++)
      annulus_map(0, i, j) = 0.0;

  tensorflow::NamedTensorList input_tensors;
  input_tensors.resize(2);
  input_tensors[0] = tensorflow::NamedTensor(input_names_.at(0), vertexTensor);
  input_tensors[1] = tensorflow::NamedTensor(input_names_.at(1), annulusTensor);
  vector<tensorflow::Tensor> outputs;
  tensorflow::run(session_, input_tensors, output_names_, &outputs);

  return (outputs.at(0).flat<float>()(1));
}

const double Distance::distance() const {
  if (entities_.first->valid() && entities_.second->valid())
    return (entities_.first->centerOfMass() - entities_.second->centerOfMass()).r();
  return numeric_limits<double>::max();
}

DEFINE_FWK_MODULE(ExoRegionProducer);
