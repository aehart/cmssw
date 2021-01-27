#include "RecoTracker/ExoRegionalTracking/interface/ExoRegion.h"

ExoRegion::ExoRegion(const edm::Handle<edm::View<reco::VertexCompositeCandidate> > &trackClusters,
                     const unsigned index,
                     const double rParam)
    : valid_(true),
      rParam_(rParam),
      sumOfCenters_(trackClusters->at(index).vx(), trackClusters->at(index).vy(), trackClusters->at(index).vz()) {
  const edm::Ref<edm::View<reco::VertexCompositeCandidate> > trackClusterRef(trackClusters, index);
  constituents_.push_back(trackClusterRef);
}

void ExoRegion::merge(const ExoRegion &other) {
  for (const auto &trackCluster : other.constituents())
    constituents_.push_back(trackCluster);
  sumOfCenters_ += other.sumOfCenters();
}

const bool ExoRegion::inRegion(const reco::TransientTrack &other, double &distance) const {
  distance = -1.0;
  const GlobalPoint point(vx(), vy(), vz());
  const auto &state = other.trajectoryStateClosestToPoint(point);
  if (!state.isValid())
    return false;
  const auto &pointOfClosestApproach = state.position();
  distance = fabs((pointOfClosestApproach - point).mag());

  return (distance < rParam());
}
