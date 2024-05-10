#include "../interface/SecondaryVertexCollectionBuilder.h"


SecondaryVertexCollectionBuilder::SecondaryVertexCollectionBuilder(const edm::ParameterSet& iConfig) {

  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  trkPtMin_ = iConfig.getUntrackedParameter<double>("trkPtMin");
  trkTimeQualityCut_ = iConfig.getUntrackedParameter<double>("trkTimeQualityCut");
  svChi2dofMax_ = iConfig.getUntrackedParameter<double>("svChi2dofMax");
}


// SecondaryVertexCollectionBuilder::~SecondaryVertexCollectionBuilder() {}


void SecondaryVertexCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesToken,
    edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesMTDTimingToken,
    const reco::Vertex& primaryVertex,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeValueMapToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeErrorMapToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeQualityMapToken) {

  cmsSecondaryVertices_ = iEvent.get(secondaryVerticesToken);
  cmsSecondaryVerticesMTDTiming_ = iEvent.get(secondaryVerticesMTDTimingToken);
  trackTimeValueMap_ = iEvent.get(trackTimeValueMapToken);
  trackTimeErrorMap_ = iEvent.get(trackTimeErrorMapToken);
  trackTimeQualityMap_ = iEvent.get(trackTimeQualityMapToken);

  secondaryVertices_.clear();
  secondaryVerticesMTDTiming_.clear();

  for (const reco::Vertex& sv : cmsSecondaryVertices_) {
    if (sv.normalizedChi2() > svChi2dofMax_) continue; // Take out poorly fitted vertices
    // Count how many good tracks and check if this is a "real" vertex
    std::vector<reco::TrackBaseRef>* goodTracks = new std::vector<reco::TrackBaseRef>;
    for (const reco::TrackBaseRef& trkRef : sv.tracks()) {
      if (goodRecoTrack(trkRef)) goodTracks->push_back(trkRef);
    }
    if (goodTracks->size() < 2) continue; // Not a vertex

    SecondaryVertex newSV(sv, goodTracks, primaryVertex,
        trackTimeValueMap_, trackTimeErrorMap_, trackTimeQualityMap_);
    secondaryVertices_.push_back(newSV);
  } // End loop over original SV collection

  for (const reco::Vertex& sv : cmsSecondaryVerticesMTDTiming_) {
    if (sv.normalizedChi2() > svChi2dofMax_) continue; // Take out poorly fitted vertices
    // Count how many good tracks and check if this is a "real" vertex
    std::vector<reco::TrackBaseRef>* goodTracks = new std::vector<reco::TrackBaseRef>;
    for (const reco::TrackBaseRef& trkRef : sv.tracks()) {
      if (goodRecoTrack(trkRef)) goodTracks->push_back(trkRef);
    }
    if (goodTracks->size() < 2) continue; // Not a vertex

    SecondaryVertex newSV(sv, goodTracks, primaryVertex,
        trackTimeValueMap_, trackTimeErrorMap_, trackTimeQualityMap_);
    secondaryVerticesMTDTiming_.push_back(newSV);
  } // End loop over original SV collection with MTD timing

  // Sort collections
  // std::sort(secondaryVertices_.begin(), secondaryVertices_.end(), compare);
  // std::sort(secondaryVerticesMTDTiming_.begin(), secondaryVerticesMTDTiming_.end(), compare);
}


template <class T>
bool SecondaryVertexCollectionBuilder::goodRecoTrack(const T& trkRef) {

    bool trkPass = true;
    if (trkRef->pt() < trkPtMin_) trkPass = false;
    if (abs(trkRef->eta()) > absEtaMax_) trkPass = false;
    if (!trackTimeValueMap_.contains(trkRef.id())) trkPass = false;
    else {
      if (trackTimeErrorMap_[trkRef] == -1.0) trkPass = false;
      // if (trackTimeQualityMap_[trkRef] < trkTimeQualityCut_) trkPass = false;
    }
    return trkPass;
}


bool SecondaryVertexCollectionBuilder::compare(const SecondaryVertex& sva, const SecondaryVertex& svb) {

  return svb.d3d() < sva.d3d();
}
