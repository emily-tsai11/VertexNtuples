#include "../interface/SecondaryVertexCollectionBuilder.h"


SecondaryVertexCollectionBuilder::SecondaryVertexCollectionBuilder(const edm::ParameterSet& iConfig) {

  trkPtMin_ = iConfig.getUntrackedParameter<double>("trkPtMin");
  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  trkTimeQualityCut_ = iConfig.getUntrackedParameter<double>("trkTimeQualityCut");
  svChi2dofMax_ = iConfig.getUntrackedParameter<double>("svChi2dofMax");
}


// SecondaryVertexCollectionBuilder::~SecondaryVertexCollectionBuilder() {}


void SecondaryVertexCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesToken,
    edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesMTDBSToken,
    edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesMTDBS4Token,
    edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesMTDPVToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeBSValueMapToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeBSErrorMapToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeBSQualityMapToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackTimePVValueMapToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackTimePVErrorMapToken,
    // edm::EDGetTokenT<edm::ValueMap<float>>& trackTimePVQualityMapToken,
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& slimmedCandSVToken,
    const reco::Vertex& primaryVertex) {

  cmsSecondaryVertices_ = iEvent.get(secondaryVerticesToken);
  cmsSecondaryVerticesMTDBS_ = iEvent.get(secondaryVerticesMTDBSToken);
  cmsSecondaryVerticesMTDBS4_ = iEvent.get(secondaryVerticesMTDBS4Token);
  cmsSecondaryVerticesMTDPV_ = iEvent.get(secondaryVerticesMTDPVToken);
  trackTimeBSValueMap_ = iEvent.get(trackTimeBSValueMapToken);
  trackTimeBSErrorMap_ = iEvent.get(trackTimeBSErrorMapToken);
  trackTimeBSQualityMap_ = iEvent.get(trackTimeBSQualityMapToken);
  trackTimePVValueMap_ = iEvent.get(trackTimePVValueMapToken);
  trackTimePVErrorMap_ = iEvent.get(trackTimePVErrorMapToken);
  // trackTimePVQualityMap_ = iEvent.get(trackTimePVQualityMapToken);
  cmsSlimmedCandSVs_ = iEvent.get(slimmedCandSVToken);

  secondaryVertices_.clear();
  secondaryVerticesMTDBS_.clear();
  secondaryVerticesMTDBS4_.clear();
  secondaryVerticesMTDPV_.clear();
  slimmedCandSVs_.clear();

  for (const reco::Vertex& sv : cmsSecondaryVertices_) {
    if (!goodRecoVertex(sv)) continue;
    SecondaryVertex newSV(sv, primaryVertex);
    secondaryVertices_.push_back(newSV);
  }

  for (const reco::Vertex& sv : cmsSecondaryVerticesMTDBS_) {
    if (!goodRecoVertex(sv)) continue;
    SecondaryVertex newSV(sv, primaryVertex, trackTimeBSValueMap_, trackTimeBSErrorMap_, trackTimeBSQualityMap_);
    secondaryVerticesMTDBS_.push_back(newSV);
  }

  for (const reco::Vertex& sv : cmsSecondaryVerticesMTDBS4_) {
    if (!goodRecoVertex(sv)) continue;
    SecondaryVertex newSV(sv, primaryVertex, trackTimeBSValueMap_, trackTimeBSErrorMap_, trackTimeBSQualityMap_);
    secondaryVerticesMTDBS4_.push_back(newSV);
  }

  for (const reco::Vertex& sv : cmsSecondaryVerticesMTDPV_) {
    if (!goodRecoVertex(sv)) continue;
    // SecondaryVertex newSV(sv, primaryVertex, trackTimePVValueMap_, trackTimePVErrorMap_, trackTimePVQualityMap_);
    SecondaryVertex newSV(sv, primaryVertex, trackTimePVValueMap_, trackTimePVErrorMap_, trackTimeBSQualityMap_); // Temporary, quality anyways not used yet
    secondaryVerticesMTDPV_.push_back(newSV);
  }

  for (const reco::VertexCompositePtrCandidate& sv : cmsSlimmedCandSVs_) {
    if (!goodRecoVertex(sv)) continue;
    SecondaryVertex newSV(sv, primaryVertex);
    slimmedCandSVs_.push_back(newSV);
  }

  // Sort collections
  // std::sort(secondaryVertices_.begin(), secondaryVertices_.end(), compare);
  // std::sort(secondaryVerticesMTDBS_.begin(), secondaryVerticesMTDBS_.end(), compare);
  // std::sort(secondaryVerticesMTDBS4_.begin(), secondaryVerticesMTDBS4_.end(), compare);
  // std::sort(secondaryVerticesMTDPV_.begin(), secondaryVerticesMTDPV_.end(), compare);
  // std::sort(slimmedCandSVs_.begin(), slimmedCandSVs_.end(), compare);
}


bool SecondaryVertexCollectionBuilder::goodRecoVertex(const reco::Vertex& v) {

  bool vtxPass = true;
  if (abs(v.p4().Eta()) > absEtaMax_) vtxPass = false;
  if (v.normalizedChi2() > svChi2dofMax_) vtxPass = false; // Take out poorly fitted vertices
  if (v.tracksSize() < 2) vtxPass = false; // Not a vertex
  return vtxPass;
}


bool SecondaryVertexCollectionBuilder::goodRecoVertex(const reco::VertexCompositePtrCandidate& v) {

  bool vtxPass = true;
  if (abs(v.eta()) > absEtaMax_) vtxPass = false;
  if (v.vertexNormalizedChi2() > svChi2dofMax_) vtxPass = false; // Take out poorly fitted vertices
  if (v.daughterPtrVector().size() < 2) vtxPass = false; // Not a vertex
  return vtxPass;
}


// bool SecondaryVertexCollectionBuilder::compare(const SecondaryVertex& sva, const SecondaryVertex& svb) {

//   return svb.d3d() < sva.d3d();
// }
