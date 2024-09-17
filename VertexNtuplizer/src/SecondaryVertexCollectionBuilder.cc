#include "../interface/SecondaryVertexCollectionBuilder.h"


SecondaryVertexCollectionBuilder::SecondaryVertexCollectionBuilder(const edm::ParameterSet& iConfig) {

  trkPtMin_ = iConfig.getUntrackedParameter<double>("trkPtMin");
  absEtaMax_ = iConfig.getUntrackedParameter<double>("absEtaMax");
  trkTimeQualityCut_ = iConfig.getUntrackedParameter<double>("trkTimeQualityCut");
  svChi2dofMax_ = iConfig.getUntrackedParameter<double>("svChi2dofMax");
}


// SecondaryVertexCollectionBuilder::~SecondaryVertexCollectionBuilder() {}


void SecondaryVertexCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& inclusiveSecondaryVerticesToken,
    // edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& inclusiveSecondaryVerticesMTDPVToken,
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& mergedSecondaryVerticesToken,
    // edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& mergedSecondaryVerticesMTDPVToken,
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& slimmedSecondaryVerticesToken,
    // edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>& slimmedSecondaryVerticesMTDPVToken,
    edm::EDGetTokenT<reco::TrackCollection>& generalTracksToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackT0FromBSToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackSigmaT0FromBSToken,
    // edm::EDGetTokenT<edm::ValueMap<float>>& trackQualityFromBSToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackT0FromPVToken,
    edm::EDGetTokenT<edm::ValueMap<float>>& trackSigmaT0FromPVToken,
    // edm::EDGetTokenT<edm::ValueMap<float>>& trackQualityFromPVToken,
    const reco::Vertex& primaryVertex) {

  iEvent.getByToken(generalTracksToken, generalTracksHandle_);

  inclusiveSecondaryVertices_ = iEvent.get(inclusiveSecondaryVerticesToken);
  // inclusiveSecondaryVerticesMTDPV_ = iEvent.get(inclusiveSecondaryVerticesMTDPVToken);
  mergedSecondaryVertices_ = iEvent.get(mergedSecondaryVerticesToken);
  // mergedSecondaryVerticesMTDPV_ = iEvent.get(mergedSecondaryVerticesMTDPVToken);
  slimmedSecondaryVertices_ = iEvent.get(slimmedSecondaryVerticesToken);
  // slimmedSecondaryVerticesMTDPV_ = iEvent.get(slimmedSecondaryVerticesMTDPVToken);
  trackT0FromBS_ = iEvent.get(trackT0FromBSToken);
  trackSigmaT0FromBS_ = iEvent.get(trackSigmaT0FromBSToken);
  // trackQualityFromBS_ = iEvent.get(trackQualityFromBSToken);
  trackT0FromPV_ = iEvent.get(trackT0FromPVToken);
  trackSigmaT0FromPV_ = iEvent.get(trackSigmaT0FromPVToken);
  // trackQualityFromPV_ = iEvent.get(trackQualityFromPVToken);

  secVerticesInclusive_.clear();
  secVerticesInclusiveMTDBS_.clear();
  secVerticesInclusiveMTDPV_.clear();
  secVerticesSlimmed_.clear();
  secVerticesSlimmedMTDBS_.clear();
  secVerticesSlimmedMTDPV_.clear();

  for (const reco::VertexCompositePtrCandidate& sv : inclusiveSecondaryVertices_) {
    if (!goodRecoVertex(sv)) continue;

    SecondaryVertex newSV(sv, primaryVertex);
    secVerticesInclusive_.push_back(newSV);

    SecondaryVertex newSVMTDBS(sv, sv, primaryVertex, generalTracksHandle_, trackT0FromBS_, trackSigmaT0FromBS_);
    // SecondaryVertex newSVMTDBS(sv, sv, primaryVertex, generalTracksHandle_, trackT0FromBS_, trackSigmaT0FromBS_, trackQualityFromBS_);
    secVerticesInclusiveMTDBS_.push_back(newSVMTDBS);

    SecondaryVertex newSVMTDPV(sv, sv, primaryVertex, generalTracksHandle_, trackT0FromPV_, trackSigmaT0FromPV_);
    // SecondaryVertex newSVMTDPV(sv, sv, primaryVertex, generalTracksHandle_, trackT0FromPV_, trackSigmaT0FromPV_, trackQualityFromPV_);
    secVerticesInclusiveMTDPV_.push_back(newSVMTDPV);
  }

  unsigned int nMergedSVs = mergedSecondaryVertices_.size();
  unsigned int nSlimmedSVs = slimmedSecondaryVertices_.size();
  if (nMergedSVs != nSlimmedSVs) std::cout << "WARNING: DIFFERENT NUMBER OF SVS!!!" << std::endl;

  for (unsigned int iSV = 0; iSV < nSlimmedSVs; iSV++) {
    const reco::VertexCompositePtrCandidate& mergedSV = mergedSecondaryVertices_.at(iSV);
    const reco::VertexCompositePtrCandidate& slimmedSV = slimmedSecondaryVertices_.at(iSV);

    unsigned int nMergedTrks = mergedSV.daughterPtrVector().size();
    unsigned int nSlimmedTrks = slimmedSV.daughterPtrVector().size();
    if (nMergedTrks != nSlimmedTrks) std::cout << "WARNING: DIFFERENT NUMBER OF TRACKS!!!" << std::endl;

    if (!goodRecoVertex(slimmedSV)) continue; // Is checking the merged SV important? 

    SecondaryVertex newSV(slimmedSV, primaryVertex);
    secVerticesSlimmed_.push_back(newSV);

    SecondaryVertex newSVMTDBS(slimmedSV, mergedSV, primaryVertex, generalTracksHandle_, trackT0FromBS_, trackSigmaT0FromBS_);
    // SecondaryVertex newSVMTDBS(sv, sv, primaryVertex, generalTracksHandle_, trackT0FromBS_, trackSigmaT0FromBS_, trackQualityFromBS_);
    secVerticesSlimmedMTDBS_.push_back(newSVMTDBS);

    SecondaryVertex newSVMTDPV(slimmedSV, mergedSV, primaryVertex, generalTracksHandle_, trackT0FromPV_, trackSigmaT0FromPV_);
    // SecondaryVertex newSVMTDPV(sv, sv, primaryVertex, generalTracksHandle_, trackT0FromPV_, trackSigmaT0FromPV_, trackQualityFromPV_);
    secVerticesSlimmedMTDPV_.push_back(newSVMTDPV);
  }

  // Sort collections
  // std::sort(secVerticesInclusive_.begin(), secVerticesInclusive_.end(), compare);
  // std::sort(secVerticesInclusiveMTDPV_.begin(), secVerticesInclusiveMTDPV_.end(), compare);
  // std::sort(secVerticesMerged_.begin(), secVerticesMerged_.end(), compare);
  // std::sort(secVerticesMergedMTDPV_.begin(), secVerticesMergedMTDPV_.end(), compare);
  // std::sort(secVerticesSlimmed_.begin(), secVerticesSlimmed_.end(), compare);
  // std::sort(secVerticesSlimmedMTDPV_.begin(), secVerticesSlimmedMTDPV_.end(), compare);
}


bool SecondaryVertexCollectionBuilder::goodRecoVertex(const reco::VertexCompositePtrCandidate& v) {

  bool vtxPass = true;
  if (abs(v.eta()) > absEtaMax_) vtxPass = false; // TODO: ETA DEFINITION?
  if (v.vertexNormalizedChi2() > svChi2dofMax_) vtxPass = false; // Take out poorly fitted vertices
  if (v.daughterPtrVector().size() < 2) vtxPass = false; // Not a vertex
  return vtxPass;
}


// bool SecondaryVertexCollectionBuilder::compare(const SecondaryVertex& sva, const SecondaryVertex& svb) {

//   return svb.d3d() < sva.d3d();
// }
