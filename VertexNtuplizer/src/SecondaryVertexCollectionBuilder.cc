#include "../interface/SecondaryVertexCollectionBuilder.h"


SecondaryVertexCollectionBuilder::SecondaryVertexCollectionBuilder(const edm::ParameterSet& iConfig) {

  trkTimeQualityCut_ = iConfig.getUntrackedParameter<double>("trkTimeQualityCut");
}


// SecondaryVertexCollectionBuilder::~SecondaryVertexCollectionBuilder() {}


void SecondaryVertexCollectionBuilder::build(const edm::Event& iEvent,
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesToken,
    edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDTimingToken) {

  cmsSecondaryVertices_ = iEvent.get(secondaryVerticesToken);
  cmsSecondaryVerticesMTDTiming_ = iEvent.get(secondaryVerticesMTDTimingToken);

  secondaryVertices_.clear();
  secondaryVerticesMTDTiming_.clear();

  //
}


bool SecondaryVertexCollectionBuilder::goodSecondaryVertex(const reco::Vertex& sv) {

  bool pass = true;
  //
  return pass;
}
