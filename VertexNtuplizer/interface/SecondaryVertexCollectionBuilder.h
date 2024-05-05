#ifndef SECONDARY_VERTEX_COLLECTION_BUILDER
#define SECONDARY_VERTEX_COLLECTION_BUILDER


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SecondaryVertex.h"


typedef std::vector<SecondaryVertex> SecondaryVertexCollection;


class SecondaryVertexCollectionBuilder {

  public:

    SecondaryVertexCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~SecondaryVertexCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesToken,
        edm::EDGetTokenT<reco::VertexCollection> secondaryVerticesMTDTimingToken);

    SecondaryVertexCollection getSecondaryVertexCollection() { return secondaryVertices_; }
    SecondaryVertexCollection getSecondaryVertexCollectionMTDTiming() { return secondaryVerticesMTDTiming_; }

  private:

    bool goodSecondaryVertex(const reco::Vertex& sv);

    double trkTimeQualityCut_;

    SecondaryVertexCollection secondaryVertices_;
    SecondaryVertexCollection secondaryVerticesMTDTiming_;

    reco::VertexCollection cmsSecondaryVertices_;
    reco::VertexCollection cmsSecondaryVerticesMTDTiming_;
};


#endif
