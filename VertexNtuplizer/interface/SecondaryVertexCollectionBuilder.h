#ifndef SECONDARY_VERTEX_COLLECTION_BUILDER
#define SECONDARY_VERTEX_COLLECTION_BUILDER


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SecondaryVertex.h"


typedef std::vector<SecondaryVertex> SecondaryVertexCollection;


class SecondaryVertexCollectionBuilder {

  public:

    SecondaryVertexCollectionBuilder(const edm::ParameterSet& iConfig);
    // ~SecondaryVertexCollectionBuilder();

    void build(const edm::Event& iEvent,
        edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesToken,
        edm::EDGetTokenT<reco::VertexCollection>& secondaryVerticesMTDTimingToken,
        const reco::Vertex& primaryVertex,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeValueMapToken,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeErrorMapToken,
        edm::EDGetTokenT<edm::ValueMap<float>>& trackTimeQualityMapToken);

    SecondaryVertexCollection getSecondaryVertexCollection() { return secondaryVertices_; }
    SecondaryVertexCollection getSecondaryVertexCollectionMTDTiming() { return secondaryVerticesMTDTiming_; }

  private:

    template <class T> bool goodRecoTrack(const T& trkRef);

    double absEtaMax_;
    double trkPtMin_;
    double trkTimeQualityCut_;
    double svChi2dofMax_;

    SecondaryVertexCollection secondaryVertices_;
    SecondaryVertexCollection secondaryVerticesMTDTiming_;

    reco::VertexCollection cmsSecondaryVertices_;
    reco::VertexCollection cmsSecondaryVerticesMTDTiming_;
    edm::ValueMap<float> trackTimeValueMap_;
    edm::ValueMap<float> trackTimeErrorMap_;
    edm::ValueMap<float> trackTimeQualityMap_;
};


#endif
